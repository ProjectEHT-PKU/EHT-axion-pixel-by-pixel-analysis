"""Explicit eight-phase performance pilot, no full batch or extra snapshots."""
import json
import os
from pathlib import Path
import sys
import threading
import time

ROOT = Path('/home/canoziia/MyProject/PHYS/eht-script-pipeline')
RUN = Path('/home/canoziia/GRMHD_DATA/production/pilot-m65-a030-t1600-8x6-pinned')
sys.path.insert(0, str(ROOT));sys.path.insert(0, str(ROOT/'scripts'))
from _producer import options, save_request
from eht_axion.producer.plan import plan_generation
from eht_axion.producer.jobs import run_jobs

args = options('generate', ['--run-root',str(RUN),'--workers','8','--threads','6',
                          '--cpu-affinity','physical','--timeout-seconds','1800'])
plan = plan_generation(args.config_data,args.fluid_root,
                       {'low':args.radial_low,'high':args.radial_high},RUN/'build')
selected = [t for t in plan['tasks'] if t['kind']=='response' and t['alpha']==.3 and
            t['label']==1600 and t['phase_index']<8]
assert len(selected)==8
plan['tasks']=selected
plan['input_roots']=[str(args.fluid_root),str(args.radial_low.parent),str(args.radial_high.parent)]
plan['source_roots']=[str(args.source_root)]
plan['expected_core_commit']=args.expected_commit
save_request(args,plan)
metrics=RUN/'benchmark';metrics.mkdir(exist_ok=False)
stop=threading.Event();samples=[]
HZ=os.sysconf('SC_CLK_TCK')
start=time.perf_counter()


def descendants(pid):
    pending=[pid];found=set()
    while pending:
        parent=pending.pop()
        # Python worker threads own their Popen children on Linux, so read all tids.
        try: tids=list(Path(f'/proc/{parent}/task').iterdir())
        except OSError: continue
        for tid in tids:
            try: children=(tid/'children').read_text().split()
            except OSError: continue
            for raw in children:
                child=int(raw)
                if child not in found:found.add(child);pending.append(child)
    return found


def sample_process(pid):
    try:
        raw=Path(f'/proc/{pid}/stat').read_text()
        if raw[raw.index('(')+1:raw.rindex(')')]!='RAPTOR':return None
        fields=raw[raw.rindex(')')+2:].split()
        status=dict(line.split(':',1) for line in Path(f'/proc/{pid}/status').read_text().splitlines())
        tids=list(Path(f'/proc/{pid}/task').iterdir())
        masks=[sorted(os.sched_getaffinity(int(tid.name))) for tid in tids]
        return {'pid':pid,'cpu_seconds':(int(fields[11])+int(fields[12]))/HZ,
                'rss_kib':int(status.get('VmRSS','0 kB').split()[0]),
                'hwm_kib':int(status.get('VmHWM','0 kB').split()[0]),
                'threads':len(tids),'thread_cpu_sets':masks}
    except (OSError,ValueError,ProcessLookupError):return None


def monitor():
    with (metrics/'samples.jsonl').open('w') as f:
        while not stop.is_set():
            row={'elapsed_seconds':time.perf_counter()-start,'processes':[]}
            for pid in descendants(os.getpid()):
                x=sample_process(pid)
                if x is not None:row['processes'].append(x)
            row['aggregate_rss_kib']=sum(x['rss_kib'] for x in row['processes'])
            samples.append(row);f.write(json.dumps(row)+'\n');f.flush()
            stop.wait(.25)


worker=threading.Thread(target=monitor,daemon=True);worker.start()
try:
    result=run_jobs(plan,RUN/'build',RUN,workers=8,threads=6,cpu_affinity='physical',
                    timeout_seconds=1800)
finally:
    elapsed=time.perf_counter()-start;stop.set();worker.join(3)
(metrics/'result.json').write_text(json.dumps(result,indent=2)+'\n')
by_pid={}
for row in samples:
    for p in row['processes']:
        by_pid.setdefault(p['pid'],[]).append(p)
records=[]
for row in result['jobs']:
    if row.get('attempt'):
        record=json.loads((Path(row['attempt'])/'result.json').read_text());pid=record.get('pid')
        observations=by_pid.get(pid,[]);requested=set(record.get('cpu_set') or [])
        observed={cpu for s in observations for mask in s['thread_cpu_sets'] for cpu in mask}
        # Samples can catch taskset-to-RAPTOR startup before OpenMP applies binding;
        # observed masks must be confined to requested worker CPU set.
        records.append({'task_id':row['id'],'pid':pid,'status':row['status'],
                        'wall_seconds':record.get('process_wall_seconds'),
                        'requested_cpu_set':sorted(requested),'observed_cpu_union':sorted(observed),
                        'affinity_verified':bool(observations) and observed.issubset(requested),
                        'peak_rss_kib':max((s['hwm_kib'] for s in observations),default=0),
                        'sampled_cpu_seconds':max((s['cpu_seconds'] for s in observations),default=0),
                        'max_threads':max((s['threads'] for s in observations),default=0)})
cpu=sum(x['sampled_cpu_seconds'] for x in records)
summary={'batch_wall_seconds':elapsed,'workers':8,'threads_per_worker':6,'jobs':records,
         'counts':result['counts'],'images_per_hour':8/elapsed*3600,
         'sampled_cpu_seconds':cpu,'average_cpu_cores':cpu/elapsed,
         'aggregate_peak_rss_kib':max((s['aggregate_rss_kib'] for s in samples),default=0),
         'observed_max_concurrent_processes':max((len(s['processes']) for s in samples),default=0),
         'all_affinity_verified':all(x['affinity_verified'] for x in records),
         'scope':'8 phase0..7 response images, alpha.3,snapshot1600,100x100,norm.1,M6.5; no400/nootheralpha',
         'cpu_accounting':'/proc sampled; per-process CPU lower estimates omit final fractional interval'}
(metrics/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
print(json.dumps(summary,indent=2),flush=True)
raise SystemExit(0 if result['complete'] else 1)
