"""Independent dense-grid check of the fixed reassessment candidates at alpha=.3."""
import argparse
import csv
import inspect
import json
from pathlib import Path

import _bootstrap  # noqa: F401
import numpy as np
from eht_axion.pipeline import common
from eht_axion.research import background_validation, joint_days
from eht_axion.research.time_sampling import H2, H3
from reassess_background_study import build_cases, fit_backgrounds, load_npz, validate_inputs


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__)
    for name in ['inputs', 'reference', 'output-root']:
        p.add_argument('--'+name, type=Path, required=True)
    a = p.parse_args(argv)
    with common.verification_session():
        common.verify_stage(a.inputs, recursive=False, allow_imported=False)
        common.verify_stage(a.reference, recursive=False, allow_imported=False)
        paths = [a.inputs/f'{kind}-{dim}.npz' for kind in ['states', 'observations']
                 for dim in ['1d', '2d']]
        common.require_declared_outputs(a.inputs, paths)
        reference = a.reference/'summary.csv'
        common.require_declared_outputs(a.reference, [reference])
        with reference.open() as f:
            records = {r['key']:r for r in csv.DictReader(f) if float(r['alpha']) == .3}
        receipt_root = a.output_root.parent/(a.output_root.name+'-source-receipts')
        receipt_root.mkdir(exist_ok=False)
        sources = paths+[reference, Path(__file__).resolve(), Path(inspect.getfile(build_cases)),
                         Path(inspect.getfile(background_validation)), Path(inspect.getfile(joint_days))]
        for label, root in [('inputs', a.inputs), ('reference', a.reference)]:
            q = receipt_root/(label+'-receipt.json')
            q.write_bytes(common.manifest_path(root).read_bytes())
            sources.append(q)
        rows = []
        identity = dict(alpha=.3, phase_count=128, coupling_points=240, no_model_changes=True)
        with common.stage_directory(a.output_root, 'research:background-reassessment-grid',
                                    identity, sources) as out:
            for dim in ['1d', '2d']:
                s, obs = [load_npz(a.inputs/f'{kind}-{dim}.npz') for kind in ['states', 'observations']]
                validate_inputs(s, obs)
                for h, cols in [(H2, [0, 2]), (H3, [0, 1, 2])]:
                    fits = fit_backgrounds(s, cols)
                    for bg, response, kwargs in build_cases(s, .3, h, fits):
                        for target in ['observed', 'matched_zero']:
                            means = s['observed_daily_local_unwrapped'] if target == 'observed' else (
                                np.zeros_like(obs['means']))
                            r = joint_days.calculate_joint_days(means, obs['covariances'], h,
                                couplings=np.logspace(np.log10(1/137), 3, 240),
                                phases=np.arange(128)*2*np.pi/128, **kwargs)
                            # Exact key constructed by the source runner is matched by fields.
                            old = next(v for v in records.values() if v['dimension'] == dim and
                                int(v['contrasts']) == len(h) and v['background'] == bg and
                                v['response'] == response and v['mode'] == target)
                            key = old['key']
                            np.savez(out/(key+'.npz'), couplings=r['couplings'], cdf=r['cdf'],
                                density=r['density'], loglike=r['loglike'], c95=r['c95'])
                            rows.append(dict(key=key, dimension=dim, contrasts=len(h), background=bg,
                                response=response, mode=target, c95=r['c95'],
                                coarse_c95=float(old['c95']),
                                delta_dex=float(np.log10(r['c95']/float(old['c95'])))))
            with (out/'comparison.csv').open('w') as f:
                w = csv.DictWriter(f, fieldnames=list(rows[0]))
                w.writeheader()
                w.writerows(rows)
            (out/'report.json').write_text(json.dumps(identity, indent=2)+'\n')


if __name__ == '__main__':
    main()
