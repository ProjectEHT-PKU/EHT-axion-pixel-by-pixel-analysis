#!/usr/bin/env python3
"""EHT imaging: summarize the D-term Monte Carlo images of eht_sample.py into EHT ensembles.

The ensembles are the EHT input of the analysis (analysis_limits.py).

pixi run -e eht-summary eht-summarize --sample-root MC --output-root OUT
MC is the output root of eht_sample.py; its batch must be finished and must have been
made with the science values of the current experiment.yaml. Draw count and sky grid
come from experiment.yaml. Only the blurred images (blurred with
shared.image.blur_fwhm_uas, like the simulated images) are summarized.
Outputs: day-D/{stokes.npy,evpa.npy,global.npz},
axes.npz, report.json (conventions, registration, scope).
"""
from pathlib import Path
import argparse
import hashlib
import json

import numpy as np

from _bootstrap import ROOT
from _data_root import add_data_root_option, fill_defaults
from eht_axion import provenance
from eht_axion.eht import ensemble as en
from eht_axion.eht.ensemble import load_stokes_fits

contained = provenance.contained_path


def read_json(path):
    return json.loads(Path(path).read_text())


def summary_settings(config):
    """Values from experiment.yaml used by this summary (recorded in its manifest)."""
    image = config['shared']['image']
    return dict(draws=config['eht_imaging']['uncertainty_draws'], shape=image['shape'],
                fov_uas=image['fov_uas'])


def verify_batch(root, settings, science):
    """Verify a finished eht_sample.py output root before any image is read.

    batch.json must carry its own identity (SHA-256 of the plan) and the science values
    ``science`` of the current experiment; finish.json must list one certified receipt
    per draw; every consumed receipt, image and certificate is checked by SHA-256.
    """
    results = Path(root).resolve(strict=True)
    count = settings['draws']
    paths = []

    def checked(path, digest=None):
        actual = provenance.sha256(path)
        if digest is not None and actual != digest:
            raise ValueError('SHA mismatch: ' + str(path))
        paths.append(path)
        return read_json(path)

    batch = checked(results / 'batch.json')
    finish = checked(results / 'finish.json')
    identity = batch.get('identity')
    unsigned = dict(batch)
    unsigned.pop('identity', None)
    encoded = json.dumps(unsigned, sort_keys=True, separators=(',', ':'), allow_nan=False).encode()
    if (batch.get('schema') != 'eht-sample-batch-v2' or hashlib.sha256(encoded).hexdigest() != identity
            or finish.get('schema') != 'eht-sample-finish-v2' or finish.get('identity') != identity
            or finish.get('numerically_certified') != count * len(en.DAYS)):
        raise ValueError('incomplete or inconsistent eht_sample batch (batch.json / finish.json)')
    if batch.get('science') != json.loads(json.dumps(science)):
        raise ValueError('eht_sample batch was made with different experiment values')
    draws = checked(results / 'all-draws.json')
    taskpath = results / 'tasks.jsonl'
    provenance.sha256(taskpath)
    paths.append(taskpath)
    tasks = [json.loads(line) for line in taskpath.read_text().splitlines()]
    expected = {f'{d}-{k:04d}' for d in en.DAYS for k in range(count)}
    if (len(tasks) != count * len(en.DAYS) or {t['id'] for t in tasks} != expected
            or set(finish['receipts']) != expected):
        raise ValueError(f'not exactly {count} unique draws/day')
    selected = {d: {} for d in en.DAYS}
    for task in tasks:
        d, k = task['day'], task['draw']
        if (d not in en.DAYS or not isinstance(k, int) or not 0 <= k < count
                or task['id'] != f'{d}-{k:04d}' or task['dterms'] != draws[d]['draws'][k]):
            raise ValueError('draw/task mismatch')
        entry = finish['receipts'][task['id']]
        directory = contained(results, entry['directory'])
        receipt = checked(directory / 'taskreceipt.json', entry['sha256'])
        if (receipt['task'] != task or receipt['identity'] != identity
                or receipt['status'] != 'certified' or receipt['schema'] != 'eht-sample-task-v2'):
            raise ValueError('selected receipt identity/status mismatch')
        products = {}
        for kind, name in [('blurred', 'best-blurred-IQU.fits'),
                           ('certificate', 'certificate.json')]:
            p = contained(directory, name)
            if provenance.sha256(p) != receipt['outputs'][name]:
                raise ValueError('consumed output SHA mismatch: ' + str(p))
            paths.append(p)
            products[kind] = p
        cert = read_json(products['certificate'])
        solver = cert['solver']
        if not (solver['numerically_certified'] and -1e-10 <= solver['best_gap'] <= cert['gap_atol']):
            raise ValueError('invalid numerical certificate')
        selected[d][k] = products
    return selected, list(dict.fromkeys(paths)), identity


def run(selected, out, settings):
    draws, n = settings['draws'], settings['shape']
    pixel = settings['fov_uas'] / n
    reference = None
    for d in en.DAYS:
        directory = out / ('day-' + d)
        directory.mkdir(parents=True)
        samples = np.lib.format.open_memmap(directory / 'stokes.npy', mode='w+',
                                           dtype='float64', shape=(draws, n, n, 4))
        for k in range(draws):
            loaded = load_stokes_fits(selected[d][k]['blurred'])
            east, north = loaded['east_uas'], loaded['north_uas']
            if (loaded['stokes'].shape != (n, n, 4)
                    or not np.allclose(np.diff(east), -pixel, rtol=0, atol=1e-9)
                    or not np.allclose(np.diff(north), pixel, rtol=0, atol=1e-9)):
                raise ValueError('FITS grid differs from shared.image (E decreasing, N increasing)')
            geometry = (east, north, loaded['frequency_hz'], loaded['bunit'])
            if reference is None:
                reference = geometry
                np.savez(out / 'axes.npz', east_uas=east, north_uas=north)
            elif (not np.array_equal(east, reference[0]) or not np.array_equal(north, reference[1])
                  or geometry[2:] != reference[2:]):
                raise ValueError('inconsistent FITS axes/frequency/units')
            samples[k] = loaded['stokes']
        samples.flush()
        # The Paper VII Sec. 5.2 alignment is a zero shift here: images used as reconstructed.
        en.summarize(samples, *reference[:2], directory)
        print(f'SUMMARIZED day-{d}', flush=True)
        del samples
    report = dict(schema='eht-ensemble-v2', days=list(en.DAYS), draws_per_day=draws,
                  map_order='draw,y,x,IQUV; native FITS E decreases along x, N increases along y',
                  angle_convention='EVPA radians modulo pi, east of north',
                  registration='none: the Paper VII Sec. 5.2 alignment is a zero shift; '
                               'images used as reconstructed',
                  resolution='blurred products of eht_sample.py (shared.image.blur_fwhm_uas); '
                             'no blur applied here',
                  evpa_mask='EVPA NaN where I <= 0 or Q+iU = 0',
                  scope='conditional fixed-I D-term MC only; no thermal bootstrap, response variance, '
                        'covariance-estimation uncertainty or scientific acceptance certification',
                  beta2='Paper VII eq. 14 with rho_min 0 and rho_max the full field: '
                        'sum((Q+iU)*exp(-2i atan2(E,N)))/sum(I) about the grid centre',
                  frequency_hz=reference[2], bunit=reference[3])
    (out / 'report.json').write_text(json.dumps(report, indent=2, allow_nan=False) + '\n')


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--sample-root', type=Path,
                   help='finished eht_sample.py output root (default: eht/sample)')
    p.add_argument('--output-root', type=Path, help='new directory (default: eht/ensemble)')
    p.add_argument('--config', type=Path, default=ROOT / 'configs' / 'experiment.yaml',
                   help='default: configs/experiment.yaml')
    add_data_root_option(p)
    args = p.parse_args(argv)
    fill_defaults(p, args, {'sample_root': ('eht_sample', {}), 'output_root': ('eht_ensemble', {})})
    from eht_axion.experiment import read_experiment
    from eht_sample import science_settings
    settings = summary_settings(read_experiment(args.config))
    sample, destination = args.sample_root.resolve(), args.output_root.resolve()
    if destination.is_relative_to(sample) or sample.is_relative_to(destination):
        p.error('sample and output roots must be separate and non-nested')
    if destination.exists():
        p.error(f'new output root required: {destination}')
    with provenance.verification_session():
        selected, inputs, identity = verify_batch(sample, settings, science_settings(args.config))
        config = dict(batch_identity=identity, draw_count=settings['draws'] * len(en.DAYS))
        with provenance.stage_directory(destination, 'eht:ensemble', config, inputs) as out:
            if out is not None:
                records = provenance.describe_inputs(inputs)
                (out / 'consumed-inputs.json').write_text(json.dumps(records, indent=2) + '\n')
                run(selected, out, settings)
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
