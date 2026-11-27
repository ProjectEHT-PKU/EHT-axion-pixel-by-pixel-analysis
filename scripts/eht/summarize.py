#!/usr/bin/env python3
"""EHT imaging: summarize the D-term Monte Carlo images of scripts/eht/sample.py into EHT ensembles.

The ensembles are the EHT input of the analysis (scripts/analysis/constraints.py).

pixi run -e eht-summary eht-summarize --sample-root MC --output-root OUT
MC is the output root of scripts/eht/sample.py; its batch must be finished and must have been
made with the science values of the current experiment.yaml. Draw count and sky grid
come from experiment.yaml. Only the blurred images (blurred with
shared.image.blur_fwhm_uas, like the simulated images) are summarized.
Outputs: stokes.npy (day, draw, y, x, IQUV), global.npz (values: day, draw, column; days;
columns), axes.npz, report.json (conventions, registration, scope).
"""

from pathlib import Path
import sys

sys.path.insert(0, str(next(p for p in Path(__file__).resolve().parents if (p / "eht_axion").is_dir())))

import argparse
import json

import numpy as np

from eht_axion.environment import add_data_root_option, fill_defaults
from eht_axion.experiment import DEFAULT_PATH
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
    """Check a finished scripts/eht/sample.py output root before any image is read.

    batch.json must record the science values ``science`` of the current experiment;
    finish.json must list one certified receipt per draw.
    """
    results = Path(root).resolve(strict=True)
    count = settings['draws']
    batch = read_json(results / 'batch.json')
    finish = read_json(results / 'finish.json')
    if (batch.get('schema') != 'eht-sample-batch' or finish.get('schema') != 'eht-sample-finish'
            or finish.get('numerically_certified') != count * len(en.DAYS)):
        raise ValueError('incomplete or inconsistent Monte Carlo batch (batch.json / finish.json)')
    if batch.get('science') != json.loads(json.dumps(science)):
        raise ValueError('Monte Carlo batch was made with different experiment values')
    from scripts.eht.sample import enumerate_tasks

    draws = read_json(results / 'all-draws.json')
    tasks = enumerate_tasks(draws)
    expected = {f'{d}-{k:04d}' for d in en.DAYS for k in range(count)}
    if (len(tasks) != count * len(en.DAYS) or {t['id'] for t in tasks} != expected
            or set(finish['receipts']) != expected):
        raise ValueError(f'not exactly {count} unique draws/day')
    paths = [results / 'batch.json', results / 'finish.json', results / 'all-draws.json']
    selected = {d: {} for d in en.DAYS}
    for task in tasks:
        d, k = task['day'], task['draw']
        if (d not in en.DAYS or not isinstance(k, int) or not 0 <= k < count
                or task['id'] != f'{d}-{k:04d}'):
            raise ValueError('draw/task mismatch')
        directory = contained(results, finish['receipts'][task['id']])
        receipt = read_json(directory / 'taskreceipt.json')
        if (receipt['task'] != task or receipt['status'] != 'certified'
                or receipt['schema'] != 'eht-sample-task'):
            raise ValueError('selected receipt task/status mismatch')
        products = {kind: contained(directory, name) for kind, name in
                    [('blurred', 'best-blurred-IQU.fits'), ('certificate', 'certificate.json')]}
        cert = read_json(products['certificate'])
        solver = cert['solver']
        if not (solver['numerically_certified'] and -1e-10 <= solver['best_gap'] <= cert['gap_atol']):
            raise ValueError('invalid numerical certificate')
        paths += [directory / 'taskreceipt.json', products['blurred'], products['certificate']]
        selected[d][k] = products
    return selected, paths


def run(selected, out, settings):
    draws, n = settings['draws'], settings['shape']
    pixel = settings['fov_uas'] / n
    reference = None
    samples = np.lib.format.open_memmap(out / 'stokes.npy', mode='w+', dtype='float64',
                                       shape=(len(en.DAYS), draws, n, n, 4))
    for i, d in enumerate(en.DAYS):
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
            samples[i, k] = loaded['stokes']
        print(f'SUMMARIZED day-{d}', flush=True)
    samples.flush()
    # The Paper VII Sec. 5.2 alignment is a zero shift here: images used as reconstructed.
    values = np.array([en.global_values(samples[i], *reference[:2]) for i in range(len(en.DAYS))])
    np.savez(out / 'global.npz', values=values, days=np.array(en.DAYS), columns=np.array(en.GLOBAL_COLUMNS))
    del samples
    report = dict(schema='eht-ensemble', days=list(en.DAYS), draws_per_day=draws,
                  map_order='day,draw,y,x,IQUV; native FITS E decreases along x, N increases along y',
                  angle_convention='EVPA radians modulo pi, east of north',
                  registration='none: the Paper VII Sec. 5.2 alignment is a zero shift; '
                               'images used as reconstructed',
                  resolution='blurred products of scripts/eht/sample.py (shared.image.blur_fwhm_uas); '
                             'no blur applied here',
                  scope='conditional fixed-I D-term MC only; no thermal bootstrap, response variance, '
                        'covariance-estimation uncertainty or scientific acceptance certification',
                  beta2='Paper VII eq. 14 with rho_min 0 and rho_max the full field: '
                        'sum((Q+iU)*exp(-2i atan2(E,N)))/sum(I) about the grid centre',
                  frequency_hz=reference[2], bunit=reference[3])
    (out / 'report.json').write_text(json.dumps(report, indent=2, allow_nan=False) + '\n')


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--sample-root', type=Path,
                   help='finished scripts/eht/sample.py output root (default: eht/sample)')
    p.add_argument('--output-root', type=Path, help='new directory (default: eht/ensemble)')
    p.add_argument('--config', type=Path, default=DEFAULT_PATH,
                   help='default: configs/experiment.yaml')
    add_data_root_option(p)
    args = p.parse_args(argv)
    fill_defaults(p, args, {'sample_root': ('eht.sample', {}), 'output_root': ('eht.ensemble', {})})
    from eht_axion.experiment import read_experiment
    from scripts.eht.sample import science_settings
    settings = summary_settings(read_experiment(args.config))
    sample, destination = args.sample_root.resolve(), args.output_root.resolve()
    if destination.is_relative_to(sample) or sample.is_relative_to(destination):
        p.error('sample and output roots must be separate and non-nested')
    if destination.exists():
        p.error(f'new output root required: {destination}')
    selected, inputs = verify_batch(sample, settings, science_settings(args.config))
    config = dict(draw_count=settings['draws'] * len(en.DAYS))
    with provenance.stage_directory(destination, 'eht:ensemble', config, inputs) as out:
        if out is not None:
            run(selected, out, settings)
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
