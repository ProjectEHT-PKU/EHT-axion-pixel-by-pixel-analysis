#!/usr/bin/env python3
"""Process 2 -> 3: summarize the copied fixed-D FITS draws of eht_sample.py into ensembles.

pixi run -e eht-research python scripts/eht_summarize.py --dataset-root ROOT --output-root OUT
Draw count, sky grid, aperture (radius, width) and the resultant threshold come
from experiment.yaml, the same values used by the analysis.
Requires ROOT/COPY_VERIFIED.json with the original finish SHA and batch identity.
Outputs: {native,icentroid}/{raw,display19}/day-D/{stokes.npy,evpa.npy,
map-summary.npz,aperture.npz,global.npz}; pair summaries under each mode's pairs/.
See report.json for estimator definitions, registration, and conditional scope.
"""
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

import argparse  # noqa: E402
import hashlib  # noqa: E402
import json  # noqa: E402

import numpy as np  # noqa: E402

from eht_axion.pipeline import common  # noqa: E402
from eht_axion.research import ensemble as en  # noqa: E402
from eht_axion.research.polarimetry import load_stokes_fits  # noqa: E402


def contained(root, relative):
    p = Path(relative)
    target = (root / p).resolve(strict=True)
    if p.is_absolute() or '..' in p.parts or not target.is_relative_to(root):
        raise ValueError('unsafe snapshot path: ' + str(relative))
    return target


def read_json(path):
    return json.loads(Path(path).read_text())


def summary_settings(config):
    """Values from experiment.yaml used by this summary (recorded in its manifest)."""
    image, aperture = config['shared']['image'], config['analysis']['aperture']
    return dict(draws=config['eht_imaging']['uncertainty_draws'], shape=image['shape'],
                fov_uas=image['fov_uas'],
                rmax_uas=aperture['radial_max_fraction'] * (image['fov_uas'] / 2),
                window_degrees=aperture['window_degrees'],
                rmin=config['analysis']['features']['resultant_min'])


def verify_copy(dataset, settings):
    """Reverify consumed bytes once/session, not old absolute producer closures."""
    dataset = Path(dataset).resolve(strict=True)
    draws = settings['draws']
    paths = []

    def checked(path, digest=None):
        actual = common.sha256(path)
        if digest is not None and actual != digest:
            raise ValueError('SHA mismatch: ' + str(path))
        paths.append(path)
        return read_json(path)

    copy = checked(dataset / 'COPY_VERIFIED.json')
    results = contained(dataset, 'results/full4000-v1')
    if (copy.get('schema') != 'eht-copyback-verified-v1' or copy.get('complete') is not True
            or copy.get('copied_root') != str(dataset) or copy.get('results_root') != str(results)
            or copy.get('certified_tasks') != draws * len(en.DAYS)
            or copy.get('days') != {d: draws for d in en.DAYS}
            or not all(copy.get(k) is True for k in ('full_attempt_outputs_SHA_verified',
                         'exported_inputs_SHA_verified', 'runtime_code_files_SHA_verified'))):
        raise ValueError('unsupported/incomplete copy verification receipt')
    finish = checked(results / 'finish.json', copy['finish_sha256'])
    batch = checked(results / 'batch.json', copy['batch_sha256'])
    identity = copy['identity']
    unsigned = dict(batch)
    unsigned.pop('identity', None)
    encoded = json.dumps(unsigned, sort_keys=True, separators=(',', ':'), allow_nan=False).encode()
    if (hashlib.sha256(encoded).hexdigest() != identity or batch['identity'] != identity
            or finish['identity'] != identity or finish['schema'] != 'eht-uncertainty-finish-v1'
            or finish['numerically_certified'] != draws * len(en.DAYS)):
        raise ValueError('batch/finish identity mismatch')
    draws = checked(results / 'all-draws.json')
    taskpath = results / 'tasks.jsonl'
    common.sha256(taskpath)
    paths.append(taskpath)
    tasks = [json.loads(line) for line in taskpath.read_text().splitlines()]
    expected = {f'{d}-{k:04d}' for d in en.DAYS for k in range(draws)}
    if (len(tasks) != draws * len(en.DAYS) or {t['id'] for t in tasks} != expected
            or set(finish['receipts']) != expected):
        raise ValueError(f'not exactly {draws} unique draws/day')
    selected = {d: {} for d in en.DAYS}
    for task in tasks:
        d, k = task['day'], task['draw']
        if (d not in en.DAYS or not isinstance(k, int) or not 0 <= k < draws
                or task['id'] != f'{d}-{k:04d}' or task['dterms'] != draws[d]['draws'][k]):
            raise ValueError('draw/task mismatch')
        entry = finish['receipts'][task['id']]
        directory = contained(results, entry['directory'])
        receipt = checked(directory / 'taskreceipt.json', entry['sha256'])
        if (receipt['task'] != task or receipt['identity'] != identity
                or receipt['status'] != 'certified' or receipt['schema'] != 'eht-uncertainty-task-v1'):
            raise ValueError('selected receipt identity/status mismatch')
        products = {}
        for kind, name in [('raw', 'best-raw-IQU.fits'),
                           ('display19', 'best-display-IQU-additional19uas.fits'),
                           ('certificate', 'certificate.json')]:
            p = contained(directory, name)
            if common.sha256(p) != receipt['outputs'][name]:
                raise ValueError('consumed output SHA mismatch: ' + str(p))
            paths.append(p)
            products[kind] = p
        cert = read_json(products['certificate'])
        solver = cert['solver']
        if not (solver['numerically_certified'] and -1e-10 <= solver['best_gap'] <= cert['gap_atol']):
            raise ValueError('invalid numerical certificate')
        selected[d][k] = products
    return selected, list(dict.fromkeys(paths)), copy


def run(selected, out, settings):
    draws, n, rmin = settings['draws'], settings['shape'], settings['rmin']
    pixel = settings['fov_uas'] / n
    reference = None
    registrations = {}
    summaries = {}
    for d in en.DAYS:
        daily = {}
        for resolution in ('raw', 'display19'):
            directory = out / 'native' / resolution / ('day-' + d)
            directory.mkdir(parents=True)
            samples = np.lib.format.open_memmap(directory / 'stokes.npy', mode='w+',
                                               dtype='float64', shape=(draws, n, n, 4))
            for k in range(draws):
                loaded = load_stokes_fits(selected[d][k][resolution])
                east, north = loaded['east_uas'], loaded['north_uas']
                if (loaded['stokes'].shape != (n, n, 4)
                        or not np.allclose(np.diff(east), -pixel, rtol=0, atol=1e-9)
                        or not np.allclose(np.diff(north), pixel, rtol=0, atol=1e-9)):
                    raise ValueError('FITS grid differs from shared.image (E decreasing, N increasing)')
                geometry = (east, north, loaded['frequency_hz'], loaded['bunit'])
                if reference is None:
                    reference = geometry
                    np.savez(out / 'axes.npz', east_uas=east, north_uas=north, angle_deg=np.arange(360.))
                elif (not np.array_equal(east, reference[0]) or not np.array_equal(north, reference[1])
                      or geometry[2:] != reference[2:]):
                    raise ValueError('inconsistent FITS axes/frequency/units')
                samples[k] = loaded['stokes']
            samples.flush()
            daily[resolution] = samples
        offset = en.centroid_offset(np.mean(daily['raw'][..., 0], axis=0), reference[0], reference[1])
        registrations[d] = dict(offset_yx_pixels=offset.tolist(), edge_I_fraction_lost={})
        for resolution, native in daily.items():
            for mode in ('native', 'icentroid'):
                directory = out / mode / resolution / ('day-' + d)
                if mode == 'native':
                    samples = native
                else:
                    directory.mkdir(parents=True)
                    samples = np.lib.format.open_memmap(directory / 'stokes.npy', mode='w+',
                                                       dtype='float64', shape=native.shape)
                    en.translate_ensemble(native, offset, samples)
                    before = native[..., 0].sum(axis=(1, 2))
                    after = samples[..., 0].sum(axis=(1, 2))
                    loss = 1 - after / before
                    registrations[d]['edge_I_fraction_lost'][resolution] = dict(
                        minimum=float(loss.min()), maximum=float(loss.max()), mean=float(loss.mean()))
                summaries[mode, resolution, d] = en.summarize(
                    samples, *reference[:2], directory, rmin, settings['rmax_uas'],
                    settings['window_degrees'])
                print(f'SUMMARIZED {mode}/{resolution}/day-{d}', flush=True)
                if mode == 'icentroid':
                    del samples
        del daily
    for mode in ('native', 'icentroid'):
        for resolution in ('raw', 'display19'):
            directory = out / mode / resolution / 'pairs'
            directory.mkdir()
            for later, earlier in en.PAIRS:
                a, ca = summaries[mode, resolution, later]
                b, cb = summaries[mode, resolution, earlier]
                np.savez(directory / f'{later}-{earlier}.npz', **en.independent_difference(a, b, ca, cb))
    report = dict(schema='eht-fixed-D-ensemble-summary-v1', days=list(en.DAYS), draws_per_day=draws,
                  registrations=registrations, rmin=rmin, rmax_uas=settings['rmax_uas'],
                  aperture_full_width_deg=float(settings['window_degrees']),
                  map_order='draw,y,x,IQUV; native FITS E decreases along x, N increases along y',
                  angle_convention='EVPA radians modulo pi, E of N; aperture angles degrees E of N; no roll',
                  registration='raw ensemble mean I centroid to sky zero; one bilinear zero-padded shift/day '
                               'shared by all channels/draws/resolutions; no angle recentering or clipping',
                  resolution='raw physical products and existing additional19uas FWHM display products; '
                             'no blur applied',
                  estimator='circular mean of per-draw EVPA, NOT EVPA of mean Stokes',
                  covariance='equal MC weights, population ddof=0, centered local-branch residuals; '
                             'independent-day covariance sum, never matched-index posterior pairs',
                  spatial_covariance='reconstruct selected columns using '
                                     'ensemble.selected_covariance(evpa,mean_rad,indices)',
                  rank_bound=draws - 1, mask='all draws finite with positive I and nonzero Q+iU; R>=rmin; '
                                      'map mean I and mean polarized amplitude positive; pair intersection',
                  scope='conditional fixed-I D-term MC only; no thermal bootstrap, response variance, '
                        'covariance-estimation uncertainty or scientific acceptance certification',
                  beta2='whole-field sum((Q+iU)*exp(-2i atan2(E,N)))/sum(I); no annular cut',
                  frequency_hz=reference[2], bunit=reference[3])
    (out / 'report.json').write_text(json.dumps(report, indent=2, allow_nan=False) + '\n')


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--dataset-root', type=Path, required=True)
    p.add_argument('--output-root', type=Path, required=True)
    p.add_argument('--config', type=Path, default=ROOT / 'experiment.yaml', help='default: experiment.yaml')
    args = p.parse_args(argv)
    from eht_axion.experiment import read_experiment
    settings = summary_settings(read_experiment(args.config))
    dataset, destination = args.dataset_root.resolve(), args.output_root.resolve()
    if destination.is_relative_to(dataset) or dataset.is_relative_to(destination):
        p.error('dataset and output roots must be separate and non-nested')
    with common.verification_session():
        selected, inputs, copy = verify_copy(dataset, settings)
        inputs += [Path(__file__).resolve(), args.config.resolve(),
                   ROOT / 'eht_axion/research/ensemble.py',
                   ROOT / 'eht_axion/research/polarimetry.py']
        config = dict(copy_identity=copy['identity'], rmin=settings['rmin'], rmax_uas=settings['rmax_uas'],
                      source_kind='verified-copy-custom-batch', draw_count=settings['draws'] * len(en.DAYS))
        with common.stage_directory(destination, 'research:eht-fixed-D-ensemble', config, inputs) as out:
            if out is not None:
                (out / 'copy-receipt.json').write_text(json.dumps(copy, indent=2) + '\n')
                records = common.describe_inputs(inputs)
                (out / 'consumed-inputs.json').write_text(json.dumps(records, indent=2) + '\n')
                run(selected, out, settings)
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
