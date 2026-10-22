"""New-root-only EHT120 theory preparation; never consumes EHT images or runs RT.

Canonical RAPTOR CLI orchestration plus mask-independent background preparation.
Selected-feature background/Cbg belongs to the parent after its EHT adapter.
"""
from pathlib import Path
import argparse
import json

import _bootstrap  # noqa: F401
from _stage import execute
from eht_axion.pipeline.common import stage_directory, verification_session, manifest_path
from eht_axion.pipeline.config import read_config
from eht_axion.pipeline.constants import analysis_geometry
from eht_axion.pipeline.layout import RunLayout, input_raptor_files
from eht_axion.pipeline.raptor import preprocess_cube, estimate_background
from eht_axion.pipeline.workflow import stage_inputs
from eht_axion.producer.delivery import generated_inputs, input_cadence_hours

PROJECT = Path(__file__).resolve().parents[1]
CONFIG = PROJECT / 'configs/eht-120uas-observed.yaml'
RAW = Path('/home/canoziia/GRMHD_DATA/production/m65-211-merged/data')
OUTPUT = Path('/home/canoziia/GRMHD_DATA/analysis/eht-120uas-observed-templates-v1')


def validate():
    config = read_config(CONFIG)
    original = read_config(PROJECT / 'configs/eht-observed-research.yaml')
    for section in ('physical', 'producer', 'data', 'cloud', 'background', 'inference', 'aperture', 'eht'):
        if config[section] != original[section]:
            raise ValueError(f'Unapproved change: {section}')
    expected = dict(original['image'], final_fov_uas=120.0)
    if config['image'] != expected:
        raise ValueError('Only theory field width may change')
    contract = config['eht_input']
    assert (contract['shape'], contract['fov_uas'], contract['pixel_size_uas'],
            contract['angle_roll']) == (64, 120.0, 1.875, -90)
    assert all(v == [0, 0] for name in ('rank_offsets', 'final_offsets')
               for v in contract[name].values())
    geometry = analysis_geometry(config)
    assert geometry['metadata']['radial_max_uas'] == 60.0
    return config, geometry


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--validate-only', action='store_true')
    args = parser.parse_args()
    config, geometry = validate()
    print(json.dumps(dict(config=str(CONFIG), raw=str(RAW), output=str(OUTPUT),
                          stages=['raptor:preprocess', 'raptor:fit', 'raptor:summarize',
                                  'raptor:background-preprocess', 'raptor:background-unmasked'],
                          geometry=geometry), indent=2), flush=True)
    if args.validate_only:
        return
    if OUTPUT.exists():
        raise FileExistsError(f'Fresh output root required: {OUTPUT}')
    with verification_session():
        execute('raptor', ['--config', str(CONFIG), '--raptor-input-root', str(RAW),
                          '--output-root', str(OUTPUT), '--no-trust-pickle', '--no-rebuild',
                          '--no-reuse-imported', '--no-env-file'])
        layout = RunLayout(OUTPUT)
        cfg = config['background']
        cubes = input_raptor_files(RAW, cfg['alpha'], cfg['source_resolution'])
        inputs, _ = generated_inputs(cubes, config, alpha=cfg['alpha'],
                                     resolution=cfg['source_resolution'], verified=set())
        cadence = input_cadence_hours(cubes, config)
        params = dict(geometry['image_preprocess'], state_indices=[0])
        with stage_directory(layout.background_processed, 'raptor:background-preprocess',
                             params, inputs) as temporary:
            if temporary is not None:
                preprocess_cube(cubes, temporary, **params)
        # Explicitly UNMASKED diagnostics: not a substitute for parent-selected Cbg.
        for dim in ('1d', '2d'):
            params = dict(lag=cfg['lag_frames'], cadence_hours=cadence, method=cfg['method'],
                          domain=dim, mask=None)
            if dim == '1d':
                params.update(geometry['aperturefit'])
            dest = OUTPUT / 'background-unmasked' / dim
            with stage_directory(dest, 'raptor:background-unmasked', params,
                                 stage_inputs(layout.background_processed)) as temporary:
                if temporary is not None:
                    estimate_background(layout.background_processed / 'stokes.npy', temporary, **params)
        record = dict(status='template-and-unmasked-background-preparation-complete',
                      response_summary_roots={dim: str(layout.templates(dim)) for dim in ('1d', '2d')},
                      background_processed_root=str(layout.background_processed),
                      background_unmasked_roots={dim: str(OUTPUT / 'background-unmasked' / dim)
                                                 for dim in ('1d', '2d')},
                      pending=['parent EHT120 native-coordinate/centroid ensemble adapter',
                               'selected-feature background statistics/Cbg and observed limits'],
                      eht_inputs_consumed=False, science_certified=False)
        with stage_directory(OUTPUT / 'template-run', 'research:eht120-template-preparation',
                             record, [CONFIG, Path(__file__).resolve()] +
                             [manifest_path(layout.templates(d)) for d in ('1d', '2d')] +
                             [manifest_path(OUTPUT / 'background-unmasked' / d) for d in ('1d', '2d')]) as out:
            if out is not None:
                (out / 'report.json').write_text(json.dumps(record, indent=2) + '\n')
    print('EXIT0: template manifests and unmasked background preparation complete; NOT science certification',
          flush=True)


if __name__ == '__main__':
    main()
