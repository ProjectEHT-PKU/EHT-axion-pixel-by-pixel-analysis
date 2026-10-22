"""Fixed-D wrapper regression: absent stations and no leakage fitting."""
from types import SimpleNamespace
from unittest.mock import Mock

import reconstruct_eht_fixed_d as w


def test_absent_station_filtering_and_fixed_schedule(monkeypatch):
    rows = {'site': ['AA', 'AZ']}
    base = SimpleNamespace(tarr=rows, copy=lambda: 'fresh')
    c = {'fsys': .005, 'maxit': 500}
    core = SimpleNamespace(add_fsys=Mock(return_value='noisy'),
        apply_fixed_d=Mock(return_value='corrected'), make_imaging_obs=Mock(return_value='imaging'))
    monkeypatch.setattr(w.driver, 'initial_image', lambda *a: 'initial')
    monkeypatch.setattr(w.driver, 'restart', Mock(return_value='blurred'))
    imaging = Mock(return_value=('image', [{'nit': 3, 'success': True}]))
    monkeypatch.setattr(w.driver, 'image_p', imaging)
    result, records = w.fixed_d_image(base, 'I', {'AA': [0j, 0j], 'AZ': [0j, 0j], 'JC': [0j, 0j]},
                                      c, core, 'eh')
    assert result == 'image' and len(records) == 6
    core.add_fsys.assert_called_once_with('fresh', .005)
    core.apply_fixed_d.assert_called_once_with('noisy', {'AA': [0j, 0j], 'AZ': [0j, 0j]})
    assert w.driver.restart.call_count == 4
    assert [r['blur_restart'] for r in records] == [False, True, True, True, True, False]
    assert all(r['unused_station_draws'] == ['JC'] for r in records)
    assert [call.args[3]['maxit'] for call in imaging.call_args_list] == [100]*5 + [500]
