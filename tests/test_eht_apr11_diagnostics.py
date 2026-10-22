"""Synthetic aggregation/axial-angle contracts; no science data or imaging."""
import json

import numpy as np
import pytest

import diagnose_eht_apr11 as d


def test_baselines_are_unordered():
    assert d.baseline_labels(['PV', 'AA', 'LM'], ['AA', 'PV', 'AZ']).tolist() == [
        'AA-PV', 'AA-PV', 'AZ-LM']


def test_marginal_chi2_and_counts():
    rows = d.aggregate(np.array([1+2j, 3+4j, 2j]), np.array([1., 2., 1.]),
                       np.array(['a', 'a', 'b']))
    assert rows[0]['count'] == 2
    assert rows[0]['chi2'] == pytest.approx((5+25/4)/4)
    assert rows[0]['chi2'] == pytest.approx(
        (rows[0]['real_chi2']+rows[0]['imag_chi2'])/2)
    assert rows[1]['chi2'] == 2
    with pytest.raises(ValueError, match='positive'):
        d.aggregate([1j], [0], ['a'])


def test_axial_evpa_wrap_and_weight():
    i = np.ones(2)
    p = np.array([.1, .5])*np.exp(2j*np.deg2rad([89., 0.]))
    q = np.array([.1, .5])*np.exp(2j*np.deg2rad([-89., 10.]))
    result = d.map_metrics(i, p, i, q, np.ones(2, bool))
    assert result['evpa_rms_deg'] == pytest.approx(np.sqrt((.1*4+.5*100)/.6))
    assert result['cos2_delta'] == pytest.approx(
        (.1*np.cos(np.deg2rad(4))+.5*np.cos(np.deg2rad(20)))/.6)
    same = d.map_metrics(i, p, i, p, np.ones(2, bool))
    assert same['qu_relative_l2'] == 0
    assert same['evpa_rms_deg'] == pytest.approx(0, abs=1e-12)
    assert same['net_m_a'] == pytest.approx(abs(p.sum())/2)
    with pytest.raises(ValueError, match='empty'):
        d.map_metrics(i, p, i, q, np.zeros(2, bool))


def test_dterms_plain_and_round_wrapper(tmp_path):
    table = {'AA': [[0., .1], [.2, 0.]]}
    path = tmp_path/'round-100.json'
    for value in (table, {'round': 100, 'dterms': table}):
        path.write_text(json.dumps(value))
        np.testing.assert_array_equal(d.read_dterms(path)['AA'], [.1j, .2])
    path.write_text(json.dumps({'AA': [[0, 1, 2]]}))
    with pytest.raises(ValueError, match='invalid D'):
        d.read_dterms(path)
