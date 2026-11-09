"""Physical constants, M87* unit scales and the axion 211 frequency."""
import pytest

from eht_axion.physics import PhysicalScale, omega


def test_physical_scale_and_validation():
    scale = PhysicalScale(mass_solar=6.5e9, distance_mpc=16.8, label_step_tg=0.1)
    assert scale.tg_hours == pytest.approx(8.89324754)
    assert scale.angular_rg_uas == pytest.approx(3.81899325)
    with pytest.raises(ValueError):
        PhysicalScale(mass_solar=float("nan"), distance_mpc=16.8, label_step_tg=0.1)


@pytest.mark.parametrize("alpha", [0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.44])
def test_omega_principal_n(alpha):
    # 211 state: radial overtone 0 + l 1 + 1 = principal N 2.
    assert omega(alpha, 2) == alpha * (1 - alpha**2 / 8)
    assert omega(alpha, 4) == alpha * (1 - alpha**2 / 32)
