"""Physical constants, M87* unit conversions (gravitational time/angle scales) and the
axion 211-state frequency omega(alpha), shared by RAPTOR production and the analysis."""

from dataclasses import asdict, dataclass
import numpy as np

# IAU 2015 Resolution B3 nominal solar GM; c and hbar use SI/CODATA values.
GM_SUN = 1.3271244e20  # m^3 s^-2
C = 299792458.0  # m s^-1
PARSEC = 3.085677581491367e16  # m
HBAR_EV_S = 6.582119569e-16


@dataclass(frozen=True)
class PhysicalScale:
    """Unit conversions for a black hole of given mass and distance (experiment.physical_scale)."""

    mass_solar: float
    distance_mpc: float
    label_step_tg: float

    def __post_init__(self):
        if not all(np.isfinite(v) and v > 0 for v in asdict(self).values()):
            raise ValueError("Physical scales must be finite and positive")

    @property
    def tg_seconds(self):
        return GM_SUN * self.mass_solar / C**3

    @property
    def tg_hours(self):
        return self.tg_seconds / 3600

    @property
    def angular_rg_uas(self):
        return (
            GM_SUN
            * self.mass_solar
            / C**2
            / (self.distance_mpc * 1e6 * PARSEC)
            * 180
            / np.pi
            * 3.6e9
        )

    def angular_to_rg(self, angle_uas):
        """Convert an angular size in microarcseconds to gravitational radii."""
        return angle_uas / self.angular_rg_uas

    def gaussian_sigma_rg(self, fwhm_uas):
        """Gaussian matching-kernel FWHM (microarcseconds) to sigma (Rg)."""
        return self.angular_to_rg(fwhm_uas) / np.sqrt(8 * np.log(2))

    def metadata(self):
        return dict(
            asdict(self),
            tg_seconds=self.tg_seconds,
            tg_hours=self.tg_hours,
            label_step_hours=self.label_step_tg * self.tg_hours,
            angular_rg_uas=self.angular_rg_uas,
            omega_conversion_per_hour=1 / self.tg_hours,
            axion_mass_per_alpha_ev=HBAR_EV_S / self.tg_seconds,
            sources=["https://arxiv.org/abs/1906.11238", "https://arxiv.org/abs/1510.07674"],
        )


def omega(alpha: float, n_principal: int) -> float:
    """Leading small-alpha M*omega_R; principal N=2 for the scalar 211 mode.

    N = radial_overtone + l + 1, not N+l+1 or n+l+m. This is not
    an exact Kerr eigenfrequency or the interpolated radial-generator table.
    """
    if not np.isfinite(alpha) or alpha <= 0:
        raise ValueError("alpha must be positive and finite")
    if isinstance(n_principal, bool) or not isinstance(n_principal, (int, np.integer)) or n_principal < 1:
        raise ValueError("n_principal must be a positive integer")
    return alpha * (1 - alpha**2 / (2 * n_principal**2))
