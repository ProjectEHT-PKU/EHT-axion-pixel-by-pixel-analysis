"""Explicit physical units versus preserved project configuration choices."""

from dataclasses import asdict, dataclass
import numpy as np

# IAU 2015 Resolution B3 nominal solar GM; c and hbar use SI/CODATA values.
GM_SUN = 1.3271244e20  # m^3 s^-2
C = 299792458.0  # m s^-1
PARSEC = 3.085677581491367e16  # m
HBAR_EV_S = 6.582119569e-16


@dataclass(frozen=True)
class PhysicalScale:
    """Physical unit conversions.

    Defaults preserve the legacy public API only. Analysis geometry always
    supplies every physical field explicitly from the validated configuration.
    """

    mass_solar: float = 6.5e9
    distance_mpc: float = 16.8
    label_step_tg: float = 0.1

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


def analysis_geometry(config) -> dict:
    """Resolve explicit analysis choices into preprocessing/aperture kwargs.

    Accepts either complete legacy Rg geometry or complete angular geometry;
    aperture extent independently uses Rg or a fraction of the final half-FOV.
    No defaults are filled and the input is not mutated. Only the analysis
    schema is checked here; read_config also validates the optional producer.
    Gaussian FWHM describes a benchmark matching kernel, not a measured beam.
    """
    from .config import _validate

    _validate(config)
    scale = PhysicalScale(**config["physical"])
    image, aperture = config["image"], config["aperture"]
    angular_mode = "final_fov_uas" in image
    final_fov = (scale.angular_to_rg(image["final_fov_uas"]) if angular_mode
                 else image["final_fov_rg"])
    sigma = (scale.gaussian_sigma_rg(image["gaussian_fwhm_uas"]) if angular_mode
             else image["gaussian_sigma_rg"])
    radial_max = (aperture["radial_max_fraction"] * (final_fov / 2)
                  if "radial_max_fraction" in aperture else aperture["radial_max_rg"])
    if "eht_input" in config:
        declared = config["eht_input"]
        if (declared["shape"] != image["final_shape"] or not np.isclose(
                final_fov * scale.angular_rg_uas, declared["fov_uas"], rtol=1e-12, atol=0)):
            raise ValueError("Analysis target grid differs from declared EHT input sky grid")
    return {
        "image_preprocess": {
            "sigma": sigma,
            "init_fov": image["initial_fov_rg"],
            "final_fov": final_fov,
            "final_shape": image["final_shape"],
            "rotation_angle": np.deg2rad(image["rotation_degrees"]),
        },
        "aperturefit": {
            "fov": final_fov,
            "weighting": aperture["weighting"],
            "radial_max": radial_max,
            "window_degrees": aperture["window_degrees"],
            "angular_bins": aperture["angular_bins"],
        },
        "metadata": dict(
            scale.metadata(),
            image_mode="angular" if angular_mode else "rg",
            aperture_mode="fraction" if "radial_max_fraction" in aperture else "rg",
            initial_fov_uas=image["initial_fov_rg"] * scale.angular_rg_uas,
            final_fov_uas=final_fov * scale.angular_rg_uas,
            gaussian_fwhm_uas=sigma * np.sqrt(8 * np.log(2)) * scale.angular_rg_uas,
            pixel_size_uas=final_fov * scale.angular_rg_uas / image["final_shape"],
            radial_max_rg=radial_max,
            radial_max_uas=radial_max * scale.angular_rg_uas,
        ),
    }
