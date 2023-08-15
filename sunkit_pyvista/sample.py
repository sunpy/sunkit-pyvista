import astropy.units as u
import sunpy.map

__all__ = ["LOW_RES_AIA_193"]

LOW_RES_AIA_193 = sunpy.map.Map(
    "https://github.com/sunpy/data/raw/main/sunpy/v1/AIA20110607_063307_0193_lowres.fits",
)
LOW_RES_AIA_193 = LOW_RES_AIA_193.resample([512, 512] * u.pix, method="spline")
