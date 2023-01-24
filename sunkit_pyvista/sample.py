import astropy.units as u
import sunpy.map

__all__ = ["low_res_aia_193"]

low_res_aia_193 = sunpy.map.Map(
    "https://github.com/sunpy/data/raw/main/sunpy/v1/AIA20110607_063307_0193_lowres.fits"
).resample([512, 512] * u.pix, method="spline")
