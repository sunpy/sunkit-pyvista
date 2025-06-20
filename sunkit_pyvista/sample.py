import warning

import astropy.units as u

import sunpy.map

__all__ = ["LOW_RES_AIA_193"]


with warning.catch_warnings():
    # Ignore the warning about the FITS file not being in the correct format
    # This is because the FITS file is a low-resolution version of an AIA image
    warning.simplefilter("ignore")
    LOW_RES_AIA_193 = sunpy.map.Map(
        "https://github.com/sunpy/data/raw/main/sunpy/v1/AIA20110607_063307_0193_lowres.fits",
    )
    LOW_RES_AIA_193 = LOW_RES_AIA_193.resample([512, 512] * u.pix, method="spline")
