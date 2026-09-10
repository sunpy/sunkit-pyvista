import warnings

import astropy.units as u
from astropy.io import fits

import sunpy.map
from sunpy.data.sample import AIA_193_IMAGE

__all__ = ["LOW_RES_AIA_193"]

with warnings.catch_warnings():
    warnings.simplefilter("ignore", fits.verify.VerifyWarning)
    LOW_RES_AIA_193 = sunpy.map.Map(AIA_193_IMAGE).resample([512, 512] * u.pix, method="spline")
