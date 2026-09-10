import astropy.units as u

import sunpy.map
from sunpy.data.sample import AIA_193_IMAGE

__all__ = ["LOW_RES_AIA_193"]


LOW_RES_AIA_193 = sunpy.map.Map(AIA_193_IMAGE).resample([512, 512] * u.pix, method="spline")
