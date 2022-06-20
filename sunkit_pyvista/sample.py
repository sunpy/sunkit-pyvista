import astropy.units as u
import sunpy.data.sample
import sunpy.map

__all__ = ["low_res_aia_193"]

low_res_aia_193 = sunpy.map.Map(sunpy.data.sample.AIA_193_IMAGE).resample([512, 512] * u.pix)
