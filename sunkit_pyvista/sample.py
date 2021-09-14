import astropy.units as u
from sunpy.data.sample import AIA_171_IMAGE
import sunpy.map


def low_res_aia_171(resolution=[512, 512] * u.pix):
    return sunpy.map.Map(AIA_171_IMAGE).resample(resolution)
