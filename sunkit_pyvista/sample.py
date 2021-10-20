import astropy.units as u
import sunpy.map
from sunpy.data.sample import AIA_171_IMAGE


def low_res_aia_171(resolution=[512, 512] * u.pix):
    return sunpy.map.Map(AIA_171_IMAGE).resample(resolution)
