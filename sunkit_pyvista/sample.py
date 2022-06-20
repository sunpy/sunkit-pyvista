import astropy.units as u
import sunpy.data.sample
import sunpy.map


def low_res_aia_193(resolution=[512, 512] * u.pix):
    return sunpy.map.Map(sunpy.data.sample.AIA_193_IMAGE).resample(resolution)
