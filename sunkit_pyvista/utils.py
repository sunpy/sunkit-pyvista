import numpy as np

import astropy.units as u
from astropy.coordinates import SkyCoord
from sunpy.coordinates import HeliographicStonyhurst
from sunpy.sun import constants

__all__ = ['get_limb_coordinates']


@u.quantity_input
def get_limb_coordinates(observer, rsun=constants.radius, resolution=1000):
    """
    Get coordinates for the solar limb as viewed by a specified observer.

    Parameters
    ----------
    observer : `~astropy.coordinates.SkyCoord`
        Observer coordinate.
    rsun : `~astropy.units.Quantity`
        Physical radius of the limb from Sun center. Defaults to the standard
        photospheric radius.
    resolution : int
        Number of coordinates to return. The coordinates are equally spaced
        around the limb as seen from the observer.
    """
    observer = observer.transform_to(
        HeliographicStonyhurst(obstime=observer.obstime))
    dsun = observer.radius
    if dsun <= rsun:
        raise ValueError('Observer distance must be greater than rsun')
    # Create the limb coordinate array using Heliocentric Radial
    limb_radial_distance = np.sqrt(dsun**2 - rsun**2)
    limb_hcr_rho = limb_radial_distance * rsun / dsun
    limb_hcr_z = dsun - np.sqrt(limb_radial_distance**2 - limb_hcr_rho**2)
    limb_hcr_psi = np.linspace(0, 2*np.pi, resolution+1)[:-1] << u.rad
    limb = SkyCoord(limb_hcr_rho, limb_hcr_psi, limb_hcr_z,
                    representation_type='cylindrical',
                    frame='heliocentric',
                    observer=observer, obstime=observer.obstime)
    return limb
