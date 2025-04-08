import warnings

import astropy.constants as const
import numpy as np
import pyvista as pv
from astropy.coordinates import SkyCoord
from sunpy.coordinates import HeliocentricInertial
from sunpy.map import all_corner_coords_from_map


def pfss_to_grid(pfss_output, coordinate_frame):
    """
    Parameters
    ----------
    pfss_output : sunkit_magex.pfss.Output
        sunkit_magex PFSS model output.
    coordinate_frame : astropy.coordinates.BaseCoordinateFrame
        Coordinate frame for the output grid.

    Returns
    -------
    pyvista.StructuredGrid
    """
    if coordinate_frame != HeliocentricInertial():
        # TODO: generalise this to any coordinate frame that has origin at the
        # solar center
        msg = "coordinate_frame must be Heliocentric"
        raise ValueError(msg)
    # Get latitude and longitude coordinate values
    corner_coords = all_corner_coords_from_map(pfss_output.input_map).T
    corner_coords = corner_coords.transform_to(coordinate_frame)
    corner_coords.representation_type = "spherical"
    lon = corner_coords.lon[:, :, np.newaxis]
    lat = corner_coords.lat[:, :, np.newaxis]
    # Get radial coordinate values
    distance = np.exp(pfss_output.grid.rg) * const.R_sun
    distance = distance[np.newaxis, np.newaxis, :]
    # Put lat/lon/distance together
    corner_coords = SkyCoord(lon=lon, lat=lat, distance=distance, frame=coordinate_frame)
    # Create grid
    corner_coords.representation_type = "cartesian"
    grid = pv.StructuredGrid(
        corner_coords.x.to_value(const.R_sun),
        corner_coords.y.to_value(const.R_sun),
        corner_coords.z.to_value(const.R_sun),
    )
    # Add data
    b_grid = pfss_output.bg
    # Internal pyvisa call is raising:
    # Call to deprecated method GetDimensions. Use GetDimensions() instead.
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        # TODO: Check these the correct order
        grid.point_data["Br"] = b_grid[..., 2].flatten(order="F")
        grid.point_data["Btheta"] = b_grid[..., 1].flatten(order="F")
        grid.point_data["Bphi"] = b_grid[..., 0].flatten(order="F")
    return grid
