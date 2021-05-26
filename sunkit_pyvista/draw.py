import matplotlib.pyplot as plt
import numpy as np

import astropy.units as u
from astropy.coordinates import Longitude, SkyCoord
from astropy.visualization.wcsaxes import Quadrangle

from sunpy.coordinates.utils import get_rectangle_coordinates
from sunpy.coordinates import HeliocentricInertial
from sunpy.map.maputils import all_corner_coords_from_map

from matplotlib.patches import Polygon
from astropy.constants import R_sun


class Draw:
    def __init__(self, coordinate_frame=None):
        if coordinate_frame is None:
            coordinate_frame = HeliocentricInertial()
        self._coordinate_frame = coordinate_frame

    @property
    def coordinate_frame(self):
        """
        Coordinate frame of the plot.
        """
        return self._coordinate_frame

    def _convert_to_Skycoord(self, m, coords, frame='heliographic_stonyhurst'):
        """
        Converts a given set of coordinates to SkyCood objects with
        """
        corners = SkyCoord(coords.xy[:,0]*u.deg, coords.xy[:,1]*u.deg , frame=frame, obstime=m.date)
        return corners

    def draw_quadrangle(self, m, bottom_left, *, width: u.deg = None, height: u.deg = None,
                        axes=None, top_right=None, **kwargs):
        
        bottom_left, top_right = get_rectangle_coordinates(
            bottom_left, top_right=top_right, width=width, height=height)

        width = Longitude(top_right.spherical.lon - bottom_left.spherical.lon)
        height = top_right.spherical.lat - bottom_left.spherical.lat
        if not axes:
            axes = plt.gca()
        kwergs = {
            "transform": axes.get_transform(bottom_left.frame.replicate_without_data()),
            "edgecolor": "white",
            "fill": False,
        }

        kwergs.update(kwargs)

        quad = Quadrangle(m._get_lon_lat(bottom_left), width, height, **kwergs)
        corners = self._convert_to_Skycoord(m, quad)
        nodes = self._coords_to_xyz(corners)
        return nodes

    def _coords_to_xyz(self, coords):
        coords = coords.transform_to(self.coordinate_frame)
        coords.representation_type = 'cartesian'
        return np.column_stack((coords.x.to_value(R_sun),
                                coords.y.to_value(R_sun),
                                coords.z.to_value(R_sun)))
    def create_map(self, m):
        """
        Create a mesh from a map.
        Parameters
        ----------
        m : sunpy.map.Map
        Returns
        -------
        pyvista.StructuredGrid
        """
        corner_coords = all_corner_coords_from_map(m)
        nodes = self._coords_to_xyz(corner_coords.ravel())
        return nodes
