import functools

import numpy as np
import pyvista as pv

import astropy.units as u
from astropy.constants import R_sun
from sunpy.coordinates import HeliocentricInertial, SkyCoord
from sunpy.map.maputils import all_corner_coords_from_map

__all__ = ['SunpyPlotter']


class SunpyPlotter:
    """
    A plotter for 3D data.

    This class wraps `pyvsita.Plotter` to provide coordinate-aware plotting.
    For now, all coordinates are converted to
    a specific frame (`~sunpy.coordinates.HeliocentricInertial` by default),
    and distance units are such that :math:`R_{sun} = 1`.

    Parameters
    ----------
    coordinate_frame : `astropy.coordinates.BaseFrame`
        Coordinate frame of the plot. The x, y, z axes of the pyvista plotter
        will be the x, y, z axes in this coordinate system.
    """
    def __init__(self, coordinate_frame=None):
        if coordinate_frame is None:
            coordinate_frame = HeliocentricInertial()
        self._coordinate_frame = coordinate_frame
        self._plotter = pv.Plotter()

    @property
    def coordinate_frame(self):
        """
        Coordinate frame of the plot.
        """
        return self._coordinate_frame

    @property
    def plotter(self):
        """
        `pyvista.Plotter`.
        """
        return self._plotter

    @functools.wraps(pv.Plotter.show)
    def show(self, *args, **kwargs):
        """
        Show the plot.
        """
        self.plotter.show(*args, **kwargs)

    def _coords_to_xyz(self, coords):
        coords = coords.transform_to(self.coordinate_frame)
        coords.representation_type = 'cartesian'
        return np.column_stack((coords.x.to_value(R_sun),
                                coords.y.to_value(R_sun),
                                coords.z.to_value(R_sun)))

    def set_camera_coordinates(self, coord):
        """
        Sets the inital camera position of the rendered plot.

        Parameters
        ----------
        coords : `astropy.coordinates.SkyCoord`
            Coordinates of the camera.
        """
        camera_position = self._coords_to_xyz(coord)
        pos = tuple(camera_position[0])
        self.plotter.camera.position = pos

    def rotate_camera(self, angle: u.deg = None):
        """
        Rotates the camera by the specified value in degrees.

        Parameters
        ----------
        angle : `astropy.units.Quantity`
            The angle of rotation.
        """
        rotation_angle = angle.to_value(u.deg)
        self.plotter.camera.roll = rotation_angle

    def _pyvista_mesh(self, m):
        """
        Create a mesh from a map.

        Parameters
        ----------
        m : `sunpy.map.Map`
            The map to use.

        Returns
        -------
        `pyvista.StructuredGrid`
        """
        corner_coords = all_corner_coords_from_map(m)
        nodes = self._coords_to_xyz(corner_coords.ravel())
        grid = pv.StructuredGrid()
        grid.points = nodes
        grid.dimensions = [m.data.shape[0] + 1,
                           m.data.shape[1] + 1,
                           1]
        data = m.data.T.reshape(-1)
        grid['data'] = m.plot_settings['norm'](data)
        return grid

    def plot_map(self, m, **kwargs):
        """
        Plot a map.

        Parameters
        ----------
        m : `sunpy.map.Map`
            Map to be plotted.
        **kwargs :
            Keyword arguments are handed to `pyvista.Plotter.add_mesh`.
        """
        cmap = kwargs.pop('cmap', m.cmap)
        mesh = self._pyvista_mesh(m)
        self.plotter.add_mesh(mesh, cmap=cmap, **kwargs)

    def plot_line(self, coords, **kwargs):
        """
        Plot a line from a set of coordinates.

        Parameters
        ----------
        coords : `astropy.coordinates.SkyCoord`
            Coordinates to plot as a line.
        **kwargs :
            Keyword arguments are passed to `pyvista.Plotter.add_mesh`.
        Notes
        -----
        This plots a `pyvista.Spline` object.
        """
        points = self._coords_to_xyz(coords)
        spline = pv.Spline(points)
        self.plotter.add_mesh(spline, **kwargs)

    def plot_solar_axis(self, length=2.5, arrow_kwargs={}, **kwargs):
        """
        Plot the solar rotation axis as an arrow.

        Parameters
        ----------
        length : float
            Length of the arrow in multiples of solar radii.
        arrow_kwargs : dict
            Keyword arguments to be handed to `pyvista.Arrow`.
            ``start``, ``direction``, and ``scale`` cannot be manually
            specified, as they are automatically set.
        **kwargs :
            Keyword arguments are handed to `pyvista.Plotter.add_mesh`.
        """
        defaults = {'shaft_radius': 0.01,
                    'tip_length': 0.05,
                    'tip_radius': 0.02}
        defaults.update(arrow_kwargs)
        arrow = pv.Arrow(start=(0, 0, -length / 2),
                         direction=(0, 0, length),
                         scale='auto',
                         **defaults)
        self.plotter.add_mesh(arrow, **kwargs)

    def plot_quadrangle(self, m, bottom_left, top_right, width, height, **kwargs):
        """
        Plots a quadrangle on the given map.
        This draws a quadrangle that has corners at ``(bottom_left, top_right)``,
        if ``width`` and ``height`` are specified, they are respectively added to the
        longitude and latitude of the ``bottom_left`` coordinate to calculate a
        ``top_right`` coordinate.

        Parameters
        ----------
        bottom_left :
        top_right :
        width :
        height :
        **kwargs :
        """
        quadrangle_patch = m.draw_quadrangle(bottom_left, width, height, top_right, resolution=500)
        quadrangle_coordinates = quadrangle_patch.get_xy()

        c = SkyCoord(quadrangle_coordinates[:, 0]*u.deg, quadrangle_coordinates[:, 1]*u.deg,
                     frame=bottom_left.frame, obstime=m.date)
        mesh = self._coords_to_xyz(c)
        self.plotter.add_mesh(mesh, **kwargs)
