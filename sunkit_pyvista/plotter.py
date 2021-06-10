import functools

import numpy as np
import pfsspy
import pyvista as pv
from pfsspy import tracing

from astropy.constants import R_sun
from astropy.coordinates import SkyCoord
from sunpy.coordinates import HeliocentricInertial
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
        self._plotter.show_axes()

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

    def plot_field_lines(self, gong_map, lon, lat, radius, nrho=25, rss=2.5, **kwargs):
        """
        Plots the field lines from `pfsspy`.

        Parameters
        ----------
        gong_map : `pfsspy.map.GongSynopticMap`
            to produce the field lines
        lon : `astropy.units.quantity`
        lat : `astropy.units.quantity`
        radius : `float`
        nrho : `int`
        rss : `int`
        **kwargs :
            Keyword arguments are handed to `pyvista.Plotter.add_mesh`.
        """
        tracer = tracing.PythonTracer()
        input_ = pfsspy.Input(gong_map, nrho, rss)
        output_ = pfsspy.pfss(input_)

        seeds = SkyCoord(lon, lat, radius*R_sun,
                         frame=gong_map.coordinate_frame)

        field_lines = tracer.trace(seeds, output_)
        grid_data = []

        for field_line in field_lines:
            grid = self._coords_to_xyz(field_line.coords.ravel())
            grid_data.append(grid)

        grid_data = np.concatenate(grid_data)
        field_line_mesh = pv.StructuredGrid(grid_data[:, 0], grid_data[:, 1], grid_data[:, 2])

        self.plotter.add_mesh(field_line_mesh)
