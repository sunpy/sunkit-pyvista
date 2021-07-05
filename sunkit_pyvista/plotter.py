import functools
from functools import partial

import numpy as np
import pyvista as pv
import pyvistaqt as pvq

import astropy.units as u
from astropy.constants import R_sun
from astropy.coordinates import Longitude, SkyCoord
from astropy.visualization import AsymmetricPercentileInterval
from sunpy.coordinates import HeliocentricInertial
from sunpy.coordinates.utils import get_rectangle_coordinates
from sunpy.map import GenericMap
from sunpy.map.maputils import all_corner_coords_from_map
from sunpy.util import expand_list
from sunpy.visualization._quadrangle import Quadrangle

from sunkit_pyvista.mapsequence_animator import SequenceAnimator

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
        self.camera = self._plotter.camera
        self.bg_plotter = pvq.BackgroundPlotter(show=True)

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

    def _toggle_animation(self, state, animate):
        animate.animation_state = state

    def _coords_to_xyz(self, coords):
        coords = coords.transform_to(self.coordinate_frame)
        coords.representation_type = 'cartesian'
        return np.column_stack((coords.x.to_value(R_sun),
                                coords.y.to_value(R_sun),
                                coords.z.to_value(R_sun)))

    def _get_clim(self, data, clip_interval):
        """
        Get vmin, vmax of a data slice when clip_interval is specified.
        """
        percent_limits = clip_interval.to('%').value
        vmin, vmax = AsymmetricPercentileInterval(*percent_limits).get_limits(data)
        return [vmin, vmax]

    def set_camera_coordinate(self, coord):
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

    @u.quantity_input
    def set_view_angle(self, angle: u.deg):
        """
        Sets the view angle of the camera to the specified value

        Parameters
        ----------
        angle : `astropy.units.Quantity`
            The viewing angle.
        """
        view_angle = angle.to_value(u.deg)
        if not (view_angle > 0 and view_angle <= 180):
            raise ValueError("specified view angle must be "
                             "0 deg < angle <= 180 deg")
        # Zoom/view_angle = current view angle (default is set to 30 degrees) / 1
        zoom_value = self.camera.view_angle / view_angle
        self.plotter.camera.zoom(zoom_value)

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

    @u.quantity_input
    def plot_map(self, m, clip_interval: u.percent = None, **kwargs):
        """
        Plot a map.

        Parameters
        ----------
        m : `sunpy.map.Map`
            Map to be plotted.
        clip_interval : two-element `~astropy.units.Quantity`, optional
            If provided, the data will be clipped to the percentile
            interval bounded by the two numbers.
        **kwargs :
            Keyword arguments are handed to `pyvista.Plotter.add_mesh`.
        """
        cmap = kwargs.pop('cmap', m.cmap)
        mesh = self._pyvista_mesh(m)

        if clip_interval is not None:
            if len(clip_interval) == 2:
                clim = self._get_clim(data=mesh['data'],
                                      clip_interval=clip_interval)
            else:
                raise ValueError("Clip percentile interval must be "
                                 "specified as two numbers.")
        else:
            clim = [0, 1]
        self.plotter.add_mesh(mesh, cmap=cmap, clim=clim, **kwargs)

    def plot_map_sequence(self, *args, interval=2, **kwargs):
        """
        Plot a sequence of maps as an animation.

        Parameters
        ----------
        m : `sunpy.map.Map`
            Map(s) to be plotted.
        **kwargs :
            Keyword arguments are handed to `pyvistaq.BackGroundplotter.add_mesh`.
        """
        map_meshes = []
        color_maps = []

        maps = expand_list(args)
        for m in maps:
            if not isinstance(m, GenericMap):
                raise ValueError(
                    'MapSequence expects pre-constructed map objects.')
        for m in maps:
            mesh = self._pyvista_mesh(m)
            map_meshes.append(mesh)
            color_maps.append(m.cmap)
        animate = SequenceAnimator(time=interval, map_meshes=map_meshes,
                                   color_maps=color_maps, **kwargs)

        self.bg_plotter.add_mesh(map_meshes[0], cmap=color_maps[0], **kwargs)
        self.bg_plotter.add_checkbox_button_widget(
            partial(self._toggle_animation, animate=animate),
            value=False, color_on='green')
        self.bg_plotter.add_callback(partial(animate,
                                     bg_plotter=self.bg_plotter,
                                     **kwargs), interval=16)
        self.bg_plotter.add_text("Play", position=(70, 10))
        self.bg_plotter.enable_anti_aliasing()
        self.bg_plotter.hide_axes()

    def plot_coordinates(self, coords, radius=0.05, **kwargs):
        """
        Plot a sphere if a single coordinate is passed and
        plots a line if multiple coordinates are passed.

        Parameters
        ----------
        coords : `astropy.coordinates.SkyCoord`
            Coordinate(s) to plot as a center of sphere or line.
        radius : `int`, optional
            Radius of the sphere times the radius of the sun
            to be plotted when a single coordinate is passed.
            Defaults to ``0.05`` times the radius of the sun.
        **kwargs :
            Keyword arguments are passed to `pyvista.Plotter.add_mesh`.

        Notes
        -----
        This plots a `pyvista.Sphere` object if a single coordinate is passed
        and plots a `pyvista.Spline` object if multiple coordinates are passed.
        ``radius`` is only considered when a sphere is plotted.
        """
        points = self._coords_to_xyz(coords)
        if points.shape[0] > 1:
            point_mesh = pv.Spline(points)
        else:
            point_mesh = pv.Sphere(radius=radius, center=points[0])
        self.plotter.add_mesh(point_mesh, smooth_shading=True, **kwargs)

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

    def plot_quadrangle(self, bottom_left, top_right=None, width: u.deg = None, height: u.deg = None, **kwargs):
        """
        Plots a quadrangle on the given map.
        This draws a quadrangle that has corners at ``(bottom_left, top_right)``,
        if ``width`` and ``height`` are specified, they are respectively added to the
        longitude and latitude of the ``bottom_left`` coordinate to calculate a
        ``top_right`` coordinate.

        Parameters
        ----------
        bottom_left :`~astropy.coordinates.SkyCoord`
            The bottom-left coordinate of the quadrangle. It can
            have shape ``(2,)`` to simultaneously define ``top_right``.
        top_right : `~astropy.coordinates.SkyCoord`
            The top-right coordinate of the quadrangle.
        width : `astropy.units.Quantity`, optional
            The width of the quadrangle. Required if ``top_right`` is omitted.
        height : `astropy.units.Quantity`
            The height of the quadrangle. Required if ``top_right`` is omitted.
        **kwargs : Keyword arguments are handed to `pyvista.Plotter.add_mesh`.
        """
        bottom_left, top_right = get_rectangle_coordinates(
            bottom_left, top_right=top_right, width=width, height=height)
        width = Longitude(top_right.spherical.lon - bottom_left.spherical.lon)
        height = top_right.spherical.lat - bottom_left.spherical.lat

        quadrangle_patch = Quadrangle((bottom_left.lon, bottom_left.lat), width, height, resolution=1000)
        quadrangle_coordinates = quadrangle_patch.get_xy()
        c = SkyCoord(quadrangle_coordinates[:, 0]*u.deg, quadrangle_coordinates[:, 1]*u.deg, frame=bottom_left.frame)
        c.transform_to(self.coordinate_frame)
        mesh = self._coords_to_xyz(c)
        self.plotter.add_mesh(mesh, **kwargs)

    def plot_field_lines(self, field_lines, **kwargs):
        """
        Plots the field lines from `pfsspy`.

        Parameters
        ----------
        field_lines : `pfsspy.fieldline.FieldLines`
            Field lines to be plotted.
        **kwargs :
            Keyword arguments are handed to `pyvista.Plotter.add_mesh`.
        """
        for field_line in field_lines:
            grid = self._coords_to_xyz(field_line.coords.ravel())
            field_line_mesh = pv.StructuredGrid(grid[:, 0], grid[:, 1], grid[:, 2])
            color = {0: 'black', -1: 'tab:blue', 1: 'tab:red'}.get(field_line.polarity)
            self.plotter.add_mesh(field_line_mesh, color=color, **kwargs)
