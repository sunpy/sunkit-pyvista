import functools
from pathlib import Path

import numpy as np
import pyvista as pv

import astropy.units as u
import sunpy.visualization.colormaps as cm
from astropy.constants import R_sun
from astropy.coordinates import Longitude, SkyCoord
from astropy.visualization import AsymmetricPercentileInterval
from sunpy.coordinates import HeliocentricInertial
from sunpy.coordinates.utils import get_rectangle_coordinates
from sunpy.map.maputils import all_corner_coords_from_map
from sunpy.visualization._quadrangle import Quadrangle

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

    Attributes
    ----------
    all_meshes : `dict`
        Stores a reference to all the plotted meshes in a dictionary.
    """
    def __init__(self, coordinate_frame=None):
        if coordinate_frame is None:
            coordinate_frame = HeliocentricInertial()
        self._coordinate_frame = coordinate_frame
        self._plotter = pv.Plotter()
        self.camera = self._plotter.camera
        self.all_meshes = {}
        self.color_maps = {v: k for k, v in cm.cmlist.items()}

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

    def _add_mesh_to_dict(self, block_name, mesh):
        """
        Adds all of the meshes to a `dict`
        that stores a reference to the meshes.
        """
        if block_name in self.all_meshes:
            self.all_meshes[block_name].append(mesh)
        else:
            self.all_meshes[block_name] = [mesh]

    def _coords_to_xyz(self, coords):
        """
        Transform coordinates to x, y, z values used internally by pyvista.

        Parameters
        ----------
        coords : astropy.coordinates.SkyCoord

        Returns
        -------
        numpy.array
            This is the same shape as the input *coords* array, with an extra
            len(3) axis at the end which corresponds to the (x, y, z) components
            of the coordinates.
        """
        coords = coords.transform_to(self.coordinate_frame)
        coords.representation_type = 'cartesian'
        return np.stack((coords.x.to_value(R_sun),
                         coords.y.to_value(R_sun),
                         coords.z.to_value(R_sun)),
                        axis=-1)

    def _get_clim(self, data, clip_interval):
        """
        Get vmin, vmax of a data slice when clip_interval is specified.
        """
        percent_limits = clip_interval.to('%').value
        vmin, vmax = AsymmetricPercentileInterval(*percent_limits).get_limits(data)
        return [vmin, vmax]

    def set_camera_coordinate(self, coord):
        """
        Set the camera position.

        Parameters
        ----------
        coords : `astropy.coordinates.SkyCoord`
            Camera coordinate.
        """
        camera_position = self._coords_to_xyz(coord)
        pos = tuple(camera_position)
        self.plotter.camera.position = pos

    @u.quantity_input
    def set_view_angle(self, angle: u.deg):
        """
        Set the camera view angle.

        Parameters
        ----------
        angle : `astropy.units.Quantity`
            The viewing angle.
        """
        view_angle = angle.to_value(u.deg)
        if not (view_angle > 0 and view_angle <= 180):
            raise ValueError("specified view angle must be "
                             "0 deg < angle <= 180 deg")
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
        nodes = self._coords_to_xyz(corner_coords)
        grid = pv.StructuredGrid(nodes[:, :, 0], nodes[:, :, 1], nodes[:, :, 2])
        grid['data'] = m.plot_settings['norm'](m.data.T.ravel())
        return grid

    @u.quantity_input
    def plot_map(self, m, clip_interval: u.percent = None, **kwargs):
        """
        Plot a sunpy map.

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
        map_mesh = self._pyvista_mesh(m)
        if clip_interval is not None:
            if len(clip_interval) == 2:
                clim = self._get_clim(data=map_mesh['data'],
                                      clip_interval=clip_interval)
            else:
                raise ValueError("Clip percentile interval must be "
                                 "specified as two numbers.")
        else:
            clim = [0, 1]

        cmap_name = self.color_maps[cmap]
        map_mesh.add_field_array([cmap_name], 'color')
        self.plotter.add_mesh(map_mesh, cmap=cmap, clim=clim, **kwargs)
        self._add_mesh_to_dict(block_name='maps', mesh=map_mesh)

    def plot_coordinates(self, coords, radius=0.05, **kwargs):
        """
        Plot a sphere if a single coordinate is passed and plots a line if
        multiple coordinates are passed.

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
        if points.ndim == 1 or (points.ndim == 2 and points.shape[0] == 1):
            # Single coordinate
            point_mesh = pv.Sphere(radius=radius, center=points)
        else:
            point_mesh = pv.Spline(points)

        color = kwargs.get('color', np.nan)
        point_mesh.add_field_array([color], 'color')
        self.plotter.add_mesh(point_mesh, smooth_shading=True, **kwargs)
        self._add_mesh_to_dict(block_name='coordinates', mesh=point_mesh)

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
        arrow_mesh = pv.Arrow(start=(0, 0, -length / 2),
                              direction=(0, 0, length),
                              scale='auto',
                              **defaults)
        color = kwargs.get('color', np.nan)
        arrow_mesh.add_field_array([color], 'color')
        self.plotter.add_mesh(arrow_mesh, **kwargs)
        self._add_mesh_to_dict(block_name='solar_axis', mesh=arrow_mesh)

    def plot_quadrangle(self, bottom_left, top_right=None, width: u.deg = None, height: u.deg = None, **kwargs):
        """
        Plot a quadrangle.

        This draws a quadrangle that has corners at ``(bottom_left, top_right)``,
        if ``width`` and ``height`` are specified, they are respectively added to
        the longitude and latitude of the ``bottom_left`` coordinate to calculate
        a ``top_right`` coordinate.

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
        quad_grid = self._coords_to_xyz(c)
        quad_block = pv.Spline(quad_grid)
        color = kwargs.get('color', np.nan)
        radius = kwargs.get('radius', 0.01)
        quad_block = quad_block.tube(radius=radius)
        quad_block.add_field_array([color], 'color')
        self.plotter.add_mesh(quad_block, **kwargs)
        self._add_mesh_to_dict(block_name='quadrangles', mesh=quad_block)

    def plot_field_lines(self, field_lines, **kwargs):
        """
        Plot magnetic field lines from `pfsspy`.

        Parameters
        ----------
        field_lines : `pfsspy.fieldline.FieldLines`
            Field lines to be plotted.
        **kwargs :
            Keyword arguments are handed to `pyvista.Plotter.add_mesh`.
        """
        field_line_meshes = pv.MultiBlock([])
        for field_line in field_lines:
            grid = self._coords_to_xyz(field_line.coords.ravel())
            field_line_mesh = pv.StructuredGrid(grid[:, 0], grid[:, 1], grid[:, 2])
            color = {0: 'black', -1: 'tab:blue', 1: 'tab:red'}.get(field_line.polarity)
            field_line_mesh.add_field_array([color], 'color')
            self.plotter.add_mesh(field_line_mesh, color=color, **kwargs)
            field_line_meshes.append(field_line_mesh)

        self._add_mesh_to_dict(block_name='field_lines', mesh=field_line_meshes)

    def save(self, filepath, overwrite=False):
        """
        Save all the meshes.

        This saves the rendered plot as a vtm extended file as well as a directory
        of the individual meshes with the specified name.

        Parameters
        ----------
        filepath : `str` or `pathlib.Path`
            Name of the file to save as, should have vtm or vtmb as an extension.

        Examples
        --------
        >>> from sunkit_pyvista import SunpyPlotter
        >>> plotter = SunpyPlotter()
        >>> plotter.plot_solar_axis()
        >>> plotter.save('./filename.vtm') # doctest: +SKIP

        """
        file_path = Path(filepath)
        directory_path = file_path.with_suffix('')

        if not overwrite:
            if file_path.is_file():
                raise ValueError(f"VTM file '{directory_path.absolute()}' already exists")
        if directory_path.exists():
            raise ValueError(f"Directory '{directory_path.absolute()}' already exists")

        mesh_block = pv.MultiBlock()
        for objects in self.all_meshes:
            for meshes in self.all_meshes[objects]:
                mesh_block.append(meshes)
        mesh_block.save(file_path)

    def _loop_through_meshes(self, mesh_block):
        """
        Recursively loop to add nested `~pyvista.core.MultiBlock` to the `pyvsita.Plotter`
        along with the color of the mesh.
        """
        for block in mesh_block:
            if isinstance(block, pv.MultiBlock):
                self._loop_through_meshes(block)
            else:
                color = dict(block.field_arrays).pop('color', [None])[0]

                if not isinstance(color, str):
                    color = None
                if color in cm.cmlist:
                    self.plotter.add_mesh(block, cmap=color)
                else:
                    self.plotter.add_mesh(block, color=color)

    def load(self, filepath):
        """
        Load saved meshes into this plotter.

        Parameters
        ----------
        filepath : `str` or `pathlib.Path`
            Name of the file to load, should have vtm or vtmb as an extension.
        """
        file_path = Path(filepath)
        mesh_block = pv.read(file_path)
        self._loop_through_meshes(mesh_block)
