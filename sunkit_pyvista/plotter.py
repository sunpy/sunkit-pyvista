import contextlib
from pathlib import Path

import numpy as np
import pyvista as pv
from matplotlib import colors

import astropy.units as u
from astropy.constants import R_sun
from astropy.coordinates import Longitude, SkyCoord
from astropy.visualization import AsymmetricPercentileInterval
from astropy.visualization.wcsaxes import Quadrangle

from streamtracer import StreamTracer, VectorGrid
from sunkit_magex.pfss.coords import strum2cart
from sunpy.coordinates import HeliocentricInertial
from sunpy.coordinates.screens import SphericalScreen
from sunpy.coordinates.utils import get_rectangle_coordinates
from sunpy.map.maputils import all_corner_coords_from_map

from sunkit_pyvista.utils import get_limb_coordinates

__all__ = ["SaveMixIn", "SunpyPlotter", "CartesianPlotter"]


class SaveMixIn:
    """
    A mixin class to add savefig used for pytest-mpl figure comparison testing.
    """

    def savefig(self, filename, **kwargs):
        """
        Save a screenshot of the current figure.

        This is a wrapper around `pyvista.Plotter.screenshot` and was written
        so we use pytest-mpl's figure comparison testing.

        Parameters
        ----------
        filename : str
            The name of the file to save the screenshot to.
        kwargs : dict
            Keyword arguments are ignored.
        """
        self.screenshot(filename)


class SunpyPlotter(SaveMixIn, pv.Plotter):
    """
    A plotter for 3D data.

    This class inherits `pyvista.Plotter` so we can provide coordinate-aware plotting.
    For now, all coordinates are converted to
    a specific frame (`~sunpy.coordinates.HeliocentricInertial` by default),
    and distance units are such that :math:`R_{sun} = 1`.

    Parameters
    ----------
    coordinate_frame : `astropy.coordinates.BaseCoordinateFrame`
        Coordinate frame of the plot. The x, y, z axes of the pyvista plotter
        will be the x, y, z axes in this coordinate system.
    obstime : `astropy.time.Time`
        The obstime to use for the default coordinate frame if
        ``coordinate_frame=`` is not specified.  Must not be specified if
        ``coordinate_frame`` is given.
    kwargs : dict
        All other keyword arguments are passed through to `pyvista.Plotter`.

    Attributes
    ----------
    all_meshes : `dict`
        Stores a reference to all the plotted meshes in a dictionary.
    """

    def __init__(self, *args, coordinate_frame=None, obstime=None, **kwargs):
        super().__init__(*args, **kwargs)
        if coordinate_frame is not None and obstime is not None:
            msg = "Only coordinate_frame or obstime can be specified, not both."
            raise ValueError(msg)
        if coordinate_frame is None:
            coordinate_frame = HeliocentricInertial(obstime=obstime)
        self._coordinate_frame = coordinate_frame
        self.all_meshes = {}

    @property
    def coordinate_frame(self):
        """
        Coordinate frame of the plot.
        """
        return self._coordinate_frame

    def _extract_color(self, mesh_kwargs):
        """
        Converts a given color string to it's equivalent rgb tuple.

        Parameters
        ----------
        mesh_kwargs : dict

        Returns
        -------
        tuple
            A tuple containing the (r, g, b) values from the strings passed
            to it. Defaults to (255, 255, 255) - white.
        """
        color_string = mesh_kwargs.pop("color", "white")
        return colors.to_rgb(color_string)

    def _add_mesh_to_dict(self, block_name, mesh):
        """
        Adds all of the meshes to a `dict` that stores a reference to the
        meshes.
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
        coords.representation_type = "cartesian"
        return np.stack(
            (
                coords.x.to_value(R_sun),
                coords.y.to_value(R_sun),
                coords.z.to_value(R_sun),
            ),
            axis=-1,
        )

    def _get_clim(self, data, clip_interval):
        """
        Get vmin, vmax of a data slice when clip_interval is specified.
        """
        percent_limits = clip_interval.to("%").value
        vmin, vmax = AsymmetricPercentileInterval(*percent_limits).get_limits(data)
        return [vmin, vmax]

    def coordinates_to_polydata(self, coords):
        """
        Convert a set of coordinates in a `~astropy.coordinates.SkyCoord` to a
        `pyvista.PolyData` mesh.

        Parameters
        ----------
        coords : `astropy.coordinates.SkyCoord`
            Coordinates to convert.

        Returns
        -------
        `pyvista.PolyData`
        """
        points = self._coords_to_xyz(coords)
        pd = pv.PolyData(points.reshape(-1, 3))
        self._add_mesh_to_dict("polydata", pd)
        return pd

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
        self.camera.position = pos

    def set_camera_focus(self, coord):
        """
        Set the camera focus.

        Parameters
        ----------
        coords : `astropy.coordinates.SkyCoord`
            Camera coordinate.
        """
        camera_position = self._coords_to_xyz(coord)
        pos = tuple(camera_position)
        self.set_focus(pos)

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
            msg = "specified view angle must be 0 deg < angle <= 180 deg"
            raise ValueError(msg)
        zoom_value = self.camera.view_angle / view_angle
        self.camera.zoom(zoom_value)

    def _map_to_mesh(self, m, *, assume_spherical=True):
        """
        Create a mesh from a map.

        Parameters
        ----------
        m : `sunpy.map.Map`
            The map to use.
        assume_spherical : bool, optional
            If `False`, only on-limb pixels are returned.

        Returns
        -------
        `pyvista.StructuredGrid`
        """
        corner_coords = all_corner_coords_from_map(m)

        if assume_spherical:
            context = SphericalScreen(
                m.observer_coordinate,
                only_off_disk=True,
            )
        else:
            context = contextlib.nullcontext()
        with context:
            # Convert SkyCoord to xyz
            verts = self._coords_to_xyz(corner_coords)
        nx, ny = verts.shape[:2]
        nverts = nx * ny
        verts = verts.reshape(nverts, 3)

        # Get vertex incices for each face
        vert_indices = np.arange(nverts).reshape(nx, ny)
        lower_left = vert_indices[:-1, :-1]
        lower_right = vert_indices[1:, :-1]
        upper_right = vert_indices[1:, 1:]
        upper_left = vert_indices[:-1, 1:]

        nfaces = (nx - 1) * (ny - 1)
        faces = np.column_stack(
            [
                np.ones(nfaces).astype(int) * 4,
                lower_left.ravel(),
                lower_right.ravel(),
                upper_right.ravel(),
                upper_left.ravel(),
            ],
        )
        # Remove faces that don't have a finite vertex
        # this can often happen with off-limb vertices)
        finite = np.sum(np.isfinite(verts[faces[:, 1], :]), axis=1) == 3
        finite = finite & (np.sum(np.isfinite(verts[faces[:, 2], :]), axis=1) == 3)
        finite = finite & (np.sum(np.isfinite(verts[faces[:, 3], :]), axis=1) == 3)
        finite = finite & (np.sum(np.isfinite(verts[faces[:, 4], :]), axis=1) == 3)
        faces = faces[finite, :]

        # Remove any non-finite vertices. This reduces the size of the mesh,
        # and is also needed for the ipygany backend to work which is used for
        # showing examples in the documentation
        vert_mask = np.isfinite(verts[:, 0])
        finite_indices = np.cumsum(vert_mask) - 1
        # Re-map face indices
        faces[:, 1:] = finite_indices[faces[:, 1:]]
        # Remove non-finite vertices
        verts = verts[vert_mask, :]
        grid = pv.PolyData(verts, faces.ravel())
        grid["data"] = m.plot_settings["norm"](m.data.ravel()[finite])
        return grid

    @u.quantity_input
    def plot_map(
        self,
        m,
        *,
        clip_interval: u.percent = None,
        assume_spherical_screen=True,
        **kwargs,
    ):
        """
        Plot a sunpy map.

        Parameters
        ----------
        m : `sunpy.map.Map`
            Map to be plotted.
        clip_interval : two-element `~astropy.units.Quantity`, optional
            If provided, the data will be clipped to the percentile
            interval bounded by the two numbers.
        assume_spherical_screen : bool, optional
            If `True` (default) then off-limb pixels are plotted using
            `sunpy.coordinates.screens.SphericalScreen`.
            If `False`, off-limb pixels are not plotted.
        **kwargs :
            Keyword arguments are handed to `pyvista.Plotter.add_mesh`.
        """
        map_mesh = self._map_to_mesh(m, assume_spherical=assume_spherical_screen)
        if clip_interval is not None:
            if len(clip_interval) == 2:
                clim = self._get_clim(
                    data=map_mesh["data"],
                    clip_interval=clip_interval,
                )
            else:
                msg = "Clip percentile interval must be specified as two numbers."
                raise ValueError(
                    msg,
                )
        else:
            clim = [0, 1]
        cmap = self._get_cmap(kwargs, m)
        kwargs.setdefault("show_scalar_bar", False)
        self.add_mesh(map_mesh, cmap=cmap, clim=clim, **kwargs)
        map_mesh.add_field_data([cmap], "cmap")
        self._add_mesh_to_dict(block_name="maps", mesh=map_mesh)

    @staticmethod
    def _get_cmap(kwargs, m):
        """
        Get the colormap as a string.

        Returns
        -------
        str
        """
        return kwargs.pop("cmap", m.plot_settings["cmap"])

    def plot_coordinates(self, coords, *, radius=0.05, **kwargs):
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

        color = self._extract_color(kwargs)
        point_mesh.add_field_data(color, "color")

        kwargs["render_lines_as_tubes"] = kwargs.pop("render_lines_as_tubes", True)
        self.add_mesh(point_mesh, color=color, smooth_shading=True, **kwargs)
        self._add_mesh_to_dict(block_name="coordinates", mesh=point_mesh)

    def plot_solar_axis(self, *, length=2.5, arrow_kwargs=None, **kwargs):
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
        if arrow_kwargs is None:
            arrow_kwargs = {}
        defaults = {"shaft_radius": 0.01, "tip_length": 0.05, "tip_radius": 0.02}
        defaults.update(arrow_kwargs)
        arrow_mesh = pv.Arrow(
            start=(0, 0, -length / 2),
            direction=(0, 0, length),
            scale="auto",
            **defaults,
        )
        color = self._extract_color(kwargs)
        arrow_mesh.add_field_data(color, "color")
        self.add_mesh(arrow_mesh, color=color, **kwargs)
        self._add_mesh_to_dict(block_name="solar_axis", mesh=arrow_mesh)

    def plot_quadrangle(
        self,
        bottom_left,
        *,
        top_right=None,
        width: u.deg = None,
        height: u.deg = None,
        radius=0.01,
        **kwargs,
    ):
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
        radius : `float`
            Radius of the `pyvista.Spline` used to create the quadrangle.
            Defaults to ``0.01`` times the radius of the sun.
        **kwargs : Keyword arguments are handed to `pyvista.Plotter.add_mesh`.
        """
        bottom_left, top_right = get_rectangle_coordinates(
            bottom_left,
            top_right=top_right,
            width=width,
            height=height,
        )
        width = Longitude(top_right.spherical.lon - bottom_left.spherical.lon)
        height = top_right.spherical.lat - bottom_left.spherical.lat

        quadrangle_patch = Quadrangle(
            (bottom_left.lon, bottom_left.lat),
            width,
            height,
            resolution=1000,
        )
        quadrangle_coordinates = quadrangle_patch.get_xy()
        c = SkyCoord(
            quadrangle_coordinates[:, 0] * u.deg,
            quadrangle_coordinates[:, 1] * u.deg,
            frame=bottom_left.frame,
        )
        c.transform_to(self.coordinate_frame)
        quad_grid = self._coords_to_xyz(c)
        quad_block = pv.Spline(quad_grid)
        radius = kwargs.get("radius", 0.01)
        quad_block = quad_block.tube(radius=radius)
        color = self._extract_color(kwargs)
        quad_block.add_field_data(color, "color")
        self.add_mesh(quad_block, color=color, **kwargs)
        self._add_mesh_to_dict(block_name="quadrangles", mesh=quad_block)

    def plot_field_lines(self, field_lines, *, color_func=None, **kwargs):
        """
        Plot magnetic field lines.

        Parameters
        ----------
        field_lines : `sunkit_magex.pfss.fieldline.FieldLines`
            Field lines to be plotted.
        color_func : function
            Function to get the color for each field line.
            If not given, defaults to showing closed field lines in black,
            and open field lines in blue (positive polarity) or red (negative polarity).
            The function must have the signature::

                def color_func(field_line: sunkit_magex.pfss.fieldline.FieldLine) -> color:

             Where ``color`` is any color that `pyvista` recognises
             (e.g. a RGBA tuple, a RGB tuple, a color string)
        **kwargs :
            Keyword arguments are handed to `pyvista.Plotter.add_mesh`.
        """
        if not color_func:

            def color_func(field_line):
                color = {0: "black", -1: "tab:blue", 1: "tab:red"}.get(
                    field_line.polarity,
                )
                return colors.to_rgb(color)

        field_line_meshes = pv.MultiBlock([])
        for field_line in field_lines:
            xyz = self._coords_to_xyz(field_line.coords.ravel())
            if not xyz.size:
                continue
            spline = pv.Spline(xyz)
            color = color_func(field_line)
            if isinstance(color, tuple):
                color = list(color)
                if len(color) == 4:
                    kwargs["opacity"] = color[3]
                    color = color[:3]

            spline.add_field_data([color], "color")

            kwargs["render_lines_as_tubes"] = kwargs.pop("render_lines_as_tubes", True)
            kwargs["line_width"] = kwargs.pop("line_width", 2)
            self.add_mesh(spline, color=color, **kwargs)
            field_line_meshes.append(spline)

        self._add_mesh_to_dict(block_name="field_lines", mesh=spline)

    def plot_current_sheet(self, pfss_out, **kwargs):
        """
        Plot current sheet, where :math:`B_r=0`.

        Parameters
        ----------
        pfss_out : `sunkit_magex.pfss.output.Output`
            Magnetic field calculated by `sunkit_magex`.
        **kwargs :
            Keyword arguments are passed to `pyvista.Plotter.add_mesh`.
        """
        sc_vect = pfss_out.grid.sc
        pc_vect = np.insert(pfss_out.grid.pc, 0, pfss_out.grid.pc[-1])
        rg_vect = pfss_out.grid.rg
        coord = SkyCoord((pfss_out.input_map.meta["crval1"] * u.deg).to(u.rad), 0.0 * u.rad, frame=HeliocentricInertial)
        coord = coord.transform_to(self.coordinate_frame)
        bc_r = np.insert(pfss_out.bc[0], 0, pfss_out.bc[0][-1], axis=0)
        [s_arr, p_arr, r_arr] = np.meshgrid(sc_vect, pc_vect, rg_vect)
        x_arr, y_arr, z_arr = strum2cart(r_arr, s_arr, p_arr + coord.lon.to(u.rad).value)
        pfss_out_pv = pv.StructuredGrid(x_arr, y_arr, z_arr)
        pfss_out_pv["Br"] = bc_r.ravel("F")
        isos_br = pfss_out_pv.contour(isosurfaces=1, rng=[0, 0])
        self.add_mesh(isos_br, **kwargs)
        self._add_mesh_to_dict(block_name="current_sheet", mesh=isos_br)

    def save(self, filepath, *, overwrite=False):
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
        >>> plotter.save("./filename.vtm")  # doctest: +SKIP
        """
        file_path = Path(filepath)
        directory_path = file_path.with_suffix("")

        if not overwrite and file_path.is_file():
            msg = f"VTM file '{directory_path.absolute()}' already exists"
            raise ValueError(
                msg,
            )
        if directory_path.exists():
            msg = f"Directory '{directory_path.absolute()}' already exists"
            raise ValueError(msg)

        mesh_block = pv.MultiBlock()
        for objects in self.all_meshes:
            for meshes in self.all_meshes[objects]:
                mesh_block.append(meshes)
        mesh_block.save(file_path)

    def _loop_through_meshes(self, mesh_block):
        """
        Recursively loop to add nested `~pyvista.core.MultiBlock` to the
        `pyvsita.Plotter` along with the color of the mesh.
        """
        for block in mesh_block:
            if isinstance(block, pv.MultiBlock):
                self._loop_through_meshes(block)
            else:
                color = dict(block.field_data).get("color", None)
                cmap = dict(block.field_data).get("cmap", [None])[0]
                self.add_mesh(block, color=color, cmap=cmap)

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

    def plot_limb(self, m, *, radius=0.02, **kwargs):
        """
        Draws the solar limb as seen by the map's observer.

        Parameters
        ----------
        m : `sunpy.map.Map`
                Map's limb to be plotted.
        radius : `float`
            Radius of the `pyvista.Spline` used to create the limb.
            Defaults to ``0.02`` times the radius of the sun.
        **kwargs : Keyword arguments are handed to `pyvista.Plotter.add_mesh`.
        """
        limb_coordinates = get_limb_coordinates(
            m.observer_coordinate,
            rsun=m.rsun_meters,
            resolution=1000,
        )
        limb_coordinates.transform_to(self.coordinate_frame)
        limb_grid = self._coords_to_xyz(limb_coordinates)
        limb_block = pv.Spline(limb_grid)
        color = self._extract_color(mesh_kwargs=kwargs)
        limb_block = limb_block.tube(radius=radius)
        limb_block.add_field_data(color, "color")
        self.add_mesh(limb_block, color=color, **kwargs)
        self._add_mesh_to_dict(block_name="limbs", mesh=limb_block)


class CartesianPlotter(SaveMixIn, pv.Plotter):
    """
    A plotter for 3D data in a Cartesian box.

    This class inherits `pyvista.Plotter`. It is used to visualize 3D vector field lines
    (e.g., magnetic field lines of solar active regions) traced by `streamtracer`.

    Parameters
    ----------
    kwargs : dict
        All other keyword arguments are passed through to `pyvista.Plotter`.
    """

    def define_vector_field(self, vectors, *args, **kwargs):
        """
        Define a 3D vector field.

        Parameters
        ----------
        vectors : numpy.ndarray
            A (nx, ny, nz, 3) array representing the 3D vector field.
        args: list
            arguments for `streamtracer.VectorGrid`.
        kwargs: dict
            keyword arguments for `streamtracer.VectorGrid`.
        """
        if kwargs.get("grid_spacing") is None:
            if kwargs.get("grid_coords") is None:
                kwargs["grid_spacing"] = [1, 1, 1]
        self._grid = VectorGrid(vectors.astype(np.float64), *args, **kwargs)
        # Create a `pyvista.StructuredGrid` mesh from the 3D vector field.
        x, y, z = np.meshgrid(self._grid.xcoords, self._grid.ycoords, self._grid.zcoords, indexing="ij")
        mesh = pv.StructuredGrid(x, y, z)
        vectors = self._grid.vectors.transpose(2, 1, 0, 3).reshape(-1, 3)
        mesh["vectors"] = vectors
        mesh.active_vectors_name = "vectors"
        magnitudes = np.linalg.norm(vectors, axis=-1)
        mesh["magnitudes"] = magnitudes
        mesh.active_scalars_name = "magnitudes"
        self._mesh = mesh

    def show_outline(self, color="black", **kwargs):
        """
        Show the outline of the 3D vector field.

        Parameters
        ----------
        color : str, optional
            The color of the outline.
            Default is 'black'.
        kwargs : dict
            Keyword arguments for `pyvista.Plotter.add_mesh`.
        """
        self.add_mesh(self._mesh.outline(), color=color, **kwargs)

    def show_boundary(self, boundary="bottom", *, component=2, cmap="gray", **kwargs):
        """
        Show the boundary of the 3D vector field.

        Parameters
        ----------
        boundary : str, optional
            The boundary to be plotted.
            'bottom', 'top', 'left', 'right', 'front', or 'back'.
            Default is 'bottom'.
        component : int, optional
            The component of the vector field to be plotted.
            0, 1, or 2 for x, y, or z component, respectively.
            Default is 2.
        cmap : str, optional
            The colormap for the boundary.
            Default is 'gray'.
        kwargs : dict
            Keyword arguments for `pyvista.Plotter.add_mesh`.
        """
        nx, ny, nz = self._mesh.dimensions
        x_min, y_min, z_min = 0, 0, 0
        x_max, y_max, z_max = nx - 1, ny - 1, nz - 1

        if boundary == "bottom":
            subset = (x_min, x_max, y_min, y_max, z_min, z_min)
        elif boundary == "top":
            subset = (x_min, x_max, y_min, y_max, z_max, z_max)
        elif boundary == "left":
            subset = (x_min, x_min, y_min, y_max, z_min, z_max)
        elif boundary == "right":
            subset = (x_max, x_max, y_min, y_max, z_min, z_max)
        elif boundary == "front":
            subset = (x_min, x_max, y_min, y_min, z_min, z_max)
        elif boundary == "back":
            subset = (x_min, x_max, y_max, y_max, z_min, z_max)
        else:
            raise ValueError("Invalid boundary. Choose from 'bottom', 'top', 'left', 'right', 'front', or 'back'.")

        surface = self._mesh.extract_subset(subset).extract_surface()
        surface.active_vectors_name = "vectors"
        surface.active_scalars_name = "magnitudes"
        self.add_mesh(surface, scalars="vectors", component=component, cmap=cmap, lighting=False, **kwargs)

    def plot_field_lines(
        self,
        seeds,
        *,
        render_lines_as_tubes=True,
        radius=1,
        max_steps=10000,
        step_size=0.1,
        seeds_config=dict(show_seeds=False, color="red", point_size=5),
        **kwargs,
    ):
        """
        Plot field lines traced from seeds using `streamtracer.StreamTracer`.

        Parameters
        ----------
        seeds : numpy.ndarray
            A (N, 3) array representing the seeds.
        render_lines_as_tubes : bool, optional
            Whether to render field lines as tubes.
            Default is True.
        radius : float, optional
            The radius of the tubes for rendering field lines.
            Default is 1.
        max_steps : int, optional
            The maximum number of steps for tracing field lines.
            Default is 10000.
        step_size : float, optional
            The step size for tracing field lines.
            Default is 0.1.
        seeds_config : dict, optional
            Configuration for plotting seeds.
            Default is dict(show_seeds=False, color='red', point_size=5).
        kwargs : dict
            Keyword arguments for `pyvista.Plotter.add_mesh`.
        """
        tracer = StreamTracer(max_steps, step_size)
        tracer.trace(seeds, self._grid)
        tracer_xs = []
        tracer_xs.append(tracer.xs)
        tracer_xs = [item for sublist in tracer_xs for item in sublist]
        for i, xl in enumerate(tracer_xs):
            assert seeds[i] in xl
            if len(xl) < 2:
                continue
            spline = pv.Spline(xl)
            if render_lines_as_tubes:
                spline = spline.tube(radius=radius)
            self.add_mesh(spline, **kwargs)
        if seeds_config.get("show_seeds"):
            seeds_config.pop("show_seeds")
            self.add_mesh(pv.PolyData(seeds), **seeds_config)
