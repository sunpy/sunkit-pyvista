0.4.0 (2025-07-01)
==================

Breaking Changes
----------------

- Instead of wrapping `pyvista.Plotter`, we now inherit from it.
  This allows us to drop the ``plotter.plotter`` lines in examples and user facing API. (`#138 <https://github.com/sunpy/sunkit-pyvista/pull/138>`__)
- Updated minimum required version of ``sunpy`` to v6.0.0. (`#165 <https://github.com/sunpy/sunkit-pyvista/pull/165>`__)
- Updated minimum required version of ``python`` to v3.11. (`#165 <https://github.com/sunpy/sunkit-pyvista/pull/165>`__)


0.3.0 (2024-06-13)
==================

Breaking Changes
----------------

- Migrated from ``pfsspy`` to ``sunkit-magex`` to handle the magnetic field fieldlines and extrapolation. (`#157 <https://github.com/sunpy/sunkit-pyvista/pull/157>`__)
- Increased minimum version of Python to 3.10 (`#157 <https://github.com/sunpy/sunkit-pyvista/pull/157>`__)


Documentation
-------------

- Added a new gallery example (:ref:`sphx_glr_generated_gallery_floating_sphere.py`) that plots a volume of points from the solar surface. (`#137 <https://github.com/sunpy/sunkit-pyvista/pull/137>`__)


0.2.1 (2023-11-27)
==================

Breaking Changes
----------------

- Due to an incompatibility between ``sunkit-magex`` and ``sunpy`` 5.1.0, ``sunkit-pyvista`` is not compatible with ``sunpy`` 5.1.0.
  The ``sunpy`` version is pinned for the time being.

0.2.0 (2023-11-17)
==================

Breaking Changes
----------------

- Increased minimum versions of dependencies.

  * Minimum Python version is now 3.9.
  * Minimum ``pyvista`` version is now 0.38.4.
  * Minimum ``sunpy`` version is now 5.0.0.
  * Minimum ``pfsspy`` version is now 1.1.2. (`#116 <https://github.com/sunpy/sunkit-pyvista/pull/116>`__)
- Default ``line_width`` for :meth:`~sunkit_pyvista.plotter.SunpyPlotter.plot_field_lines` has been reduced to 2, from 5. (`#117 <https://github.com/sunpy/sunkit-pyvista/pull/117>`__)

Internal Changes
----------------

- Fixed the online documentation to display static images of the examples instead of interactive ones. (`#117 <https://github.com/sunpy/sunkit-pyvista/pull/117>`__)
- Fixed a failing unit test due to upstream deprecation. (`#118 <https://github.com/sunpy/sunkit-pyvista/pull/118>`__)
- Fixed the opacity keyword argument in :meth:`~sunkit_pyvista.plotter.SunpyPlotter.plot_field_lines`. (`#150 <https://github.com/sunpy/sunkit-pyvista/pull/150>`__)


0.1.0 (2022-08-04)
==================

Features
--------

- Creation of :class:`~sunkit_pyvista.plotter.SunpyPlotter` class which allows for plotting of a :class:`~sunpy.map.GenericMap` through pyvista. (`#4 <https://github.com/sunpy/sunkit-pyvista/pull/4>`__)
- Adds :meth:`~sunkit_pyvista.plotter.SunpyPlotter.set_camera_coordinate` which allows for specifying the initial camera coordinates. (`#10 <https://github.com/sunpy/sunkit-pyvista/pull/10>`__)
- Added :meth:`~sunkit_pyvista.plotter.SunpyPlotter.plot_field_lines` method which allows for plotting of magnetic field lines from ``pfsspy`` through pyvista. (`#12 <https://github.com/sunpy/sunkit-pyvista/pull/12>`__)
- Adds :meth:`~sunkit_pyvista.plotter.SunpyPlotter.set_view_angle` which allows for setting of the camera's width/view angle. (`#16 <https://github.com/sunpy/sunkit-pyvista/pull/16>`__)
- Adds :meth:`~sunkit_pyvista.plotter.SunpyPlotter.plot_quadrangle` which allows for drawing a quadrangle defined
  by coordinates passed to it. (`#17 <https://github.com/sunpy/sunkit-pyvista/pull/17>`__)
- Adds :meth:`~sunkit_pyvista.plotter.SunpyPlotter._add_mesh_to_dict` stores each mesh in a dictionary for later access. (`#24 <https://github.com/sunpy/sunkit-pyvista/pull/24>`__)
- :meth:`~sunkit_pyvista.plotter.SunpyPlotter.plot_map` accepts the keyword ``clip_interval`` as argument, which clips the data
  according to the percentile interval bounded by the two numbers. (`#26 <https://github.com/sunpy/sunkit-pyvista/pull/26>`__)
- ``plot_line()`` is renamed to :meth:`~sunkit_pyvista.plotter.SunpyPlotter.plot_coordinates`
  which also allows for a sphere with given ``radius`` to be plotted if a single coordinate is passed to it. (`#29 <https://github.com/sunpy/sunkit-pyvista/pull/29>`__)
- Added :meth:`~sunkit_pyvista.plotter.SunpyPlotter.save` method which allows for saving of plots to a vtm file. (`#37 <https://github.com/sunpy/sunkit-pyvista/pull/37>`__)
- Allows for figure tests to be performed with ``pytest``. (`#38 <https://github.com/sunpy/sunkit-pyvista/pull/38>`__)
- Adds :meth:`~sunkit_pyvista.plotter.SunpyPlotter.plot_limb` which allows for drawing a limb as seen by the map's observer. (`#59 <https://github.com/sunpy/sunkit-pyvista/pull/59>`__)
- Allows user to specify color via a color function to :meth:`~sunkit_pyvista.plotter.SunpyPlotter.plot_field_lines`. (`#70 <https://github.com/sunpy/sunkit-pyvista/pull/70>`__)
- Added narrative documentation for sunkit-pyvista. (`#84 <https://github.com/sunpy/sunkit-pyvista/pull/84>`__)


Bug Fixes
---------

- Adds parenthesis to fix check in :meth:`~sunkit_pyvista.plotter.SunpyPlotter.set_view_angle`. (`#34 <https://github.com/sunpy/sunkit-pyvista/pull/34>`__)
- Fixes error while loading color map in :meth:`~sunkit_pyvista.plotter.SunpyPlotter.load`. (`#55 <https://github.com/sunpy/sunkit-pyvista/pull/55>`__)


Internal Changes
----------------

- Increases test coverage for :class:`~sunkit_pyvista.plotter.SunpyPlotter`. (`#23 <https://github.com/sunpy/sunkit-pyvista/pull/23>`__)
- Rearranged existing examples and added an example brightest pixel with :meth:`~sunkit_pyvista.plotter.SunpyPlotter.plot_coordinates`. (`#30 <https://github.com/sunpy/sunkit-pyvista/pull/30>`__)
- :meth:`~sunkit_pyvista.plotter.SunpyPlotter.plot_quadrangle` uses a :meth`~pyvista.utilities.Spline` for combining the individual points. (`#52 <https://github.com/sunpy/sunkit-pyvista/pull/52>`__)
- Adds an example using :meth:`~sunpy.coordinates.frames.Helioprojective.assume_spherical_screen`. (`#69 <https://github.com/sunpy/sunkit-pyvista/pull/69>`__)
- Changed the manner that colors or colormaps are saved.
  Changed default of meshes to be white. (`#73 <https://github.com/sunpy/sunkit-pyvista/pull/73>`__)
- Removes colorbars when displaying plots. (`#79 <https://github.com/sunpy/sunkit-pyvista/pull/79>`__)
