0.1.0 (2022-08-04)
==================

Features
--------

- Creation of :class:`~sunkit_pyvista.plotter.SunpyPlotter` class which allows for plotting of a :class:`~sunpy.map.GenericMap` through pyvista. (`#4 <https://github.com/sunpy/sunkit-pyvista/pull/4>`__)
- Adds :meth:`~sunkit_pyvista.plotter.SunpyPlotter.set_camera_coordinate` which allows for specifying the initial camera coordinates. (`#10 <https://github.com/sunpy/sunkit-pyvista/pull/10>`__)
- Added :meth:`~sunkit_pyvista.plotter.SunpyPlotter.plot_field_lines` method which allows for plotting of magnetic field lines from `pfsspy` through pyvista. (`#12 <https://github.com/sunpy/sunkit-pyvista/pull/12>`__)
- Adds :meth:`~sunkit_pyvista.plotter.SunpyPlotter.set_view_angle` which allows for setting of the camera's width/view angle. (`#16 <https://github.com/sunpy/sunkit-pyvista/pull/16>`__)
- Adds :meth:`~sunkit_pyvista.plotter.SunpyPlotter.plot_quadrangle` which allows for drawing a quadrangle defined
  by coordinates passed to it. (`#17 <https://github.com/sunpy/sunkit-pyvista/pull/17>`__)
- Adds :meth:`~sunkit_pyvista.plotter.SunpyPlotter._add_mesh_to_dict` stores each mesh in a dictionary for later access. (`#24 <https://github.com/sunpy/sunkit-pyvista/pull/24>`__)
- :meth:`~sunkit_pyvista.plotter.SunpyPlotter.plot_map` accepts the keyword ``clip_interval`` as argument, which clips the data
  according to the percentile interval bounded by the two numbers. (`#26 <https://github.com/sunpy/sunkit-pyvista/pull/26>`__)
- ``plot_line()`` is renamed to :meth:`~sunkit_pyvista.plotter.SunpyPlotter.plot_coordinates`
  which also allows for a sphere with given ``radius`` to be plotted if a single coordinate is passed to it. (`#29 <https://github.com/sunpy/sunkit-pyvista/pull/29>`__)
- Added :meth:`~sunkit_pyvista.plotter.SunpyPlotter.save` method which allows for saving of plots to a vtm file. (`#37 <https://github.com/sunpy/sunkit-pyvista/pull/37>`__)
- Allows for figure tests to be performed with `pytest`. (`#38 <https://github.com/sunpy/sunkit-pyvista/pull/38>`__)
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
