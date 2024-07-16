Getting started
===============

Once you have everything installed, you can start off by importing the main plotting class::

    from sunkit_pyvista import SunpyPlotter

From here, you will need to initialize a plotter object::

    plotter = SunpyPlotter()

This has the methods required to plot data.
For example, :meth:`~sunkit_pyvista.SunpyPlotter.plot_map` will plot a sunpy Map::

    plotter.plot_map(map)

You can plot a fixed coordinate as a sphere using :meth:`~sunkit_pyvista.SunpyPlotter.plot_coordinates`::

    plotter.plot_coordinates(sky_coord)

You can plot the solar rotation axis using :meth:`~sunkit_pyvista.SunpyPlotter.plot_solar_axis`::

    plotter.plot_solar_axis()

You can plot a quadrangle using  :meth:`~sunkit_pyvista.SunpyPlotter.plot_quadrangle`::

    plotter.plot_quadrangle(bottom_left, width, height)

You can find an example of this in ``Extending functionality from sunpy``.

You can plot a field lines from `sunkit_magex` using :meth:`~sunkit_pyvista.SunpyPlotter.plot_field_lines`::

    plotter.plot_field_lines(field_lines)

You can find an example of this in ``Plotting Field Lines from pfsspy``.

You can plot the solar limb using :meth:`~sunkit_pyvista.SunpyPlotter.plot_limb`::

    plotter.plot_limb(sunpy_map)

Changing the camera coordinate, you can use :meth:`~sunkit_pyvista.SunpyPlotter.set_camera_coordinate`::

    plotter.set_camera_coordinate(sky_coord)

Finally you will want to show the plot::

    plotter.show()

In addition, saving and loading of the plotted meshes is supported::

    plotter.save(filepath)

and::

    plotter.load(filepath)

Hopefully, this gives you a good idea of how to use the plotting class.
