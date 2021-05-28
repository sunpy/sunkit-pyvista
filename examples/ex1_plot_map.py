"""
=======================
Plot a GenericMap in 3D
=======================

This example shows how to plot a GenericMap in 3D.
"""
from sunpy.data.sample import AIA_193_IMAGE
from sunpy.map import Map
from sunkit_pyvista import SunpyPlotter
import pyvista as pv

m = Map(AIA_193_IMAGE)
m.plot()

pv.start_xvfb()
plotter = SunpyPlotter()
map_mesh = plotter.plot_map(m)
line_mesh = plotter.plot_solar_axis()

plotter.show(cpos=(-100,0,0), screenshot="../docs/generated/images/ex1_plot_map.png")

# %%
# .. image:: ../images/ex1_plot_map.png
