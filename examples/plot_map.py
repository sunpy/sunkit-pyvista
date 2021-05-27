"""
======================
Plot a GenericMap in 3D
======================

This example shows how to plot a GenericMap in 3D.
"""
from sunpy.data.sample import AIA_193_IMAGE
from sunpy.map import Map
from sunkit_pyvista import SunpyPlotter
from pyvirtualdisplay import 

display = Display(visible=0, size=(1200, 1000))
display.start()

m = Map(AIA_193_IMAGE)
m.plot()

plotter = SunpyPlotter()
map_mesh = plotter.plot_map(m)
line_mesh = plotter.plot_solar_axis()

display.stop()
# plotter.show(cpos=(-100,0,0))
