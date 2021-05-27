"""
=====================
Plot GenericMap in 3
=====================

This example shows how to plot map in 3D.
"""
from sunpy.data.sample import AIA_193_IMAGE
from sunpy.map import Map
from sunkit_pyvista import SunpyPlotter
from xvfbwrapper import Xvfb

vdisplay = Xvfb()
vdisplay.start()

m = Map(AIA_193_IMAGE)
m.plot()

plotter = SunpyPlotter()
map_mesh = plotter.plot_map(m)
line_mesh = plotter.plot_solar_axis()

vdisplay.stop()
# plotter.show(cpos=(-100,0,0))
