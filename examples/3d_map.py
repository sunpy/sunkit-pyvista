"""
=======================================
Three dimensional plots with sunpy Maps
=======================================

Using sunkit-pyvista, one can interface with the `pyvista` package to
produce interactive 3D plots for sunpy Maps.
"""
from sunkit_pyvista import SunpyPlotter
from sunkit_pyvista.sample import low_res_aia_171

###############################################################################
# We will firstly use an AIA 193 image from the sunpy sample data as the base image.
m = low_res_aia_171()

###############################################################################
# 3D plots are done on "plotter" objects, which are similar to matplotlib axes.
# sunkit-pyvista has a built in `SunpyPlotter` class that can be used to plot maps
# and coordinate aware objects.

# Start by creating a plotter
plotter = SunpyPlotter()

# Plot the map
plotter.plot_map(m)
# Add an arrow to show the solar rotation axis
plotter.plot_solar_axis()
plotter.show()
