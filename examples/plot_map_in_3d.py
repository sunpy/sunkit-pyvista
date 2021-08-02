"""
=======================================
Three dimensional plots with sunpy Maps
=======================================

Using sunkit-pyvista, one can interface with the `pyvista` package to
produce interactive 3D plots for sunpy Maps.
"""
import astropy.constants as const
import astropy.units as u
from astropy.coordinates import SkyCoord
from sunpy.data.sample import AIA_193_IMAGE
from sunpy.map import Map

from sunkit_pyvista import SunpyPlotter

###############################################################################
# We will firstly use an AIA 193 image from the sunpy sample data as the base image.
m = Map(AIA_193_IMAGE)

###############################################################################
# 3D plots are done on "plotter" objects, which are similar to matplotlib axes.
# sunkit-pyvista has a built in `SunpyPlotter` class that can be used to plot maps
# and coordinate aware objects.

# Start by creating a plotter
plotter = SunpyPlotter()

# Plot a map by specifying a clip interval to clip data
# according to the percentile interval specified.
plotter.plot_map(m, clip_interval=(0, 99)*u.percent)

# Add an arrow to show the solar rotation axis
plotter.plot_solar_axis()

# Plot an arbitrary line
line = SkyCoord(lon=[180, 190, 200] * u.deg,
                lat=[0, 10, 20] * u.deg,
                distance=[1, 2, 3] * const.R_sun,
                frame='heliocentricinertial')

plotter.plot_coordinates(line)
plotter.show(cpos=(0, 1, 1))
