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
from sunpy.coordinates import frames
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
# Plot a map
plotter.plot_map(m, clip_interval=(0, 99)*u.percent)
# Add an arrow to show the solar rotation axis
plotter.plot_solar_axis()
# Plot an arbitrary line
line = SkyCoord(lon=[180, 190, 200] * u.deg,
                lat=[0, 10, 20] * u.deg,
                distance=[1, 2, 3] * const.R_sun,
                frame='heliocentricinertial')
plotter.plot_coordinates(line)

# Plot an arbitrary point
point = SkyCoord(lon=30 * u.deg,
                 lat=-10 * u.deg,
                 obstime=m.date,
                 frame=frames.HeliographicStonyhurst)
plotter.plot_coordinates(point, color='yellow')

plotter.set_camera_coordinate(m.observer_coordinate)

# Set the view angle of the plot
plotter.set_view_angle(m.rsun_obs)

# Plot a quadrangle with width of 20 degrees and a height of 60 degrees
bottom_left = SkyCoord(30*u.deg, -10*u.deg,
                       frame=frames.HeliographicStonyhurst,
                       obstime=m.date)
plotter.plot_quadrangle(bottom_left=bottom_left, width=20*u.deg, height=60*u.deg, color='blue')
plotter.show()
