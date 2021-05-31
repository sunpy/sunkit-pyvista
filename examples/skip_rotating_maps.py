"""
=======================================
Three dimensional plots with sunpy Maps
=======================================

sunkit-pyvista also allows for rotation of maps to render the initial plot with the
specified angle.
"""

import astropy.constants as const
import astropy.units as u
from astropy.coordinates import SkyCoord
from sunpy.data.sample import AIA_193_IMAGE
from sunpy.map import Map

from sunkit_pyvista import SunpyPlotter

# pv.start_xvfb()

###############################################################################
# Import some sample data
m = Map(AIA_193_IMAGE)

# Start by creating a plotter
plotter = SunpyPlotter()
# Plot a map
plotter.plot_map(m)
# Add an arrow to show the solar rotation axis
plotter.plot_solar_axis()
# Plot an arbitrary line
line = SkyCoord(lon=[180, 190, 200] * u.deg,
                lat=[0, 10, 20] * u.deg,
                distance=[1, 2, 3] * const.R_sun,
                frame='heliocentricinertial')
plotter.plot_line(line)

# `SunpyPlotter` provides a rotate method which accepts an angle. This rotates
# all the rendered meshes to be rotated by that particular angle.
plotter.rotate(angle=30 * u.deg)

plotter.show(cpos=(-100, 0, 0))
