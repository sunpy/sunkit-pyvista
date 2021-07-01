"""
=====================
Adding extra features
=====================

sunkit-pyvista also extends features from `sunpy`, producing
them in 3D.
"""

import astropy.constants as const
import astropy.units as u
from astropy.constants import R_sun
from astropy.coordinates import SkyCoord
from sunpy.coordinates import frames
from sunpy.data.sample import AIA_171_IMAGE
from sunpy.map import Map

from sunkit_pyvista import SunpyPlotter

###############################################################################
# We will firstly use an AIA 171 image from the sunpy sample data as the base image.
m = Map(AIA_171_IMAGE)

# Start by creating a plotter
plotter = SunpyPlotter()
# Plot a map
plotter.plot_map(m)
# Add an arrow to show the solar rotation axis
plotter.plot_solar_axis()

# We can also plot an arbitrary point by passing a single coordinate to
# :meth:`~sunkit_pyvista.plotter.SunpyPlotter.plot_coordinates`.
point = SkyCoord(lon=30 * u.deg,
                 lat=-10 * u.deg,
                 obstime=m.date,
                 frame=frames.HeliographicStonyhurst)
plotter.plot_coordinates(point, color='blue')

# Plot a quadrangle with width of 20 degrees and a height of 60 degrees
bottom_left = SkyCoord(30*u.deg, -10*u.deg,
                       frame=frames.HeliographicStonyhurst,
                       obstime=m.date)
plotter.plot_quadrangle(bottom_left=bottom_left, width=20*u.deg,
                        height=60*u.deg, color='blue')

# Set the camera coordinate
bottom_left = SkyCoord(30*u.deg, -10*u.deg, 7*R_sun,
                       frame=frames.HeliographicStonyhurst,
                       obstime=m.date)
plotter.set_camera_coordinate(bottom_left)
plotter.show()
