"""
=====================
Adding extra features
=====================

sunkit-pyvista also extends ``draw_quadrangle``, and ``plot_coord`` from `sunpy`, producing them in 3D.
"""

import numpy as np

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

###############################################################################
# Start by creating a plotter
plotter = SunpyPlotter()
# Plot a map
plotter.plot_map(m, clip_interval=(1, 99.9)*u.percent)
# Add an arrow to show the solar rotation axis
plotter.plot_solar_axis()

###############################################################################
# We can also plot an arbitrary point by passing a single coordinate to
# :meth:`~sunkit_pyvista.plotter.SunpyPlotter.plot_coordinates`.
# Plotting a point on the brightest pixel in the map
pixel_pos = np.argwhere(m.data == m.data.max()) * u.pixel
hpc_max = m.pixel_to_world(pixel_pos[:, 1], pixel_pos[:, 0])
plotter.plot_coordinates(hpc_max, color='blue')

# Plot a quadrangle with width of 20 degrees and a height of 60 degrees
bottom_left = SkyCoord(30*u.deg, -10*u.deg,
                       frame=frames.HeliographicStonyhurst,
                       obstime=m.date)
plotter.plot_quadrangle(bottom_left=bottom_left, width=20*u.deg,
                        height=60*u.deg, color='blue')

# Set the camera coordinate to view the plot correctly
camera_coord = SkyCoord(30*u.deg, -10*u.deg, 6*R_sun,
                        frame=frames.HeliographicStonyhurst,
                        obstime=m.date)
plotter.set_camera_coordinate(camera_coord)
plotter.show()
