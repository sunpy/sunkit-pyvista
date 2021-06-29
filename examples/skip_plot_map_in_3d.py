"""
===========================
Finding the brightest pixel
===========================

sunkit-pyvista can be used to overplot the location of pixels
on the map.
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

# Start by creating a plotter
plotter = SunpyPlotter()

# Plot a map
plotter.plot_map(m, clip_interval=(1, 99)*u.percent)
# Add an arrow to show the solar rotation axis
plotter.plot_solar_axis()

# To find the brightest pixel, we find the maximum in the AIA image data then
# transform that pixel coordinate to a map coordinate.
pixel_pos = np.argwhere(m.data == m.data.max()) * u.pixel
hpc_max = m.pixel_to_world(pixel_pos[:, 1], pixel_pos[:, 0])

# To plot a overplot particular coordinate, we simply
# use the :meth:`~sunkit_pyvista.plotter.SunpyPlotter.plot_coordinates`
# with a single coordinate as an argument.
plotter.plot_coordinates(hpc_max, radius=0.04, color='blue')

# Setting the camera coordinate to view the point.
camera_position = SkyCoord(30*u.deg, -10*u.deg, 6*R_sun,
                           frame=frames.HeliographicStonyhurst,
                           obstime=m.date)
plotter.set_camera_coordinate(camera_position)
plotter.show()
