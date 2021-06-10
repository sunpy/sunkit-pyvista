"""
================================
Plotting Field Lines from pfsspy
================================

sunkit-pyvista and can also be used to plot field lines from `pfsspy`.
"""

import numpy as np
from pfsspy.sample_data import get_gong_map

import astropy.units as u
from sunpy.data.sample import AIA_193_IMAGE
from sunpy.map import Map

from sunkit_pyvista import SunpyPlotter

###############################################################################
# We will firstly use an AIA 193 image from the sunpy sample data as the base image.
m = Map(AIA_193_IMAGE)

# Start by creating a plotter
plotter = SunpyPlotter()

# Plot a map
plotter.plot_map(m)
# Add an arrow to show the solar rotation axis
plotter.plot_solar_axis()

# We load a gong_map from pfsspy
gong_fname = get_gong_map()
gong_map = Map(gong_fname)

# Define the number of grid points in rho and solar surface rarius
nrho = 35
rss = 2.5

# Create 5 points spaced between lat={-90, 90} degrees
lat = np.linspace(-np.pi / 2, np.pi / 2, 8, endpoint=False)
# Create 5 points spaced between long={0, 180} degrees
lon = np.linspace(0, 2 * np.pi, 8, endpoint=False)
# Make a 2D grid from these 1D points
lat, lon = np.meshgrid(lat, lon, indexing='ij')

# Create lon, lat and radial coordinate values by using a pfsspy map.
lat, lon = lat.ravel() * u.rad, lon.ravel() * u.rad
radius = 1.2

# We plot the field lines
plotter.plot_field_lines(gong_map, lon, lat, radius, nrho, rss)

plotter.show()
