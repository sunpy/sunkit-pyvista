"""
=========================================================
Plotting field lines and current sheets from sunkit-magex
=========================================================

``sunkit-pyvista`` can be used to plot field lines from ``sunkit-magex``.

This example requires the `streamtracer <https://streamtracer.readthedocs.io/en/stable/>`__ package.

This example uses specific sample data from ``sunkit-magex`` to demonstrate the plotting of field lines and current sheets.
If you want to adapt this, you have to make sure that the input magnetic field data has the same date as the AIA image.
"""

import matplotlib.pyplot as plt
import numpy as np
from matplotlib import colors

import astropy.units as u
from astropy.constants import R_sun
from astropy.coordinates import SkyCoord

import sunpy.map
from sunkit_magex import pfss
from sunkit_magex.pfss import tracing
from sunkit_magex.pfss.sample_data import get_gong_map
from sunpy.coordinates import frames

from sunkit_pyvista import SunpyPlotter
from sunkit_pyvista.sample import LOW_RES_AIA_193

###############################################################################
# We will start by doing the magnetic field extrapolation using ``sunkit-magex``.

# We load a gong_map from sunkit-magex
gong_fname = get_gong_map()
gong_map = sunpy.map.Map(gong_fname)

# Create points spaced between lat={-89, 89} degrees
lat = np.linspace(-np.pi / 2, np.pi / 2, 32, endpoint=False)
# Create 32 points spaced between long={0, 360} degrees
lon = np.linspace(1, 2 * np.pi, 32, endpoint=False)
# Make a 2D grid from these 1D points
lat, lon = np.meshgrid(lat, lon, indexing="ij")
# Create lon, lat and radial coordinate values by using a sunkit-magex
# and trace them using tracer
lat, lon = lat.ravel() * u.rad, lon.ravel() * u.rad
# Define the number of grid points in rho and solar surface radius
nrho = 30
rss = 1.5
input_ = pfss.Input(gong_map, nrho, rss)
pfss_output = pfss.pfss(input_)
seeds = SkyCoord(lon, lat, 1.2 * R_sun, frame=gong_map.coordinate_frame)
tracer = tracing.PythonTracer()
field_lines = tracer.trace(seeds, pfss_output)

###############################################################################
# We will be using an AIA 193 image from the sunpy sample data as the base image.

# sphinx_gallery_defer_figures

# Start by creating a plotter
plotter = SunpyPlotter()

# Plot the AIA  map
plotter.plot_map(LOW_RES_AIA_193, clip_interval=[1, 99] * u.percent, assume_spherical_screen=False)
# Add an arrow to show the solar rotation axis
plotter.plot_solar_axis()
# Now we plot the Gong Map to fill in the farside.
plotter.plot_map(gong_map, cmap="hmimag", clip_interval=[1, 99] * u.percent)

###############################################################################
# We can also specify a color function while plotting the field lines.
# This function takes a single field line, and returns a color either
# in the form of a string, (r,g,b) or (r,g,b,a) tuple.
# In this case we use a Matplotlib norm and colormap to return a tuple of RGBA values.


def my_fline_color_func(field_line):
    norm = colors.LogNorm(vmin=1, vmax=1000)
    cmap = plt.get_cmap("viridis")
    return cmap(norm(np.abs(field_line.expansion_factor)))


# Plotting the field lines
plotter.plot_field_lines(field_lines, color_func=my_fline_color_func)
# Plotting the current sheet
plotter.plot_current_sheet(pfss_output, opacity=0.5)

# Set the camera coordinate to view the plot correctly
camera_coord = SkyCoord(
    0 * u.deg,
    0 * u.deg,
    6 * R_sun,
    frame=frames.HeliographicStonyhurst,
    obstime=LOW_RES_AIA_193.date,
)
plotter.set_camera_coordinate(camera_coord)

plotter.show()
