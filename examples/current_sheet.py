"""
================================
Plotting Current Sheets from pfsspy
================================

``sunkit-pyvista`` can be used to extract and plot the current sheet from ``pfsspy``.

This example requires the pre-compiled FORTRAN
code from `streamtracer <https://streamtracer.readthedocs.io/en/stable/>`__.
Wheels are provided for all major platforms, so this example should run without issue.
"""
import astropy.units as u
import matplotlib.pyplot as plt
import numpy as np
import pfsspy
import sunpy.map
from astropy.constants import R_sun
from astropy.coordinates import SkyCoord

from matplotlib import colors
from pfsspy import tracing
from pfsspy.sample_data import get_gong_map
from sunpy.coordinates import frames

from sunkit_pyvista import SunpyPlotter

###############################################################################
# We now need to do the magnetic field extrapolation using ``pfsspy``.

# sphinx_gallery_defer_figures

# We load a gong_map from pfsspy
gong_fname = get_gong_map()
gong_map = sunpy.map.Map(gong_fname)
gong_map.plot()

plotter = SunpyPlotter()
# Now we plot the Gong Map to fill in the farside.
plotter.plot_map(gong_map, cmap="hmimag", clip_interval=[1, 99] * u.percent)

# Create points spaced between lat={-90, 90} degrees
lat = np.linspace(-np.pi / 2, np.pi / 2, 16, endpoint=False)
# Create 32 points spaced between long={0, 360} degrees
lon = np.linspace(0, np.pi*2, 16, endpoint=False)
# Make a 2D grid from these 1D points
lat, lon = np.meshgrid(lat, lon, indexing="ij")
# Create lon, lat and radial coordinate values by using a pfsspy
# and trace them using tracer
lat, lon = lat.ravel() * u.rad, lon.ravel() * u.rad
# Define the number of grid points in rho and solar surface radius
nrho = 30
rss = 2.5
input_ = pfsspy.Input(gong_map, nrho, rss)
output_ = pfsspy.pfss(input_)
seeds = SkyCoord(lon, lat, 2.4 * R_sun, frame=gong_map.coordinate_frame)
tracer = tracing.PythonTracer()
field_lines = tracer.trace(seeds, output_)

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
plotter.plot_current_sheet(output_)

plotter.show()
