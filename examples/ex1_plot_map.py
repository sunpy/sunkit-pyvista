"""
=======================
Three dimensional plots
=======================
sunpy can interface with the `pyvista` package to produce interactive 3D plots.
"""
###############################################################################
# Start by importing the required modules
import astropy.constants as const
import astropy.units as u
from astropy.coordinates import SkyCoord

from sunpy.data.sample import AIA_193_IMAGE
from sunpy.map import Map
from sunkit_pyvista import SunpyPlotter

###############################################################################
# Import some sample data
m = Map(AIA_193_IMAGE)

###############################################################################
# 3D plots are done on "plotter" objects, which are similar to matplotlib axes.
# sunpy has a built in `PyVistaPlotter` class that can be used to plot maps
# and coordinate aware objects.

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

# Uncomment this line to show the plot when running locally
# plotter.show()

plotter.show(cpos=(-100,0,0), screenshot=os.path.join(ROOT, "docs/generated/ex1_plot_map.png"))

# %%
# .. image:: ../ex1_plot_map.png