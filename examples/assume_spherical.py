"""
========================
Ignoring off-limb pixels
========================

How to not plot off-limb pixels.

By default sunkit-pyvista plots all pixels in a map, using
:meth:`sunpy.coordinates.Helioprojective.assume_spherical_screen` to
project off-limb pixels on to a spherical screen.

This example shows how off-limb pixels can be ignored altogether.
"""

import astropy.units as u

from sunkit_pyvista import SunpyPlotter
from sunkit_pyvista.sample import low_res_aia_193

###############################################################################
# We will use an AIA 193 image from the sunpy sample data as the base image.

# Start by creating a plotter
plotter = SunpyPlotter()
# Plot a map setting the `assume_spherical_screen` to False
plotter.plot_map(
    low_res_aia_193, clip_interval=[1, 99] * u.percent,
    assume_spherical_screen=False
)

plotter.show()
