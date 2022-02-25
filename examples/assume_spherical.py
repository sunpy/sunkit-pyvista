"""
========================
Ignoring off-limb pixels
========================

How to not plot off-limb pixels.

By default sunkit-pyvista plots all pixels in a map, using
:meth:`sunpy.coordinates.Helioprojective.assume_spherical_screen` to project
off-limb pixels on to a spherical screen. This example shows how off-limb
pixels can be ingored altogether.
"""

import astropy.units as u
from sunpy.data.sample import AIA_171_IMAGE
from sunpy.map import Map

from sunkit_pyvista import SunpyPlotter

m = Map(AIA_171_IMAGE)

plotter = SunpyPlotter()
plotter.plot_map(m, clip_interval=[1, 99] * u.percent,
                 assume_spherical_screen=False)
plotter.show()
