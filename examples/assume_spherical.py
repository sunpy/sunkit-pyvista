"""
===========================
Assuming a Spherical Screen
===========================

We can plot the Using the assumption that the image lies on the surface of a spherical screen centered at
AIA with a radius equal to the Sun-AIA distance.
"""
from sunpy.coordinates import Helioprojective

from sunkit_pyvista import SunpyPlotter
from sunkit_pyvista.sample import low_res_aia_171

m = low_res_aia_171()

plotter = SunpyPlotter()

with Helioprojective.assume_spherical_screen(m.observer_coordinate):
    plotter.plot_map(m)
    plotter.show()
