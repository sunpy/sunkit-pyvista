"""
===========================
Assuming a Spherical Screen
===========================

We can plot the Using the assumption that the image lies on the surface of a spherical screen centered at
AIA with a radius equal to the Sun-AIA distance.
"""


from sunpy.coordinates import Helioprojective
from sunpy.data.sample import AIA_171_IMAGE
from sunpy.map import Map

from sunkit_pyvista import SunpyPlotter

m = Map(AIA_171_IMAGE)

plotter = SunpyPlotter()

with Helioprojective.assume_spherical_screen(m.observer_coordinate):
    plotter.plot_map(m)
    plotter.show()
