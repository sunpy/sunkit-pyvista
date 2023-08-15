"""
==================================
Extending functionality from sunpy
==================================

sunkit-pyvista also extends :meth:`~sunpy.map.GenericMap.draw_quadrangle` from
`sunpy` and :meth:`~astropy.visualization.wcsaxes.WCSAxes.plot_coord` from `astropy`
to produce them in 3D.
"""
import astropy.units as u
import numpy as np
from astropy.constants import R_sun
from astropy.coordinates import SkyCoord
from sunpy.coordinates import frames

from sunkit_pyvista import SunpyPlotter
from sunkit_pyvista.sample import LOW_RES_AIA_193

###############################################################################
# We will use an AIA 193 image from the sunpy sample data as the base image.

# sphinx_gallery_defer_figures

# Start by creating a plotter
plotter = SunpyPlotter()
# Plot a map
plotter.plot_map(LOW_RES_AIA_193, clip_interval=(1, 99.9) * u.percent)
# Add an arrow to show the solar rotation axis
plotter.plot_solar_axis()

###############################################################################
# We can also plot an arbitrary point by passing a single coordinate to
# :meth:`~sunkit_pyvista.plotter.SunpyPlotter.plot_coordinates`.
# Plotting a point on the brightest pixel in the map

pixel_pos = np.argwhere(LOW_RES_AIA_193.data == LOW_RES_AIA_193.data.max()) * u.pixel
hpc_max = LOW_RES_AIA_193.pixel_to_world(pixel_pos[:, 1], pixel_pos[:, 0])
plotter.plot_coordinates(hpc_max, color="blue")

# Plot a quadrangle with width of 20 degrees and a height of 60 degrees
bottom_left = SkyCoord(
    30 * u.deg,
    -10 * u.deg,
    frame=frames.HeliographicStonyhurst,
    obstime=LOW_RES_AIA_193.date,
)
plotter.plot_quadrangle(
    bottom_left=bottom_left,
    width=20 * u.deg,
    height=60 * u.deg,
    color="blue",
)

# Set the camera coordinate to view the plot correctly
camera_coord = SkyCoord(
    30 * u.deg,
    -10 * u.deg,
    6 * R_sun,
    frame=frames.HeliographicStonyhurst,
    obstime=LOW_RES_AIA_193.date,
)
plotter.set_camera_coordinate(camera_coord)

plotter.show()
