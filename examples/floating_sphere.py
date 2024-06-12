"""
===========================================
Plotting a volume from a cloud of SkyCoords
===========================================

This example demonstrates how to take a collection of points which you want to
visualize as either a cloud of points or a surface.
"""

import astropy.units as u
import numpy as np
import sunpy.map
from astropy.coordinates import SkyCoord, SphericalRepresentation
from sunpy.sun.constants import radius as Rsun

from sunkit_pyvista import SunpyPlotter
from sunkit_pyvista.sample import LOW_RES_AIA_193

################################################################################
# Let's start by loading a map to give the visualization some context.

aia = sunpy.map.Map(LOW_RES_AIA_193)

################################################################################
# Next let's generate some point data.
#
# For this example we are going to generate a could which is roughly spherical.
# We start by defining the center of our sphere.

sphere_center = SkyCoord(10 * u.deg, 10 * u.deg, 1.2 * Rsun, frame="heliographic_stonyhurst")

################################################################################
# Let's generate some spherical coordinates for theta and phi, with 50 points in
# each coordinate. Then we generate a sphere of radius 0.05 RSun but with some
# random noise.

theta, phi = np.mgrid[0:360:50j, -90:90:50j] * u.deg
rng = np.random.default_rng()
radius = 0.05 * Rsun + rng.random(theta.shape) * 0.01 * Rsun

################################################################################
# Next we use these vectors to make a SkyCoord around the defined sphere center.

vectors = SphericalRepresentation(theta, phi, radius)
new_points = sphere_center.cartesian + vectors.to_cartesian()
point_cloud = sphere_center.frame.realize_frame(new_points)

################################################################################
# Now we setup the pyvista visualization.
# For the coordinate conversions to work we need to set the ``obstime`` correctly
# for the visualization frame.

# sphinx_gallery_defer_figures
plotter = SunpyPlotter(obstime=aia.date)

# Draw the map for context
plotter.plot_map(aia)

################################################################################
# We now convert our cloud of coordinates to the native coordinate system of the
# visualization and then build a `pyvista.PolyData` object (which needs a 1D
# list of vectors).
# Then we plot these points

# sphinx_gallery_defer_figures

cloud = plotter.coordinates_to_polydata(point_cloud)
_ = plotter.plotter.add_points(cloud, point_size=0.7, color="cyan", style="points_gaussian")

################################################################################
# Next we want to build a surface from these points.
# We use the Delauny triangulation method as the results are better than other
# options.

# sphinx_gallery_defer_figures

surf = cloud.delaunay_3d()
_ = plotter.plotter.add_mesh(surf)

################################################################################
# Finally set up the camera position and focus.

plotter.set_camera_focus(sphere_center)
cam_coord = SkyCoord(0 * u.deg, 0 * u.deg, 5 * Rsun, frame="heliographic_stonyhurst")
plotter.set_camera_coordinate(cam_coord)

plotter.show()
