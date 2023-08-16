"""
===========================================
Plotting A Volume from a cloud of SkyCoords
===========================================

This example demonstrates how to take a collection of points which you want to
visualise as either a cloud of points or a surface.
"""
import astropy.units as u
import numpy as np
import pyvista as pv
import sunpy.map
from astropy.coordinates import SkyCoord, SphericalRepresentation
from astropy.time import Time
from sunpy.coordinates import HeliocentricInertial
from sunpy.net import Fido
from sunpy.net import attrs as a
from sunpy.sun.constants import radius as Rsun

from sunkit_pyvista import SunpyPlotter
from sunkit_pyvista.sample import low_res_aia_193

################################################################################
# Let's start by loading a map to give the visualization some context.
aia = sunpy.map.Map(low_res_aia_193)

################################################################################
# Next let's generate some point data.
# For this example we are going to generate a could which is roughly spherical.
# We start by defining the center of our sphere.
sphere_center = SkyCoord(10*u.deg, 10*u.deg, 1.2 * Rsun,
                         frame="heliographic_stonyhurst", obstime=aia.date)

################################################################################
# Let's generate some spherical coordinates for theta and phi, with 50 points in
# each coordinate.  Then we generate a sphere of radius 0.05 RSun but with some
# random noise.
theta, phi = np.mgrid[0:360:50j, -90:90:50j] * u.deg
radius = 0.05 * Rsun + np.random.random(theta.shape) * 0.01 * Rsun

################################################################################
# Next we use these vectors to make a SkyCoord around the defined sphere center.
vectors = SphericalRepresentation(theta, phi, radius)
point_cloud = sphere_center.frame.realize_frame(sphere_center.cartesian + vectors.to_cartesian())

################################################################################
# Now we setup the pyvista visualization.
# For the coordinate conversions to work we need to set the obstime correctly
# for the visualization frame.
plotter = SunpyPlotter(HeliocentricInertial(obstime=aia.date))

# Draw the map for context
plotter.plot_map(aia)

################################################################################
# We now convert our cloud of coordinates to the native coordinate system of the
# visualization and then build a `pyvista.PolyData` object (which needs a 1D
# list of vectors).
# Then we plot these points
surface_points = plotter._coords_to_xyz(point_cloud)
cloud = pv.PolyData(surface_points.reshape(-1, 3))
plotter.plotter.add_points(cloud, point_size=0.7, color="cyan", style="points_gaussian")

################################################################################
# Next we want to build a surface from these points.
# We use the Delauny triangulation method as the results are better than other
# options.
surf = cloud.delaunay_3d()
plotter.plotter.add_mesh(surf)

################################################################################
# Finally set up the camera position and focus.
plotter.plotter.set_focus(plotter._coords_to_xyz(sphere_center))
plotter.set_camera_coordinate(SkyCoord(0*u.deg, 0*u.deg, 5*Rsun, frame="heliographic_stonyhurst"))

plotter.show()
