"""
=========================================================
Plotting field lines in a Cartesian box
=========================================================

This example requires the `streamtracer <https://streamtracer.readthedocs.io/en/stable/>`__ package.

This example uses an analytical magnetic field for a magnetic flux rope in the solar corona from `Hesse et al. (2005) <https://ui.adsabs.harvard.edu/abs/2005ApJ...631.1227H/abstract>`__ to demonstrate how to plot field lines in a Cartesian box.
"""
import numpy as np
from sunkit_pyvista import CartesianPlotter

###############################################################################
# We begin with the magnetic flux rope model from `Hesse et al. (2005) <https://ui.adsabs.harvard.edu/abs/2005ApJ...631.1227H/abstract>`__. The implementation is taken from the `flhtools <https://github.com/antyeates1983/flhtools>`__ package, which is accompanied by the paper by `Yeates and Page (2018) <https://ui.adsabs.harvard.edu/abs/2018JPlPh..84f7702Y/abstract>`__.


def Bx(x, y, z):
    return x * 0 - 2

def By(x, y, z, t=2):
    return -z - t * (1 - z**2) / (1 + z**2 / 25)**2 / (1 + x**2 / 25)

def Bz(x, y, z):
    return y

nx, ny, nz = 64, 64, 64
x1 = np.linspace(-20, 20, nx)
y1 = np.linspace(-20, 20, ny)
z1 = np.linspace(0, 40, nz)
x, y, z = np.meshgrid(x1, y1, z1, indexing='ij')

bx, by, bz = Bx(x, y, z), By(x, y, z), Bz(x, y, z)
b = np.stack([bx, by, bz], axis=-1)
print(b.shape)

#############################################################################
# We define seed points for tracing field lines.
x_seed = np.linspace(-20, 20, 8)
y_seed = np.linspace(-20, 20, 8)
seeds = np.array([[x, y, 0] for x in x_seed for y in y_seed])
print(seeds.shape)

###############################################################################
# We use the `~sunkit_pyvista.plotter.CartesianPlotter` to plot the field lines in a Cartesian box.
plotter = CartesianPlotter()
plotter.set_background('antiquewhite')

# Define the vector field using the magnetic field data and the grid coordinates.
plotter.define_vector_field(b, grid_coords=(x1, y1, z1))

# Show the axes of the box
plotter.show_bounds()

# Show the outline of the box
plotter.show_outline(color='black')

# Show Bz component (0 is x, 1 is y, 2 is z) on the bottom boundary of the box with a gray colormap.
plotter.show_boundary(
    'bottom',
    component=2,
    cmap='gray',
    clim=[-10, 10],
    show_scalar_bar=True,
    scalar_bar_args=dict(
        title='Bz',
        vertical=True,
    )
)

# Plot field lines from the defined seeds.
plotter.plot_field_lines(
    seeds, 
    color='cyan', 
    radius=0.1,
    seeds_config=dict(
        show_seeds=True, 
        color='red', 
        point_size=10
    )
)

# Set the camera position for a better view
plotter.camera.azimuth = 200
plotter.camera.elevation = -5
plotter.camera.zoom(0.8)

plotter.show()