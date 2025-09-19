"""
This file contains tests for the main plotting routines.
"""

import numpy as np
import pytest
import pyvista as pv

import astropy.constants as const
import astropy.units as u
from astropy.coordinates import SkyCoord

import sunpy.data.test as test
import sunpy.map as smap
from sunpy.coordinates import frames

from sunkit_pyvista import CartesianPlotter, SunpyPlotter

pv.OFF_SCREEN = True


@pytest.fixture()
def aia171_test_map():
    return smap.Map(test.get_test_filepath("aia_171_level1.fits"))


@pytest.fixture()
def plotter():
    return SunpyPlotter()


@pytest.fixture()
def cartesian_plotter():
    return CartesianPlotter()


def test_plot_map_with_functionality(
    aia171_test_map,
    plotter,
    verify_cache_image,
    tmp_path,
):
    plotter.plot_map(aia171_test_map, clip_interval=(0, 99) * u.percent)
    plotter.plot_solar_axis()

    bottom_left = SkyCoord(
        30 * u.deg,
        -10 * u.deg,
        frame=frames.HeliographicStonyhurst,
        obstime=aia171_test_map.date,
    )
    plotter.plot_quadrangle(
        bottom_left=bottom_left,
        width=20 * u.deg,
        height=60 * u.deg,
        color="blue",
    )
    plotter.plot_limb(aia171_test_map)

    line = SkyCoord(
        lon=[90, 2200, 2200] * u.deg,
        lat=[10, 10, 20] * u.deg,
        distance=[1, 2, 3] * const.R_sun,
        frame="heliocentricinertial",
    )
    plotter.plot_coordinates(line)

    coordinate = SkyCoord(
        30 * u.deg,
        -10 * u.deg,
        frame=frames.HeliographicStonyhurst,
        obstime=aia171_test_map.date,
    )
    plotter.plot_coordinates(coordinate, color="blue")

    filepath = tmp_path / "save_data.vtm"
    plotter.save(filepath=filepath)

    plotter = SunpyPlotter()
    plotter.load(filepath)

    plotter.show(cpos=(0, 1, 0), before_close_callback=verify_cache_image)


def test_cartesian_plotter(cartesian_plotter, verify_cache_image):
    def Bx(x, y, z):
        return x * 0 - 2

    def By(x, y, z, t=2):
        return -z - t * (1 - z**2) / (1 + z**2 / 25) ** 2 / (1 + x**2 / 25)

    def Bz(x, y, z):
        return y

    nx, ny, nz = 64, 64, 64
    x1 = np.linspace(-20, 20, nx)
    y1 = np.linspace(-20, 20, ny)
    z1 = np.linspace(0, 40, nz)
    x, y, z = np.meshgrid(x1, y1, z1, indexing="ij")
    bx, by, bz = Bx(x, y, z), By(x, y, z), Bz(x, y, z)
    b = np.stack([bx, by, bz], axis=-1)  # (64, 64, 64, 3)

    x_seed = np.linspace(-20, 20, 8)
    y_seed = np.linspace(-20, 20, 8)
    seeds = np.array([[x, y, 0] for x in x_seed for y in y_seed])  # (64, 3)

    plotter = cartesian_plotter
    plotter.set_background("antiquewhite")
    plotter.define_vector_field(b, grid_coords=(x1, y1, z1))
    plotter.show_bounds()
    plotter.show_outline(color="black")
    plotter.show_boundary(
        "bottom",
        component=2,
        cmap="gray",
        clim=[-10, 10],
        show_scalar_bar=True,
        scalar_bar_args=dict(
            title="Bz",
            vertical=True,
        ),
    )
    plotter.plot_field_lines(
        seeds, color="cyan", radius=0.1, seeds_config=dict(show_seeds=True, color="red", point_size=10)
    )
    plotter.camera.azimuth = 200
    plotter.camera.elevation = -5
    plotter.camera.zoom(0.8)
    plotter.show(before_close_callback=verify_cache_image)
