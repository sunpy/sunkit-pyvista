"""
This file contains figure comparison tests.
"""

import astropy.constants as const
import astropy.units as u
import pytest
import pyvista
import sunpy.data.test as test
import sunpy.map as smap
from astropy.coordinates import SkyCoord
from sunpy.coordinates import frames

from sunkit_pyvista import SunpyPlotter

pyvista.OFF_SCREEN = True


@pytest.fixture()
def aia171_test_map():
    return smap.Map(test.get_test_filepath("aia_171_level1.fits"))


@pytest.fixture()
def plotter():
    return SunpyPlotter()


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
