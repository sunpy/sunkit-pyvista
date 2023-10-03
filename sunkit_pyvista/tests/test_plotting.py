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


def test_current_sheet_figure(plotter, verify_cache_image):
    gong_fname = get_gong_map()
    gong_map = smap.Map(gong_fname)
    nrho = 35
    rss = 2.5
    input_ = pfsspy.Input(gong_map, nrho, rss)
    output_ = pfsspy.pfss(input_)
    plotter.plot_current_sheet(output_)
    plotter.show(before_close_callback=verify_cache_image)


def test_field_lines_figure(aia171_test_map, plotter, verify_cache_image):
    gong_fname = get_gong_map()
    gong_map = smap.Map(gong_fname)
    nrho = 35
    rss = 2.5
    lat = np.linspace(-np.pi / 2, np.pi / 2, 8, endpoint=False)
    lon = np.linspace(0, 2 * np.pi, 8, endpoint=False)
    lat, lon = np.meshgrid(lat, lon, indexing="ij")
    lat, lon = lat.ravel() * u.rad, lon.ravel() * u.rad
    radius = 1.2
    tracer = tracing.PythonTracer()
    input_ = pfsspy.Input(gong_map, nrho, rss)
    output_ = pfsspy.pfss(input_)
    seeds = SkyCoord(lon, lat, radius * const.R_sun, frame=gong_map.coordinate_frame)
    field_lines = tracer.trace(seeds, output_)

    def color_function(field_line):
        norm = colors.LogNorm(vmin=1, vmax=1000)
        cmap = plt.get_cmap("magma")
        return cmap(norm(np.abs(field_line.expansion_factor)))

    plotter.plot_map(aia171_test_map)
    plotter.plot_field_lines(field_lines, color_func=color_function)
    plotter.show(cpos=(0, 1, 0), before_close_callback=verify_cache_image)
