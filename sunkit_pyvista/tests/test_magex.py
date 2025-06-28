import matplotlib.pyplot as plt
import numpy as np
import pyvista as pv
from matplotlib import colors

import astropy.constants as const
import astropy.units as u
from astropy.coordinates import SkyCoord

import sunpy.map as smap
from sunkit_magex import pfss
from sunkit_magex.pfss import tracing
from sunkit_magex.pfss.sample_data import get_gong_map

from sunkit_pyvista import SunpyPlotter


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
    input_ = pfss.Input(gong_map, nrho, rss)
    output_ = pfss.pfss(input_)
    seeds = SkyCoord(lon, lat, radius * const.R_sun, frame=gong_map.coordinate_frame)
    field_lines = tracer.trace(seeds, output_)

    def color_function(field_line):
        norm = colors.LogNorm(vmin=1, vmax=1000)
        cmap = plt.get_cmap("magma")
        return cmap(norm(np.abs(field_line.expansion_factor)))

    plotter.plot_map(aia171_test_map)
    plotter.plot_field_lines(field_lines, color_func=color_function)
    plotter.show(cpos=(0, 1, 0), before_close_callback=verify_cache_image)


def test_field_lines_and_color_func(plotter):
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
    input_ = pfss.Input(gong_map, nrho, rss)
    output_ = pfss.pfss(input_)
    seeds = SkyCoord(lon, lat, radius * const.R_sun, frame=gong_map.coordinate_frame)
    field_lines = tracer.trace(seeds, output_)
    plotter.plot_field_lines(field_lines)
    assert isinstance(plotter.all_meshes["field_lines"][0], pv.PolyData)

    def color_func(field_line):
        norm = colors.LogNorm(vmin=1, vmax=1000)
        cmap = plt.get_cmap("viridis")
        return cmap(norm(np.abs(field_line.expansion_factor)))

    plotter = SunpyPlotter()
    plotter.plot_field_lines(field_lines, color_func=color_func)
