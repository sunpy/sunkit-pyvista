import numpy as np
import pfsspy
import pytest
import pyvista as pv
from pfsspy import tracing
from pfsspy.sample_data import get_gong_map

import astropy.constants as const
import astropy.units as u
import sunpy.data.test as test
import sunpy.map as smap
from astropy.coordinates import SkyCoord
from sunpy.coordinates import HeliocentricInertial, HeliographicStonyhurst

from sunkit_pyvista import SunpyPlotter


@pytest.fixture
def aia171_test_map():
    return smap.Map(test.get_test_filepath('aia_171_level1.fits'))


@pytest.fixture
def plotter():
    return SunpyPlotter()


@pytest.mark.display_server
def test_basic(plotter):
    assert isinstance(plotter.plotter, pv.Plotter)
    plotter.show()


def test_coordinate_frame(plotter):
    assert isinstance(plotter.coordinate_frame, HeliocentricInertial)


def test_coord_to_xyz(aia171_test_map, plotter):
    coordinate = aia171_test_map.observer_coordinate
    plot_coord = plotter._coords_to_xyz(coordinate)
    assert plot_coord.shape[1] == 3


def test_camera_position(aia171_test_map, plotter):
    coord = aia171_test_map.observer_coordinate
    camera_position = plotter._coords_to_xyz(coord)
    plotter.set_camera_coordinate(coord)
    assert (plotter.camera.position == camera_position[0]).all()


def test_set_view_angle(plotter):
    plotter.set_view_angle(45*u.deg)
    assert plotter.camera.view_angle == 45
    with pytest.raises(ValueError, match=r"specified view angle must be "
                       r"0 deg < angle <= 180 deg"):
        plotter.set_view_angle(190*u.deg)


def test_plot_map(aia171_test_map, plotter):
    plotter.plot_map(aia171_test_map)
    assert plotter.plotter.mesh.n_cells == 128**2
    assert plotter.plotter.mesh.n_points == 129**2


def test_plot_solar_axis(plotter):
    plotter.plot_solar_axis()
    assert plotter.plotter.mesh.n_cells == 43
    assert plotter.plotter.mesh.n_points == 101


def test_plot_quadrangle(aia171_test_map, plotter):
    bottom_left = SkyCoord(30*u.deg, -10*u.deg,
                           frame=HeliographicStonyhurst,
                           obstime=aia171_test_map.date)
    plotter.plot_quadrangle(bottom_left=bottom_left, width=20*u.deg,
                            height=60*u.deg, color='blue')
    assert plotter.plotter.mesh.n_cells == 4000
    assert plotter.plotter.mesh.n_points == 4001


def test_plot_coordinates(plotter):
    # Tests the plot for a line
    line = SkyCoord(lon=[180, 190, 200] * u.deg,
                    lat=[0, 10, 20] * u.deg,
                    distance=[1, 2, 3] * const.R_sun,
                    frame='heliocentricinertial')
    plotter.plot_coordinates(line)
    assert plotter.plotter.mesh.n_cells == 1
    assert plotter.plotter.mesh.n_points == 3

    # Tests plotting of a small sphere
    sphere = SkyCoord(lon=225*u.deg,
                      lat=45*u.deg,
                      distance=1*const.R_sun,
                      frame='heliocentricinertial')
    plotter.plot_coordinates(sphere)
    assert plotter.plotter.mesh.n_cells == 1680
    assert plotter.plotter.mesh.n_points == 842
    expected_center = [-0.5000000149011612, -0.5, 0.7071067690849304]
    assert np.allclose(plotter.plotter.mesh.center, expected_center)


def test_clip_interval(aia171_test_map, plotter):
    plotter.plot_map(aia171_test_map, clip_interval=(1, 99)*u.percent)
    clim = plotter._get_clim(data=plotter.plotter.mesh['data'],
                             clip_interval=(1, 99)*u.percent)
    expected_clim = [0.006716044038535769, 0.8024368512284383]
    assert np.allclose(clim, expected_clim)

    expected_clim = [0, 1]
    clim = plotter._get_clim(data=plotter.plotter.mesh['data'],
                             clip_interval=(0, 100)*u.percent)
    assert np.allclose(clim, expected_clim)

    with pytest.raises(ValueError, match=r"Clip percentile interval must be "
                       r"specified as two numbers."):
        plotter.plot_map(aia171_test_map, clip_interval=(1, 50, 99)*u.percent)


def test_multi_block(plotter):
    plotter.plot_solar_axis()
    assert plotter.all_meshes['solar_axis'][0].n_cells == 43
    assert plotter.all_meshes['solar_axis'][0].n_points == 101


def test_field_lines(aia171_test_map, plotter):
    gong_fname = get_gong_map()
    gong_map = smap.Map(gong_fname)
    nrho = 35
    rss = 2.5
    lat = np.linspace(-np.pi / 2, np.pi / 2, 8, endpoint=False)
    lon = np.linspace(0, 2 * np.pi, 8, endpoint=False)
    lat, lon = np.meshgrid(lat, lon, indexing='ij')
    lat, lon = lat.ravel() * u.rad, lon.ravel() * u.rad
    radius = 1.2
    tracer = tracing.PythonTracer()
    input_ = pfsspy.Input(gong_map, nrho, rss)
    output_ = pfsspy.pfss(input_)
    seeds = SkyCoord(lon, lat, radius*const.R_sun,
                     frame=gong_map.coordinate_frame)
    field_lines = tracer.trace(seeds, output_)
    plotter.plot_field_lines(field_lines)
    assert len(plotter.all_meshes['field_lines'][0]) == 64
    assert isinstance(plotter.all_meshes['field_lines'][0],
                      pv.core.composite.MultiBlock)


def test_save_and_load(aia171_test_map, plotter, tmp_path):
    plotter.plot_solar_axis()

    filepath = (tmp_path / "save_data.vtm")
    plotter.save(filepath=filepath)

    plotter.plotter.clear()
    plotter.load(filepath)

    assert plotter.plotter.mesh.n_cells == 43
    assert plotter.plotter.mesh.n_points == 101

    directory_path = filepath.with_suffix('')
    with pytest.raises(ValueError, match=f"Directory '{str(directory_path.absolute())}' already exists"):
        plotter.save(filepath=filepath)


def test_loop_through_meshes(plotter):
    sphere = pv.Cube()
    sphere2 = pv.Cube(center=(0, 1, 1))
    inner_block = pv.MultiBlock([sphere])
    outer_block = pv.MultiBlock([inner_block, sphere2])

    plotter._loop_through_meshes(outer_block)

    assert plotter.plotter.mesh.center == [0, 1, 1]
