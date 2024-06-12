import pathlib

import astropy.constants as const
import astropy.units as u
import numpy as np
import pytest
import pyvista as pv
from astropy.coordinates import SkyCoord
from sunpy.coordinates import HeliocentricInertial, HeliographicStonyhurst


@pytest.mark.display_server()
def test_basic(plotter):
    assert isinstance(plotter.plotter, pv.Plotter)
    plotter.show()


def test_coordinate_frame(plotter):
    assert isinstance(plotter.coordinate_frame, HeliocentricInertial)


def test_coord_to_xyz(aia171_test_map, plotter):
    coordinate = aia171_test_map.observer_coordinate
    plot_coord = plotter._coords_to_xyz(coordinate)
    assert plot_coord.shape == (3,)


def test_camera_position(aia171_test_map, plotter):
    coord = aia171_test_map.observer_coordinate
    camera_position = plotter._coords_to_xyz(coord)
    plotter.set_camera_coordinate(coord)
    assert (plotter.camera.position == camera_position).all()


def test_set_view_angle(plotter):
    plotter.set_view_angle(45 * u.deg)
    assert plotter.camera.view_angle == 45
    with pytest.raises(
        ValueError,
        match=r"specified view angle must be " r"0 deg < angle <= 180 deg",
    ):
        plotter.set_view_angle(190 * u.deg)


def test_plot_map(aia171_test_map, plotter):
    plotter.plot_map(aia171_test_map)
    assert plotter.plotter.mesh.n_cells == 16384
    assert plotter.plotter.mesh.n_points == 16641


def test_plot_solar_axis(plotter):
    plotter.plot_solar_axis()
    assert plotter.plotter.mesh.n_cells == 43
    assert plotter.plotter.mesh.n_points == 101


def test_plot_quadrangle(aia171_test_map, plotter):
    bottom_left = SkyCoord(
        30 * u.deg,
        -10 * u.deg,
        frame=HeliographicStonyhurst,
        obstime=aia171_test_map.date,
    )
    plotter.plot_quadrangle(
        bottom_left=bottom_left,
        width=20 * u.deg,
        height=60 * u.deg,
        color="blue",
    )
    assert plotter.plotter.mesh.n_cells == 22
    assert plotter.plotter.mesh.n_points == 80060


def test_plot_coordinates(aia171_test_map, plotter):
    # Tests the plot for a line
    line = SkyCoord(
        lon=[180, 190, 200] * u.deg,
        lat=[0, 10, 20] * u.deg,
        distance=[1, 2, 3] * const.R_sun,
        frame="heliocentricinertial",
    )
    plotter.plot_coordinates(line)
    assert plotter.plotter.mesh.n_cells == 1
    assert plotter.plotter.mesh.n_points == 3

    # Tests plotting of a small sphere
    sphere = SkyCoord(
        lon=225 * u.deg,
        lat=45 * u.deg,
        distance=1 * const.R_sun,
        frame="heliocentricinertial",
    )
    plotter.plot_coordinates(sphere)
    assert plotter.plotter.mesh.n_cells == 1680
    assert plotter.plotter.mesh.n_points == 842
    expected_center = [-0.5000000149011612, -0.5, 0.7071067690849304]
    assert np.allclose(plotter.plotter.mesh.center, expected_center)

    pixel_pos = np.argwhere(aia171_test_map.data == aia171_test_map.data.max()) * u.pixel
    hpc_max = aia171_test_map.pixel_to_world(pixel_pos[:, 1], pixel_pos[:, 0])
    plotter.plot_coordinates(hpc_max, color="blue")
    assert plotter.plotter.mesh.n_cells == 1680
    assert plotter.plotter.mesh.n_points == 842


def test_clip_interval(aia171_test_map, plotter):
    plotter.plot_map(aia171_test_map, clip_interval=(1, 99) * u.percent)
    clim = plotter._get_clim(
        data=plotter.plotter.mesh["data"],
        clip_interval=(1, 99) * u.percent,
    )
    expected_clim = [0.006716044038535769, 0.8024368512284383]
    assert np.allclose(clim, expected_clim)

    expected_clim = [0, 1]
    clim = plotter._get_clim(
        data=plotter.plotter.mesh["data"],
        clip_interval=(0, 100) * u.percent,
    )
    assert np.allclose(clim, expected_clim)

    with pytest.raises(
        ValueError,
        match=r"Clip percentile interval must be " r"specified as two numbers.",
    ):
        plotter.plot_map(aia171_test_map, clip_interval=(1, 50, 99) * u.percent)


def test_multi_block(plotter):
    plotter.plot_solar_axis()
    assert plotter.all_meshes["solar_axis"][0].n_cells == 43
    assert plotter.all_meshes["solar_axis"][0].n_points == 101


def test_save_and_load(aia171_test_map, plotter, tmp_path):
    plotter.plot_map(aia171_test_map)

    filepath = tmp_path / "save_data.vtm"
    plotter.save(filepath=filepath)

    plotter.plotter.clear()
    plotter.load(filepath)

    assert plotter.plotter.mesh.n_cells == 16384
    assert plotter.plotter.mesh.n_points == 16641
    assert dict(plotter.plotter.mesh.field_data)["cmap"][0] == "sdoaia171"

    with pytest.raises(ValueError, match="VTM file"):
        plotter.save(filepath=filepath)
    pathlib.Path(tmp_path / "save_data_dir").mkdir(parents=True, exist_ok=True)
    filepath = tmp_path / "save_data_dir.vtm"
    with pytest.raises(ValueError, match="already exists"):
        plotter.save(filepath=filepath)


def test_loop_through_meshes(plotter):
    sphere = pv.Cube()
    sphere2 = pv.Cube(center=(0, 1, 1))
    inner_block = pv.MultiBlock([sphere])
    outer_block = pv.MultiBlock([inner_block, sphere2])
    plotter._loop_through_meshes(outer_block)

    assert plotter.plotter.mesh.center == [0, 1, 1]


def test_plot_limb(aia171_test_map, plotter):
    plotter.plot_limb(aia171_test_map)
    assert plotter.plotter.mesh.n_cells == 22
    assert plotter.plotter.mesh.n_points == 20040
