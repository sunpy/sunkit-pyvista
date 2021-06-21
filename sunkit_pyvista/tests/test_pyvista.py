import pytest
import pyvista as pv

import astropy.units as u
import sunpy.data.test as test
import sunpy.map as smap
from sunpy.coordinates import HeliocentricInertial

from sunkit_pyvista import SunpyPlotter


@pytest.fixture
def aia171_test_map():
    return smap.Map(test.get_test_filepath('aia_171_level1.fits'))


@pytest.fixture
def plotter():
    return SunpyPlotter()


def test_basic(plotter):
    assert isinstance(plotter.plotter, pv.Plotter)


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
