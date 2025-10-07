import logging

import pytest
import pyvista
import pyvista as pv

import sunpy.data.test as test
import sunpy.map as smap

from sunkit_pyvista import CartesianPlotter, SunpyPlotter

pv.OFF_SCREEN = True

try:
    pyvista.start_xvfb()
except Exception as e:  # NOQA: BLE001
    msg = f"Could not start xvfb server:\n{e}"
    logging.info(msg)


@pytest.fixture()
def aia171_test_map():
    return smap.Map(test.get_test_filepath("aia_171_level1.fits"))


@pytest.fixture()
def plotter():
    return SunpyPlotter()


@pytest.fixture()
def cartesian_plotter():
    return CartesianPlotter()
