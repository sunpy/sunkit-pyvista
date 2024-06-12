import pytest
import sunpy.map as smap
from sunkit_pyvista import SunpyPlotter
import sunpy.data.test as test

@pytest.fixture()
def aia171_test_map():
    return smap.Map(test.get_test_filepath("aia_171_level1.fits"))


@pytest.fixture()
def plotter():
    return SunpyPlotter()
