import warnings
from pathlib import Path
from functools import wraps

import matplotlib.pyplot as plt
import pytest
import pyvista

import astropy

import sunpy

__all__ = ["get_hash_library_name", "figure_test"]


@pytest.fixture
def warnings_as_errors():
    warnings.simplefilter("error")
    yield
    warnings.resetwarnings()


def get_hash_library_name():
    """
    Generate the hash library name for this env.
    """
    sunpy_version = (
        "dev" if (("dev" in sunpy.__version__) or ("rc" in sunpy.__version__)) else sunpy.__version__.replace(".", "")
    )
    pyvista_version = (
        "dev"
        if (("dev" in pyvista.__version__) or ("rc" in pyvista.__version__))
        else pyvista.__version__.replace(".", "")
    )
    astropy_version = (
        "dev"
        if (("dev" in astropy.__version__) or ("rc" in astropy.__version__))
        else astropy.__version__.replace(".", "")
    )
    return f"figure_hashes_pyvista_{pyvista_version}_astropy_{astropy_version}_sunpy_{sunpy_version}.json"


def figure_test(test_function):
    """
    A decorator for a test that verifies the hash of the current figure or the
    returned figure, with the name of the test function as the hash identifier
    in the library. A PNG is also created in the 'result_image' directory,
    which is created on the current path.

    All such decorated tests are marked with `pytest.mark.mpl_image` for convenient filtering.

    Examples
    --------
    @figure_test
    def test_simple_plot():
        plt.plot([0,1])
    """
    hash_library_name = get_hash_library_name()
    hash_library_file = Path(__file__).parent / hash_library_name

    @pytest.mark.mpl_image_compare(
        hash_library=hash_library_file, savefig_kwargs={"metadata": {"Software": None}}, style="default"
    )
    @wraps(test_function)
    def test_wrapper(*args, **kwargs):
        ret = test_function(*args, **kwargs)
        if ret is None:
            ret = plt.gcf()
        return ret

    return test_wrapper
