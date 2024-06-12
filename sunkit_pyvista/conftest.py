import inspect
import logging
import platform
import warnings
from pathlib import Path

import pytest
import pyvista

try:
    pyvista.start_xvfb()
except Exception as e:  # NOQA: BLE001
    msg = f"Could not start xvfb server:\n{e}"
    logging.info(msg)

IMAGE_CACHE_DIR = Path(__file__).parent.absolute() / "tests" / "image_cache"
if not IMAGE_CACHE_DIR.is_dir():
    IMAGE_CACHE_DIR.mkdir()
# Normal image warning/error thresholds (assumes using use_vtk)
IMAGE_REGRESSION_ERROR = 400  # major differences
IMAGE_REGRESSION_WARNING = 200  # minor differences


@pytest.fixture(scope="session", autouse=True)
def _get_cmd_opt(pytestconfig):
    global update_image_cache, ignore_image_cache
    update_image_cache = pytestconfig.getoption("update_image_cache")
    ignore_image_cache = pytestconfig.getoption("ignore_image_cache")


def pytest_addoption(parser):
    parser.addoption("--update_image_cache", action="store_true", default=False)
    parser.addoption("--ignore_image_cache", action="store_true", default=False)


def verify_cache_images(plotter):
    """
    Either store or validate an image.

    This is function should only be called within a pytest environment.
    Pass it to either the ``Plotter.show()`` or the ``pyvista.plot()``
    functions as the before_close_callback keyword arg.

    Assign this only once for each test you'd like to validate the
    previous image of. This will not work with parameterized tests.
    """
    import vtk

    if not vtk.__version__ >= "9.2.0" or platform.system() != "Linux":
        pytest.skip("VTK 9.2 or above and linux required for figure tests.")

    # Since each test must contain a unique name,
    # we can use the function test to name the image.
    stack = inspect.stack()
    test_name = None
    for item in stack:
        if item.function == "check_gc":
            return None
        if item.function[:5] == "test_":
            test_name = item.function
            break

    allowed_error = IMAGE_REGRESSION_ERROR
    allowed_warning = IMAGE_REGRESSION_WARNING

    if test_name is None:
        msg = (
            "Unable to identify calling test function. This function "
            "should only be used within a pytest environment."
        )
        raise RuntimeError(
            msg,
        )

    image_filename = IMAGE_CACHE_DIR / (test_name[5:] + ".png")

    if update_image_cache or not image_filename.is_file():
        logging.info("Image doesn't exist, saving file in image_cache")
        return plotter.screenshot(str(image_filename))

    if ignore_image_cache:
        return None

    error = pyvista.compare_images(str(image_filename), plotter)
    if error > allowed_error:
        msg = "Exceeded image regression error of " f"{IMAGE_REGRESSION_ERROR} with an image error of " f"{error}"
        raise RuntimeError(
            msg,
        )
    if error > allowed_warning:
        warnings.warn(
            "Exceeded image regression warning of " f"{IMAGE_REGRESSION_WARNING} with an image error of " f"{error}",
            stacklevel=2,
        )
    return None


@pytest.fixture()
def verify_cache_image():
    return verify_cache_images
