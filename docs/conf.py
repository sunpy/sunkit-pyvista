"""
Configuration file for the Sphinx documentation builder.
"""
import warnings
import numpy as np
import pyvista

# Use the sunpy theme
from sunpy_sphinx_theme.conf import *
from packaging.version import Version
from sunkit_pyvista import __version__
import os
from datetime import datetime

# -- Project information -----------------------------------------------------
on_rtd = os.environ.get("READTHEDOCS", None) == "True"
os.environ["HIDE_PARFIVE_PROGESS"] = "True"
project = "sunkit-pyvista"
author = "SunPy Community"
copyright = "{}, {}".format(datetime.now().year, author)

# The full version, including alpha/beta/rc tags
release = __version__
sunkit_pyvista_version = Version(__version__)
is_release = not (
    sunkit_pyvista_version.is_prerelease or sunkit_pyvista_version.is_devrelease
)

# -- General configuration ---------------------------------------------------
extensions = [
    "sphinx_automodapi.automodapi",
    "sphinx_automodapi.smart_resolver",
    "sphinx_changelog",
    "sphinx.ext.autodoc",
    "sphinx.ext.coverage",
    "sphinx.ext.doctest",
    "sphinx.ext.inheritance_diagram",
    "sphinx.ext.intersphinx",
    "sphinx.ext.mathjax",
    "sphinx.ext.napoleon",
    "sphinx.ext.todo",
    "sphinx.ext.viewcode",
    "jupyter_sphinx",
]

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]

# The suffix(es) of source filenames.
# You can specify multiple suffix as a list of string:
source_suffix = ".rst"

# The master toctree document.
master_doc = "index"

# Enable nitpicky mode, which forces links to be non-broken
nitpicky = True
nitpick_ignore = [
    # Prevents sphinx nitpicky mode picking up on optional
    # (see https://github.com/sphinx-doc/sphinx/issues/6861)
    ("py:class", "optional"),
    # See https://github.com/numpy/numpy/issues/10039
    ("py:obj", "numpy.datetime64"),
    # There's no specific file or function classes to link to
    ("py:class", "file object"),
    ("py:class", "function"),
    ("py:obj", "function"),
    ("py:class", "any type"),
    ("py:class", "Unit('pix')"),
    ("py:class", "Unit('deg')"),
    ("py:class", "Unit('arcsec')"),
    ("py:class", "Unit('%')"),
    ("py:class", "Unit('s')"),
    ("py:class", "Unit('Angstrom')"),
    ("py:class", "Unit('arcsec / pix')"),
    ("py:class", "Unit('W / m2')"),
    ("py:class", "array-like"),
    ("py:obj", "parfive"),
    ("py:class", "string"),
    ("py:class", "floats"),
]

# -- Options for intersphinx extension ---------------------------------------
intersphinx_mapping = {
    "python": (
        "https://docs.python.org/3/",
        (None, "http://www.astropy.org/astropy-data/intersphinx/python3.inv"),
    ),
    "numpy": (
        "https://numpy.org/doc/stable/",
        (None, "http://www.astropy.org/astropy-data/intersphinx/numpy.inv"),
    ),
    "matplotlib": (
        "https://matplotlib.org/",
        (None, "http://www.astropy.org/astropy-data/intersphinx/matplotlib.inv"),
    ),
    "astropy": ("https://docs.astropy.org/en/stable/", None),
    "sunpy": ("https://docs.sunpy.org/en/stable", None),
    "pyvista": ("https://docs.pyvista.org/", None),
}

# -- pyvista configuration ---------------------------------------------------
# Manage errors
pyvista.set_error_output_file("errors.txt")
pyvista.OFF_SCREEN = True  # Not necessary - simply an insurance policy
# Preferred plotting style for documentation
pyvista.set_plot_theme("document")
pyvista.global_theme.window_size = np.array([512, 512]) * 2
