"""
Configuration file for the Sphinx documentation builder.
"""

import os
from datetime import datetime
from pathlib import Path

import pyvista
from pyvista.plotting.utilities.sphinx_gallery import DynamicScraper
from sunpy_sphinx_theme import PNG_ICON

from sunkit_pyvista import __version__

# -- Read the Docs Specific Configuration --------------------------------------
# This needs to be done before sunpy is imported
on_rtd = os.environ.get("READTHEDOCS", None) == "True"
if on_rtd:
    os.environ["SUNPY_CONFIGDIR"] = "/home/docs/"
    os.environ["HOME"] = "/home/docs/"
    os.environ["LANG"] = "C"
    os.environ["LC_ALL"] = "C"
    os.environ["PARFIVE_HIDE_PROGRESS"] = "True"
    os.environ["PYDEVD_DISABLE_FILE_VALIDATION"] = "1"

# -- Project information -----------------------------------------------------
project = "sunkit-pyvista"
author = "SunPy Community"
copyright = f"{datetime.now().year}, {author}"  # NOQA: A001, DTZ005
release = __version__

# -- General configuration ---------------------------------------------------
extensions = [
    "sphinx_automodapi.automodapi",
    "sphinx_automodapi.smart_resolver",
    "sphinx_changelog",
    "sphinx_gallery.gen_gallery",
    "pyvista.ext.viewer_directive",
    "sphinx.ext.autodoc",
    "sphinx.ext.coverage",
    "sphinx.ext.doctest",
    "sphinx.ext.inheritance_diagram",
    "sphinx.ext.intersphinx",
    "sphinx.ext.mathjax",
    "sphinx.ext.napoleon",
    "sphinx.ext.todo",
    "sphinx.ext.viewcode",
    "sphinx_design",
]

html_theme = "sunpy"

# For the linkcheck
linkcheck_ignore = [
    r"https://doi.org/\d+",
    r"https://element.io/\d+",
    r"https://github.com/\d+",
    r"https://docs.sunpy.org/\d+",
]
linkcheck_anchors = False

exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]
source_suffix = ".rst"
master_doc = "index"
nitpicky = True
nitpick_ignore = []
with Path("nitpick-exceptions.txt").open() as f:
    for line in f:
        if line.strip() == "" or line.startswith("#"):
            continue
        dtype, target = line.split(None, 1)
        target = target.strip()
        nitpick_ignore.append((dtype, target))

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
    "sunkit_magex": ("https://docs.sunpy.org/projects/sunkit-magex/en/stable/", None),
}

# -- pyvista configuration ---------------------------------------------------
os.environ["PYVISTA_BUILDING_GALLERY"] = "True"
pyvista.BUILDING_GALLERY = True
pyvista.global_theme.font.label_size = 18
pyvista.global_theme.font.size = 18
pyvista.global_theme.font.title_size = 18
pyvista.global_theme.window_size = [512, 512]
pyvista.OFF_SCREEN = True
pyvista.set_error_output_file("errors.txt")
pyvista.set_plot_theme("document")
# We also need to start this on CI services and GitHub Actions has a CI env var
if on_rtd or os.environ.get("CI"):
    pyvista.start_xvfb()

# -- Sphinx Gallery ------------------------------------------------------------
sphinx_gallery_conf = {
    "backreferences_dir": Path("generated") / "modules",
    "filename_pattern": "^((?!skip_).)*$",
    "examples_dirs": Path("..") / "examples",
    "gallery_dirs": Path("generated") / "gallery",
    "matplotlib_animations": True,
    "default_thumb_file": PNG_ICON,
    "abort_on_example_error": False,
    "plot_gallery": "True",
    "remove_config_comments": True,
    "doc_module": ("sunpy"),
    "only_warn_on_example_error": True,
    "image_scrapers": (DynamicScraper(), "matplotlib"),
}
