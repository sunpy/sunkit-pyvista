"""
Configuration file for the Sphinx documentation builder.
"""

import os
from datetime import datetime
from pathlib import Path

import pyvista
from packaging.version import Version
from pyvista.plotting.utilities.sphinx_gallery import DynamicScraper
from sunpy_sphinx_theme.conf import PNG_ICON

from sunkit_pyvista import __version__

# -- Project information -----------------------------------------------------
on_rtd = os.environ.get("READTHEDOCS")
os.environ["HIDE_PARFIVE_PROGESS"] = "True"
os.environ["PYDEVD_DISABLE_FILE_VALIDATION"] = "1"
project = "sunkit-pyvista"
author = "SunPy Community"
copyright = f"{datetime.now(datetime.timezone.utc).year}, {author}"  # NOQA: A001
release = __version__
sunkit_pyvista_version = Version(__version__)
is_release = not (sunkit_pyvista_version.is_prerelease or sunkit_pyvista_version.is_devrelease)

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
