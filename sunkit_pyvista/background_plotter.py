from sunkit_pyvista.mapsequence_animator import SequenceAnimator
from sunkit_pyvista import SunpyPlotter
from functools import partial
from sunpy.util import expand_list
import pyvistaqt as pvq
import functools
from functools import partial
from pathlib import Path

import numpy as np
import pyvista as pv


import astropy.units as u
from astropy.constants import R_sun
from astropy.coordinates import Longitude, SkyCoord
from astropy.visualization import AsymmetricPercentileInterval
from sunpy.coordinates import HeliocentricInertial
from sunpy.coordinates.utils import get_rectangle_coordinates
from sunpy.map import GenericMap
from sunpy.map.maputils import all_corner_coords_from_map
from sunpy.util import expand_list
from sunpy.visualization._quadrangle import Quadrangle

from sunkit_pyvista.mapsequence_animator import SequenceAnimator

class SunpyBackgroundPlotter(SunpyPlotter):
    def __init__(self, coordinate_frame=None):
        super().__init__(coordinate_frame=coordinate_frame)
        self._plotter = pvq.BackgroundPlotter()

    @property
    def coordinate_frame(self):
        """
        Coordinate frame of the plot.
        """
        return self._coordinate_frame

    @property
    def plotter(self):
        """
        `pyvistaqt.BackgroundPlotter`.
        """
        return self._plotter

    def _toggle_animation(self, state, animate):
        animate.animation_state = state

    def plot_map_sequence(self, *args, interval=2, **kwargs):
        """
        Plot a sequence of maps as an animation.

        Parameters
        ----------
        m : `sunpy.map.Map`
            Map(s) to be plotted.
        **kwargs :
            Keyword arguments are handed to `pyvistaq.BackGroundplotter.add_mesh`.
        """
        map_meshes = []
        color_maps = []

        maps = expand_list(args)
        for m in maps:
            if not isinstance(m, GenericMap):
                raise ValueError(
                    'MapSequence expects pre-constructed map objects.')
        for m in maps:
            mesh = self._pyvista_mesh(m)
            map_meshes.append(mesh)
            color_maps.append(m.cmap)
        animate = SequenceAnimator(time=interval, map_meshes=map_meshes,
                                   color_maps=color_maps, **kwargs)

        self.plotter.add_mesh(map_meshes[0], cmap=color_maps[0], **kwargs)
        self.plotter.add_checkbox_button_widget(
            partial(self._toggle_animation, animate=animate),
            value=False, color_on='green')
        self.plotter.add_callback(partial(animate,
                                     bg_plotter=self.plotter,
                                     **kwargs), interval=16)
        self.plotter.add_text("Play", position=(70, 10))
        self.plotter.enable_anti_aliasing()
        self.plotter.hide_axes()
