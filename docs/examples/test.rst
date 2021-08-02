=======================================
Three dimensional plots with sunpy Maps
=======================================

Using sunkit-pyvista, one can interface with the `pyvista` package to
produce interactive 3D plots for sunpy Maps.

.. jupyter-execute::
    :hide-code:

    import pyvista
    pyvista.set_jupyter_backend('ipygany')

.. jupyter-execute::

    import astropy.constants as const
    import astropy.units as u
    from astropy.coordinates import SkyCoord
    from sunpy.data.sample import AIA_193_IMAGE
    from sunpy.map import Map

    from sunkit_pyvista import SunpyPlotter

    m = Map(AIA_193_IMAGE).resample([512, 512] * u.pixel)

    p = SunpyPlotter()
    p.plot_map(m)
    p.show(jupyter_backend='ipygany')
