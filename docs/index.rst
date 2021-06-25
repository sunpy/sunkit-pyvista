****************************
sunkit-pyvista Documentation
****************************

.. jupyter-execute::
   :hide-code:

   # using ipyvtk as it loads faster
   import pyvista
   pyvista.set_jupyter_backend('ipygany')

`sunkit-pyvista` is a python package for visualizing solar physics data in 3D.
As a quick example:

.. jupyter-execute::

    grid = pyvista.Sphere()
    plotter = pyvista.Plotter()
    plotter.add_mesh(grid)
    plotter.show()

.. toctree::
   :maxdepth: 2

   api
   generated/gallery/index
   changelog

Indexes
=======

* :ref:`genindex`
* :ref:`modindex`
