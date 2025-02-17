**************
sunkit-pyvista
**************

|Latest Version| |Latest Documentation| |codecov| |matrix| |Powered by NumFOCUS|

.. |Latest Documentation| image:: https://readthedocs.org/projects/sunkit-pyvista/badge/?version=latest
   :target: https://docs.sunpy.org/projects/sunkit-pyvista/en/latest/?badge=latest
   :alt: Documentation Status
.. |Latest Version| image:: https://img.shields.io/pypi/v/sunkit-pyvista.svg
   :target: https://pypi.python.org/pypi/sunkit-pyvista/
.. |matrix| image:: https://img.shields.io/matrix/sunpy:openastronomy.org.svg?colorB=%23FE7900&label=Chat&logo=matrix&server_fqdn=openastronomy.modular.im
   :target: https://openastronomy.element.io/#/room/#sunpy:openastronomy.org
.. |codecov| image:: https://codecov.io/gh/sunpy/sunkit-pyvista/branch/main/graph/badge.svg
   :target: https://codecov.io/gh/sunpy/sunkit-pyvista
.. |Powered by NumFOCUS| image:: https://img.shields.io/badge/powered%20by-NumFOCUS-orange.svg?style=flat&colorA=E1523D&colorB=007D8A
   :target: https://numfocus.org

The goal of sunkit-pyvista is to allow 3D visualization of solar physics data with pyvista.

Installation
============

See the `installation instructions <https://docs.sunpy.org/projects/sunkit-pyvista/en/latest/installing.html>`__ for more information.

Developing
==========

If you want to develop sunkit-pyvista you will need to install it from GitHub.
For detailed installation instructions, see `Development installation`_ in the sunpy docs.
These are for "sunpy" but the same works for sunkt-pyvista with a quick name substitution.

.. _Development installation:  https://docs.sunpy.org/en/latest/dev_guide/contents/newcomers.html#setting-up-a-development-environment

Running Figure tests
--------------------

You need to use your local environment and make sure its up to date.
This does not with tox for now.

2 command line arguments can be passed to pytest:

* ``update_image_cache`` : Resets the ``image_cache`` directory with new figures generated from your local environment.
* ``ignore_image_cache`` : Ignores the image cache and will not perform any image comparison.

Example usage:

.. code-block:: bash

   pytest -s -vvv --update_image_cache sunkit_pyvista/

Using tox for this does not work for now.

Getting Help
============

For more information or to ask questions about sunkit-pyvista, check out:

-  `sunkit-pyvista Documentation`_
-  `SunPy Element Channel`_

.. _sunkit-pyvista Documentation: https://docs.sunpy.org/projects/sunkit-pyvista/en/latest/
.. _SunPy Element Channel: https://app.element.io/#/room/#sunpy:openastronomy.org

Usage of Generative AI
----------------------

We expect authentic engagement in our community.
Be wary of posting output from Large Language Models or similar generative AI as comments on GitHub or any other platform, as such comments tend to be formulaic and low quality content.
If you use generative AI tools as an aid in developing code or documentation changes, ensure that you fully understand the proposed changes and can explain why they are the correct approach and an improvement to the current state.

Contributing
============

If you would like to get involved, check out the `Developers Guide`_ section of the SunPy docs.
Stop by our chat room `#sunpy:openastronomy.org`_ if you have any questions.
Help is always welcome so let us know what you like to work on, or check out the `issues page`_ for the list of known outstanding items.

For more information on contributing to SunPy, please read our `Newcomers' guide`_.

.. _Developers Guide: https://docs.sunpy.org/en/latest/dev_guide/index.html
.. _`#sunpy:openastronomy.org`: https://app.element.io/#/room/#sunpy:openastronomy.org
.. _issues page: https://github.com/sunpy/sunkit-pyvista/issues
.. _Newcomers' guide: https://docs.sunpy.org/en/latest/dev_guide/contents/newcomers.html

Code of Conduct
===============

When you are interacting with the SunPy community you are asked to follow our `Code of Conduct`_.

.. _Code of Conduct: https://sunpy.org/coc
