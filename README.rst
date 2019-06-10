=====================
🐛 agr_file_generator
=====================

A package for generating VCF files for the Alliance of Genome Resources.


Developemnt
===========

Setup
-----
Use virtualenv_ for developing.

.. note::

   You may find that virtualenvwrapper_ is a more convenient tool.

Create a virtualenv:

.. code-block:: bash

   VENV_HOME="~/.virtualenv"
   mkdir -p "${VENV_HOME}"
   VENV="${VENV_HOME}/azanium"
   python3 -m venv "${VENV}"
   # Activate the virtualenv (type deactivate to return to normal shell)
   source "${VENV}/bin/activate"


Install the source code in "editable" mode.

.. note::

   This means that we don't have to re-install our package
   each time we make a change.

.. code-block:: bash

   pip install --editable ".[dev]"


Running tests
-------------

.. code-block:: bash

   pytest


Running the file generator
==========================

.. code-block:: bash

   agr_file_generator


.. _virtualenv: http://docs.python-guide.org/en/latest/dev/virtualenvs/
.. _virtualenvwrapper: https://virtualenvwrapper.readthedocs.io/en/latest/

