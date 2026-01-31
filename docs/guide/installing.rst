Installing
==========

.. important::

    IGUA is a Python library but relies on `MMseqs2 <https://mmseqs.com/>`_ to
    handle the different clustering stages. Most installation modes will also 
    ensure the ``mmseqs`` binary is installed on your machine or bundle a 
    working MMseqs2 install (if using a container). However, for direct Python 
    installs with ``pip``, you will need to ensure yourself that MMseqs2 is 
    available.

Conda
^^^^^

IGUA and all of its dependencies are available via `Bioconda <https://anaconda.org/channels/bioconda/packages/igua/overview>`_
and can be installed using e.g., ``conda`` or ``pixi``:

1. First, `set up Bioconda with Pixi or Conda. <https://bioconda.github.io/>`_.
2. Then, install IGUA using the appropriate method:

   - With ``conda``::
     
       $ conda install igua

   - With ``pixi``::
        
       $ pixi add igua


Apptainer, Docker, and Singularity
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

IGUA (and all of its dependencies) can be run using e.g., Docker, Apptainer, 
and Singularity, using images available on `quay.io <https://quay.io/repository/biocontainers/igua?tab=tags&tag=latest>`_.

For instance, to install IGUA 0.1.0 with AppTainer, use::

    $ apptainer pull docker://quay.io/biocontainers/igua:0.1.0--py39h5b94c0b_0


PyPi
^^^^

IGUA can be downloaded directly from PyPI, which hosts pre-compiled distributions
for Linux, MacOS and Windows. To install IGUA with ``pip``, use:

.. code:: console

    $ pip install --user igua

**Note that you need to have a working install of MMSeqs2 on your system, 
as MMseqs2 cannot be installed directly from PyPI.**


Arch User Repository
^^^^^^^^^^^^^^^^^^^^

A package recipe for Arch Linux can be found in the Arch User Repository
under the name `igua <https://aur.archlinux.org/packages/igua>`_.
It will always match the latest release from PyPI.

Steps to install on ArchLinux depend on your `AUR helper <https://wiki.archlinux.org/title/AUR_helpers>`_
(``yaourt``, ``aura``, ``yay``, etc.). For ``aura``, you'll need to run:

.. code:: console

    $ aura -A igua


GitHub + ``pip``
^^^^^^^^^^^^^^^^

If, for any reason, you prefer to download the library from GitHub, you can clone
the repository and install the repository by running (with the admin rights):

.. code:: console

    $ pip install -U git+https://github.com/zellerlab/IGUA

.. caution::

    Keep in mind this will install always try to install the latest commit,
    which may not even build, so consider using a versioned release instead.


GitHub + ``build``
^^^^^^^^^^^^^^^^^^

If you do not want to use ``pip``, you can still clone the repository and
use ``build`` and ``installer`` manually:

.. code:: console

    $ git clone https://github.com/zellerlab/IGUA
    $ cd IGUA
    $ python -m build .
    # python -m installer dist/*.whl

.. Danger::

    Installing packages without ``pip`` is strongly discouraged, as they can
    only be uninstalled manually, and may damage your system.
