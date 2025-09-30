===================
CAMB
===================
:CAMB: Code for Anisotropies in the Microwave Background
:Author: Antony Lewis and Anthony Challinor
:Homepage: https://camb.info/

.. image:: https://img.shields.io/pypi/v/camb.svg?style=flat
   :target: https://pypi.python.org/pypi/camb/
.. image:: https://img.shields.io/conda/vn/conda-forge/camb.svg
   :target: https://anaconda.org/conda-forge/camb
.. image:: https://readthedocs.org/projects/camb/badge/?version=latest
   :target: https://camb.readthedocs.io/en/latest
.. image:: https://github.com/cmbant/camb/actions/workflows/tests.yml/badge.svg?branch=master
  :target: https://github.com/cmbant/CAMB/actions
.. image:: https://mybinder.org/badge_logo.svg
  :target: https://mybinder.org/v2/gh/cmbant/CAMB/HEAD?filepath=docs%2FCAMBdemo.ipynb

Installation of current gCAMB version
=====================================================
Clode the repository:

```
$ git clone https://github.com/lstorchi/CAMB.git
Cloning into 'CAMB'...
remote: Enumerating objects: 8834, done.
remote: Counting objects: 100% (3096/3096), done.
remote: Compressing objects: 100% (375/375), done.
remote: Total 8834 (delta 2956), reused 2721 (delta 2721), pack-reused 5738 (from 2)
Receiving objects: 100% (8834/8834), 112.37 MiB | 35.22 MiB/s, done.
Resolving deltas: 100% (6557/6557), done.
```

checkout the proper branch: 

```
$ cd CAMB/
$ git checkout gpuport
branch 'gpuport' set up to track 'origin/gpuport'.
Switched to a new branch 'gpuport'
```

Clone forutils:

```
$ git clone  https://github.com/lstorchi/forutils.git
Cloning into 'forutils'...
remote: Enumerating objects: 451, done.
remote: Counting objects: 100% (91/91), done.
remote: Compressing objects: 100% (68/68), done.
remote: Total 451 (delta 54), reused 56 (delta 23), pack-reused 360 (from 1)
Receiving objects: 100% (451/451), 138.92 KiB | 2.67 MiB/s, done.
Resolving deltas: 100% (288/288), done.
```

and checkout the proper branch: 

```
$ cd forutils/
$ git checkout gpuport
branch 'gpuport' set up to track 'origin/gpuport'.
Switched to a new branch 'gpuport'
```

Now we can compile the code using the nvfortran compiler:

```
$ cd ../fortran/
$ make 
make -C Release --no-print-directory -f../Makefile FORUTILS_SRC_DIR=.. libforutils.a
nvfortran -DUSEACC -cpp -Mextend -acc=gpu -gpu=deepcopy -Minfo=accel  -openmp -O3 -o MiscUtils.o -c ../Mi
scUtils.f90
.........
```

test the code:

```
$ export OMP_NUM_THREADS=32
$ export NVCOMPILER_ACC_CUDA_HEAPSIZE=3G
$ time ./camb  params_high.ini
 Start reading params_high.ini
 Read number_of_threads            0
 Read DebugParam    0.000
........
```

you can check the exectution of the code on the GPU via nvidia-smi -l 1

Description and installation
=============================

CAMB is a cosmology code for calculating cosmological observables, including
CMB, lensing, source count and 21cm angular power spectra, matter power spectra, transfer functions
and background evolution. The code is in Python, with numerical code implemented in fast modern Fortran.

See the `CAMB python example notebook <https://camb.readthedocs.io/en/latest/CAMBdemo.html>`_ for a
quick introduction to how to use the CAMB Python package.

For a standard non-editable installation use::

    pip install camb [--user]

The --user is optional and only required if you don't have write permission to your main python installation.
To install from source, clone from github using::

    git clone --recursive https://github.com/cmbant/CAMB

Then install using::

    pip install -e ./CAMB [--user]

You will need gfortran 6 or higher installed to compile (usually included with gcc by default).
If you have gfortran installed, "python setup.py make" (and other standard setup commands) will build the Fortran
library on all systems (including Windows without directly using a Makefile).

The python wrapper provides a module called "camb" documented in the Python `CAMB documentation <https://camb.readthedocs.io/en/latest/>`_.

After installation you can also run CAMB from the command line reading parameters from a .ini file, e.g.::

  camb inifiles/planck_2018.ini

To compile the Fortran command-line code run "make camb" in the fortran directory. For full details
see the  `ReadMe <https://camb.info/readme.html>`_.

Branches
=============================

The master branch contains latest changes to the main release version.

There is a test suite, which runs automatically on GitHub actions for new commits and pull requests.
Reference results and test outputs are stored in the `test outputs repository <https://github.com/cmbant/CAMB_test_outputs/>`_. Tests can also be run locally.

To reproduce legacy results, see these branches:

 - *CAMB_sources* is the old public `CAMB Sources <https://camb.info/sources/>`_ code.
 - *CAMB_v0* is the old Fortran-oriented (gfortran 4.8-compatible) version as used by the Planck 2018 analysis.
 - *rayleigh* includes frequency-dependent Rayleigh scattering
 - *python2* is the last Python 2 compatible version

===================

.. raw:: html

    <a href="https://www.sussex.ac.uk/astronomy/"><img src="https://cdn.cosmologist.info/antony/Sussex_white.svg" style="height:200px" height="200px"></a>
    <a href="https://erc.europa.eu/"><img src="https://cdn.cosmologist.info/antony/ERC_white.svg" style="height:200px" height="200px"></a>
    <a href="https://stfc.ukri.org/"><img src="https://cdn.cosmologist.info/antony/STFC_white.svg" style="height:200px" height="200px"></a>
