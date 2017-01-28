NBLists: A Neighbor List Handler for Molecular Simulation
=========================================================

NBLists is a shared library that can manage neighbor lists in molecular simulations. Its main
features are:

1. Serves for both Molecular Dynamics (with Newton's 3rd Law) or Monte Carlo simulations.
2. Uses OpenMP for improved performance in shared-memory, parallel processors.
3. Includes both a C header file and a Fortran-90 module, so that the library can be easily employed
in simulation codes written in various programming languages.

Standard compilation
--------------------

### Dependencies (considering a Ubuntu 16.04.1 LTS fresh install):

#### Basic for compiling the library

* gfortran

> **Tested with:**
> - GNU Fortran (Ubuntu 5.4.0-6ubuntu1~16.04.2) 5.4.0 20160609
>
> **Can usually be installed via apt:**
>
>      sudo apt-get install gfortran

### Compiling the library:

* Clone the repository

        git clone https://github.com/atoms-ufrj/NBLists

  This will create a local copy of the repository on your device.

* Execute the Makefile in the root directory of the repository tree

        make

  This will build the shared library (lib/libnblists.so).

## Installing the library in the system path:

* Use building option **install**

        sudo make install

  This will copy the static library (lib/libnblists.so) to `/usr/local/lib` and the C header file
(include/nblists.h) and fortran module (include/nblists.f90) to `/usr/local/include`.

## Uninstalling the library from the system path:

* Use building option **uninstall**

        sudo make uninstall

