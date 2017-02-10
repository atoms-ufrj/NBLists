NBLists: Neighbor Lists for Molecular Simulation
================================================

NBLists is a library of procedures used to manage neighbor lists in molecular simulation codes. Its
main features are:

1. Works with both Molecular Dynamics (considering Newton's 3rd Law) or Monte Carlo simulations.
2. Uses OpenMP for improved performance in shared-memory, parallel processors.
3. Includes a C header and a Fortran-90 module, thus being easily integrable to simulation codes
written in various programming languages.

Standard compilation
--------------------

### Dependencies (considering a Ubuntu 16.04.1 LTS fresh install):

* gfortran

> **Tested with:**
> - GNU Fortran (Ubuntu 5.4.0-6ubuntu1~16.04.2) 5.4.0 20160609
>
> **Can usually be installed via apt:**
>
>      sudo apt-get install gfortran

### Downloading the source code:

* Clone the repository:

        git clone https://github.com/atoms-ufrj/NBLists

  This will create a local copy of the repository on your device.

### Updating the repository:

* Execute the following command in the root directory of the repository tree:

        git pull

  This will synchronize the local copy with the online repository.

### Compiling the library:

* Execute the Makefile in the root directory of the repository tree:

        make

  This will build the shared library (lib/libnblists.so).

### Testing the compiled library:

* Use building option **test**

        make test

  This will compile and run some testing codes in the `test` directory.

### Installing the library in the system path:

* Use building option **install**

        sudo make install

  This will copy the static library (lib/libnblists.so) to `/usr/local/lib` and the C header file
(include/nblists.h) and fortran module (include/nblists.f90) to `/usr/local/include`.

### Uninstalling the library from the system path:

* Use building option **uninstall**

        sudo make uninstall


Using the Library
-----------------

```c++
    #include <nblists.h>

    int threads; // Number of parallel threads to be used
    double rc;   // Cutoff distance for pair interactions
    double skin; // Extra distance for the Verlet-type neighbor list
    int N;       // Number of atoms
    double **R;  // Pointers to a contiguously stored array of atomic coordinates
    double L;    // Side length of a cubic box with periodic boundary condition
    int *group;  // Either NULL or a pointer to an array of group indices
    ...

    // Initializing the list:
    nbList list = neighbor_list( threads, rc, skin, N, group );
    /* If group = NULL, all possible atom pairs will be condidered.
       Otherwise, if each atom i belongs to a group whose index is stored in group[i],
       same-group atom pairs will ignored during the neighbor list building. */

    // Determining how coordinates are stored in memory:
    list.options.jointXYZ = 1; // Default
    /* Coordinates are stored as [x1,y1,z1, ..., xN,yN,zN]. If they are stored
       as [x1,...,xN, y1,...,yN, z1,...,zN] instead, then make jointXYZ = 0. */

    // Determining whether Newton's Third Law is taken into account:
    list.options.thirdLaw = 1; // Default
    /* Whenever distance(i,j) < Rc + Rs, either j will be listed as a neighbor of i
       or vice-versa. If both pairs must be listed, then make thirdLaw = 0. */

    // Determining which base is used for array indexing:
    list.options.zeroBase = 1; // Default
    /* A C-style array indexing (0-based) is employed. If a Fortran-style
       indexing (1-based) is to be used, then make zeroBase = 0. */

    // Using the library to allocate a contiguous array of coordinates (optional):
    _Bool jointXYZ = 1;
    neighbor_allocate_2d_array( &R, N, jointXYZ );
    ...

    // Basic neighbor list usage in Molecular Dynamics:
    if (neighbor_list_outdated( list, R[0] ))
      neighbor_list_build( &list, L, R[0] );
    for (int i = 0; i < N; i++)
      for (int k = list.start[i]; k <= list.end[i]; k++) {
        int j = list.item[k];
        // Compute interaction of pair i-j
        ...
      }
```
