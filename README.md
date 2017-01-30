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

### Downloading the library:

* Clone the repository

        git clone https://github.com/atoms-ufrj/NBLists

  This will create a local copy of the repository on your device.

### Compiling the library:

* Execute the Makefile in the root directory of the repository tree

        make

  This will build the shared library (lib/libnblists.so).

### Installing the library in the system path:

* Use building option **install**

        sudo make install

  This will copy the static library (lib/libnblists.so) to `/usr/local/lib` and the C header file
(include/nblists.h) and fortran module (include/nblists.f90) to `/usr/local/include`.

## Uninstalling the library from the system path:

* Use building option **uninstall**

        sudo make uninstall


Using the C Binding
-------------------

```c++
    #include <nblists.h>

    int threads; // Number of parallel threads to be used
    double rc;   // Cutoff distance for pair interactions
    double skin; // Extra distance for the Verlet-type neighbor list
    int N;       // Number of atoms
    double **R;  // Array of atom coordinates
    double L;    // Side length of a cubic box with periodic boundary condition
    ...

    // If groups are defined so that same-group atom pairs are ignored:
    int *group = new int[N]; // group[i] is the one that contains atom i
    ...
    nbList list1 = neighbor_list( threads, rc, skin, N, group );

    // Otherwise, in order to consider all possible atom pairs:
    nbList list2 = neighbor_list( threads, rc, skin, N, NULL );

    // If either pair i-j or j-i must appear in the neighbor list:
    list1.Options.thirdLaw = True; // Default

    // If both pairs i-j and j-i must appear in the neighbor list
    list2.Options.thirdLaw = False;

    // If coordinates are stored as R = [x1,y1,z1, x2,y2,z2, ..., xN,yN,zN]:
    list1.Options.jointXYZ = True; // Default

    // If coordinates are stored as R = [x1,x2,...,xN, y1,y2,...,yN, z1,z2,...,zN]:
    list2.Options.jointXYZ = False;

    // If C-style array indexing (zero-based) is employed:
    list1.Options.zeroBase = True; // Default

    // Basic neighbor list usage in Molecular Dynamics:
    if (neighbor_list_outdated( list1, R, 0, NULL ))
      neighbor_list_build( &list1, L, R );
    for (int i = 0; i < N; i++)
      for (int k = list1.first[i]-1; k < list1.last[i]; k++) {
        int j = list1.item[k] - 1;
        // Compute interaction of pair i-j
        ...
      }

    // Basic neighbor list usage in Monte Carlo:
    int M;          // Number of atoms being moved
    int *moved;     // Index of each atom being moved
    double **Rnew;  // New coordinates of the M atoms

    for (int i = 0; i < N; i++)
      for (int k = list1.first[i]-1; k < list1.last[i]; k++) {
        int j = list1.item[k] - 1;
        // Compute interaction of pair i-j
        ...
      }
    

    if (atom_i_move_accepted) {
      double d[3];
      for (int k=0; k < 3; k++) d[k] = R[i][k] - list2.R0[i][k];
      if (d[0]*d[0] + d[1]*d[1] + d[2]*d[2] < list2.maxSD)
        
    }

    if ()

    if neighbor_list_outdated( list1, R )
      neighbor_list_build( &list1, L, R );

    for (int i = 0; i < N; i++)
      for (int k = list1.first[i]-1; k < list1.last[i]; k++) {
        int j = list1.item[k] - 1;
        // Compute interaction of pair i-j
        ...
      }
```
