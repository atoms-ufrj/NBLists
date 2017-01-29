typedef struct {
  int*     first;          // Location of the first neighbor of each atom
  int*     last;           // Location of the last neighbor of each atom
  int*     item;           // Indices of neighbor atoms
  double** R0;             // Atom positions at latest neighbor list build
  int      nitems;         // Number of atoms in neighbor list
  int      nmax;           // Maximum number of atoms in neighbor list
  int      builds;         // Number of neighbor-list builds
  double   time;           // Time taken in neighbor list handling
  void*    Data;           // Pointer to system data
  struct {
    _Bool thirdLaw;
    _Bool jointXYZ;
  } Options;
} nbList;

nbList neighbor_list( int threads, double rc, double skin, int N, int* body );
_Bool  neighbor_list_outdated( nbList* list, double* positions );
void   neighbor_list_build( nbList* list, double Lbox, double* positions );

