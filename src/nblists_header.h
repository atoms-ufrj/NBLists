typedef struct {
  int*   first;          // Location of the first neighbor of each atom
  int*   last;           // Location of the last neighbor of each atom
  int*   item;           // Indices of neighbor atoms
  int    nitems;         // Number of atoms in neighbor list
  int    nmax;           // Maximum number of atoms in neighbor list
  int    builds;         // Number of neighbor-list builds
  double time;           // Time taken in neighbor list handling
  void*  Data;           // Pointer to system data
} nbList;

nbList neighbor_list( int threads, double rc, double skin, int N, int* body, int MC );
void neighbor_handle( nbList* list, double Lbox, double* positions );

