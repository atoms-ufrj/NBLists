typedef struct {
  int*     start;          // Location of the first neighbor of each atom
  int*     end;            // Location of the last neighbor of each atom
  int*     item;           // Indices of neighbor atoms
  int      nitems;         // Number of atoms in neighbor list
  int      nmax;           // Maximum number of atoms in neighbor list
  int      builds;         // Number of neighbor-list builds
  double   time;           // Time taken in neighbor list handling
  struct {
    _Bool thirdLaw;
    _Bool jointXYZ;
    _Bool zeroBase;
  } options;
  void*    Data;           // Pointer to system data
} nbList;

nbList neighbor_list( int threads, double rc, double skin, int N, int* group );
void   neighbor_allocate_array( nbList*, double*** array, _Bool jointXYZ );
_Bool  neighbor_list_outdated( nbList* list, double** coords, int N, int* atoms );
void   neighbor_list_build( nbList* list, double L, double** coords );

