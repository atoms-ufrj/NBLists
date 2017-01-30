#include <nblists.h>
#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <time.h>

typedef struct {
  int N, Nx3, seed, Nsteps, Nprop;
  double rho, L, Rc, Rs, Rc2, Temp, Press, Dt, InvL, Ws, SixWs, Ec, Dt_2, Potential, Virial;
  double **R, **V, **F;
} simpar;

#define M_PI 3.14159265358979323846
#define TRUE 1
#define FALSE 0

//--------------------------------------------------------------------------------------------------

double drand()   /* uniform distribution, (0..1] */
{
  return (rand() + 1.0)/(RAND_MAX + 1.0);
}

//--------------------------------------------------------------------------------------------------

double random_normal()  /* normal distribution, centered on 0, std dev 1 */
{
  return sqrt(-2*log(drand())) * cos(2*M_PI*drand());
}

//--------------------------------------------------------------------------------------------------

void read_data( simpar *par, char *filename )
{
  FILE *file;
  file = fopen(filename,"r");
  #define readline \
    if (fgets(line, sizeof(line), file)) \
    if (!fgets(line, sizeof(line), file)) { \
      printf("ERROR: could not read data."); \
      exit(0); \
    }
  double InvRc2, InvRc6, InvRc12;
  char line[256];
  if (file != NULL) {
    readline; sscanf( line, "%d",  &par->N );
    readline; sscanf( line, "%lf", &par->Rc );
    readline; sscanf( line, "%lf", &par->Rs );
    readline; sscanf( line, "%d",  &par->seed );
    readline; sscanf( line, "%lf", &par->Dt );
    readline; sscanf( line, "%d",  &par->Nsteps );
    readline; sscanf( line, "%d",  &par->Nprop );
    readline; sscanf( line, "%lf", &par->rho );
    readline; sscanf( line, "%lf", &par->Temp );
    readline; sscanf( line, "%lf", &par->Press );
  }
  #undef readline
  par->Rc2 = par->Rc*par->Rc;
  par->L = pow(par->N/par->rho,1.0/3.0);
  par->InvL = 1.0/par->L;
  InvRc2 = 1.0/par->Rc2;
  InvRc6 = InvRc2*InvRc2*InvRc2;
  InvRc12 = InvRc6*InvRc6;
  par->Ws = 2.0*InvRc12 - InvRc6;
  par->SixWs = 6.0*par->Ws;
  par->Ec = InvRc12 - InvRc6;
  par->Dt_2 = 0.5*par->Dt;
  par->Nx3 = 3*par->N;
}

//--------------------------------------------------------------------------------------------------

void create_configuration( simpar *par )
{
  double Vcm[3] = { 0.0, 0.0, 0.0 };
  int Nd = ceil(pow(par->N,1.0/3.0));
  int a[3];
  for (int i = 0; i < par->N; i++) {
    a[2] = i/(Nd*Nd);
    a[1] = (i - a[2]*Nd*Nd)/Nd;
    a[0] = i - a[1]*Nd - a[2]*Nd*Nd;
    for (int j = 0; j < 3; j++) {
      par->R[i][j] = (par->L/Nd)*a[j] + 0.5;
      par->V[i][j] = random_normal();
      Vcm[j] += par->V[i][j];
    }
  }
  for (int j = 0; j < 3; j++)
    Vcm[j] /= par->N;
  double factor = 0.0;
  for (int i = 0; i < par->N; i++)
    for (int j = 0; j < 3; j++) {
      par->V[i][j] -= Vcm[j];
      factor += par->V[i][j]*par->V[i][j];
    }
  factor = sqrt(par->Temp*(3*par->N - 3)/factor);
  for (int i = 0; i < par->N; i++)
    for (int j = 0; j < 3; j++)
      par->V[i][j] *= factor;
}

//--------------------------------------------------------------------------------------------------

void compute_pairs( simpar *par, nbList *list )
{
  int i, j, k, m;
  double **R, **F, E, W, vec[3], r2, invR2, invR6, invR12, Eij, Wij, Fijk;

  if (neighbor_list_outdated( list, par->R[0]))
    neighbor_list_build( list, par->L, par->R[0] );

  E = 0.0;
  W = 0.0;
  R = par->R;
  F = par->F;
  for (i = 0; i < par->N; i++)
    F[i][0] = F[i][1] = F[i][2] = 0.0;

  for (i = 0; i < par->N; i++) {
    for (m = list->start[i]; m <= list->end[i]; m++) {
      j = list->item[m];
      for (k = 0; k < 3; k++) {
        vec[k] = R[i][k] - R[j][k];
        vec[k] = vec[k] - par->L*round(par->InvL*vec[k]);
      }
      r2 = vec[0]*vec[0] + vec[1]*vec[1] + vec[2]*vec[2];
      if (r2 < par->Rc2) {
        invR2 = 1.0/r2;
        invR6 = invR2*invR2*invR2;
        invR12 = invR6*invR6;
        Eij = invR12 - invR6;
        E += Eij;
        Wij = invR12 + Eij;
        W += Wij;
        Wij *= 24.0*invR2;
        for (k = 0; k < 3; k++) {
          Fijk = Wij*vec[k];
          F[i][k] += Fijk;
          F[j][k] -= Fijk;
        }
      }
    }
  }
  par->Potential = 4.0*E;
  par->Virial = 8.0*W;
}

//--------------------------------------------------------------------------------------------------

void compute_all_pairs( simpar *par )
{
  int i, j, k;
  double **R, **F, E, W, vec[3], r2, invR2, invR6, invR12, Eij, Wij, Fijk;

  E = 0.0;
  W = 0.0;
  R = par->R;
  F = par->F;
  for (i = 0; i < par->N; i++)
    F[i][0] = F[i][1] = F[i][2] = 0.0;

  for (i = 0; i < par->N-1; i++) {
    for (j = i+1; j < par->N; j++) {
      for (k = 0; k < 3; k++) {
        vec[k] = R[i][k] - R[j][k];
        vec[k] = vec[k] - par->L*round(par->InvL*vec[k]);
      }
      r2 = vec[0]*vec[0] + vec[1]*vec[1] + vec[2]*vec[2];
      if (r2 < par->Rc2) {
        invR2 = 1.0/r2;
        invR6 = invR2*invR2*invR2;
        invR12 = invR6*invR6;
        Eij = invR12 - invR6;
        E += Eij;
        Wij = invR12 + Eij;
        W += Wij;
        Wij *= 24.0*invR2;
        for (k = 0; k < 3; k++) {
          Fijk = Wij*vec[k];
          F[i][k] += Fijk;
          F[j][k] -= Fijk;
        }
      }
    }
  }
  par->Potential = 4.0*E;
  par->Virial = 8.0*W;
}

//--------------------------------------------------------------------------------------------------

double energy( simpar par )
{
  double K = 0.0;
  double *V = par.V[0];
  for (int i = 0; i < par.Nx3; i++)
    K += V[i]*V[i];
  return par.Potential + 0.5*K;
}

//--------------------------------------------------------------------------------------------------

int main( int argc, char *argv[] )  {
  int threads;
  char *filename;
  if (argc == 2) {
    threads = 1;
    filename = argv[1];
  }
  else if (argc == 3) {
    threads = atoi(argv[1]);
    filename = argv[2];
  }
  else {
    printf("Usage: testc [number-of-threads] input-file\n");
    exit(EXIT_FAILURE);
  }
  simpar par;
  read_data( &par, filename );
  nbList list = neighbor_list( threads, par.Rc, par.Rs, par.N, NULL );
  neighbor_allocate_2d_array( &list, &par.R, TRUE );
  neighbor_allocate_2d_array( &list, &par.V, TRUE );
  neighbor_allocate_2d_array( &list, &par.F, TRUE );
  create_configuration( &par );

  compute_pairs( &par, &list );
  printf("Step Potential Virial Total\n");
  printf("%d %lf %lf %lf\n", 0, par.Potential, par.Virial, energy( par ));
  for (int step = 1; step <= par.Nsteps; step++) {
    if (step % par.Nprop == 0) printf("%d %lf %lf %lf\n", step, par.Potential, par.Virial, energy( par ));
      for (int i = 0; i < par.N; i++)
        for (int j = 0; j < 3; j++) {
          par.V[i][j] += par.F[i][j]*par.Dt_2;
          par.R[i][j] += par.V[i][j]*par.Dt;
        }
/*      compute_all_pairs( &par );*/
      compute_pairs( &par, &list );
      for (int i = 0; i < par.N; i++)
        for (int j = 0; j < 3; j++)
          par.V[i][j] += par.F[i][j]*par.Dt_2;
  }
  printf("neighbor list builds = %d\n", list.builds);
  printf("list building time = %f s.\n", list.time);

  return EXIT_SUCCESS;
}

