#include <nblists.h>
#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <time.h>

typedef struct {
  int N, Nx3, seed, Npassos, Nprop;
  double rho, L, Rc, Rs, Rc2, Temp, Press, Dt, InvL, Ws, SixWs, Ec, Dt_2, Energy, Virial;
  double *R, *V, *F;
} simpar;

#define M_PI 3.14159265358979323846

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
    readline; sscanf( line, "%d",  &par->Npassos );
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
  par->R = malloc( par->Nx3*sizeof(double) );
  par->V = malloc( par->Nx3*sizeof(double) );
  par->F = malloc( par->Nx3*sizeof(double) );
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
      int k = 3*i + j;
      par->R[k] = (par->L/Nd)*a[j] + 0.5;
      par->V[k] = random_normal();
      Vcm[j] += par->V[k];
    }
  }
  for (int j = 0; j < 3; j++)
    Vcm[j] /= par->N;
  double factor = 0.0;
  for (int i = 0; i < par->N; i++)
    for (int j = 0; j < 3; j++) {
      int k = 3*i+j;
      par->V[k] -= Vcm[j];
      factor += par->V[k]*par->V[k];
    }
  factor = sqrt(par->Temp*(3*par->N - 3)/factor);
  for (int i = 0; i < 3*par->N; i++)
    par->V[i] *= factor;
}

//--------------------------------------------------------------------------------------------------

void compute( simpar *par, nbList *list )
{
  double *R = par->R;
  double *F = par->F;

  double E = 0.0;
  double W = 0.0;
  for (int i = 0; i < par->Nx3; i++)
    F[i] = 0.0;

  for (int i = 0; i < par->N; i++) {
    int ix = 3*i;
    int iy = ix + 1;
    int iz = iy + 1;
    for (int k = list->start[i]; k <= list->end[i]; k++) {
      int j = list->item[k];
      int jx = 3*j;
      int jy = jx + 1;
      int jz = jy + 1;
      double rx = R[ix] - R[jx];
      double ry = R[iy] - R[jy];
      double rz = R[iz] - R[jz];
      double r2 = rx*rx + ry*ry + rz*rz;
      if (r2 < par->Rc2) {
        double invR2 = 1.0/r2;
        double invR6 = invR2*invR2*invR2;
        double invR12 = invR6*invR6;
        double Eij = invR12 - invR6;
        E += Eij;
        double Wij = invR12 + Eij;
        W += Wij;
        Wij *= invR2;
        double Fijx = Wij*rx;
        double Fijy = Wij*ry;
        double Fijz = Wij*rz;
        F[ix] += Fijx;
        F[iy] += Fijy;
        F[iz] += Fijz;
        F[jx] -= Fijx;
        F[jy] -= Fijy;
        F[jz] -= Fijz;
      }
    }
  }
  par->Energy = 4.0*E;
  par->Virial = 8.0*W;
}

//--------------------------------------------------------------------------------------------------

void brute_force_compute( simpar *par )
{
  double *R = par->R;
  double *F = par->F;

  double E = 0.0;
  double W = 0.0;
  for (int i = 0; i < par->Nx3; i++)
    F[i] = 0.0;

  for (int i = 0; i < par->N-1; i++) {
    int ix = 3*i;
    int iy = ix + 1;
    int iz = iy + 1;
    for (int j = i+1; j < par->N; j++) {
      int jx = 3*j;
      int jy = jx + 1;
      int jz = jy + 1;
      double rx = R[ix] - R[jx];
      double ry = R[iy] - R[jy];
      double rz = R[iz] - R[jz];
      double r2 = rx*rx + ry*ry + rz*rz;
      if (r2 < par->Rc2) {
        double invR2 = 1.0/r2;
        double invR6 = invR2*invR2*invR2;
        double invR12 = invR6*invR6;
        double Eij = invR12 - invR6;
        E += Eij;
        double Wij = invR12 + Eij;
        W += Wij;
        Wij *= invR2;
        double Fijx = Wij*rx;
        double Fijy = Wij*ry;
        double Fijz = Wij*rz;
        F[ix] += Fijx;
        F[iy] += Fijy;
        F[iz] += Fijz;
        F[jx] -= Fijx;
        F[jy] -= Fijy;
        F[jz] -= Fijz;
      }
    }
  }
  par->Energy = 4.0*E;
  par->Virial = 8.0*W;
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
  create_configuration( &par );

  nbList list = neighbor_list( threads, par.Rc, par.Rs, par.N, NULL );

  neighbor_list_build( &list, par.L, par.R );
  compute( &par, &list );
  printf("Energy, Virial = %f %f\n",par.Energy,par.Virial);

  brute_force_compute( &par );
  printf("Energy, Virial = %f %f\n",par.Energy,par.Virial);

/*  if (neighbor_list_outdated( &list, par.R, par.N, NULL ))*/
    
/*  printf("%d\n", test );*/

/*
  tEmDee md = EmDee_system( threads, 1, par.Rc, par.Rs, par.N, NULL, NULL );
  void* lj = EmDee_pair_lj_cut( 1.0, 1.0 );
  EmDee_set_pair_model( md, 1, 1, lj );

  EmDee_upload( &md, "box", &par.L );
  EmDee_upload( &md, "coordinates", par.R );
  EmDee_upload( &md, "momenta", par.V );
  EmDee_random_momenta( &md, par.Temp, 1, par.seed );

  printf("%d %lf %lf\n", 0, md.Potential, md.Virial);
  for (int passo = 1; passo <= par.Npassos; passo++) {
    if (passo % par.Nprop == 0) printf("%d %lf %lf\n", passo, md.Potential, md.Virial);
    EmDee_boost( &md, 1.0, 0.0, par.Dt_2 );
    EmDee_move( &md, 1.0, 0.0, par.Dt );
    EmDee_boost( &md, 1.0, 0.0, par.Dt_2 );
  }
  printf("neighbor list builds = %d\n", md.builds);
  printf("pair time = %f s.\n", md.pairTime);
  printf("excecution time = %f s.\n", md.totalTime);
*/
  return EXIT_SUCCESS;
}

