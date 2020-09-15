/* *********************************************************
   ****  atomistic spin dynamics, Heisenberg-Modell ********
   ****  version from april 2020 ********************
   ********************************************************* */


/***************************** Include libraries **************************************************/
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include "input.c"
#include "initial.c"
#include "solver-steps.c"
#include "solver.c"


/* ************************** main program *************************** */
int main() {
  static double sx[SIZE_X][SIZE_Y][SIZE_Z]={0.0}, sy[SIZE_X][SIZE_Y][SIZE_Z]={0.0}, sz[SIZE_X][SIZE_Y][SIZE_Z]={0.0};      /* spin components*/
  int spin = 0;									/* number of spins */
  int find_index[LX*LY*LZ][3]={0};  /* lattice positions */

  /*Initialize the system*/
  #ifdef KUBUS_SC
  spin = Kubus(sx, sy, sz, find_index);             /* initials system, define lattice */
  #endif

  /* variables for integration loop */
  int n, zaehler=0;  								/* counter for number of steps */
  double temp;
  double rauschen;
  double ceta_x[spin], ceta_y[spin], ceta_z[spin]; /* thermal noise */
  int temp_i, i,j,k,l;
  double dt=DT;


  /* variables for output */
  double magx[LZ], magy[LZ], magz[LZ];
  double mx, my,mz;
  char output[30];
  FILE *test;

  const gsl_rng_type *rng_type;
  gsl_rng *rng;                     /* Pointer to rng */

  /* Initialisierung des RNG */
   rng_type = gsl_rng_taus2;  /* Type of generator, default mt19937 */
   rng = gsl_rng_alloc(rng_type);
   gsl_rng_set(rng,SEED);           /* initialize rng with seed */
   temp=sqrt(2.0*ALPH*TEMP*MU_S*dt/GAMMA);

  /* time-integration loop */
  for(n = 1; n <= NMAX; n++){
     for(l = 0; l < spin; l++){
         i = find_index[l][0]; j = find_index[l][1]; k = find_index[l][2];
         ceta_x[l] = gsl_ran_gaussian_ziggurat(rng, temp);
         ceta_y[l] = gsl_ran_gaussian_ziggurat(rng, temp);
         ceta_z[l] = gsl_ran_gaussian_ziggurat(rng, temp);
         }
         /*time evolution*/
         heun(spin, n, find_index, sx,sy,sz, ceta_x, ceta_y, ceta_z);
         if(n%100000==0){
           mx=my=mz=0.0;
           for(k=RAND_PLUS;k<LZ+RAND_PLUS;k++){
             magx[k-RAND_PLUS]=magy[k-RAND_PLUS]=magz[k-RAND_PLUS]=0.0;
             for(i=RAND_PLUS;i<LX+RAND_PLUS;i++){
               for(j=RAND_PLUS;j<LY+RAND_PLUS;j++){
                 magx[k-RAND_PLUS]+=sx[i][j][k];
                 magy[k-RAND_PLUS]+=sy[i][j][k];
                 magz[k-RAND_PLUS]+=sz[i][j][k];
               }
             }
             mx+=magx[k];
             my+=magy[k];
             mz+=magz[k];
             magx[k-RAND_PLUS]=magx[k-RAND_PLUS]/(LX*LY);
             magy[k-RAND_PLUS]=magy[k-RAND_PLUS]/(LX*LY);
             magz[k-RAND_PLUS]=magz[k-RAND_PLUS]/(LX*LY);
           }
           if(n%100==0){
             test=fopen("mag.dat","a");
             fprintf(test, "%d \t %le \t %le \t %le \n", n, mx/LZ, my/LZ, mz/LZ);
             fclose(test);
           }
           if(n%100000==0){
           sprintf(output,"mag_%d.dat",n );
           test=fopen(output,"a");
           for(k=0;k<LZ;k++){
             fprintf(test, "%d \t %le \t %le \t %le \n", k, magx[k], magy[k], magz[k]);
           }
           fclose(test);
         }
     }
  exit(0);
}
