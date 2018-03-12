#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "timing.h"
#define N 1500                 /* matrix size */
#include "omp.h"
int main (int argc, char *argv[]) 
{
int i, j, k;
timing_t tstart, tend;
double	a[N][N],           /* matrix A to be multiplied */
	b[N][N],           /* matrix B to be multiplied */
	c[N][N];           /* result matrix C */


  for (i=0; i<N; i++)
    for (j=0; j<N; j++)
      a[i][j]= i+j;

  for (i=0; i<N; i++)
    for (j=0; j<N; j++)
      b[i][j]= i*j;

  for (i=0; i<N; i++)
    for (j=0; j<N; j++)
      c[i][j]= 0.0;
 int bsize=25;
 int en = bsize*(N/bsize);
 int kk, jj;
 double sum;

get_time(&tstart);
#pragma omp parallel for collapse(4)
 for (kk=0;kk<en;kk+=bsize){
   for (jj=0;jj<en;jj+=bsize){
       for (i=0;i<N;i++){
	 for (j=jj;j<jj+bsize;j++){
	   sum=c[i][j];
#pragma omp for reduction (+:sum)
	   for (k=kk;k<kk+bsize;k++){
	     sum+=a[i][k]*b[k][j];
	   }
	   c[i][j]=sum;
	 }
     }
   }
 }
get_time(&tend);


   printf("Elapsed time: %g s\n", timespec_diff(tstart, tend));

  printf("*****************************************************\n");

  printf("Result Matrix:\n");
  for (i=0; i<N; i++)
  {
    for (j=0; j<N; j++) 
      printf("%6.2f   ", c[i][j]);
    printf("\n"); 
  }
  printf("******************************************************\n");
  printf ("Done.\n");
printf("Elapsed time: %g s\n", timespec_diff(tstart, tend));
}

