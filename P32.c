#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "timing.h"
#include "openacc.h"
#define N 1500                 /* matrix size */

int main (int argc, char *argv[]) 
{
int	i, j, k;
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
      c[i][j]= 0;

get_time(&tstart);
#pragma acc data copyin(a,b,c) copyout(c)
    {
#pragma acc parallel loop
  for (i=0; i<N; i++)    
    {
    for(j=0; j<N; j++)
      for (k=0; k<N; k++)
        c[i][j] += a[i][k] * b[k][j];
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

