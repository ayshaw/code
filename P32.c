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
	b[N][N],/* matrix B to be multiplied */
    tmp,
	c[N][N];           /* result matrix C */
#pragma acc data copy(a)
#pragma acc kernels
  for (i=0; i<N; i++)
   
    for (j=0; j<N; j++)
      a[i][j]= i+j;
#pragma acc data copy(b)
#pragma acc kernels
  for (i=0; i<N; i++)
  
    for (j=0; j<N; j++)
      b[i][j]= i*j;
#pragma acc data copy(c)
#pragma acc kernels
  for (i=0; i<N; i++)
   
    for (j=0; j<N; j++)
      c[i][j]= 0;

get_time(&tstart);

#pragma acc data copyin(a,b) copy(c)
#pragma acc kernels loop independent vector(32)
  for (i=0; i<N; i++)    {

      for(j=0; j<N; j++) {
          tmp=0.0f;
#pragma acc loop reduction(+:tmp)
          for (k=0; k<N; k++){
        tmp += a[i][k] * b[k][j];
          }
          c[i][j]=tmp;
    }
  }
get_time(&tend);


/*

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
 */
printf("Elapsed time: %g s\n", timespec_diff(tstart, tend));
}

