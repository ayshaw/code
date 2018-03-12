#include <omp.h>
#include <time.h>
#include <stdio.h>
#include <stdlib.h>

#define N 1500                 /* matrix size */

int main (int argc, char *argv[]) 
{
int	i, j, k;
double	a[N][N],           /* matrix A to be multiplied */
	b[N][N],           /* matrix B to be multiplied */
	c[N][N];           /* result matrix C */



  /** Spawn a parallel region explicitly scoping all variables ***/
#pragma omp parallel shared(a,b,c) private(i,j,k)
  {

  #pragma omp for schedule (static)
  for (i=0; i<N; i++)
    for (j=0; j<N; j++)
      a[i][j]= i+j;

  #pragma omp for schedule (static)
  for (i=0; i<N; i++)
    for (j=0; j<N; j++)
      b[i][j]= i*j;

  #pragma omp for schedule (static)
  for (i=0; i<N; i++)
    for (j=0; j<N; j++)
     c[i][j]= 0;

  #pragma omp for schedule (static)
  for (i=0; i<N; i++)    
    {
    for(j=0; j<N; j++)       
      for (k=0; k<N; k++)
        c[i][j] += a[i][k] * b[k][j];
    }
}


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
}

