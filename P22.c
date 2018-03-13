#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "timing.h"
#define N 1500                /* matrix size */
#define block 25

int main (int argc, char *argv[])
{
    int    i, j, k, jj, kk, c1, c2, c3, c4;
    int en = block*(N/block);
    timing_t tstart, tend;
    double    a[N][N],           /* matrix A to be multiplied */
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
    
    for (kk = 0; kk < en; kk += block) {
        for (jj = 0; jj < en; jj += block) {
            for (i = 0; i < N; i+=2) {
                for (j = jj; j < jj + block; j+=2) {
                    c1 = c2 = c3 = c4 = 0;
                    for (k=kk; k<kk+block; k++) {
                        c1 += a[i][k] * b[k][j];
                        c2 += a[i][k] * b[k][j+1];
                        c3 += a[i+1][k] * b[k][j];
                        c4 += a[i+1][k] * b[k][j+1];
                    }
                    c[i][j] += c1;
                    c[i][j+1] += c2;
                    c[i+1][j] += c3;
                    c[i+1][j+1] += c4;
                }
            }
        }
    }
    printf("%g",c[0][0]);
    get_time(&tend);
    
    printf("Elapsed time: %g s\n", timespec_diff(tstart, tend));
}

