#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "timing.h"
#define N 1500                /* matrix size */
#define block 25

int main (int argc, char *argv[])
{
    int    i, j, k, jj, kk, c00, c01, c10, c11;
    int n = block*(N/block);
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
    
    for (kk = 0; kk < n; kk += block) {
        for (jj = 0; jj < n; jj += block) {
            for (i = 0; i < N; i+=2) {
                for (j = jj; j < jj + block; j+=2) {
                    c00 = c01 = c10 = c11 = 0;
                    for (k=kk; k<kk+block; k++) {
                        c00 += a[i][k] * b[k][j];
                        c01 += a[i][k] * b[k][j+1];
                        c10 += a[i+1][k] * b[k][j];
                        c11 += a[i+1][k] * b[k][j+1];
                    }
                    c[i][j] += c00;
                    c[i][j+1] += c01;
                    c[i+1][j] += c10;
                    c[i+1][j+1] += c11;
                }
            }
        }
    }
    
    get_time(&tend);
    
    printf("Elapsed time: %g s\n", timespec_diff(tstart, tend));
}

