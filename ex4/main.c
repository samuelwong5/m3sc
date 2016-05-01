#include <time.h>

#include "dst.h"
#include "debug_utils.h"

#define _STATIC_SF

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

void
sn_time(int N)
{
    double *x = malloc(sizeof(*x)*(N+1));
    double *y = malloc(sizeof(*x)*(N+1));
    double *w = malloc(sizeof(*x)*(N+1));
    double *s = SFactors(8*N);
    
    double *y2 = malloc(sizeof(*y)*(N+1));
    double *x2 = malloc(sizeof(*x)*(N+1));
    
    // Work
    static int ITERATIONS = 10000;

    for (int i = 1; i <= N; i++) {
        y[i] = (double) i;
        y2[i] = (double) i;
    }
    
    /* S flag is set to indicate an 
       offset for the matrix size. */
    clock_t fsn_start = clock();
    //for (int i = 0; i < ITERATIONS; i++)
    FastSN(x, y, w, s, N, 1);
    clock_t fsn_end = clock();
    //for (int i = 0; i < ITERATIONS; i++)
    //    SN(x2, y2, N); 
    clock_t sn_end = clock();
    
    printf("FastSN(%d): %fs\n", N, (double)(fsn_end - fsn_start) / CLOCKS_PER_SEC);
    printf("SN(%d)    : %fs\n", N, (double)(sn_end - fsn_end) / CLOCKS_PER_SEC);
    //test(x2, x, N - 1);
    
    // Cleanup
    free(x);
    free(y);
    free(w);
    free(s);
    
    free(x2);
    free(y2);
}




int poisson(int N)
{
    // Init
    double *x = malloc(sizeof(*x)*(N));
    double *y = malloc(sizeof(*x)*(N));
    double *w = malloc(sizeof(*x)*(N));
    double *s = SFactors(4*N);
    
    double max = -1.f;
    int max_x, max_y;
    
    // Work
    clock_t start = clock();
    int max_i = 0, max_j = 0; 
    double max_p = -1.f;
    for (int j = 1; j < N; j++) {
        for (int i = 1; i < N; i++) {
            if (j >= N / 4 && j <= N / 2 && i >= N / 4 && i <= N / 2) { 
                double psi = 100.f;
                psi /= (N * N);
                if (j == N/4 || j == N / 2)
                    psi *= .5f;
                if (i == N/4 || i == N / 2)
                    psi *= .5f;
                psi /= (4 - 2 * cos(i*M_PI/N) - 2 * cos(j*M_PI/N));
                if (psi > max_p) {
                    max_i = i;
                    max_j = j;
                    max_p = psi;
                }
            } 
        }
    }
    /*
    for (int j = 1; j < N; j++) {
        for (int i = 1; i < N; i++) {
            if (j >= N / 4 && j <= N / 2 && i >= N / 4 && i <= N / 2) 
                y[i] = 100.f;
            else
                y[i] = 0.f;
        }
        // Transform p_ij into p^ij
        FastSN(x, y, w, s, N, 1);
        for (int i = 1; i < N; i++)
            y[i] = x[i];
        FastSN(x, y, w, s, N, 1);
        for (int i = 1; i < N; i++)
            y[i] = x[i] * 4.f / N / N;
        // Transform psi^kl into psi_kl
        FastSN(x, y, w, s, N, 1);
        for (int i = 1; i < N; i++)
            y[i] = x[i];
        FastSN(x, y, w, s, N, 1);
        for (int i = 1; i < N; i++) {
            x[i] = x[i] / (M_PI * M_PI * (i * i + j * j));
            if (x[i] > max) {
                max = x[i];
                max_x = j;
                max_y = i;
            }
        }
    }*/
    clock_t end = clock();
    printf("%d,%d,%d,%g,%g\n", N, max_i, max_j, max_p, (double)(end - start) / CLOCKS_PER_SEC);
   
    // Cleanup
    free(x);
    free(y);
    free(w);
    free(s);
}



void 
debug(int N, int func)
{
    // Init
    double *x = malloc(sizeof(*x)*(N+1));
    double *y = malloc(sizeof(*x)*(N+1));
    double *w = malloc(sizeof(*x)*(N+1));
    double *s = SFactors(8*N);
    
    double *y2 = malloc(sizeof(*y)*(N+1));
    double *x2 = malloc(sizeof(*x)*(N+1));
    
    // Work
    for (int i = 1; i <= N; i++) {
        y[i] = (double) i;
        y2[i] = (double) i;
    }
    
    /* S flag is set to indicate an 
       offset for the matrix size. */
    int S = 0;
    
    if (func == 1) {
        S = 1;
        printf("Testing S(%d): ", N);
        FastSN(x, y, w, s, N, 1);
        SN(x2, y2, N, s); 
    } else if (func == 2) {    
        printf("Testing T(%d): ", N);
        FastTN(x, y, w, s, N, 1);
        TN(x2, y2, N, s);     
    } else if (func == 3) {
        printf("Testing U(%d): ", N);
        FastUN(x, y, w, s, N, 1);
        UN(x2, y2, N, s);     
    } else {
       printf("Unknown function id: %d!\n", func);
       exit(-1);
    }
    test(x2, x, N - S);
    
    // Cleanup
    free(x);
    free(y);
    free(w);
    free(s);
    
    free(x2);
    free(y2);
}


int main()
{
#ifdef _STATIC_SF
    printf("Using static version of SFactors...\n");
#endif

    int it = 20; //50;
    int M = 2;
    int N = 3;
    int O = 5;
    while (it --> 0) {
        //poisson(M);
        for (int i = 1; i <= 3; i++) {
            debug(M, i);
            debug(N, i);  
            debug(O, i);              
        }
        //sn_time(O);
        // }
        M *= 2;
        N *= 2; O *= 2;
    }
     
#ifdef _STATIC_SF
    // Cleanup for SFactors
    sf_release();
#endif
}