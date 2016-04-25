#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

int FastSN(double *, double *, double *, double *, int, int);
int FastTN(double *, double *, double *, double *, int, int);
int FastUN(double *, double *, double *, double *, int, int);
void dbl_arr_print(double *, int);


static double *sf; /* Pre-calculated sine values */


double 
sf_index(int i, int j, double *S) 
{
    //printf("    - Called sf_index(%d, %d)\n", i, j);
    
    i = i % (2 * j);
    int k = 1;
    if (i > j) {
        k = -1;
        i -= j;
    }
    if (i * 2 > j) { i = j - i; }
    if (i * 2 == j) { return k; }
    if (i == 0) { return 0; }
    while (i % 2 == 0) { i /= 2; j /= 2; }
    
    double res;
    if (j <= 3) {
        res = S[j - 1] * k;
    } else if (j == 6) {
        res = S[3] * k;
    } else if (j % 3 == 0) {
        res = S[j / 4 + i / 2 + 1] * k;
    } else {
        res = S[(j * 3) / 4 + (i * 3) / 2 + 1] * k;    
    }
    //printf("[%dpi/%d] %g \n", i, j, res);
    return res;
}


void 
sf_init(void)
{
    sf = malloc(sizeof(*sf) * 3);
    sf[0] = 1;
    sf[1] = 1;
    sf[2] = sin(M_PI / 3);
}


void 
sf_release(void)
{
    free(sf);
}


double * 
SFactors(int n)
{
    /* Initialize array of sine values for base case. */
    if (!sf) { sf_init(); }
 
    /* Maximum N for which values are pre-calculated */
    static int MAX_N = 3;
 
    /* Calculate new values if n > MAX_N */
    if (n > MAX_N) {
        /* Sine values previously calculated 
           up to index equal to the old MAX_N
           divided by 2. */
        int index = MAX_N == 3 ? 2 : MAX_N / 2;
        
        int NEW_MAX_N = MAX_N;
        /* Lazy evaluation for sine values */
        while (NEW_MAX_N < n) { NEW_MAX_N *= 2; } 
        
        /* Allocate memory for new sine array.
           The memory needed is NEW_MAX_N / 2, 
           with an extra value from 1-indexing. */
        int sf_size = NEW_MAX_N / 2;
        sf = realloc(sf, sizeof(*sf) * (sf_size + 1));
        
        /* Calculate new sine values. */
        for (int i = MAX_N * 2; i <= NEW_MAX_N; i *= 2) {
            int j = 1;
            do {
                //printf("i: %d    j: %d\n", j, i);
                sf[++index] = sin(M_PI * j / i);
                j += 2;
            } while (j < i / 2);            
        }
        
        MAX_N = NEW_MAX_N;
    }
    
    /* Copy sine values into new array and return the results. */
    double *v = malloc(sizeof(*v) * (MAX_N / 2 + 1));
    memcpy(v, sf, sizeof(*v) * (MAX_N / 2 + 1));
    return v;
}


void 
debug_sin(int index)
{
    int k = 2;
    int d = 1;
    while (index > d) {
        d *= 2;
        k *= 2;
    }
    printf("sin(%d*pi/%d)", -1 + 2 * index - d, k);
}


void
test(double *exp, double *got, int N)
{
    const double EPSILON = 1.e-5;
    int j = 1;
    for (int i = 1; i <= N; i++) 
        if (fabs(exp[i] - got[i]) > EPSILON)
            j = 0;
    if (j) {
        printf("Match!\n");
        return;
    }
    printf("Mismatch!\n");
    printf("Expected: ");
    dbl_arr_print(exp, N);
    printf("Got:      ");
    dbl_arr_print(got, N);
    return;
}



int SN(double *x, double *y, int N) 
{
    for (int i = 1; i <= N; i++) {
        x[i] = 0;
        for (int j = 1; j <= N; j++) {
            x[i] += sin(M_PI*i*j/N) * y[j];
        }
    }
    return 1;
}
   

int 
FastSN(double *x, double *y, double *w, double *S, int N, int skip)
{
    //printf("  - Called S(%d)\n", N);
    /* i,j th element of S(N) = sin(i*j*pi/N) */
    if (N == 1) {
    
    } else if (N == 2) {
        /* S(2) = [[ 1 ]] */
        x[skip] = y[skip];
    } else if (N == 3) {
        x[skip] = y[skip] * sf_index(1, 3, S) + y[skip * 2] * sf_index(2, 3, S);
        x[skip * 2] = y[skip] * sf_index(2, 3, S) + y[skip * 2] * sf_index(4, 3, S);  
    } else if (N % 2 == 0 && N > 2) {
        /* Separate S into symmetric and antisymmetric parts */
        for (int j = 1; j < N; j++)
            w[j * skip] = y[j * skip];
        for (int j = 1; j < N / 2; j++)
            y[j * skip] = w[j * skip] - w[(N-j) * skip];
        for (int j = 1; j <= N / 2; j++)
            y[(N / 2 + j - 1) * skip] = w[j * skip] + w[(N-j) * skip];
        y[(N - 1) * skip] = w[N / 2 * skip];
        
        /* Recursive step */
        int a = FastSN(x, y, w, S, N / 2, skip);
        int b = FastTN(x + N / 2 - 1, y + N / 2 - 1, w + N / 2 - 1, S, N / 2, skip);
        
        /* Return -1 if one of the recursive calls failed. */
        if (a || b != 0) { return -1; }
        
        for (int j = 1; j < N; j++)
            y[j * skip] = x[j * skip];

        /* Permutation matrix:
           For odd i, x[i] = y[N/2 + (i - 1)/2]
           For even j, x[j] = y[j/2]      */
        for (int j = 1; j <= N; j += 2)
            x[j * skip] = y[((j - 1)/ 2 + N / 2) * skip];
        for (int j = 2; j <= N; j += 2)
            x[j * skip] = y[j / 2 * skip];
    } else {
        return -1;
    }
    //printf("  - Exit S(%d)\n", N);
    return 0;
}

int TN(double *x, double *y, int N)
{
    for (int i = 1; i <= N; i++) {
        x[i] = 0;
        for (int j = 1; j <= N; j++) {
            x[i] += sin(M_PI*(2*i-1)*j/N/2) * y[j];
        }
    }
    return 1;
}

int 
FastTN(double *x, double *y, double *w, double *S, int N, int skip)
{
    //printf("  - Called T(%d)\n", N);
    /* i,j th element of T(N) = sin((2*i-1)*j*pi/2N) */
    if (N == 1) {
        /* T(1) = [[ 0 ]] */
        x[skip] = 0;  
    } else if (N == 2) {
        /* T(2) = [ [sin(pi/4), sin(2pi/4)], 
                    [sin(3pi/4), sin(6pi/4)] ] */
        x[skip] = y[skip] * sf_index(1, 4, S) + y[skip * 2] * sf_index(2, 4, S);
        x[skip * 2] = y[skip] * sf_index(3, 4, S) - y[skip * 2] * sf_index(2, 4, S);
    } else if (N == 3) {
        x[skip] = y[skip] * sf_index(1, 6, S) + y[skip * 2] * sf_index(2, 6, S) 
                + y[skip * 3] * sf_index(3, 6, S);
        x[skip * 2] = y[skip] * sf_index(3, 6, S)
                    + y[skip * 3] * sf_index(9, 6, S);
        x[skip * 3] = y[skip] * sf_index(5, 6, S) + y[skip * 2] * sf_index(10, 6, S) 
                    + y[skip * 3] * sf_index(15, 6, S);        
    } else if (N % 2 == 0 || N % 3 == 0) {
        /* Permutation matrix:
           Map [1,2,3,4,...] -> [2,4,6,...,1,3,5,...] */
        for (int j = 1; j <= N; j++)
            w[j * skip] = y[j * skip];
        for (int j = 1; j <= N / 2; j++)
            y[j * skip] = w[2 * j * skip];
        for (int j = 1; j <= N / 2; j++)
            y[(j + N / 2) * skip] = w[(2 * j - 1) * skip];

        /* Recursive step */
        int a = FastTN(x, y, w, S, N / 2, skip);
        int b = FastUN(x + N / 2, y + N / 2, w + N / 2, S, N / 2, skip);
        
        /* Return -1 if one of the recursive calls failed. */
        if (a || b != 0) { return -1; }
        
        for (int j = 1; j <= N; j++)
            w[j * skip] = x[j * skip];
        
        /* Find solution from symmetric and antisymmetric parts */
        for (int j = 1; j <= N / 2; j++) {
            x[j * skip] = w[j * skip] + w[(j + N / 2) * skip];
            x[(N + 1 - j) * skip] = - w[j * skip] + w[(j + N / 2) * skip];
        }        
    } else {
        //printf("  - T(%d) not supported!\n", N);
        return -1;
    }
    //printf("  - Exit T(%d)\n", N);
    return 0;
}


int UN(double *x, double *y, int N)
{
    //printf("  - Called U(%d)\n", N);
    for (int i = 1; i <= N; i++) {
        x[i] = 0;
        for (int j = 1; j <= N; j++) {
            x[i] += sin(M_PI*(2*i-1)*(2*j-1)/N/4) * y[j];
        }
    }
    return 1;
}


int 
FastUN(double *x, double *y, double *w, double *S, int N, int skip)
{
    /* i,j th element of U(N) = sin((2*i-1)*(2*j-1)*pi/4N) */
    if (N == 1) {
        /* U(1) = [ [ sin(pi) ] ] = [ [ 0 ] ]  */
        x[skip] = 0;
    } else if (N == 2) {
        /* U(2) = [ [sin(pi/8), sin(3pi/8)], 
                    [sin(3pi/8), sin(9pi/8)] ] */
        x[skip] = y[skip] * sf_index(1, 8, S) + y[skip * 2] * sf_index(3, 8, S);
        x[skip * 2] = y[skip] * sf_index(3, 8, S) - y[skip * 2] * sf_index(1, 8, S);
    } else if (N == 3) {
        /* U(3) = [ [sin(pi/12), sin(3pi/12), sin(5pi/12) ],
                    [sin(3pi/12), sin(9pi/12), sin(15pi/12) ],
                    [sin(5pi/12), sin(15pi/12), sin(25pi/12) ] ] */
        x[skip] = y[skip] * sf_index(1, 12, S) + y[skip * 2] * sf_index(3, 12, S)
                + y[skip * 3] * sf_index(5, 12, S);
        x[skip * 2] = y[skip] * sf_index(3, 12, S) + y[skip * 2] * sf_index(9, 12, S)
                    + y[skip * 3] * sf_index(15, 12, S);
        x[skip * 3] = y[skip] * sf_index(5, 12, S) + y[skip * 2] * sf_index(15, 12, S)
                    + y[skip * 3] * sf_index(25, 12, S);
    } else if (N % 2 == 0 || N % 3 == 0) {
        /* Decompose U(N) = X(N) * T(N/2)T(N/2) * Y(N) */
        
        /* 
            Create u[i] and v[i] from y[i] using Y(N):
            Col  : 1   2   3  ...  N-4 N-3 N-2 N-1   N
            1    : 0,  0,  0, ... ,  0,  0, -1,  1,  0
            2    : 0,  0,  0, ... , -1,  1,  0,  0,  0
                    .................................
            N/2-1: 0, -1,  1, ... ,  0,  0,  0,  0,  0
            N/2  : 1,  0,  0, ... ,  0,  0,  0,  0,  0
            N/2+1: 0,  1,  1, ... ,  0,  0,  0,  0,  0
                    ................................. 
            N-2  : 0,  0,  0, ... ,  1,  1,  0,  0,  0
            N-1  : 0,  0,  0, ... ,  0,  0,  1,  1,  0
            N    : 0,  0,  0, ... ,  0,  0,  0,  0,  1
         */   
        for (int j = 1; j <= N; j++)
            w[j * skip] = y[j * skip];
        for (int j = 1; j < N / 2; j++)
            y[j * skip] = - w[(N / 2 - j) * 2 * skip] + w[((N / 2 - j) * 2 + 1) * skip];
        y[N / 2] = w[skip];
        for (int j = 1; j < N / 2; j++)
            y[(j + N / 2) * skip] = w[j * 2 * skip] + w[(j * 2 + 1) * skip];
        y[N] = w[N * skip];
        
        /* Recursive step */
        int a = FastTN(x, y, w, S, N / 2, skip);
        int b = FastTN(x + N / 2, y + N / 2, w + N / 2, S, N / 2, skip);
       
        //TN(x, y, N / 2);
        //TN(x + N / 2, y + N / 2, N / 2);
        
        /* Return -1 if one of the recursive calls failed. */
        if (a || b != 0) { return -1; }

        for (int j = 1; j <= N; j++)
            y[j * skip] = x[j * skip];
            
        /*
            Matrix X(N):
            sin(*/
           
        for (int j = 1; j <= N; j++)
            w[j * skip] = y[j * skip];
        
        for (int j = 1; j <= N / 2; j++) {
            int k = j % 2 == 0 ? -1 : 1; 
            x[j * skip] = sf_index(-1 + 2 * j, 4 * N, S) * w[j * skip] * k
                        + sf_index(2 * N + 1 - 2 * j, 4 * N, S) * w[(j + N / 2) * skip];
            x[(N + 1 - j) * skip] = sf_index(2 * N + 1 - 2 * j, 4 * N, S) * w[j * skip] * k
                                  - sf_index(-1 + 2 * j, 4 * N, S) * w[(j + N / 2) * skip];
        }
    } else {
        return -1;
    }
    return 0;
}


void dbl_arr_print(double *arr, int N)
{
    for (int i = 1; i <= N; i++)
        printf("%g ", arr[i]);
    printf("\n");
}


void debug(int N, int func)
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
    
    int S = 0;
    if (func == 1) {
        S = 1;
        printf("Testing S(%d): ", N);
        FastSN(x, y, w, s, N, 1);
        //printf("FastSN Complete\n");
        SN(x2, y2, N); 
        //printf("SN Complete\n");
    } else if (func == 2) {    
        printf("Testing T(%d): ", N);
        FastTN(x, y, w, s, N, 1);
        TN(x2, y2, N);     
    } else if (func == 3) {
        printf("Testing U(%d): ", N);
        FastUN(x, y, w, s, N, 1);
        UN(x2, y2, N);     
    } else {
       printf("Unknown function id: %d!\n", func);
       exit(-1);
    }
    //printf("Comparing: ");
    test(x2, x, N - S);
    //dbl_arr_print(x2, N - S);
    //dbl_arr_print(x, N - S);
    
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
    int it = 10;
    int M = 2;
    int N = 3;
    while (it --> 0) {
        for (int i = 1; i <= 3; i++) {
            debug(M, i);
            debug(N, i);    
        }
        M *= 2;
        N *= 2;
    }
     
    // Cleanup for SFactors
    sf_release();
}
