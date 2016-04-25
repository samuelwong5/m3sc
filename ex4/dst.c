#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

static double *sf; /* Pre-calculated sine values */

int FastSN(double *, double *, double *, double *, int, int);
int FastTN(double *, double *, double *, double *, int, int);
int FastUN(double *, double *, double *, double *, int, int);
void dbl_arr_print(double *, int);

void 
release(void)
{
    free(sf);
}

void
test(double *exp, double *got, int N)
{
    const double EPSILON = 1.e-5;
    int i = 1;
    for (int i = 1; i <= N; i++) 
        if (fabs(exp[i] - got[i]) > EPSILON)
            i = 0;
    if (i) {
        printf("Match!\n");
        return;
    }
    printf("Mismatch!\n");
    printf("Expected:\n");
    dbl_arr_print(exp, N);
    printf("Got:\n");
    dbl_arr_print(got, N);
}

double * 
SFactors(int n)
{
    /* Maximum N for which values are pre-calculated */
    static int MAX_N = 1; 
    
    /* Calculate new values if n > MAX_N */
    if (n > MAX_N) {
        sf = realloc(sf, sizeof(double) * (n / 2 + 1));
        int index = 1;
        for (int i = MAX_N + 1; i <= n; i *= 2) {
            int j = 1;
            do {
                sf[index++] = sin(M_PI * j / i);
                //printf("i: %d   j: %d\n", j, i);
                j += 2;
                
            } while (j < i / 2);            
        }
    }
    
    /* Copy sine values and return */
    double *v = malloc(sizeof(double) * (n / 2 + 1));
    memcpy(v, sf, sizeof(double) * (n / 2 + 1));
    return v;
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
    /* i,j th element of S(N) = sin(i*j*pi/N) */
    if (N == 1) {
        return 0;   
    } else if (N == 2) {
        /* S(2) = [[ 1 ]] */
        x[skip] = y[skip];
        return 0;
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
        return 0;
    } else {
        return -1;
    }
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
    /* i,j th element of T(N) = sin((2*i-1)*j*pi/2N) */
    if (N == 1) {
        /* T(1) = [[ 0 ]] */
        x[skip] = 0;  
    } else if (N == 2) {
        /* T(2) = [ [sin(pi/4), sin(2pi/4)], 
                    [sin(3pi/4), sin(6pi/4)] ] */
        x[skip] = y[skip] * S[2] + y[skip * 2] * S[1];
        x[skip * 2] = y[skip] * S[2] - y[skip * 2] * S[1];
    } else if (N % 2 == 0) {
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
        return -1;
    }
    return 0;
}


int UN(double *x, double *y, int N)
{
    for (int i = 1; i <= N; i++) {
        x[i] = 0;
        for (int j = 1; j <= N; j++) {
            x[i] += sin(M_PI*(2*i-1)*(2*j-1)/N/4) * y[j];
        }
    }
    return 1;
}

void debug_sin(int index)
{
    int k = 2;
    int d = 1;
    while (index > d) {
        d *= 2;
        k *= 2;
    }
    printf("sin(%d*pi/%d)", -1 + 2 * index - d, k);
}

int 
FastUN(double *x, double *y, double *w, double *S, int N, int skip)
{
    /* i,j th element of U(N) = sin((2*i-1)*(2*j-1)*pi/4N) */
    if (N == 1) {
        /* U(1) = [[ sin(pi) ]] = [[ 0 ]]  */
        x[skip] = 0;
    } else if (N == 2) {
        /* U(2) = [ [sin(pi/8), sin(3pi/8)], 
                    [sin(3pi/8), sin(9pi/8)] ] */
        x[skip] = y[skip] * S[3] + y[skip * 2] * S[4];
        x[skip * 2] = y[skip] * S[4] - y[skip * 2] * S[3];
    } else if (N % 2 == 0) {
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
        
        // Finding index of precomputed sine array
        int base_index =  N;
        
        for (int j = 1; j <= N; j++)
            w[j * skip] = y[j * skip];
        
        for (int j = 1; j <= N / 2; j++) {
            int k = j % 2 == 0 ? -1 : 1;    
            /*printf("%d * ", k);
            debug_sin(base_index + j);    
            printf(" + ");
            debug_sin(base_index + N + 1 - j);
            printf("\n");*/
            x[j * skip] = S[base_index + j] * w[j * skip] * k
                           + S[base_index + N + 1 - j] * w[(j + N / 2) * skip];
            x[(N + 1 - j) * skip] = S[base_index + N + 1 - j] * w[j * skip] * k
                                - S[base_index + j] * w[(j + N / 2) * skip];
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
    double *s = SFactors(4*N);
    
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
        SN(x2, y2, N); 
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
    int it = 10;
    int N = 2;
    while (it --> 0) {
        for (int i = 1; i <= 3; i++)
            debug(N, i);
        N *= 2;
    }
     
    // Cleanup for SFactors
    release();
}
