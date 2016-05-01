#include "dst.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif


#ifdef _STATIC_SF
static double *sf; /* Static pre-calculated sine values */


/*  ---------------------------------------------
 *  sf_init:
 *  ---------------------------------------------
 *  Initializes the basic precomputed sine array.
 */
void 
sf_init(int base)
{
    if (base == 5) {
        sf = malloc(sizeof(*sf) * 6);
        sf[0] = 1.f;
        sf[1] = 1.f;
        sf[2] = sin(M_PI/5.f);
        sf[3] = sin(2.f*M_PI/5.f);
        sf[4] = sin(M_PI/10.f);
        sf[5] = sin(3.f*M_PI/10.f);
    } else if (base == 3) {
        sf = malloc(sizeof(*sf) * 4);
        sf[0] = 1;
        sf[1] = 1;
        sf[2] = sin(M_PI / 3);
        sf[3] = sin(M_PI / 6);
    } else if (base == 2) {
        sf = malloc(sizeof(*sf) * 2);
        sf[0] = 1.f;
        sf[1] = 1.f;
    } else {
        printf("Unsupported base: %d!\n", base);
        exit(-1);
    }
}


/*  ---------------------------------------------
 *  sf_release:
 *  ---------------------------------------------
 *  Frees the precomputed sine array.
 */
void 
sf_release(void)
{
    free(sf);
}
#endif


/*  ---------------------------------------------
 *  sf_index:
 *  ---------------------------------------------
 *  Finds the corresponding sine value of 
 *  sin(i*PI/j) from the precomputed arrays.
 *
 *  Parameters:
 *    int i     - numerator of fraction 
 *    int j     - denominator of fraction
 *    double *S - array of sine values
 *
 *  Returns:
 *    double    - the value of sin(i*PI/j)
 */
double 
sf_index(int i, int j, double *S) 
{
    int i2 = i, j2 = j;
    if (i == 0 || i == j) { return 0; }
    i = i % (2 * j);
    int k = 1;
    if (i > j) {
        k = -1;
        i -= j;
    }
    if (i * 2 > j) { i = j - i; }
    if (i * 2 == j) { return k; }
    
    while (i % 2 == 0 && j % 2 == 0) { i /= 2; j /= 2; }
    
    double res;
    int index;
    if (j % 5 == 0 && j <= 10) {
        index = j / 3 + i / 2 + 1;
        res = S[j / 3 + i / 2 + 1];
    } else if (j % 3 == 0 && j <= 6) {
        index = j / 3 + 1;
        res = S[j / 3 + 1];
    } else {
        index = j / 4 + i / 2 + 1;
        res = S[j / 4 + i / 2 + 1];
    }
    //printf("index: [%d]\n", index);
    //printf("\n%d %d -> %d %d {%d} = %g: ", i2, j2, i, j, index, res * k); // = {%d} %g\n", 
    //cnt++, i2, j2, i, j, index, res * k);
    return res * k;
}


/*  ---------------------------------------------
 *  SFactors:
 *  ---------------------------------------------
 *  Returns an array of precomputed sine values. 
 *  The sine values are computed lazily - i.e. 
 *  the computation happens when it is explicitly 
 *  required. 
 *
 *  Parameters:
 *    int n    - the denominator needed
 *
 *  Returns:
 *    double * - the array of sine values
 */
double * 
SFactors(int n)
{
    /* Sine value array to be returned. */
    double *v = malloc(sizeof(*v) * (n / 2 + 1));
    
#ifdef _STATIC_SF
    /* Maximum N for which values are pre-calculated */
    static int MAX_N = 0;
    static int BASE_N = 0;
    
    /* Check if new base N. */
    if (n % 5 == 0 && BASE_N != 5) {
        sf_release();
        sf_init(5);
        MAX_N = 10;
        BASE_N = 5;
    } else if (n % 3 == 0 && BASE_N != 3) {
        sf_release();
        sf_init(3);
        MAX_N = 6;
        BASE_N = 3;
    } else if (n % 5 != 0 && n % 3 != 0 && BASE_N != 2) {
        sf_release();
        sf_init(2);
        MAX_N = 2;
        BASE_N = 1;
    }
    
    /* Calculate new values if n > MAX_N */
    if (n > MAX_N) {
        /* Sine values previously calculated 
           up to index equal to the old MAX_N
           divided by 2. */
        int index = MAX_N / 2;
        
        /* Lazy evaluation for sine values. */
        int NEW_MAX_N = n;
        
        /* Allocate memory for new sine array.
           The memory needed is NEW_MAX_N / 2, 
           with an extra value for 1-indexing. */
        int sf_size = NEW_MAX_N / 2;
        sf = realloc(sf, sizeof(*sf) * (sf_size + 1));
        
        /* Calculate new sine values. */
        for (int i = MAX_N * 2; i <= NEW_MAX_N; i *= 2) {
            int j = 1;
            do {
                sf[++index] = sin(M_PI * j / i);
                //printf("[%d] %d %d = %g\n",index, j, i, sf[index]);
                j += 2;   
            } while (j < i / 2);            
        }
        
        MAX_N = NEW_MAX_N;
    }
   
    /* Copy sine values into new array and return the results. */
    memcpy(v, sf, sizeof(*v) * (MAX_N / 2 + 1));
#else
    int index, i;
    v[0] = 1.f;
    /* Check if new base N. */
    if (n % 5 == 0) {
        v[0] = 1.f;
        v[1] = 1.f;
        v[2] = sin(M_PI/5.f);
        v[3] = sin(2.f*M_PI/5.f);
        v[4] = sin(M_PI/10.f);
        v[5] = sin(3.f*M_PI/10.f);
        i = 20;
        index = 5;
    } else if (n % 3 == 0) {
        v[0] = 1;
        v[1] = 1;
        v[2] = sin(M_PI / 3);
        v[3] = sin(M_PI / 6);
        i = 12;
        index = 3;
    } else {
        v[0] = 1.f;
        v[1] = 1.f;
        i = 4;
        index = 1;
    }
    
    /* Calculate new sine values. */
    for (; i <= n; i *= 2) {
        int j = 1;
        do {
            v[++index] = sin(M_PI * j / i);
            j += 2;
        } while (j < i / 2);            
    }
#endif  
    return v;
}



/*  ---------------------------------------------
 *  FastSN:
 *  ---------------------------------------------
 *  Computes the matrix-vector multiply:
 *                x = S(N) * y
 *  where x and y are vectors of size N-1 and 
 *  S(N) is a matrix of size N-1 by N-1 whose 
 *  j, kth element is sin(jk*PI/N).
 *
 *  Parameters:
 *    double *x - the result of the matrix-vector
 *                multiplication
 *    double *y - the input vector
 *    double *w - an array to perform intermediate
 *                calculations
 *    double *S - an array of precomputed sine values
 *    int N     - the size of the matrix and vectors
 *    int skip  - the inputs and outputs should be
 *                in y[i*skip] and x[i*skip] for 
 *                1 <= i < N
 *
 *  Returns:
 *    int       - 0 if function implemented for N, 
 *                -1 otherwise
 */
int 
FastSN(double *x, double *y, double *w, double *S, int N, int skip)
{
    //printf("  - S(%d) called.\n", N);
    if (N == 1) {
        /* Nothing to do here. */
    } else if (N == 2) {
        /* S(2) = [[ 1 ]] */
        x[skip] = y[skip];
    } else if (N == 3) {
        x[skip] = y[skip] * sf_index(1, 3, S) + y[skip * 2] * sf_index(2, 3, S);
        x[skip * 2] = y[skip] * sf_index(2, 3, S) + y[skip * 2] * sf_index(4, 3, S);  
    } else if (N == 5) {
        x[skip] = y[skip] * sf_index(1, 5, S) + y[skip * 2] * sf_index(2, 5, S)
                + y[skip * 3] * sf_index(2, 5, S) + y[skip * 4] * sf_index(1, 5, S);
        x[skip * 2] = y[skip] * sf_index(2, 5, S) + y[skip * 2] * sf_index(4, 5, S)
                    + y[skip * 3] * sf_index(6, 5, S) + y[skip * 4] * sf_index(8, 5, S);
        x[skip * 3] = y[skip] * sf_index(3, 5, S) + y[skip * 2] * sf_index(6, 5, S)
                    + y[skip * 3] * sf_index(9, 5, S) + y[skip * 4] * sf_index(2, 5, S);
        x[skip * 4] = y[skip] * sf_index(1, 5, S) + y[skip * 2] * sf_index(8, 5, S)
                    + y[skip * 3] * sf_index(2, 5, S) + y[skip * 4] * sf_index(6, 5, S);
    } else if (N % 2 == 0 && N > 2) {
        /* Separate S into symmetric and antisymmetric parts */
        for (int j = 1; j < N; j++)
            w[j * skip] = y[j * skip];
        for (int j = 1; j < N / 2; j++)
            y[2 * j * skip] = w[j * skip] - w[(N-j) * skip];
        for (int j = 1; j < N / 2; j++)
            y[(2 * j - 1) * skip] = w[j * skip] + w[(N-j) * skip];
        y[(N - 1) * skip] = w[N / 2 * skip];
        /* Recursive step */
        int a = FastSN(x, y, w, S, N / 2, 2 * skip);
        int b = FastTN(x - skip, y - skip, w - skip, S, N / 2, 2 * skip);
        /* Return -1 if one of the recursive calls failed. */
        if (a || b != 0) { return -1; }
    } else {
        return -1;
    }
    //printf("  - S(%d) exited.\n", N);
    return 0;
}


/*  ---------------------------------------------
 *  FastTN:
 *  ---------------------------------------------
 *  Computes the matrix-vector multiply:
 *                x = T(N) * y
 *  where x and y are vectors of size N and 
 *  T(N) is a matrix of size N by N whose 
 *  j, kth element is sin(j*(2*k-1)*PI/2N).
 *
 *  Parameters:
 *    double *x - the result of the matrix-vector
 *                multiplication
 *    double *y - the input vector
 *    double *w - an array to perform intermediate
 *                calculations
 *    double *S - an array of precomputed sine values
 *    int N     - the size of the matrix and vectors
 *    int skip  - the inputs and outputs should be
 *                in y[i*skip] and x[i*skip] for 
 *                1 <= i <= N
 *
 *  Returns:
 *    int       - 0 if function implemented for N, 
 *                -1 otherwise
 */
int 
FastTN(double *x, double *y, double *w, double *S, int N, int skip)
{
    //printf("  - T(%d) called.\n", N);
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
    } else if (N == 5) {
        x[skip] = y[skip] * sf_index(1, 10, S) + y[skip * 2] * sf_index(2, 10, S)
                + y[skip * 3] * sf_index(3, 10, S) + y[skip * 4] * sf_index(4, 10, S)
                + y[skip * 5] * sf_index(1, 2, S);
        x[skip * 2] = y[skip] * sf_index(3, 10, S) + y[skip * 2] * sf_index(2, 5, S)
                    + y[skip * 3] * sf_index(1, 10, S) + y[skip * 4] * sf_index(6, 5, S)
                    + y[skip * 5] * sf_index(3, 2, S);
        x[skip * 3] = y[skip] * sf_index(1, 2, S) + y[skip * 2] * sf_index(1, 1, S)
                    + y[skip * 3] * sf_index(3, 2, S) + y[skip * 5] * sf_index(1, 2, S);
        x[skip * 4] = y[skip] * sf_index(3, 10, S) + y[skip * 2] * sf_index(7, 5, S)
                    + y[skip * 3] * sf_index(1, 10, S) + y[skip * 4] * sf_index(1, 5, S)
                    + y[skip * 5] * sf_index(3, 2, S);
        x[skip * 5] = y[skip] * sf_index(1, 10, S) + y[skip * 2] * sf_index(9, 5, S)
                    + y[skip * 3] * sf_index(3, 10, S) + y[skip * 4] * sf_index(8, 5, S)
                    + y[skip * 5] * sf_index(1, 2, S);
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
        int b = FastUN(x + skip * N / 2, y + skip * N / 2, w + skip * N / 2, S, N / 2, skip);
        
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
    //printf("  - T(%d) exited.\n", N);
    return 0;
}


/*  ---------------------------------------------
 *  FastUN:
 *  ---------------------------------------------
 *  Computes the matrix-vector multiply:
 *                x = U(N) * y
 *  where x and y are vectors of size N and 
 *  U(N) is a matrix of size N by N whose 
 *  j, kth element is sin((2*j-1)*(2*k-1)*PI/4N).
 *
 *  Parameters:
 *    double *x - the result of the matrix-vector
 *                multiplication
 *    double *y - the input vector
 *    double *w - an array to perform intermediate
 *                calculations
 *    double *S - an array of precomputed sine values
 *    int N     - the size of the matrix and vectors
 *    int skip  - the inputs and outputs should be
 *                in y[i*skip] and x[i*skip] for 
 *                1 <= i <= N
 *
 *  Returns:
 *    int       - 0 if function implemented for N, 
 *                -1 otherwise
 */
int 
FastUN(double *x, double *y, double *w, double *S, int N, int skip)
{
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
    } else if (N == 5) {
        x[skip] = y[skip] * sf_index(1, 20, S) + y[skip * 2] * sf_index(3, 20, S)
                + y[skip * 3] * sf_index(5, 20, S) + y[skip * 4] * sf_index(7, 20, S)
                + y[skip * 5] * sf_index(9, 20, S);
        x[skip * 2] = y[skip] * sf_index(3, 20, S) + y[skip * 2] * sf_index(9, 20, S)
                    + y[skip * 3] * sf_index(15, 20, S) + y[skip * 4] * sf_index(21, 20, S)
                    + y[skip * 5] * sf_index(27, 20, S);
        x[skip * 3] = y[skip] * sf_index(5, 20, S) + y[skip * 2] * sf_index(15, 20, S)
                    + y[skip * 3] * sf_index(25, 20, S) + y[skip * 4] * sf_index(35, 20, S)
                    + y[skip * 5] * sf_index(5, 20, S);
        x[skip * 4] = y[skip] * sf_index(7, 20, S) + y[skip * 2] * sf_index(21, 20, S)
                    + y[skip * 3] * sf_index(35, 20, S) + y[skip * 4] * sf_index(9, 20, S)
                    + y[skip * 5] * sf_index(23, 20, S);
        x[skip * 5] = y[skip] * sf_index(9, 20, S) + y[skip * 2] * sf_index(27, 20, S)
                    + y[skip * 3] * sf_index(5, 20, S) + y[skip * 4] * sf_index(23, 20, S)
                    + y[skip * 5] * sf_index(1, 20, S);
    } else if (N % 2 == 0 && N > 2) {
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
        y[N / 2 * skip] = w[skip];
        for (int j = 1; j < N / 2; j++)
            y[(j + N / 2) * skip] = w[j * 2 * skip] + w[(j * 2 + 1) * skip];
        y[N * skip] = w[N * skip];
        
        /* Recursive step */
        int a = FastTN(x, y, w, S, N / 2, skip);
        int b = FastTN(x + skip * N / 2, y + skip * N / 2, w + skip * N / 2, S, N / 2, skip);
       
        /* Return -1 if one of the recursive calls failed. */
        if (a || b != 0) { return -1; }

        /* Compute results from decomposed function calls. */          
        for (int j = 1; j <= N; j++)
            w[j * skip] = x[j * skip];
        
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



/******************************************
 * DIRECT MATRIX-VECTOR MULTIPLY VERSIONS *
 ******************************************/
int 
SN(double *x, double *y, int N, double *S) 
{
    for (int i = 1; i <= N; i++) {
        x[i] = 0.f;
        for (int j = 1; j <= N; j++) {
            uintmax_t ui = (uintmax_t) i;
            uintmax_t uj = (uintmax_t) j;
            uintmax_t unum = ((ui%(2*N))*(uj%(2*N)))%(2*N);
            int num = (int) unum;
            x[i] += sf_index(num,N,S) * y[j];
        }
    }
    return 1;
}


int 
TN(double *x, double *y, int N, double *S)
{
    for (int i = 1; i <= N; i++) {
        x[i] = 0;
        for (int j = 1; j <= N; j++) {
            
            uintmax_t ui = (uintmax_t) 2 * i - 1;
            uintmax_t uj = (uintmax_t) j;
            int denom = 2*N;
            uintmax_t unum = ((ui%(denom*2))*(uj%(denom*2)))%(denom*2);
            int num = (int) unum;
            double val = sf_index(num,denom,S);
            //printf("i: %llu j: %llu num: %d dem: %d  = %g\n", ui, uj, num, 2*N, val);
            x[i] += val * y[j];
        }
    }
    return 1;
}


int 
UN(double *x, double *y, int N, double *S)
{
    for (int i = 1; i <= N; i++) {
        x[i] = 0;
        for (int j = 1; j <= N; j++) {
            uintmax_t ui = (uintmax_t) 2 * i - 1;
            uintmax_t uj = (uintmax_t) 2 * j - 1;
            int denom = 4 * N;
            uintmax_t unum = ((ui%(denom*2))*(uj%(denom*2)))%(denom*2);
            int num = (int) unum;
            double val = sf_index(num,denom,S);
            
            x[i] += val * y[j];
        }
    }
    return 1;
}
