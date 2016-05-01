#include "debug_utils.h"


void
test(double *exp, double *got, int N)
{
    const double EPSILON = 1.e-5;
    int j = 1;
    int inc = 0;
    for (int i = 1; i <= N; i++) {
        if (fabs(exp[i] - got[i]) > EPSILON) {
            inc++;
            if (j == 1) { printf("Mismatch!\nIndex    Expected    Got         Difference\n"); }
            j = 0;
            printf("%7d  %10g  %10g  %10g\n", i, exp[i], got[i], fabs(exp[i] - got[i]));
        }
    }
    if (j) {
        printf("Match!\n");
    } else {
        printf("        %d/%d incorrect (%5f%%)\n", inc, N, 100.f * inc / N);
    }
}


void 
dbl_arr_print(double *arr, int N)
{
    dbl_arr_print2(arr, N, 1);
}


void 
dbl_arr_print2(double *arr, int N, int skip)
{
    for (int i = 1; i <= N; i++)
        printf("%g ", arr[i * skip]);
    printf("\n");
}
