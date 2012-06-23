/*
 * Evan Lezar
 * 18 November 2010
 *
 * Common interface for various routines driving the ARPACK eigenvalue solver
 */
#ifndef GEV_INTERFACE_H
#define GEV_INTERFACE_H
#include <math.h>
#include <sys/time.h>

#ifdef __cplusplus
    #define EXT extern "C"
#else
    #define EXT extern
#endif

//#define DEBUG_OUTPUT

EXT
// definition of the ARPACK routine
void arpack_ssev ( int N, void* DATA, int NEV, int NCV, float* eigenvalues, float* eigenvectors, float* residuals, char* which );

EXT
int dense_gev ( int N, float* A, float* B, int LDMAT,
        int NEV, float shift, float* eigenvalues, float* eigenvectors,
        double* timing_data_10, int* int_data_10 );

EXT
int dense_propagation ( int N, float* A, float* B, int LDMAT,
        int NEV, float shift, float* eigenvalues, float* eigenvectors,
        double* timing_data_10, int* int_data_10, int LR_or_SR );

struct checkpoint {
    int sec;
    int usec;
};
typedef struct checkpoint checkpoint;

EXT
checkpoint tic ( )
{
    struct timeval t;
    t.tv_sec = 0;
    t.tv_usec = 0;
    gettimeofday(&t,NULL);

    checkpoint s_now;
    s_now.sec = (int)t.tv_sec;
    s_now.usec = (int)t.tv_usec;

    return s_now;
}

EXT
// get the elapsed time in seconds
double toc ( checkpoint t_start )
{
    checkpoint elapsed = tic();

    elapsed.sec -= t_start.sec;
    elapsed.usec -= t_start.usec;

    return (double)elapsed.sec + (double)elapsed.usec*1.0e-6;
}

#endif // #ifndef GEV_INTERFACE_H
