/*
 * Evan Lezar
 * 18 November 2010
 *
 * Implementation of dense generalized eigenvalue solver using arpack
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <omp.h>
#include <acml.h>

#include "interface.h"
#include "common_functions.h"

#define ARPACK_CALL 0
#define SHIFT 1
#define SGETRF_S 2
#define SGEMV 9

struct data_struct {
    float* A;
    double sgemv_time;
    int sgemv_calls;
    int LDA;
};
typedef struct data_struct data_struct;

void print_data ( const char* desc, data_struct DATA )
{
    printf ( "%s: A = %p LDA = %d\n", desc, DATA.A, DATA.LDA );
}

void sgemv_wrapper ( int N, void** DATA, float* x, float* y )
{
    checkpoint t0 = tic();
    data_struct* pDATA = (data_struct*)DATA;
    // read the relevant data from the struct
    float* A = pDATA->A;
    int LDA = pDATA->LDA;

    float f_one = 1.0f;
    float f_zero = 0.0f;
    int i_one = 1;
    // calculate y <-- Ax
    sgemv_ ( "N", N, N, &f_one, A, &LDA, x, &i_one, &f_zero, y, &i_one, 1 );
    pDATA->sgemv_time += toc ( t0 );
    pDATA->sgemv_calls += 1;
}

void setup_S_and_T ( int N, float* S, float* T, int LDMAT, float shift, double* timing_data )
{
    // first perform the shift-invert process using LAPACK
    // T = T - shift*S
    checkpoint t0 = tic();
    apply_shift ( N, S, T, LDMAT, shift );
    timing_data[SHIFT] = toc(t0);
    t0 = tic();
    // perform the LU decomposition of T
    int IPIV[N];
    int INFO;
    // solve using the LU decomposition
    sgesv_( (int*)&N, (int*)&N, T, (int*)&LDMAT, IPIV, S, (int*)&LDMAT, &INFO );
    timing_data[SGETRF_S] = toc(t0);
}

void init_data ( data_struct* DATA, float* S, int LDMAT )
{
    DATA->A = S;
    DATA->LDA = LDMAT;
    DATA->sgemv_calls = 0;
    DATA->sgemv_time = 0.0;
}

// solve the eigen system Sx = lTx
int dense_sgev ( int N, float* S, float* T, int LDMAT, int NEV, float shift, float* eigenvalues, float* eigenvectors, double* timing_data_10, int* int_data_10 )
{
    checkpoint t0;
    printf ( "dense_gev: %d %d %f\n", N, LDMAT, shift );
    int result = 0;

    // apply the shift to the matrices
    setup_S_and_T ( N, S, T, LDMAT, shift, timing_data_10 );
    data_struct DATA;
    init_data( &DATA, S, LDMAT );

    t0 = tic();
    calculate_eigen_values ( N, &DATA, NEV, shift, eigenvalues, eigenvectors, "SR" );
    timing_data_10[ARPACK_CALL] = toc( t0 );

    timing_data_10[SGEMV] = DATA.sgemv_time;
    int_data_10[SGEMV] = DATA.sgemv_calls;
    return result;
}

