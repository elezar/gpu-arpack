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

void sgemv_wrapper ( int N, data_struct* DATA, float* x, float* y )
{
    checkpoint t0 = tic();
    // read the relevant data from the struct
    float* A = DATA->A;
    int LDA = DATA->LDA;
    float f_one = 1.0f;
    float f_zero = 0.0f;
    int i_one = 1;
    // calculate y <-- Ax
    sgemv_ ( "N", &N, &N, &f_one, A, &LDA, x, &i_one, &f_zero, y, &i_one, 1 );
    DATA->sgemv_time += toc ( t0 );
    DATA->sgemv_calls += 1;
}

void init_data ( data_struct* DATA, float* A, int LDMAT )
{
    DATA->A = A;
    DATA->LDA = LDMAT;
    DATA->sgemv_calls = 0;
    DATA->sgemv_time = 0.0;
}


void calculate_eigen_values ( int N, void* DATA, int NEV, float* eigenvalues, float* eigenvectors, char* which )
{
    int use_N_ev = 2*NEV;
    if ( use_N_ev > ( N/2 - 1 ) )
        use_N_ev = N/2 - 1;
    // select the number of Arnoldi vectors to generate
    int NCV = 2*use_N_ev + 1;
    if ( NCV > N )
        NCV = N;

    // allocate temporary storage for the vectors
    float* temp_ev = (float*)malloc ( NCV*2*sizeof(float) );
    float* temp_vectors = (float*) malloc (N*NCV*sizeof(float));
    float* temp_residuals = (float*)malloc ( (NCV )*sizeof(float));

    // solve the eigenvalue problem using ARPACK
    arpack_ssev(N, (void*)DATA, use_N_ev, NCV, temp_ev, temp_vectors, temp_residuals, which );
    
    // Copy the resultant eigenvalues to the previously allocated space.
    memcpy(eigenvalues, temp_ev, NEV*2*sizeof(float));
    memcpy(eigenvectors, temp_vectors, NEV*N*sizeof(float));

    // free the temporary storage
    free ( temp_ev );
    free ( temp_vectors );
    free ( temp_residuals );
}


// solve the standard eigensystem Ax = lx
int dense_seev ( int N, float* A, int LDMAT, int NEV, float* eigenvalues, float* eigenvectors, double* timing_data_10, int* int_data_10 )
{
    checkpoint t0;
    int result = 0;

    // Initialise the data structure that is passed to the ARPACK routines.
    data_struct DATA;
    init_data( &DATA, A, LDMAT );

    t0 = tic();
    // Call a C wrapper function that allows for the calculation of the NEV largest eigenvalues.
    calculate_eigen_values ( N, &DATA, NEV, eigenvalues, eigenvectors, "LM" );
    timing_data_10[ARPACK_CALL] = toc( t0 );

    timing_data_10[SGEMV] = DATA.sgemv_time;
    int_data_10[SGEMV] = DATA.sgemv_calls;
    return result;
}


