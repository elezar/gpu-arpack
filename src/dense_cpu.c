/*
 * Evan Lezar
 * An implementation of a dense eigensolver that uses the ARPACK Fortran backend.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <omp.h>
#include <cblas.h>
// #include <acml.h>

// Include the interfaces to the Fortran routines.
#include "interface.h"

#define ARPACK_CALL 0
#define SHIFT 1
#define SGETRF_S 2
#define SGEMV 9

//!
//! A data structure used to pass the eigenproblem configuration to Fortran and 
//! back to C again.
struct data_struct {
    float* A;           //! A pointer to the matrix representing the eigenproblem.
    double sgemv_time;  //! The total time required to calculate the matrix-vector products.
    int sgemv_calls;    //! The number of times the matrix-vector product was called.
    int LDA;            //! The leading dimension of the matrix.
};
typedef struct data_struct data_struct;

//! 
//! A utility function to print the structure representing the eigenproblem.
void print_data ( const char* desc, data_struct DATA )
{
    printf ( "%s: A = %p LDA = %d\n", desc, DATA.A, DATA.LDA );
}

//!
//! Calculate the matrix-vector product  y <- Ax for the eigensystem defined by DATA. This
//! routine is called from the Fortran backend whenever ARPACK requires a matrix-vector
//! product to be calculated.
//! \param[in] N The dimension of the eigensystem.
//! \param[in,out] DATA The structure representing the eigenproblem.
//! \param[in] x The vector that must be multiplied by A.
//! \param[out] y The vector that must store the result.
void sgemv_wrapper ( int N, data_struct* DATA, float* x, float* y )
{
    checkpoint t0 = tic();
    // read the relevant data from the struct
    float* A = DATA->A;
    int LDA = DATA->LDA;
    // calculate y <-- Ax
    cblas_sgemv ( CblasColMajor, CblasNoTrans, N, N, 1.0, A, LDA, x, 1, 0.0, y, 1);
    DATA->sgemv_time += toc ( t0 );
    DATA->sgemv_calls += 1;
}

//!
//! Initialise the struct representing the eigensystem.
void init_data ( data_struct* DATA, float* A, int LDMAT )
{
    DATA->A = A;
    DATA->LDA = LDMAT;
    DATA->sgemv_calls = 0;
    DATA->sgemv_time = 0.0;
}

//!
//! Solve the actual eigensystem. This takes the eigensystem defined in the structure DATA and allocates
//! the required temporary workspaces before calling the Fortran backend.
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


//!
//! Solve the standard eigensystem Ax = lx.eigen
//! \param[in] N The number of columns in the matrix A
//! \param[in] A The matrix representing the eigensystem to be solved.
//! \param[in] NEV The number of eigenvalues to calculate.
//! \param[out] eigenvalues A vector of the NEV eigenvalues. Note that the eigenvectors are complex.
//! \param[out] eigenvectors An NxNEV matrix with the eigenvectors as columns.
//! \param[out] timing_data_10 A 10-vector representing the timing data for various phases of the process.
//! \param[out] int_data_10 A 10-vector containing some run information.
//! \return A non-zero error code if an error occured.
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
