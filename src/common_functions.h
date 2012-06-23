/*
 * Evan Lezar
 * 18 November 2010
 *
 * Post processing of the eigenvalues calcuated using the ARPACK-based solvers routines
 *
 */

void apply_shift ( int N, float* S, float* T, int LDMAT, float shift )
{
    // T = T - shift*S
    checkpoint t0 = tic();
    int row, col;
#pragma omp parallel default(shared) private(row, col)
    {
#pragma omp for schedule(runtime)
    for ( col = 0; col < N; ++col )
        for ( row = 0; row < N; ++row )
            T[col*LDMAT + row] = T[col*LDMAT + row] - shift*S[col*LDMAT + row];
    }
}

struct float_complex {
    float x, y;
};
typedef struct float_complex float_complex;

inline float_complex make_float_complex ( float x, float y )
{
    float_complex c = {x, y};
    return c;

}

float_complex invert ( float_complex v )
{
    float_complex l;
    l.x = v.x/(v.x*v.x + v.y*v.y);
    l.y = -v.y/(v.x*v.x + v.y*v.y);

    return l;
}

// return v / ( 1.0 + v*shift );
float_complex unshift ( float_complex v, float shift )
{
    float_complex l = v;
    // calculate ( 1.0 + v*shift )
    l.x = l.x*shift + 1.0;
    l.y *= shift;

    float_complex result;

    // calculate v/l
    result.x = ( v.x*l.x + v.y*l.y ) / ( l.x*l.x + l.y*l.y );
    result.y = ( v.y*l.x - v.x*l.y ) / ( l.x*l.x + l.y*l.y );

    return result;
}

float eigenvalue_magnitude ( const float real, const float imag )
{
    return sqrt ( real*real + imag*imag );
}

int insert_smallest_into_list ( const float* value_list, int* index_list, const int N, const float value, const int index )
{
    if ( isnan(value) | isinf(value) )
        {
            return -1;
        }
        if ( index == 0 )
        { // first element to insert
            index_list[0] = 0;
            return 0;
        }

        int i, j, insert_at = -1;
        int limit = index;
        if ( limit > (N-1) )
            limit = (N-1);

        for ( i=0; i < limit; ++i )
        {
            if ( value_list[index_list[i]] > value )
            {
                insert_at = i;
                int last = index_list[limit];
                for ( j = limit; j > i; j-- )
                {
                    index_list[j] = index_list[j-1];
                }
                if ( (limit+1) < N )
                    index_list[limit+1] = last;
                break;
            }
        }

        if ( insert_at < 0 )
        {
            insert_at = limit;
        }

        index_list[insert_at] = index;

        return insert_at;
}


int a_is_larger_than_b ( float a, float b )
{
    if ( isnan(a) | isinf(a) )
        return 1;

    if ( isnan(b) | isinf(b) )
        return 0;

    if ( a > b )
        return 1;
    else
        return 0;
}

int insert_largest_into_list ( const float* value_list, int* index_list,
                                const int N, const float value, const int index )
{
    if ( index == 0 )
    { // first element to insert
        index_list[0] = 0;
        return 0;
    }

    int i, j, insert_at = -1;
    int limit = index;
    if ( limit > (N-1) )
        limit = (N-1);

    for ( i=0; i < limit; ++i )
    {
        if ( a_is_larger_than_b(value, value_list[index_list[i]]) )
        {
            insert_at = i;
            int last = index_list[limit];
            for ( j = limit; j > i; j-- )
            {
                index_list[j] = index_list[j-1];
            }
            if ( (limit+1) < N )
                index_list[limit+1] = last;
            break;
        }
    }

    if ( insert_at < 0 )
    {
        insert_at = limit;
    }

    index_list[insert_at] = index;

    return insert_at;
}


void get_smallest_magnitude_eigenvalues ( int N, int NEV, int num_calculated, float shift,
                                          float_complex* eigenvalues, float* eigenvectors, float* residuals,
                                          float_complex* temp_eigenvalues, const float* temp_eigenvectors, const float* temp_residuals )
{
    printf("pointers: %p %p %p %p %p %p\n", eigenvalues, eigenvectors, residuals, temp_eigenvalues, temp_eigenvectors, temp_residuals );
    int i;
    int index;
    int smallest_index[num_calculated];
    float magnitude[num_calculated];

    float_complex v;
    float_complex l;

    for ( i = 0; i < num_calculated; ++i )
    {
        v = temp_eigenvalues[i];
        l = unshift(v, shift);

        printf("%d: %f +i%f : %f ::: %f +i%f\n", i, v.x, v.y, temp_residuals[i], l.x, l.y );

        magnitude[i] = eigenvalue_magnitude ( l.x, l.y );
        insert_smallest_into_list( magnitude, smallest_index, num_calculated, magnitude[i], i );
    }

    int j = 0;
    i = 0;

    while ( i < NEV && j < num_calculated )
    {
        index = smallest_index[j];
        if ( 0 <= index  && index < num_calculated )
        {
            v = temp_eigenvalues[index];
            l = unshift(v, shift);
            if (  v.x < 0 && l.x > 0 && temp_residuals[index] < 0.5 )
            {
                eigenvalues[i] = l;

                residuals[i] = temp_residuals[index];
                memcpy ( eigenvectors + N*i, temp_eigenvectors + N*index, N*sizeof(float) );
                printf("using: %d : %f+j%f : %e\n", index, v.x, v.y, residuals[i] );
                i++;
            }
        }
        j++;
    }
}

void calculate_eigen_values ( int N, void* DATA, int NEV, float shift, float* eigenvalues, float* eigenvectors, char* which )
{
    int i;
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
    float* residuals = (float*)malloc ( (use_N_ev)*sizeof(float));

    // solve the eigenvalue problem using ARPACK
    arpack_ssev(N, (void*)DATA, use_N_ev, NCV, temp_ev, temp_vectors, temp_residuals, which );

    
    // Copy the resultant eigenvalues
    memcpy(eigenvalues, temp_ev, NEV*2*sizeof(float));
    memcpy(eigenvectors, temp_vectors, NEV*N*sizeof(float));

    for (i=0; i < NEV; ++i)
    {
        printf("%d: %f\n", i, temp_residuals[i] );
    }

    // free the temporary storage
    free ( temp_ev );
    printf("1:\n");
    free ( temp_vectors );
    printf("2:\n");
    free ( temp_residuals );
    printf("3:\n");
    free ( residuals );
    printf("4:\n");
}

void calculate_desired_eigenvalues ( int N, void* DATA, int NEV, float shift, float* eigenvalues, float* eigenvectors  )
{
    // solve the eigenvalue problem using ARPACK
    int use_N_ev = 2*NEV;
    if ( use_N_ev > ( N/2 - 1 ) )
        use_N_ev = N/2 - 1;
    if ( use_N_ev <= 0 )
        use_N_ev = 1;
    // select the number of Arnoldi vectors to generate
    int NCV = 2*use_N_ev + 1;
    if ( NCV > N )
        NCV = N;

    // allocate temporary storage for the vectors
    float_complex* temp_ev = (float_complex*)malloc ( NCV*sizeof(float_complex) );
    float* temp_vectors = (float*) malloc (N*NCV*sizeof(float));
    float* temp_residuals = (float*)malloc ( (NCV )*sizeof(float));
    float* residuals = (float*)malloc ( (NCV)*sizeof(float));
#ifdef DEBUG_OUTPUT
    printf( "N=%d, use_N_ev = %d, NCV = %d, NEV = %d\n", N, use_N_ev, NCV, NEV );
#endif
    arpack_ssev( (int*)&N, (void*)DATA, &use_N_ev, &NCV, (float*)temp_ev, temp_vectors, temp_residuals, "LR" );
    // free the temporary storage



    int largest_index[NCV];
    float real[NCV];
    float_complex v;
    float_complex l;

    int i;
    for ( i = 0; i < NCV; i++ )
    {
        v = temp_ev[i];
        l = invert (v);
#ifdef DEBUG_OUTPUT
        printf("%d: %f+j%f :: %f(%f)+j%f :: %e\n", i, v.x, v.y, l.x, l.x + shift, l.y, temp_residuals[i]);
#endif
        real[i] = v.x;
        insert_largest_into_list ( real, largest_index, NCV, real[i], i );
    }
    int index, j = 0;
    i = 0;

    while ( i < NEV && j < NCV )
    {
        index = largest_index[j];
        if ( 0 <= index  && index < NCV )
        {
            v = temp_ev[index];
            l = invert(v);
            float_complex t = invert(l);

            if (  temp_residuals[index] > 1e-9 && temp_residuals[index] < 0.5 && t.x > 1e-5 )
            {
                eigenvalues[2*i] = l.x + shift;
                eigenvalues[2*i+1] = l.y;

                residuals[i] = temp_residuals[index];
                memcpy ( eigenvectors + N*i, temp_vectors + N*index, N*sizeof(float) );
#ifdef DEBUG_OUTPUT
                printf("using: %d : %f+j%f : %e+j%e : %e (%f +j%f)\n", index, v.x, v.y, l.x, l.y, residuals[i], t.x, t.y );
#endif
                i++;
            }
        }
        j++;
    }



    free ( temp_ev );
    free ( temp_vectors );
    free ( temp_residuals );
    free ( residuals );
}

