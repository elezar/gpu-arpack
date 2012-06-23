! Evan Lezar
! 18 November 2010
! Wrapper for the ARPACK standard eigenval problems
! including GPU acceleration

module arpack_driver
    use iso_c_binding

    implicit none

    private

    interface 
        subroutine sgemv_wrapper ( N, DATA, x, y ) &
        bind(c, name='sgemv_wrapper')
        import
            integer(c_int), value :: N
            type(c_ptr), value :: DATA
            real(c_float), dimension(N), intent(in) :: x
            real(c_float), dimension(N), intent(out) :: y
        end subroutine sgemv_wrapper

    end interface

    ! External LAPACK FORTRAN function.
    real(4), external :: snrm2, slapy2


    contains
! single precision standard eigenvalue problem
subroutine arpack_ssev ( n, data, nev, ncv, eigenvalues, eigenvectors, residuals, which ) &
bind(c, name='arpack_ssev')
    
    implicit none
    character(kind=c_char) :: which(2)
    ! setup the parameters
    integer(c_int), value :: n, nev, ncv ! matrix size, number of eigenvalues to compute, number of arnoldi vectors
    ! matrices
    real(c_float) :: &
        eigenvectors(n, ncv), &
        residuals (ncv)
    complex(c_float) :: &
        eigenvalues(ncv)
    
    type(c_ptr), value :: data

    ! local variables
    ! ARPACK arrays
    integer iparam(11), ipntr(14)
    logical select(ncv)
    real :: &
        ax(n), e_real(ncv), e_imag(ncv), resid(n), &
        workd(n*3), &
        workev(3*ncv), &
        workl(3*ncv*ncv + 6*ncv)

    ! ARPACK scalars
    integer ido, lworkl, info, j, &
            ierr, nconv, maxitr, ishifts, mode, lda

    real    tol, one, zero, sigma, temp_res
    logical first
    
    ! initialise constants
    one = (1.0)
    zero = (0.0)
    sigma = (0.0)
    ! initialise ARPACK parameters
    lworkl = 3*ncv*ncv + 6*ncv
    tol    = 0.0
    ido    = 0
    INFO   = 0
    ishifts = 1
    maxitr = 1000
    mode   = 1 ! Standard eigenvalue problem
    ! Setup ARPACK parameter array
    iparam(1) = ishifts
    iparam(3) = maxitr
    iparam(7) = mode

    10 CONTINUE

    ! call the ARPACK routines snaupd
        call SNAUPD  ( ido, 'I', N, WHICH, NEV, tol, resid, NCV, &
                       eigenvectors, N, iparam, ipntr, WORKD, &
                       WORKL, LWORKL, INFO )
    ! depending on the value of ido, various operations must be performed
        IF ( ido .EQ. -1 .OR. ido .EQ. 1) THEN
    ! Perform matrix vector multiplication
    !       y <--- A*x
    !       x = workd(ipntr(1))
    !       y = workd(ipntr(2))
    !       the wrapper is provided externally and takes the parameters N, DATA, x, y
    !           DATA contains implemetnation-specific information such as pointers to A or its factors
            call sgemv_wrapper ( N, DATA, workd(ipntr(1):), workd(ipntr(2):) )
            GO TO 10
        END IF


    IF ( INFO .LT. 0 ) THEN
        print *, ' Error with _naupd, info = ', info
        print *, ' Check the documentation of _naupd'
    ELSE
! Post process the result to get the eigenvalues and eigenvectors
        call SNEUPD  ( .TRUE., 'A', select, e_real, e_imag, &
                       eigenvectors, N, sigma, sigma,  &
                        WORKEV, 'I', N, WHICH, NEV, tol, resid, NCV, &
                        eigenvectors, N, iparam, ipntr, WORKD, &
                        WORKL, LWORKL, IERR )
         IF ( IERR .NE. 0) THEN
             print *, ' Error with _neupd, info = ', IERR
             print *, ' Check the documentation of _neupd. '
         ELSE

             first = .TRUE.
             nconv = iparam(5)
             LDA = N
             DO 20 J=1, nconv
!               Compute the residual norm
!               ||  A*x - lambda*x ||

                IF ( e_imag(j) .EQ. zero ) THEN
                    !       the wrapper is provided externally and takes the parameters N, DATA, x, y
                    !           DATA contains implemetnation-specific information such as pointers to A or its factors
                    call sgemv_wrapper ( N, DATA, eigenvectors(1,j), ax )
                    call saxpy (N, -e_real(j), eigenvectors(1,j), 1, ax, 1)
                    temp_res = snrm2 (N, ax, 1)

                    temp_res = temp_res / abs( e_real(j) )

                ELSE IF ( first ) THEN
                    !       the wrapper is provided externally and takes the parameters N, DATA, x, y
                    !           DATA contains implemetnation-specific information such as pointers to A or its factors
                    call sgemv_wrapper ( N, DATA, eigenvectors(1,j), ax )
                    call saxpy (N, -e_real(j), eigenvectors(1,j), 1, ax, 1)
                    call saxpy (N, e_imag(j), eigenvectors(1,j+1), 1, ax, 1)
                    temp_res = snrm2(N, ax, 1)

                    !       the wrapper is provided externally and takes the parameters N, DATA, x, y
                    !           DATA contains implemetnation-specific information such as pointers to A or its factors
                    call sgemv_wrapper ( N, DATA, eigenvectors(1,j+1), ax )
                    call saxpy (N, -e_real(j), eigenvectors(1,j+1), 1, ax, 1)
                    call saxpy (N, -e_imag(j), eigenvectors(1,j), 1, ax, 1)
                    temp_res = slapy2( temp_res, snrm2(n, ax, 1) )

                    temp_res = temp_res / slapy2( e_real(j), e_imag(j) )
                    first = .FALSE.
                ELSE
                  first = .TRUE.
                END IF

                eigenvalues(j) = CMPLX( e_real(j), e_imag(j) )
                residuals(j) = temp_res

 20          CONTINUE
         END IF

         IF ( INFO .EQ. 1) THEN
             print *, ' Maximum number of iterations reached.'
         ELSE IF ( INFO .EQ. 3) THEN
             print *, ' No shifts could be applied during implicit', &
                      ' Arnoldi update, try increasing NCV.'
         END IF
    END IF

    end subroutine arpack_ssev

end module arpack_driver
