! Evan Lezar
! 18 November 2010
! Wrapper for the ARPACK standard eigenval problems
! including GPU acceleration

module arpack_driver
    use iso_c_binding

    implicit none

    interface 
        subroutine SGEMV_wrapper ( N, DATA, x, y ) &
        bind(c, name='SGEMV_wrapper')
        import
            integer, intent(in) :: N
            type(c_ptr) :: DATA
            real(4), intent(in) :: x
            real(4), intent(out) :: y
        end subroutine SGEMV_wrapper    

    end interface

    contains
! single precision standard eigenvalue problem
subroutine arpack_ssev ( n, data, nev, ncv, eigenvalues, eigenvectors, residuals, which )
    
    implicit none
    CHARACTER(len=2) WHICH
    ! setup the parameters
    INTEGER N, NEV, NCV ! matrix size, number of eigenvalues to compute, number of arnoldi vectors
    ! matrices
    REAL :: &
        eigenvectors(N, NCV), &
        residuals (NCV)
    COMPLEX*8 :: &
        eigenvalues(NCV)
    
    type(c_ptr) :: DATA

    ! local variables
    ! ARPACK arrays
    integer iparam(11), ipntr(14)
    logical select(NCV)
    REAL :: &
        ax(N), e_real(NCV), e_imag(NCV), resid(N), &
        WORKD(N*3), &
        WORKEV(3*NCV), &
        WORKL(3*NCV*NCV + 6*NCV)

    ! ARPACK scalars
    INTEGER ido, LWORKL, INFO, J, &
            IERR, nconv, maxitr, ishifts, mode, LDA

    REAL    tol, one, zero, sigma, temp_res
    LOGICAL first
    ! declare the SGEMV wrapper as external
    !EXTERNAL SGEMV_wrapper

    ! initialise constants
    one = (1.0)
    zero = (0.0)
    sigma = (0.0)
    ! initialise ARPACK parameters
    LWORKL = 3*NCV*NCV + 6*NCV
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
            call SGEMV_wrapper ( N, DATA, workd(ipntr(1)), workd(ipntr(2)) )
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
                    call SGEMV_wrapper ( N, DATA, eigenvectors(1,j), ax )
                    call saxpy (N, -e_real(j), eigenvectors(1,j), 1, ax, 1)
                    temp_res = snrm2 (N, ax, 1)

                    temp_res = temp_res / abs( e_real(j) )

                ELSE IF ( first ) THEN
                    !       the wrapper is provided externally and takes the parameters N, DATA, x, y
                    !           DATA contains implemetnation-specific information such as pointers to A or its factors
                    call SGEMV_wrapper ( N, DATA, eigenvectors(1,j), ax )
                    call saxpy (N, -e_real(j), eigenvectors(1,j), 1, ax, 1)
                    call saxpy (N, e_imag(j), eigenvectors(1,j+1), 1, ax, 1)
                    temp_res = snrm2(N, ax, 1)

                    !       the wrapper is provided externally and takes the parameters N, DATA, x, y
                    !           DATA contains implemetnation-specific information such as pointers to A or its factors
                    call SGEMV_wrapper ( N, DATA, eigenvectors(1,j+1), ax )
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
