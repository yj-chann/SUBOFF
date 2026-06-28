!--------------------------------------------------------------
!     GRID GENERATION IN CHANNEL FLOW
!     X2 Clustering only
!     JICHOI
!--------------------------------------------------------------

      SUBROUTINE GET_GRID(COUNT1)
        INCLUDE 'head.fi'
        INTEGER I,J,K
        REAL Y(0:NI)
        REAL RECT_X(NL+1,NL+1),RECT_Y(NL+1,NL+1),RECT_Z(NL+1,NL+1)
        REAL CY_X(NK,NJ+1),CY_Y(NK,NJ+1),CY_Z(NK,NJ+1)
        INTEGER*8 COUNT1, COUNT2, COUNT3, COUNT_RATE,COUNT_MAX
      
    IF (iMPI_MyID.EQ.0) THEN

    ! 读取壁面网格CY XYZ坐标数据
    OPEN(20,FILE='ReadData/CY_surface_XYZ.plt',STATUS='OLD')
    DO K=1,NK
    DO J=1,NJ+1
        READ(20,*) CY_X(K,J),CY_Y(K,J),CY_Z(K,J)
    ENDDO
    ENDDO
    CLOSE(20)

    OPEN(20, FILE='ReadData/RECT_surface_XYZ.plt', STATUS='OLD', FORM='FORMATTED')
    DO J = 1, NL+1
        DO K = 1, NL+1
        READ(20, *) RECT_X(J, K), RECT_Y(J, K), RECT_Z(J, K)
        END DO
    END DO
    CLOSE(20)

    OPEN(1,FILE='ReadData/Y.txt',STATUS='UNKNOWN')
    DO I=0,NI
        READ(1,*) Y(I)
    ENDDO
    CLOSE(1)
      
    
      
    ENDIF      
      

      call mpi_bcast( CY_X,NK*(NJ+1), mpi_real8, 0, MPI_COMM_WORLD, iMPI_ErrorInfo )
      call mpi_bcast( CY_Y,NK*(NJ+1), mpi_real8, 0, MPI_COMM_WORLD, iMPI_ErrorInfo )
      call mpi_bcast( CY_Z,NK*(NJ+1), mpi_real8, 0, MPI_COMM_WORLD, iMPI_ErrorInfo )
      call mpi_bcast( RECT_X,(NL+1)*(NL+1), mpi_real8, 0, MPI_COMM_WORLD, iMPI_ErrorInfo )
      call mpi_bcast( RECT_Y,(NL+1)*(NL+1), mpi_real8, 0, MPI_COMM_WORLD, iMPI_ErrorInfo )
      call mpi_bcast( RECT_Z,(NL+1)*(NL+1), mpi_real8, 0, MPI_COMM_WORLD, iMPI_ErrorInfo )
      call mpi_bcast( Y,NI+1, mpi_real8, 0, MPI_COMM_WORLD, iMPI_ErrorInfo )
      CALL MPI_Barrier( MPI_COMM_WORLD, iMPI_ErrorInfo )
      
      If ( iMPI_MyID.EQ.0 ) THEN
      CALL SYSTEM_CLOCK(COUNT2,COUNT_RATE,COUNT_MAX)
      WRITE(*,*) 'GRID_READ',(COUNT2-COUNT1)/DBLE(COUNT_RATE)
      ENDIF
      CALL GET_POINTS_INI_UP(CY_X,CY_Y,CY_Z,RECT_X,RECT_Y,RECT_Z,Y,COUNT2)
      
      If ( iMPI_MyID.EQ.0 ) THEN
      CALL SYSTEM_CLOCK(COUNT3,COUNT_RATE,COUNT_MAX)
      WRITE(*,*) 'GRID_INI_UP',(COUNT3-COUNT2)/DBLE(COUNT_RATE)
      ENDIF


    END SUBROUTINE
    
    
    
SUBROUTINE GET_P(XH,YH,ZH,XIELV,DY,XOUT,YOUT,ZOUT)
IMPLICIT NONE
REAL*8 XH,YH,ZH,XIELV,DY,XOUT,YOUT,ZOUT
REAL*8 X,Y

X=XH
Y=SQRT(YH*YH+ZH*ZH)
XOUT=X-DY*XIELV/SQRT(1.+XIELV*XIELV)
YOUT=( Y+DY/SQRT(1.+XIELV*XIELV) )*YH/Y
ZOUT=( Y+DY/SQRT(1.+XIELV*XIELV) )*ZH/Y

END SUBROUTINE GET_P



subroutine GET_K(x_in, k_out)
    implicit none
    ! Input and output variables
    real(8), intent(in)  :: x_in
    real(8), intent(out) :: k_out
    
    ! Constants
    real(8), parameter :: A = 1.126395101d0
    real(8), parameter :: B = 0.442874707d0
    real(8), parameter :: RATIO = 0.2d0
    real(8), parameter :: R_MAX = 5.0d0 / 6.0d0
    real(8), parameter :: ALPHA_MINUS_ONE = (1.0d0 / 2.1d0) - 1.0d0
    
    ! Local variables
    real(8) :: x
    real(8) :: px, px2, px3, px4, p12x
    real(8) :: term_base, term_mult
    
    ! Apply the initial scaling
    x = x_in / (RATIO * 0.3048d0)
    
    ! Check condition (equivalent to x <= 10/3)
    if (x <= (10.0d0 / 3.0d0)) then
        ! Precompute repeated terms for cleaner, more efficient math
        px = 0.3d0 * x - 1.0d0
        px2 = px * px
        px3 = px2 * px
        px4 = px3 * px
        p12x = 1.2d0 * x + 1.0d0
    
        ! First part of the complex expression (the base being raised to the power)
        term_base = A * x * px4 + B * (x**2) * px3 + 1.0d0 - px4 * p12x
    
        ! Second part of the complex expression (the multiplier)
        term_mult = 4.0d0 * A * x * px3 * 0.3d0 &
                    + A * px4 &
                    + 3.0d0 * B * (x**2) * px2 * 0.3d0 &
                    + 2.0d0 * B * x * px3 &
                    - 4.0d0 * px3 * 0.3d0 * p12x &
                    - 1.2d0 * px4
    
        ! Combine everything together
        k_out = (1.0d0 / 2.1d0) * R_MAX * (term_base**ALPHA_MINUS_ONE) * term_mult
    else
        k_out = 0.0d0
    end if
    
end subroutine GET_K
    
    