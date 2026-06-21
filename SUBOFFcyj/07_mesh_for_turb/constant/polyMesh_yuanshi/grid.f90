PROGRAM MAIN
IMPLICIT NONE
INCLUDE 'parameter.h'
INTEGER :: I, J, K, P_NUM
REAL :: Y(0:NI)
REAL :: RECT_X(NL+1, NL+1), RECT_Y(NL+1, NL+1), RECT_Z(NL+1, NL+1)
REAL :: XP, YP, ZP, XIELV
REAL :: TEMP


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


TEMP=Y(NI)-Y(NI-1)
WRITE(*,'(A,E22.15)') 'Translating dx = ', TEMP




!write points file
P_NUM=(NL+1)**2*(3+1)
OPEN(2,FILE='../polyMesh/points',STATUS='UNKNOWN')
WRITE(2,*) "/*--------------------------------*- C++ -*----------------------------------*\"
WRITE(2,*) "  =========                 |"
WRITE(2,*) "  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox"
WRITE(2,*) "   \\    /   O peration     | Website:  https://openfoam.org"
WRITE(2,*) "    \\  /    A nd           | Version:  v2312"
WRITE(2,*) "     \\/     M anipulation  |"
WRITE(2,*) "\*---------------------------------------------------------------------------*/"
WRITE(2,*) "FoamFile"
WRITE(2,*) "{"
WRITE(2,*) "    version     2.0;"
WRITE(2,*) "    format      ascii;"
WRITE(2,*) "    class       vectorField;"
WRITE(2,*) '    location    "constant/polyMesh";'
WRITE(2,*) "    object      points;"
WRITE(2,*) "}"
WRITE(2,*) "// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //"
WRITE(2,*) " "
WRITE(2,*) " "
WRITE(2,231) P_NUM
WRITE(2,*) "("
DO J=NL+1,1,-1 !YMAX->YMIN
    DO K=1,NL+1 ! ZMIN->ZMAX
        IF (J.EQ.(NL/2+1) .AND. K.EQ.(NL/2+1)) THEN
                XP=-Y(NI)
                YP=0.0
                ZP=0.0
                XP=XP
                YP=YP
                ZP=ZP
                WRITE(2,230) '( ',XP,' ',YP,' ',ZP,' )'
            DO I=1,3
                XP=XP+TEMP
                WRITE(2,230) '( ',XP,' ',YP,' ',ZP,' )'
            ENDDO
        ELSE
                CALL GET_K(RECT_X(J,K),XIELV)
                CALL GET_P(RECT_X(J,K),RECT_Y(J,K),RECT_Z(J,K),XIELV,Y(NI),XP,YP,ZP)
                XP=XP
                YP=YP
                ZP=ZP
                WRITE(2,230) '( ',XP,' ',YP,' ',ZP,' )'                 
            DO I=1,3
                XP=XP+TEMP
                WRITE(2,230) '( ',XP,' ',YP,' ',ZP,' )'
            ENDDO
        ENDIF
    ENDDO
ENDDO
WRITE(2,*) ")"
WRITE(2,*) " "
WRITE(2,*) " "
WRITE(2,*) "// ************************************************************************* //"
CLOSE(2)
230 FORMAT(A2,E22.15,A1,E22.15,A1,E22.15,A2)
231 FORMAT(I0)        

END PROGRAM



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