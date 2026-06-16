!--------------------------------------------------------------
!     GRID GENERATION IN CHANNEL FLOW
!     X2 Clustering only
!     JICHOI
!--------------------------------------------------------------

      PROGRAM MAIN
      IMPLICIT NONE
      INCLUDE 'parameter.h'

      REAL*8 PAI
      REAL*8 R,XP,YP,ZP
      REAL*8 CY_X(NJ+1,NK),CY_Y(NJ+1,NK),CY_Z(NJ+1,NK)
      REAL*8 CY_X0(NJ+1,NK),CY_Y0(NJ+1,NK),CY_Z0(NJ+1,NK)
      REAL*8 Y0(0:3),Y(0:NI)
      INTEGER ITER,I,J,K,KK
      REAL *8 XIELV
      !write points
      INTEGER P_NUM
      
      
      
      PAI=ACOS(-1.0d0)
      
   
       ! 读取壁面网格XYZ坐标数据
        OPEN(20,FILE='ReadData/CY_surface_XYZ.plt',STATUS='OLD')
        DO K=1,NK
        DO J=1,NJ+1
            READ(20,*) CY_X0(J,K),CY_Y0(J,K),CY_Z0(J,K)
        ENDDO
        ENDDO
        CLOSE(20)

        OPEN(1,FILE='ReadData/Y.txt',STATUS='UNKNOWN')
        DO I=0,NI
           READ(1,*) Y(I)
        ENDDO
        CLOSE(1)

        ! 计算层流入口的XYZ坐标
        DO K=1,NK
            DO J=1,NJ+1
            DO I=NI,NI
                CALL GET_K(CY_X0(J,K),XIELV)
                CALL GET_P(CY_X0(J,K),CY_Y0(J,K),CY_Z0(J,K),XIELV,Y(I),CY_X(J,K),CY_Y(J,K),CY_Z(J,K))
            ENDDO
            ENDDO
        ENDDO
      
      
      ! 随机厚度
      Y0(0)=0.0
      DO J=1,3
          Y0(J)=DBLE(J)*ABS(CY_X(NJ+1,1)-CY_X(NJ,1))
      ENDDO
      
      
      
      !write points file
      P_NUM=(NK+1)*(NJ+1)*(3+1) ! 圆柱双连通拓扑改为单连通拓扑
      OPEN(2,FILE='../polyMesh/points',STATUS='UNKNOWN')
    WRITE(2,*) "/*--------------------------------*- C++ -*----------------------------------*\"
    WRITE(2,*) "  =========                 |"
    WRITE(2,*) "  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox"
    WRITE(2,*) "   \\    /   O peration     | Website:  https://openfoam.org"
    WRITE(2,*) "    \\  /    A nd           | Version:  8"
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
      DO K=1,NK
      DO J=1,NJ+1
      DO I=0,3
          XP=CY_X(J,K)
          R=SQRT(CY_Y(J,K)**2+CY_Z(J,K)**2)
          YP=CY_Y(J,K)*(R-Y0(I))/R
          ZP=CY_Z(J,K)*(R-Y0(I))/R
          WRITE(2,230) '( ',XP,' ',YP,' ',ZP,' )' 
      ENDDO
      ENDDO
      ENDDO
      DO K=1,1
      DO J=1,NJ+1
      DO I=0,3
          XP=CY_X(J,K)
          R=SQRT(CY_Y(J,K)**2+CY_Z(J,K)**2)
          YP=CY_Y(J,K)*(R-Y0(I))/R
          ZP=CY_Z(J,K)*(R-Y0(I))/R
          WRITE(2,230) '( ',XP,' ',YP,' ',ZP,' )' 
      ENDDO
      ENDDO
      ENDDO

    
    WRITE(2,*) ")"
    WRITE(2,*) " "
    WRITE(2,*) " "
    WRITE(2,*) "// ************************************************************************* //"
    CLOSE(2)

    
230 FORMAT(A2,E22.15,A1,E22.15,A1,E22.15,A2)
231 FORMAT(I0)        
    END
    
    
    
    SUBROUTINE GET_P(XH,YH,ZH,XIELV,DY,XOUT,YOUT,ZOUT)
    IMPLICIT NONE
    REAL*8 XH,YH,ZH,XIELV,DY,XOUT,YOUT,ZOUT
    REAL*8 X,Y
    
    X=XH
    Y=SQRT(YH*YH+ZH*ZH)
    XOUT=X-DY*XIELV/SQRT(1.+XIELV*XIELV)
    YOUT=( Y+DY/SQRT(1.+XIELV*XIELV) )*YH/Y
    ZOUT=( Y+DY/SQRT(1.+XIELV*XIELV) )*ZH/Y
  
    END SUBROUTINE
    
    

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
    
