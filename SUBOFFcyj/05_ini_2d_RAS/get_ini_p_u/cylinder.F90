PROGRAM cylinder_main
IMPLICIT NONE
INCLUDE 'cylinder.h'

REAL*8 X0(M1,M2),Y0(M1,M2),THETA1,THETA2,UREFX
REAL*8 X_CY(NJ+1,NK),Y_CY(NJ+1,NK),Z_CY(NJ+1,NK)
REAL*8 X_RECT(NL+1,NL+1),Y_RECT(NL+1,NL+1),Z_RECT(NL+1,NL+1)
REAL*8 PAI
REAL*8 X,Y,Z
REAL*8 PREF(N1,N2),UREF(2,N1,N2),U1,U2,U3
REAL*8 PREF_IN(N2),UREF_IN(2,N2)
REAL*8 X_IN(N2),Y_IN(N2)
REAL*8 P_RECT(NL,NL),U_RECT0(2,NL,NL)
REAL*8 P_CY(NJ,NK),U_CY0(2,NJ,NK)
INTEGER*8 I,J,K,IP1,II,J1,K1
CHARACTER*400 DIR
CHARACTER*400 CHAR3

PAI=ACOS(-1.0)
UREFX=1.649194

DIR='../50000/'
OPEN(1,FILE=TRIM(ADJUSTL(DIR))//'p',STATUS='OLD')
OPEN(2,FILE=TRIM(ADJUSTL(DIR))//'U',STATUS='OLD')

DO I=1,23
    READ(1,*)
    READ(2,*)
ENDDO
DO J=1,N2
DO I=1,N1 ! 法向自外向内
    READ(1,*) PREF(I,J)


    READ(2,'(A)') CHAR3
    IP1=INDEX(CHAR3,'(') 
    CHAR3(1:IP1)=''
    IP1=INDEX(CHAR3,')')
    CHAR3(IP1:IP1)=''
    READ(CHAR3,*) UREF(1,I,J),UREF(2,I,J),U3
ENDDO
ENDDO
CLOSE(1)
CLOSE(2)

OPEN(1,FILE='ReadData/suboff_mesh_2d_ADD.plt',STATUS='OLD')   
READ(1,*)
DO I=1,M1  
DO J=1,M2
    READ(1,*) X0(I,J),Y0(I,J)
ENDDO
ENDDO  
CLOSE(1)



I=N1-NI+1
! fields distributed on the face center of the outemost grid points of suboff_mesh_2d
DO J=1,N2
    X_IN(J)=0.5*(X0(I,J)+X0(I,J+1))
    Y_IN(J)=0.5*(Y0(I,J)+Y0(I,J+1))
    PREF_IN(J)=PREF(I,J)
    UREF_IN(:,J)=UREF(:,I,J)
ENDDO

! zf 论文的page34 图3.1(b)
OPEN(1,FILE='PU_REF_XLINE.plt',STATUS='UNKNOWN')
DO J=1,N2
    WRITE(1,236) X_IN(J),PREF_IN(J),UREF_IN(1,J),UREF_IN(2,J)
ENDDO
CLOSE(1) 


! =====================================interpolation================================
OPEN(1,FILE='ReadData/INLET2.plt',STATUS='OLD')
DO K=1,NK
DO J=1,NJ
    READ(1,*) X_CY(J,K),Y_CY(J,K),Z_CY(J,K)
ENDDO
ENDDO
CLOSE(1)
    
OPEN(1,FILE='ReadData/INLET1.plt',STATUS='OLD')
DO K=1,NL
DO J=1,NL
    READ(1,*) X_RECT(J,K),Y_RECT(J,K),Z_RECT(J,K)
ENDDO
ENDDO
CLOSE(1)




! interpolation from 2d data to square-part1（INLET1) data according to x coordinate
DO K=1,NL
DO J=1,NL
    DO II=1,N2
    IF (X_RECT(J,K).GE.X_IN(II) .AND. X_RECT(J,K).LT.X_IN(II+1)) GOTO 101
    ENDDO
101     P_RECT(J,K)=( PREF_IN(II)*(X_IN(II+1)-X_RECT(J,K))+PREF_IN(II+1)*(X_RECT(J,K)-X_IN(II)) )/(X_IN(II+1)-X_IN(II))  
    U_RECT0(:,J,K)=( UREF_IN(:,II)*(X_IN(II+1)-X_RECT(J,K))+UREF_IN(:,II+1)*(X_RECT(J,K)-X_IN(II)) )/(X_IN(II+1)-X_IN(II))
ENDDO
ENDDO

! interpolation from 2d data to square-part2&3(INLET2) data according to x coordinate
DO K=1,NK
DO J=1,NJ
    DO II=1,N2
    IF (X_CY(J,K).GE.X_IN(II) .AND. X_CY(J,K).LT.X_IN(II+1)) GOTO 102
    ENDDO
102     P_CY(J,K)=( PREF_IN(II)*(X_IN(II+1)-X_CY(J,K))+PREF_IN(II+1)*(X_CY(J,K)-X_IN(II)) )/(X_IN(II+1)-X_IN(II))  
    U_CY0(:,J,K)=( UREF_IN(:,II)*(X_IN(II+1)-X_CY(J,K))+UREF_IN(:,II+1)*(X_CY(J,K)-X_IN(II)) )/(X_IN(II+1)-X_IN(II))
ENDDO
ENDDO
    




OPEN(1,FILE='PU_RECT_SUR.plt',STATUS='UNKNOWN')
WRITE(1,237) NL,NL
DO K=1,NL
DO J=1,NL
    U1=U_RECT0(1,J,K)
    U2=U_RECT0(2,J,K)*Y_RECT(J,K)/SQRT(Y_RECT(J,K)**2+Z_RECT(J,K)**2)
    U3=U_RECT0(2,J,K)*Z_RECT(J,K)/SQRT(Y_RECT(J,K)**2+Z_RECT(J,K)**2)
    WRITE(1,238) X_RECT(J,K),Y_RECT(J,K),Z_RECT(J,K),P_RECT(J,K),U1,U2,U3
ENDDO
ENDDO
CLOSE(1) 

OPEN(20,FILE='Tecplot_InputFiles/PU_RECT_Inlet1.plt',STATUS='UNKNOWN')
WRITE(20, '(A)') 'TITLE = "PU_RECT_Inlet1"'
WRITE(20, '(A)') 'VARIABLES = "X", "Y", "Z" ,"p" , "Ux" ,"Uy", "Uz"'
WRITE(20, '(A, I0, A, I0, A)') 'ZONE T="PU_RECT_Inlet1", I=', NL, ', J=', NL, ', F=POINT'
DO K=1,NL
DO J=1,NL
    U1=U_RECT0(1,J,K)
    U2=U_RECT0(2,J,K)*Y_RECT(J,K)/SQRT(Y_RECT(J,K)**2+Z_RECT(J,K)**2)
    U3=U_RECT0(2,J,K)*Z_RECT(J,K)/SQRT(Y_RECT(J,K)**2+Z_RECT(J,K)**2)
    WRITE(20,'(7ES18.10)') X_RECT(J,K),Y_RECT(J,K),Z_RECT(J,K),P_RECT(J,K),U1,U2,U3
ENDDO
ENDDO
CLOSE(20) 



OPEN(1,FILE='PU_CY_SUR.plt',STATUS='UNKNOWN')
WRITE(1,237) NJ,NK
DO K=1,NK
DO J=1,NJ
    U1=U_CY0(1,J,K)
    U2=U_CY0(2,J,K)*Y_CY(J,K)/SQRT(Y_CY(J,K)**2+Z_CY(J,K)**2)
    U3=U_CY0(2,J,K)*Z_CY(J,K)/SQRT(Y_CY(J,K)**2+Z_CY(J,K)**2)
    WRITE(1,238) X_CY(J,K),Y_CY(J,K),Z_CY(J,K),P_CY(J,K),U1,U2,U3
ENDDO
ENDDO
CLOSE(1)


OPEN(1,FILE='Tecplot_InputFiles/PU_CY_SUR_Inlet2.plt',STATUS='UNKNOWN')
WRITE(1, '(A)') 'TITLE = "PU_CY_Inlet2"'
WRITE(1,'(A)') 'VARIABLES = "X", "Y", "Z" ,"p" , "Ux" ,"Uy", "Uz"'
WRITE(1, '(A, I0, A, I0, A)') 'ZONE T="PU_CY_Inlet2", I=', NJ+1, ', J=', NK+1, ', F=POINT'
DO K=1,NK
    IF (K.LE.NL) THEN
        J1=K
        K1=1
    ELSEIF (K.LE.(2*NL)) THEN
        J1=NL
        K1=K-NL
    ELSEIF (K.LE.(3*NL)) THEN
        J1=3*NL-K+1
        K1=NL
    ELSE
        J1=1
        K1=4*NL-K+1
    ENDIF
    U1=U_RECT0(1,J1,K1)
    U2=U_RECT0(2,J1,K1)*Y_RECT(J1,K1)/SQRT(Y_RECT(J1,K1)**2+Z_RECT(J1,K1)**2)
    U3=U_RECT0(2,J1,K1)*Z_RECT(J1,K1)/SQRT(Y_RECT(J1,K1)**2+Z_RECT(J1,K1)**2)
    WRITE(1,'(7ES18.10)') X_RECT(J1,K1),Y_RECT(J1,K1),Z_RECT(J1,K1),P_RECT(J1,K1),U1,U2,U3 
DO J=1,NJ
    U1=U_CY0(1,J,K)
    U2=U_CY0(2,J,K)*Y_CY(J,K)/SQRT(Y_CY(J,K)**2+Z_CY(J,K)**2)
    U3=U_CY0(2,J,K)*Z_CY(J,K)/SQRT(Y_CY(J,K)**2+Z_CY(J,K)**2)
    WRITE(1,'(7ES18.10)') X_CY(J,K),Y_CY(J,K),Z_CY(J,K),P_CY(J,K),U1,U2,U3
ENDDO
ENDDO
DO K=1,1
        J1=1
        K1=1
    U1=U_RECT0(1,J1,K1)
    U2=U_RECT0(2,J1,K1)*Y_RECT(J1,K1)/SQRT(Y_RECT(J1,K1)**2+Z_RECT(J1,K1)**2)
    U3=U_RECT0(2,J1,K1)*Z_RECT(J1,K1)/SQRT(Y_RECT(J1,K1)**2+Z_RECT(J1,K1)**2)
    WRITE(1,'(7ES18.10)') X_RECT(J1,K1),Y_RECT(J1,K1),Z_RECT(J1,K1),P_RECT(J1,K1),U1,U2,U3  
DO J=1,NJ
    U1=U_CY0(1,J,K)
    U2=U_CY0(2,J,K)*Y_CY(J,K)/SQRT(Y_CY(J,K)**2+Z_CY(J,K)**2)
    U3=U_CY0(2,J,K)*Z_CY(J,K)/SQRT(Y_CY(J,K)**2+Z_CY(J,K)**2)
    WRITE(1,'(7ES18.10)') X_CY(J,K),Y_CY(J,K),Z_CY(J,K),P_CY(J,K),U1,U2,U3
ENDDO
ENDDO
CLOSE(1) 

! OPEN(1,FILE='Tecplot_InputFiles/PU_CY_SUR_Inlet2_Quarter.plt',STATUS='UNKNOWN')
! WRITE(1, '(A)') 'TITLE = "PU_CY_Inlet2_Quarter"'
! WRITE(1,'(A)') 'VARIABLES = "X", "Y", "Z" ,"p" , "u" ,"v", "w"'
! WRITE(1, '(A, I0, A, I0, A)') 'ZONE T="PU_CY_Inlet2_Quarter", I=', NJ+1, ', J=', NL, ', F=POINT'
! DO K=1,NL
!     J1=K
!     K1=1
!     U1=U_RECT0(1,J1,K1)
!     U2=U_RECT0(2,J1,K1)*Y_RECT(J1,K1)/SQRT(Y_RECT(J1,K1)**2+Z_RECT(J1,K1)**2)
!     U3=U_RECT0(2,J1,K1)*Z_RECT(J1,K1)/SQRT(Y_RECT(J1,K1)**2+Z_RECT(J1,K1)**2)
!     WRITE(1,'(7ES18.10)') X_RECT(J1,K1),Y_RECT(J1,K1),Z_RECT(J1,K1),P_RECT(J1,K1),U1,U2,U3 
! DO J=1,NJ
!     U1=U_CY0(1,J,K)
!     U2=U_CY0(2,J,K)*Y_CY(J,K)/SQRT(Y_CY(J,K)**2+Z_CY(J,K)**2)
!     U3=U_CY0(2,J,K)*Z_CY(J,K)/SQRT(Y_CY(J,K)**2+Z_CY(J,K)**2)
!     WRITE(1,'(7ES18.10)') X_CY(J,K),Y_CY(J,K),Z_CY(J,K),P_CY(J,K),U1,U2,U3
! ENDDO
! ENDDO
! CLOSE(1) 





236 FORMAT(4(E21.14,2X))
237 FORMAT('Zone i= ',I0,' j= ',I0,' f=point')
238 FORMAT(7(E22.15,2X))

END PROGRAM


