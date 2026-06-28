
    PROGRAM cylinder_main
    IMPLICIT NONE
    INCLUDE 'cylinder.h'

    REAL*8 X0(M1,M2),Y0(M1,M2),THETA1,THETA2
    REAL*8 PAI
    REAL*8 X,Y,Z
    INTEGER*8 I,J,K,P_NUM
 
    PAI=ACOS(-1.0)
    
    
    P_NUM=M1*M2+M1*N2
    OPEN(1,FILE='ReadData/suboff_mesh_2d_ADD.plt',STATUS='OLD')   
    READ(1,*)
    DO I=1,M1  
    DO J=1,M2
        READ(1,*) X0(I,J),Y0(I,J)
    ENDDO
    ENDDO  
    CLOSE(1)
      
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
    
    THETA1=0.5*2.0*PAI/REAL(M3)
    THETA2=(REAL(M3)-0.5)*2.0*PAI/REAL(M3)
    
    DO J=1,1  
    DO I=1,M1
        X=X0(I,J)
        Y=0.0
        Z=0.0
        WRITE(2,230) '( ',X,' ',Y,' ',Z,' )'
    ENDDO
    ENDDO
    
    DO J=2,M2  
    DO I=1,M1
        X=X0(I,J)
        Y=Y0(I,J)*COS(THETA1)
        Z=-Y0(I,J)*SIN(THETA1)
        WRITE(2,230) '( ',X,' ',Y,' ',Z,' )'
    ENDDO
    ENDDO
    
    DO J=2,M2  
    DO I=1,M1
        X=X0(I,J)
        Y=Y0(I,J)*COS(THETA2)
        Z=-Y0(I,J)*SIN(THETA2)
        WRITE(2,230) '( ',X,' ',Y,' ',Z,' )'
    ENDDO
    ENDDO
    
    WRITE(2,*) ")"
    WRITE(2,*) " "
    WRITE(2,*) " "
    WRITE(2,*) "// ************************************************************************* //"
    CLOSE(2)
    
    
230 FORMAT(A2,E21.14,A1,E21.14,A1,E21.14,A2)     
231 FORMAT(I0)   
232 FORMAT(A24,I0,A1) 
233 FORMAT(A26,I0,A9,I0,A9,I0,A17,I0,A2)   
234 FORMAT(A2,I0,A1,I0,A1,I0,A1) 
235 FORMAT(A2,I0,A1,I0,A1,I0,A1,I0,A1)

    END PROGRAM


