
    PROGRAM cylinder_main
    IMPLICIT NONE
    INCLUDE 'parameter.h'


    REAL*8 X0(M1,M2), Y0(M1,M2), THETA(M3), DTHETA
    REAL*8 PAI, RATIO_TAR
    REAL*8 X, Y, Z
    INTEGER*8 I, J, K, L, K_START, K_END,L_ALL
    CHARACTER(LEN=64) :: filename

    


    PAI=ACOS(-1.0)

   
    
    
    OPEN(1,FILE='ReadData/suboff_mesh_2d.plt',STATUS='OLD')   
    READ(1,*)
    DO I=1,M1  
    DO J=1,M2
        READ(1,*) X0(I,J),Y0(I,J)
    ENDDO
    ENDDO  
    CLOSE(1)
    
    DTHETA = 2 * PAI / DBLE(NK)
    THETA(1)= -PAI/4.0 + DTHETA/2
    DO K = 2, NK
        THETA(K) = THETA(K-1) + DTHETA
    END DO



    L_ALL=4

    DO L = 0, L_ALL-1
        ! Dynamically create the filename (e.g., '../sampleBL0.c', '../sampleBL1.c', etc.)
        WRITE(filename, '(A,I1,A)') '../sampleBow', L+1, '.c'
        
        ! Open the file
        OPEN(2, FILE=TRIM(filename), STATUS='REPLACE',ACTION='WRITE')
        
        ! Calculate the K range for the current file
        K_START = L * (N3 / L_ALL) + 1
        K_END   = (L + 1) * (N3 / L_ALL)
        
        DO K = K_START, K_END
            DO J = 1, NL/2 + NJ1
            ! Use the dynamically calculated K_START and K_END      
                DO I = 1, N1      
                    X=0.25*( X0(I,J) + X0(I,J+1) + X0(I+1,J) + X0(I+1,J+1) )
                    Y=0.25*( Y0(I,J) + Y0(I,J+1) + Y0(I+1,J) + Y0(I+1,J+1) )*COS(THETA(K))
                    Z=0.25*( Y0(I,J) + Y0(I,J+1) + Y0(I+1,J) + Y0(I+1,J+1) )*SIN(THETA(K))
                    WRITE(2,*) "(", X, " ", Y, " ", Z, ")"                 
                ENDDO       
            ENDDO       
        ENDDO       
        CLOSE(2)

        WRITE(filename, '(A,I1)') '../sampleBow', L+1
        OPEN(3, FILE=TRIM(ADJUSTL(filename)), STATUS='REPLACE', ACTION='WRITE')
        
        WRITE(3,'(A)') "/*--------------------------------*- C++ -*----------------------------------*\"
        WRITE(3,'(A)') "| =========                 |                                                 |"
        WRITE(3,'(A)') "| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |"
        WRITE(3,'(A)') "|  \\    /   O peration     | Version:  v2312                                 |"
        WRITE(3,'(A)') "|   \\  /    A nd           | Website:  www.openfoam.com                      |"
        WRITE(3,'(A)') "|    \\/     M anipulation  |                                                 |"
        WRITE(3,'(A)') "\*---------------------------------------------------------------------------*/"
        WRITE(3,'(A)') "FoamFile"
        WRITE(3,'(A)') "{"
        WRITE(3,'(A)') "    version     2.0;"
        WRITE(3,'(A)') "    format      ascii;"
        WRITE(3,'(A)') "    class       dictionary;"
        WRITE(3,'(A)') "    object      sampleBow;"
        WRITE(3,'(A)') "}"
        WRITE(3,'(A)') "// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //"
        WRITE(3,'(A)') ""
        WRITE(3,'(A)') "type            sets;"
        WRITE(3,'(A)') ""
        WRITE(3,'(A)') "libs            (sampling);"
        WRITE(3,'(A)') ""
        WRITE(3,'(A)') "interpolationScheme cellPoint;"
        WRITE(3,'(A)') ""
        WRITE(3,'(A)') "setFormat       raw;"
        WRITE(3,'(A)') ""
        WRITE(3,'(A)') "fields          (U p);"
        WRITE(3,'(A)') ""
        WRITE(3,'(A)') "sets"
        WRITE(3,'(A)') "("
        WRITE(3,'(A)') "    cloud"
        WRITE(3,'(A)') "    {"
        WRITE(3,'(A)') "        type    cloud;"
        WRITE(3,'(A)') "        axis    xyz;"
        WRITE(3,'(A)') "        points"
        WRITE(3,'(A)') "        ("
        ! Insert the dynamically generated #include statement
        WRITE(3,'(A,I1,A)') '            #include     "sampleBow', L+1, '.c"'
        WRITE(3,'(A)') "        );"
        WRITE(3,'(A)') "    }"
        WRITE(3,'(A)') ");"
        
        CLOSE(3)
    ENDDO

    
    END PROGRAM


