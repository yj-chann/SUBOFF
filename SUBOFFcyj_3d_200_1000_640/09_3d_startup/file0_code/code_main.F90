!-------------------------------------------------------------------------
!  MAIN PROGRAM
!-------------------------------------------------------------------------

    SUBROUTINE MAIN_CODE
    
!----------------------------------------
    
    INCLUDE 'head.fi'
    

    REAL PAI
    REAL UU1,UU2,UU3
    REAL X,Y,Z
    INTEGER NIF,NIF1,NIF2
    INTEGER I,J,K
    INTEGER PN,MI,MJ,MJ1,MJ2
    INTEGER P_NUM,F_NUM,F_INTERNAL_NUM,F_INTERNAL_NUM_ALL
    INTEGER F11,F12,F13,F14,F21,F22,F23,F24
    INTEGER FAD1,FAD2,FAD3,OWN,NEIGH,PRO1,PRO2
    INTEGER IP,JP,KP,N
    INTEGER JISHU,NCENTER_ALL,CELL_ID_ALL,NCENTER,I1
    INTEGER CELL_ID_CY(CNI,CNJ,0:NK+1),CELL_ID2_CY(0:C_NUM_CY-1)
    INTEGER CELL_ID_RECT(CNI,CNL,CNL),CELL_ID2_RECT(0:C_NUM_RECT-1)
    INTEGER NUM(7),INDEX_S(7)
    INTEGER FACE_ID,IF_OWN,IF_NEIGH
    INTEGER I_NEIGH,J_NEIGH,K_NEIGH,IP_NEIGH,JP_NEIGH,KP_NEIGH
    INTEGER F_LEFT_CY(CNJ,NK),F_RIGHT_CY(CNJ,NK),F_UP_CY(CNI,NK),F_DOWN_CY(CNI,NK)
    INTEGER F_LEFT_RECT(CNL,CNL),F_RIGHT_RECT(CNL,CNL),F_UP_RECT(CNI,CNL),F_DOWN_RECT(CNI,CNL),F_FRONT_RECT(CNI,CNL),F_BACK_RECT(CNI,CNL)
    INTEGER IF_CY,NN
    CHARACTER*20 CHAR1
    CHARACTER*90000 DIR,TIME,CHAR3
 
    IF_CY=1 !INITIALIZE THE VALUE
    MJ1=0 !INITIALIZE THE VALUE
    MJ2=0 !INITIALIZE THE VALUE
    PAI=ACOS(-1.0)
    !PROCESSOR NUMBER
    MI=COORD(0)+1
    MJ=COORD(1)+1
    IF (MJ.LE.PNJ) THEN
    PN=PNI*(MJ-1)+MI-1
    ELSE
    MJ1=( (COORD(1)-PNJ)-MOD(COORD(1)-PNJ,PNL) )/PNL+1
    MJ2=MJ-PNJ-(MJ1-1)*PNL
    PN=PNI*PNJ+PNL*PNI*(MJ2-1)+PNI*(MJ1-1)+MI-1
    ENDIF
    !POINTS NUMBER
    IF (MJ.LE.PNJ) THEN
    F_NUM=3*C_NUM_CY+CNJ*NK+CNI*NK
    F_INTERNAL_NUM=3*C_NUM_CY-CNJ*NK-CNI*NK
    ELSE
    F_NUM=3*C_NUM_RECT+CNL*CNL+CNL*CNI+CNL*CNI
    F_INTERNAL_NUM=3*C_NUM_RECT-CNL*CNL-CNL*CNI-CNL*CNI
    ENDIF
    F_INTERNAL_NUM_ALL=3*NI*NJ*NK-NJ*NK-NI*NK + 3*NI*NL*NL-NL*NL-NL*NI-NL*NI + NK*NI
    
    !WRITE FIELDS LOCATION
    LOCATION='/0/'
    
    WRITE(CHAR1,231) PN    
    !cellProcAddressing
    OPEN(2,FILE='../processor'//TRIM(ADJUSTL(CHAR1))//'/constant/polyMesh/cellProcAddressing',STATUS='UNKNOWN')
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
    WRITE(2,*) "    class       labelList;"
    WRITE(2,*) '    location    "constant/polyMesh";'
    WRITE(2,*) "    object      cellProcAddressing;"
    WRITE(2,*) "}"
    WRITE(2,*) "// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //"
    WRITE(2,*) " "
    WRITE(2,*) " "
    IF (MJ.LE.PNJ) THEN
    WRITE(2,231) C_NUM_CY
    ELSE
    WRITE(2,231) C_NUM_RECT    
    ENDIF  
    WRITE(2,*) "("
    
    IF (MJ.LE.PNJ) THEN
        
    JISHU=0
    DO NCENTER_ALL=0,C_NUM_ALL-1
    !CALL GET_INDEX(RENUM2(NCENTER_ALL),I,J,K,IF_CY)
    CALL GET_INDEX(NCENTER_ALL,I,J,K,IF_CY)
    IF (IF_CY.EQ.1) THEN
    IP=I-CNI*(MI-1)
    JP=J-CNJ*(MJ-1)
    KP=K
    IF (IP.LE.CNI.AND.IP.GE.1) THEN
    IF (JP.LE.CNJ.AND.JP.GE.1) THEN
        CELL_ID_ALL=NJ*NI*(K-1)+NI*(J-1)+I-1
        CELL_ID2_CY(JISHU)=CELL_ID_ALL
        CELL_ID_CY(IP,JP,KP)=JISHU
        JISHU=JISHU+1
        WRITE(2,231) NCENTER_ALL
        !IF (PN.EQ.0) THEN 
        !    WRITE(*,*) 'TEST01  ',I,J,K
        !    WRITE(*,*) 'TEST01  ',IP,JP,KP
        !    WRITE(*,*) 'TEST01  ',CELL_ID_ALL
        !ENDIF  
    ENDIF
    ENDIF
    ENDIF
    ENDDO
    CELL_ID_CY(:,:,NK+1)=CELL_ID_CY(:,:,1)
    CELL_ID_CY(:,:,0)=CELL_ID_CY(:,:,NK)
    
    ELSE
    
    JISHU=0
    DO NCENTER_ALL=0,C_NUM_ALL-1
    !CALL GET_INDEX(RENUM2(NCENTER_ALL),I,J,K,IF_CY)
    CALL GET_INDEX(NCENTER_ALL,I,J,K,IF_CY)
    IF (IF_CY.EQ.0) THEN
    IP=I-CNI*(MI-1)
    JP=J-CNL*(MJ1-1)
    KP=K-CNL*(MJ2-1)
    IF (IP.LE.CNI.AND.IP.GE.1) THEN
    IF (JP.LE.CNL.AND.JP.GE.1) THEN
    IF (KP.LE.CNL.AND.KP.GE.1) THEN
        CELL_ID_ALL=NI*NL*(K-1)+NI*(J-1)+I-1 +NI*NJ*NK
        CELL_ID2_RECT(JISHU)=CELL_ID_ALL
        CELL_ID_RECT(IP,JP,KP)=JISHU
        JISHU=JISHU+1
        WRITE(2,231) NCENTER_ALL
        !IF (PN.EQ.0) THEN 
        !    WRITE(*,*) 'TEST01  ',I,J,K
        !    WRITE(*,*) 'TEST01  ',IP,JP,KP
        !    WRITE(*,*) 'TEST01  ',CELL_ID_ALL
        !ENDIF  
    ENDIF
    ENDIF
    ENDIF
    ENDIF
    ENDDO    
        
        
    ENDIF
    
    !!TESTINGɾ��
    !IF (PN.EQ.40) THEN
    !    WRITE(*,*) CELL_ID_RECT
    !ENDIF
    
    WRITE(2,*) ")"
    WRITE(2,*) " "
    WRITE(2,*) " "
    WRITE(2,*) "// ************************************************************************* //"
    CLOSE(2)


    !faces
    OPEN(2,FILE='../processor'//TRIM(ADJUSTL(CHAR1))//'/constant/polyMesh/faces',STATUS='UNKNOWN')
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
    WRITE(2,*) "    class       faceList;"
    WRITE(2,*) '    location    "constant/polyMesh";'
    WRITE(2,*) "    object      faces;"
    WRITE(2,*) "}"
    WRITE(2,*) "// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //"
    WRITE(2,*) " "
    WRITE(2,*) " "
    WRITE(2,231) F_NUM
    WRITE(2,*) "("
    !faceProcAddressing
    OPEN(3,FILE='../processor'//TRIM(ADJUSTL(CHAR1))//'/constant/polyMesh/faceProcAddressing',STATUS='UNKNOWN')
    WRITE(3,*) "/*--------------------------------*- C++ -*----------------------------------*\"
    WRITE(3,*) "  =========                 |"
    WRITE(3,*) "  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox"
    WRITE(3,*) "   \\    /   O peration     | Website:  https://openfoam.org"
    WRITE(3,*) "    \\  /    A nd           | Version:  v2312"
    WRITE(3,*) "     \\/     M anipulation  |"
    WRITE(3,*) "\*---------------------------------------------------------------------------*/"
    WRITE(3,*) "FoamFile"
    WRITE(3,*) "{"
    WRITE(3,*) "    version     2.0;"
    WRITE(3,*) "    format      ascii;"
    WRITE(3,*) "    class       labelList;"
    WRITE(3,*) '    location    "constant/polyMesh";'
    WRITE(3,*) "    object      faceProcAddressing;"
    WRITE(3,*) "}"
    WRITE(3,*) "// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //"
    WRITE(3,*) " "      
    WRITE(3,*) " "      
    WRITE(3,231) F_NUM  
    WRITE(3,*) "("      
    !owner              
    OPEN(4,FILE='../processor'//TRIM(ADJUSTL(CHAR1))//'/constant/polyMesh/owner',STATUS='UNKNOWN')
    WRITE(4,*) "/*--------------------------------*- C++ -*----------------------------------*\"
    WRITE(4,*) "  =========                 |"
    WRITE(4,*) "  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox"
    WRITE(4,*) "   \\    /   O peration     | Website:  https://openfoam.org"
    WRITE(4,*) "    \\  /    A nd           | Version:  v2312"
    WRITE(4,*) "     \\/     M anipulation  |"
    WRITE(4,*) "\*---------------------------------------------------------------------------*/"
    WRITE(4,*) "FoamFile"
    WRITE(4,*) "{"
    WRITE(4,*) "    version     2.0;"
    WRITE(4,*) "    format      ascii;"
    WRITE(4,*) "    class       labelList;"
    IF (MJ.LE.PNJ) THEN
    WRITE(4,233) '    note        "nPoints: ',P_NUM,' nCells: ',C_NUM_CY,' nFaces: ',F_NUM,' nInternalFaces: ',F_INTERNAL_NUM,'";'
    ELSE
    WRITE(4,233) '    note        "nPoints: ',P_NUM,' nCells: ',C_NUM_RECT,' nFaces: ',F_NUM,' nInternalFaces: ',F_INTERNAL_NUM,'";'
    ENDIF
    WRITE(4,*) '    location    "constant/polyMesh";'
    WRITE(4,*) "    object      owner;"
    WRITE(4,*) "}"
    WRITE(4,*) "// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //"
    WRITE(4,*) " "
    WRITE(4,*) " "
    WRITE(4,231) F_NUM
    WRITE(4,*) "("
    !neighbour
    OPEN(5,FILE='../processor'//TRIM(ADJUSTL(CHAR1))//'/constant/polyMesh/neighbour',STATUS='UNKNOWN')
    WRITE(5,*) "/*--------------------------------*- C++ -*----------------------------------*\"
    WRITE(5,*) "  =========                 |"
    WRITE(5,*) "  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox"
    WRITE(5,*) "   \\    /   O peration     | Website:  https://openfoam.org"
    WRITE(5,*) "    \\  /    A nd           | Version:  v2312"
    WRITE(5,*) "     \\/     M anipulation  |"
    WRITE(5,*) "\*---------------------------------------------------------------------------*/"
    WRITE(5,*) "FoamFile"
    WRITE(5,*) "{"
    WRITE(5,*) "    version     2.0;"
    WRITE(5,*) "    format      ascii;"
    WRITE(5,*) "    class       labelList;"
    IF (MJ.LE.PNJ) THEN
    WRITE(5,233) '    note        "nPoints: ',P_NUM,' nCells: ',C_NUM_CY,' nFaces: ',F_NUM,' nInternalFaces: ',F_INTERNAL_NUM,'";'
    ELSE
    WRITE(5,233) '    note        "nPoints: ',P_NUM,' nCells: ',C_NUM_RECT,' nFaces: ',F_NUM,' nInternalFaces: ',F_INTERNAL_NUM,'";'
    ENDIF
    WRITE(5,*) '    location    "constant/polyMesh";'
    WRITE(5,*) "    object      neighbour;"
    WRITE(5,*) "}"
    WRITE(5,*) "// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //"
    WRITE(5,*) " "
    WRITE(5,*) " "
    WRITE(5,231) F_INTERNAL_NUM
    WRITE(5,*) "("
    
    
    
    !boundary
    OPEN(21,FILE='../processor'//TRIM(ADJUSTL(CHAR1))//'/constant/polyMesh/boundary',STATUS='UNKNOWN')
    WRITE(21,*) "/*--------------------------------*- C++ -*----------------------------------*\"
    WRITE(21,*) "  =========                 |"
    WRITE(21,*) "  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox"
    WRITE(21,*) "   \\    /   O peration     | Website:  https://openfoam.org"
    WRITE(21,*) "    \\  /    A nd           | Version:  v2312"
    WRITE(21,*) "     \\/     M anipulation  |"
    WRITE(21,*) "\*---------------------------------------------------------------------------*/"
    WRITE(21,*) "FoamFile"
    WRITE(21,*) "{"
    WRITE(21,*) "    version     2.0;"
    WRITE(21,*) "    format      ascii;"
    WRITE(21,*) "    class       polyBoundaryMesh;"
    WRITE(21,*) '    location    "constant/polyMesh";'
    WRITE(21,*) "    object      boundary;"
    WRITE(21,*) "}"
    WRITE(21,*) "// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //"
    WRITE(21,*) " "
    !boundaryProcAddressing
    OPEN(31,FILE='../processor'//TRIM(ADJUSTL(CHAR1))//'/constant/polyMesh/boundaryProcAddressing',STATUS='UNKNOWN')
    WRITE(31,*) "/*--------------------------------*- C++ -*----------------------------------*\"
    WRITE(31,*) "  =========                 |"
    WRITE(31,*) "  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox"
    WRITE(31,*) "   \\    /   O peration     | Website:  https://openfoam.org"
    WRITE(31,*) "    \\  /    A nd           | Version:  v2312"
    WRITE(31,*) "     \\/     M anipulation  |"
    WRITE(31,*) "\*---------------------------------------------------------------------------*/"
    WRITE(31,*) "FoamFile"
    WRITE(31,*) "{"
    WRITE(31,*) "    version     2.0;"
    WRITE(31,*) "    format      ascii;"
    WRITE(31,*) "    class       labelList;"
    WRITE(31,*) '    location    "constant/polyMesh";'
    WRITE(31,*) "    object      boundaryProcAddressing;"
    WRITE(31,*) "}"
    WRITE(31,*) "// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //"
    WRITE(31,*) " "
    !p
    OPEN(41,FILE='../processor'//TRIM(ADJUSTL(CHAR1))//TRIM(ADJUSTL(LOCATION))//'p',STATUS='UNKNOWN')
    WRITE(41,*) "/*--------------------------------*- C++ -*----------------------------------*\"
    WRITE(41,*) "  =========                 |"
    WRITE(41,*) "  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox"
    WRITE(41,*) "   \\    /   O peration     | Website:  https://openfoam.org"
    WRITE(41,*) "    \\  /    A nd           | Version:  v2312"
    WRITE(41,*) "     \\/     M anipulation  |"
    WRITE(41,*) "\*---------------------------------------------------------------------------*/"
    WRITE(41,*) "FoamFile"
    WRITE(41,*) "{"
    WRITE(41,*) "    version     2.0;"
    WRITE(41,*) "    format      ascii;"
    WRITE(41,*) "    class       volScalarField;"
    WRITE(41,*) '    location    "0";'
    WRITE(41,*) "    object      p;"
    WRITE(41,*) "}"
    WRITE(41,*) "// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //"
    WRITE(41,*) " "
    WRITE(41,*) "dimensions      [0 2 -2 0 0 0 0];"
    WRITE(41,*) " "
    WRITE(41,*) "internalField   nonuniform List<scalar> "
    IF (MJ.LE.PNJ) THEN
    WRITE(41,231) C_NUM_CY
    ELSE
    WRITE(41,231) C_NUM_RECT  
    ENDIF
    WRITE(41,*) "("
    IF (MJ.LE.PNJ) THEN
    DO I1=0,C_NUM_CY-1
    CALL GET_INDEX(CELL_ID2_CY(I1),I,J,K,IF_CY)
    IP=I-CNI*(MI-1)
    JP=J-CNJ*(MJ-1)
    KP=K
    WRITE(41,'(E21.14)') PREF_CY(IP,JP,KP) 
    ENDDO
    ELSE
    DO I1=0,C_NUM_RECT-1
    CALL GET_INDEX(CELL_ID2_RECT(I1),I,J,K,IF_CY)
    IP=I-CNI*(MI-1)
    JP=J-CNL*(MJ1-1)
    KP=K-CNL*(MJ2-1)
    WRITE(41,'(E21.14)') PREF_RECT(IP,JP,KP) 
    ENDDO    
    ENDIF
    WRITE(41,*) ")"
    WRITE(41,*) ";"
    WRITE(41,*) " "
    WRITE(41,*) "boundaryField"
    WRITE(41,*) "{"
    WRITE(41,*) "    INLET1"
    WRITE(41,*) "    {"
    WRITE(41,*) "        type            zeroGradient;"
    WRITE(41,*) "    }"
    WRITE(41,*) "    INLET2"
    WRITE(41,*) "    {"
    WRITE(41,*) "        type            zeroGradient;"
    WRITE(41,*) "    }"
    WRITE(41,*) "    OUTLET"
    WRITE(41,*) "    {"
    WRITE(41,*) "        type            fixedValue;"
    IF (MJ.EQ.PNJ) THEN
    WRITE(41,*) "        value           uniform 0;"
    ELSE
    WRITE(41,*) "        value           nonuniform List<scalar> 0();"    
    ENDIF
    WRITE(41,*) "    }"
    WRITE(41,*) "    HULL"
    WRITE(41,*) "    {"
    WRITE(41,*) "        type            zeroGradient;"
    WRITE(41,*) "    }"
    
    !U
    OPEN(51,FILE='../processor'//TRIM(ADJUSTL(CHAR1))//TRIM(ADJUSTL(LOCATION))//'U',STATUS='UNKNOWN')
    WRITE(51,*) "/*--------------------------------*- C++ -*----------------------------------*\"
    WRITE(51,*) "  =========                 |"
    WRITE(51,*) "  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox"
    WRITE(51,*) "   \\    /   O peration     | Website:  https://openfoam.org"
    WRITE(51,*) "    \\  /    A nd           | Version:  v2312"
    WRITE(51,*) "     \\/     M anipulation  |"
    WRITE(51,*) "\*---------------------------------------------------------------------------*/"
    WRITE(51,*) "FoamFile"
    WRITE(51,*) "{"
    WRITE(51,*) "    version     2.0;"
    WRITE(51,*) "    format      ascii;"
    WRITE(51,*) "    class       volVectorField;"
    WRITE(51,*) '    location    "0";'
    WRITE(51,*) "    object      U;"
    WRITE(51,*) "}"
    WRITE(51,*) "// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //"
    WRITE(51,*) " "
    WRITE(51,*) "dimensions      [0 1 -1 0 0 0 0];"
    WRITE(51,*) " "
    WRITE(51,*) "internalField   nonuniform List<vector> "
    IF (MJ.LE.PNJ) THEN
    WRITE(51,231) C_NUM_CY
    ELSE
    WRITE(51,231) C_NUM_RECT  
    ENDIF
    WRITE(51,*) "("
    IF (MJ.LE.PNJ) THEN
    DO I1=0,C_NUM_CY-1
    CALL GET_INDEX(CELL_ID2_CY(I1),I,J,K,IF_CY)
    IP=I-CNI*(MI-1)
    JP=J-CNJ*(MJ-1)
    KP=K
    WRITE(51,230) '( ',UREF_CY(1,IP,JP,KP),' ',UREF_CY(2,IP,JP,KP),' ',UREF_CY(3,IP,JP,KP),' )'
    ENDDO
    ELSE
    DO I1=0,C_NUM_RECT-1
    CALL GET_INDEX(CELL_ID2_RECT(I1),I,J,K,IF_CY)
    IP=I-CNI*(MI-1)
    JP=J-CNL*(MJ1-1)
    KP=K-CNL*(MJ2-1)
    WRITE(51,230) '( ',UREF_RECT(1,IP,JP,KP),' ',UREF_RECT(2,IP,JP,KP),' ',UREF_RECT(3,IP,JP,KP),' )'
    ENDDO       
    ENDIF       
    WRITE(51,*) ")"
    WRITE(51,*) ";"
    WRITE(51,*) " "
    WRITE(51,*) "boundaryField"
    WRITE(51,*) "{"
    WRITE(51,*) "    INLET1"
    WRITE(51,*) "    {"
    WRITE(51,*) "        type            timeVaryingMappedFixedValue;"
    WRITE(51,*) "        fieldTable      U;"
    WRITE(51,*) "        offset          constant ( 0 0 0 );"
    !WRITE(51,*) "        type            fixedValue;"
    IF (MI.EQ.1 .AND. MJ.GT.PNJ) THEN
    N=CNL*CNL
    WRITE(51,*) "        value           nonuniform List<vector> "  
    WRITE(51,'(I0)') N
    WRITE(51,*) '('
    DO KP=1,CNL
    DO JP=1,CNL
    !UU1=U_RIGHT_RECT(1,JP,KP)
    !UU2=U_RIGHT_RECT(2,JP,KP)
    !UU3=U_RIGHT_RECT(3,JP,KP)
    UU1=U_INLET_RECT(1,JP,KP)
    UU2=U_INLET_RECT(2,JP,KP)
    UU3=U_INLET_RECT(3,JP,KP)
    WRITE(51,230) '( ',UU1,' ',UU2,' ',UU3,' )'
    ENDDO
    ENDDO
    WRITE(51,*) ')'
    WRITE(51,*) ';'
    ELSE 
    WRITE(51,*) "        value           nonuniform List<vector> 0();" 
    ENDIF
    WRITE(51,*) '    }'
    WRITE(51,*) "    INLET2"
    WRITE(51,*) "    {"
    WRITE(51,*) "        type            fixedValue;"
    IF (MI.EQ.1 .AND. MJ.LE.PNJ) THEN
    N=CNJ*NK
    WRITE(51,*) "        value           nonuniform List<vector> "  
    WRITE(51,'(I0)') N
    WRITE(51,*) '('
    DO KP=1,NK
    DO JP=1,CNJ
    !UU1=U_RIGHT_CY(1,JP,KP)
    !UU2=U_RIGHT_CY(2,JP,KP)
    !UU3=U_RIGHT_CY(3,JP,KP)
    UU1=U_INLET_CY(1,JP,KP)
    UU2=U_INLET_CY(2,JP,KP)
    UU3=U_INLET_CY(3,JP,KP)
    WRITE(51,230) '( ',UU1,' ',UU2,' ',UU3,' )'
    ENDDO
    ENDDO
    WRITE(51,*) ')'
    WRITE(51,*) ';'
    ELSE 
    WRITE(51,*) "        value           nonuniform List<vector> 0();" 
    ENDIF
    WRITE(51,*) '    }'
    WRITE(51,*) "    OUTLET"
    WRITE(51,*) "    {"
    WRITE(51,*) "        type            inletOutlet;"
    IF (MJ.EQ.PNJ) THEN
    WRITE(51,*) "        inletValue      uniform (0 0 0);"
    WRITE(51,*) "        value           nonuniform List<vector> "
    N=CNI*NK
    WRITE(51,'(I0)') N
    WRITE(51,*) '('
    DO KP=1,NK
    K=KP
    DO IP=1,CNI
    UU1=U_UP_CY(1,IP,KP)
    UU2=U_UP_CY(2,IP,KP)
    UU3=U_UP_CY(3,IP,KP)
    WRITE(51,230) '( ',UU1,' ',UU2,' ',UU3,' )'
    ENDDO
    ENDDO
    WRITE(51,*) ')'
    WRITE(51,*) ';'
    ELSE
    WRITE(51,*) "        inletValue      nonuniform List<vector> 0();"
    WRITE(51,*) "        value           nonuniform List<vector> 0();"   
    ENDIF
    WRITE(51,*) '    }'
    WRITE(51,*) "    HULL"
    WRITE(51,*) "    {"
    WRITE(51,*) "        type            fixedValue;"
    IF (MI.EQ.PNI) THEN
    WRITE(51,*) "        value           uniform (0 0 0);"
    ELSE
    WRITE(51,*) "        value           nonuniform List<vector> 0();"    
    ENDIF
    WRITE(51,*) "    }"
    
    !nut
    OPEN(61,FILE='../processor'//TRIM(ADJUSTL(CHAR1))//TRIM(ADJUSTL(LOCATION))//'nut',STATUS='UNKNOWN')
    WRITE(61,*) "/*--------------------------------*- C++ -*----------------------------------*\"
    WRITE(61,*) "  =========                 |"
    WRITE(61,*) "  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox"
    WRITE(61,*) "   \\    /   O peration     | Website:  https://openfoam.org"
    WRITE(61,*) "    \\  /    A nd           | Version:  v2312"
    WRITE(61,*) "     \\/     M anipulation  |"
    WRITE(61,*) "\*---------------------------------------------------------------------------*/"
    WRITE(61,*) "FoamFile"
    WRITE(61,*) "{"
    WRITE(61,*) "    version     2.0;"
    WRITE(61,*) "    format      ascii;"
    WRITE(61,*) "    class       volScalarField;"
    WRITE(61,*) '    location    "0";'
    WRITE(61,*) "    object      nut;"
    WRITE(61,*) "}"
    WRITE(61,*) "// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //"
    WRITE(61,*) " "
    WRITE(61,*) "dimensions      [0 2 -1 0 0 0 0];"
    WRITE(61,*) " "
    WRITE(61,*) "internalField   uniform 0;"
    WRITE(61,*) " "
    WRITE(61,*) "boundaryField"
    WRITE(61,*) "{"
    WRITE(61,*) "    INLET1"
    WRITE(61,*) "    {"
    WRITE(61,*) "        type            calculated;"
    IF (MI.EQ.1 .AND. MJ.GT.PNJ) THEN
    WRITE(61,*) "        value           uniform 0;"
    ELSE
    WRITE(61,*) "        value           nonuniform List<scalar> 0();"
    ENDIF
    WRITE(61,*) "    }" 
    WRITE(61,*) "    INLET2"
    WRITE(61,*) "    {"
    WRITE(61,*) "        type            calculated;"
    IF (MI.EQ.1 .AND. MJ.LE.PNJ) THEN
    WRITE(61,*) "        value           uniform 0;"
    ELSE
    WRITE(61,*) "        value           nonuniform List<scalar> 0();"
    ENDIF
    WRITE(61,*) "    }" 
    WRITE(61,*) "    OUTLET"
    WRITE(61,*) "    {"
    WRITE(61,*) "        type            calculated;"
    IF (MJ.EQ.PNJ) THEN
    WRITE(61,*) "        value           uniform 0;"
    ELSE
    WRITE(61,*) "        value           nonuniform List<scalar> 0();"
    ENDIF
    WRITE(61,*) "    }" 
    WRITE(61,*) "    HULL"
    WRITE(61,*) "    {"
    WRITE(61,*) "        type            calculated;"
    IF (MI.EQ.PNI) THEN
    WRITE(61,*) "        value           uniform 0;"
    ELSE
    WRITE(61,*) "        value           nonuniform List<scalar> 0();"
    ENDIF
    WRITE(61,*) "    }" 
    
    
    
    
    
    
    
    !BACK FACE
    !F23**********F24
    ! *            *
    ! *            *
    ! *            *
    ! *            *
    !F22**********F21
    
    !FRONT FACE
    !F13**********F14
    ! *            *
    ! *            *
    ! *            *
    ! *            *
    !F12**********F11
    
    !block direction: i-F11-F12,j-F11-F14,k-F11-F21
    
    !FRONT: 4(F11 F14 F13 F12)
    !BACK:  4(F21 F22 F23 F24)
    !LEFT:  4(F12 F13 F23 F22)
    !RIGHT: 4(F11 F21 F24 F14)
    !UP:    4(F14 F24 F23 F13)
    !DOWN:  4(F11 F12 F22 F21)
    
    !NUM(7)=CELL_ID(IP,JP,KP)
    !NUM(1)=CELL_ID(IP,JP,KP-1) !FRONT FACE
    !NUM(2)=CELL_ID(IP,JP,KP+1) !BACK  FACE
    !NUM(3)=CELL_ID(IP+1,JP,KP) !LEFT  FACE
    !NUM(4)=CELL_ID(IP-1,JP,KP) !RIGHT FACE
    !NUM(5)=CELL_ID(IP,JP+1,KP) !UP    FACE
    !NUM(6)=CELL_ID(IP,JP-1,KP) !DOWN  FACE
    
    
    !!FACE'S NUMBER FOR CYLINDER (ALL SITUATIONS)
    !F11=(CNJ+1)*(CNI+1)*(KP-1)+(CNI+1)*(JP-1)+IP-1
    !F12=F11+1
    !F13=F12+CNI+1
    !F14=F13-1
    !F21=F11+(CNJ+1)*(CNI+1)
    !IF (KP.EQ.NK) F21=(CNI+1)*(JP-1)+IP-1
    !F22=F21+1
    !F23=F22+CNI+1
    !F24=F23-1
    
    !!FACE'S NUMBER FOR RECTANGLE (ALL SITUATIONS)
    !F11=(CNL+1)*(CNI+1)*(KP-1)+(CNI+1)*(JP-1)+IP-1
    !F12=F11+1
    !F13=F12+CNI+1
    !F14=F13-1
    !F21=F11+(CNL+1)*(CNI+1)
    !F22=F21+1
    !F23=F22+CNI+1
    !F24=F23-1
    
    !**********************************************************
    
    !* * * * * * * * * *     *
    !* 3-3 * 2-3 * 1-3 *   * * * j direction
    !* * * * * * * * * *     *
    !* 3-2 * 2-2 * 1-2 *     *
    !* * * * * * * * * *     *
    !* 3-1 * 2-1 * 1-1 *     *
    !* * * * * * * * * *     *
    !                        *
    !  *                     *
    !* * * * * * * * * * * * *
    !  *                     
    !  i direction
    
    !RECTANGLE:
    !* * * * * * * * * *     *
    !* 3-3 * 3-2 * 3-1 *   * * * j direction
    !* * * * * * * * * *     *
    !* 2-3 * 2-2 * 2-1 *     *
    !* * * * * * * * * *     *
    !* 1-3 * 1-2 * 1-1 *     *
    !* * * * * * * * * *     *
    !                        *
    !  *                     *
    !* * * * * * * * * * * * *
    !  *                     
    !  k direction
    
    If ( iMPI_MyID.EQ.0 ) THEN
    write(*,*) 'test 02'
    ENDIF


    DO NCENTER_ALL=0,C_NUM_ALL-1
    !CALL GET_INDEX(RENUM2(NCENTER_ALL),I,J,K,IF_CY)
    CALL GET_INDEX(NCENTER_ALL,I,J,K,IF_CY)
    IF (IF_CY.EQ.1) THEN
!*******************CYLINDER***************
    IP=I-CNI*(MI-1)
    JP=J-CNJ*(MJ-1)
    KP=K
    IF (IP.LE.CNI.AND.IP.GE.1) THEN
    IF (JP.LE.CNJ.AND.JP.GE.1) THEN
    NCENTER=CELL_ID_CY(IP,JP,KP)   
    IF (IP.EQ.1) THEN 
    IF (JP.EQ.1) THEN 
    !CASE 1-1  
    F11=(CNJ+1)*(CNI+1)*(KP-1)+(CNI+1)*(JP-1)+IP-1
    F12=F11+1
    F13=F12+CNI+1
    F14=F13-1
    F21=F11+(CNJ+1)*(CNI+1)
    IF (KP.EQ.NK) F21=(CNI+1)*(JP-1)+IP-1
    F22=F21+1
    F23=F22+CNI+1
    F24=F23-1
    CALL CAL_FIJ_CY(F11,F12,F13,F14,F21,F22,F23,F24)
    NUM(7)=CELL_ID_CY(IP,JP,KP)
    NUM(1)=CELL_ID_CY(IP,JP,KP-1) !FRONT FACE
    NUM(2)=CELL_ID_CY(IP,JP,KP+1) !BACK  FACE
    NUM(3)=CELL_ID_CY(IP+1,JP,KP) !LEFT  FACE
    NUM(4)=-1 !RIGHT FACE
    NUM(5)=CELL_ID_CY(IP,JP+1,KP) !UP    FACE
    NUM(6)=-1 !DOWN  FACE
    CALL BUBBLE_SORT(NUM,INDEX_S)
    DO I1=1,7
    IF (NUM(I1).GT.NCENTER) THEN   
    IF (INDEX_S(I1).EQ.1) THEN
        WRITE(2,235) '4(',F11,' ',F14,' ',F13,' ',F12,')'
        WRITE(4,231) NCENTER
        WRITE(5,231) NUM(I1)
    ENDIF
    IF (INDEX_S(I1).EQ.2) THEN
        WRITE(2,235) '4(',F21,' ',F22,' ',F23,' ',F24,')'
        WRITE(4,231) NCENTER
        WRITE(5,231) NUM(I1)
    ENDIF
    IF (INDEX_S(I1).EQ.3) THEN
        WRITE(2,235) '4(',F12,' ',F13,' ',F23,' ',F22,')'
        WRITE(4,231) NCENTER
        WRITE(5,231) NUM(I1)
    ENDIF  
    IF (INDEX_S(I1).EQ.4) THEN
        WRITE(2,235) '4(',F11,' ',F21,' ',F24,' ',F14,')'
        WRITE(4,231) NCENTER
        WRITE(5,231) NUM(I1)
    ENDIF
    IF (INDEX_S(I1).EQ.5) THEN
        WRITE(2,235) '4(',F14,' ',F24,' ',F23,' ',F13,')'
        WRITE(4,231) NCENTER
        WRITE(5,231) NUM(I1)
    ENDIF
    IF (INDEX_S(I1).EQ.6) THEN
        WRITE(2,235) '4(',F11,' ',F12,' ',F22,' ',F21,')'
        WRITE(4,231) NCENTER
        WRITE(5,231) NUM(I1)
    ENDIF       
    ENDIF  
    ENDDO
    
    ELSEIF (JP.EQ.CNJ) THEN
    !CASE 1-3
    F11=(CNJ+1)*(CNI+1)*(KP-1)+(CNI+1)*(JP-1)+IP-1
    F12=F11+1
    F13=F12+CNI+1
    F14=F13-1
    F21=F11+(CNJ+1)*(CNI+1)
    IF (KP.EQ.NK) F21=(CNI+1)*(JP-1)+IP-1
    F22=F21+1
    F23=F22+CNI+1
    F24=F23-1
    CALL CAL_FIJ_CY(F11,F12,F13,F14,F21,F22,F23,F24)
    NUM(7)=CELL_ID_CY(IP,JP,KP)
    NUM(1)=CELL_ID_CY(IP,JP,KP-1) !FRONT FACE
    NUM(2)=CELL_ID_CY(IP,JP,KP+1) !BACK  FACE
    NUM(3)=CELL_ID_CY(IP+1,JP,KP) !LEFT  FACE
    NUM(4)=-1 !RIGHT FACE
    NUM(5)=-1 !UP    FACE
    NUM(6)=CELL_ID_CY(IP,JP-1,KP) !DOWN  FACE
    CALL BUBBLE_SORT(NUM,INDEX_S)
    DO I1=1,7
    IF (NUM(I1).GT.NCENTER) THEN   
    IF (INDEX_S(I1).EQ.1) THEN
        WRITE(2,235) '4(',F11,' ',F14,' ',F13,' ',F12,')'
        WRITE(4,231) NCENTER
        WRITE(5,231) NUM(I1)
    ENDIF
    IF (INDEX_S(I1).EQ.2) THEN
        WRITE(2,235) '4(',F21,' ',F22,' ',F23,' ',F24,')'
        WRITE(4,231) NCENTER
        WRITE(5,231) NUM(I1)
    ENDIF
    IF (INDEX_S(I1).EQ.3) THEN
        WRITE(2,235) '4(',F12,' ',F13,' ',F23,' ',F22,')'
        WRITE(4,231) NCENTER
        WRITE(5,231) NUM(I1)
    ENDIF  
    IF (INDEX_S(I1).EQ.4) THEN
        WRITE(2,235) '4(',F11,' ',F21,' ',F24,' ',F14,')'
        WRITE(4,231) NCENTER
        WRITE(5,231) NUM(I1)
    ENDIF
    IF (INDEX_S(I1).EQ.5) THEN
        WRITE(2,235) '4(',F14,' ',F24,' ',F23,' ',F13,')'
        WRITE(4,231) NCENTER
        WRITE(5,231) NUM(I1)
    ENDIF
    IF (INDEX_S(I1).EQ.6) THEN
        WRITE(2,235) '4(',F11,' ',F12,' ',F22,' ',F21,')'
        WRITE(4,231) NCENTER
        WRITE(5,231) NUM(I1)
    ENDIF       
    ENDIF  
    ENDDO
        
    ELSE
    !CASE 1-2 
    F11=(CNJ+1)*(CNI+1)*(KP-1)+(CNI+1)*(JP-1)+IP-1
    F12=F11+1
    F13=F12+CNI+1
    F14=F13-1
    F21=F11+(CNJ+1)*(CNI+1)
    IF (KP.EQ.NK) F21=(CNI+1)*(JP-1)+IP-1
    F22=F21+1
    F23=F22+CNI+1
    F24=F23-1
    CALL CAL_FIJ_CY(F11,F12,F13,F14,F21,F22,F23,F24)
    NUM(7)=CELL_ID_CY(IP,JP,KP)
    NUM(1)=CELL_ID_CY(IP,JP,KP-1) !FRONT FACE
    NUM(2)=CELL_ID_CY(IP,JP,KP+1) !BACK  FACE
    NUM(3)=CELL_ID_CY(IP+1,JP,KP) !LEFT  FACE
    NUM(4)=-1 !RIGHT FACE
    NUM(5)=CELL_ID_CY(IP,JP+1,KP) !UP    FACE
    NUM(6)=CELL_ID_CY(IP,JP-1,KP) !DOWN  FACE
    CALL BUBBLE_SORT(NUM,INDEX_S)
    DO I1=1,7
    IF (NUM(I1).GT.NCENTER) THEN   
    IF (INDEX_S(I1).EQ.1) THEN
        WRITE(2,235) '4(',F11,' ',F14,' ',F13,' ',F12,')'
        WRITE(4,231) NCENTER
        WRITE(5,231) NUM(I1)
    ENDIF
    IF (INDEX_S(I1).EQ.2) THEN
        WRITE(2,235) '4(',F21,' ',F22,' ',F23,' ',F24,')'
        WRITE(4,231) NCENTER
        WRITE(5,231) NUM(I1)
    ENDIF
    IF (INDEX_S(I1).EQ.3) THEN
        WRITE(2,235) '4(',F12,' ',F13,' ',F23,' ',F22,')'
        WRITE(4,231) NCENTER
        WRITE(5,231) NUM(I1)
    ENDIF  
    IF (INDEX_S(I1).EQ.4) THEN
        WRITE(2,235) '4(',F11,' ',F21,' ',F24,' ',F14,')'
        WRITE(4,231) NCENTER
        WRITE(5,231) NUM(I1)
    ENDIF
    IF (INDEX_S(I1).EQ.5) THEN
        WRITE(2,235) '4(',F14,' ',F24,' ',F23,' ',F13,')'
        WRITE(4,231) NCENTER
        WRITE(5,231) NUM(I1)
    ENDIF
    IF (INDEX_S(I1).EQ.6) THEN
        WRITE(2,235) '4(',F11,' ',F12,' ',F22,' ',F21,')'
        WRITE(4,231) NCENTER
        WRITE(5,231) NUM(I1)
    ENDIF       
    ENDIF  
    ENDDO
    
    ENDIF    
    ELSEIF (IP.EQ.CNI) THEN
    IF (JP.EQ.1) THEN 
    !CASE 3-1   
    F11=(CNJ+1)*(CNI+1)*(KP-1)+(CNI+1)*(JP-1)+IP-1
    F12=F11+1
    F13=F12+CNI+1
    F14=F13-1
    F21=F11+(CNJ+1)*(CNI+1)
    IF (KP.EQ.NK) F21=(CNI+1)*(JP-1)+IP-1
    F22=F21+1
    F23=F22+CNI+1
    F24=F23-1
    CALL CAL_FIJ_CY(F11,F12,F13,F14,F21,F22,F23,F24)
    NUM(7)=CELL_ID_CY(IP,JP,KP)
    NUM(1)=CELL_ID_CY(IP,JP,KP-1) !FRONT FACE
    NUM(2)=CELL_ID_CY(IP,JP,KP+1) !BACK  FACE
    NUM(3)=-1 !LEFT  FACE
    NUM(4)=CELL_ID_CY(IP-1,JP,KP) !RIGHT FACE
    NUM(5)=CELL_ID_CY(IP,JP+1,KP) !UP    FACE
    NUM(6)=-1 !DOWN  FACE
    CALL BUBBLE_SORT(NUM,INDEX_S)
    DO I1=1,7
    IF (NUM(I1).GT.NCENTER) THEN   
    IF (INDEX_S(I1).EQ.1) THEN
        WRITE(2,235) '4(',F11,' ',F14,' ',F13,' ',F12,')'
        WRITE(4,231) NCENTER
        WRITE(5,231) NUM(I1)
    ENDIF
    IF (INDEX_S(I1).EQ.2) THEN
        WRITE(2,235) '4(',F21,' ',F22,' ',F23,' ',F24,')'
        WRITE(4,231) NCENTER
        WRITE(5,231) NUM(I1)
    ENDIF
    IF (INDEX_S(I1).EQ.3) THEN
        WRITE(2,235) '4(',F12,' ',F13,' ',F23,' ',F22,')'
        WRITE(4,231) NCENTER
        WRITE(5,231) NUM(I1)
    ENDIF  
    IF (INDEX_S(I1).EQ.4) THEN
        WRITE(2,235) '4(',F11,' ',F21,' ',F24,' ',F14,')'
        WRITE(4,231) NCENTER
        WRITE(5,231) NUM(I1)
    ENDIF
    IF (INDEX_S(I1).EQ.5) THEN
        WRITE(2,235) '4(',F14,' ',F24,' ',F23,' ',F13,')'
        WRITE(4,231) NCENTER
        WRITE(5,231) NUM(I1)
    ENDIF
    IF (INDEX_S(I1).EQ.6) THEN
        WRITE(2,235) '4(',F11,' ',F12,' ',F22,' ',F21,')'
        WRITE(4,231) NCENTER
        WRITE(5,231) NUM(I1)
    ENDIF       
    ENDIF  
    ENDDO
        
    ELSEIF (JP.EQ.CNJ) THEN
    !CASE 3-3
    F11=(CNJ+1)*(CNI+1)*(KP-1)+(CNI+1)*(JP-1)+IP-1
    F12=F11+1
    F13=F12+CNI+1
    F14=F13-1
    F21=F11+(CNJ+1)*(CNI+1)
    IF (KP.EQ.NK) F21=(CNI+1)*(JP-1)+IP-1
    F22=F21+1
    F23=F22+CNI+1
    F24=F23-1
    CALL CAL_FIJ_CY(F11,F12,F13,F14,F21,F22,F23,F24)
    NUM(7)=CELL_ID_CY(IP,JP,KP)
    NUM(1)=CELL_ID_CY(IP,JP,KP-1) !FRONT FACE
    NUM(2)=CELL_ID_CY(IP,JP,KP+1) !BACK  FACE
    NUM(3)=-1 !LEFT  FACE
    NUM(4)=CELL_ID_CY(IP-1,JP,KP) !RIGHT FACE
    NUM(5)=-1 !UP    FACE
    NUM(6)=CELL_ID_CY(IP,JP-1,KP) !DOWN  FACE
    CALL BUBBLE_SORT(NUM,INDEX_S)
    DO I1=1,7
    IF (NUM(I1).GT.NCENTER) THEN   
    IF (INDEX_S(I1).EQ.1) THEN
        WRITE(2,235) '4(',F11,' ',F14,' ',F13,' ',F12,')'
        WRITE(4,231) NCENTER
        WRITE(5,231) NUM(I1)
    ENDIF
    IF (INDEX_S(I1).EQ.2) THEN
        WRITE(2,235) '4(',F21,' ',F22,' ',F23,' ',F24,')'
        WRITE(4,231) NCENTER
        WRITE(5,231) NUM(I1)
    ENDIF
    IF (INDEX_S(I1).EQ.3) THEN
        WRITE(2,235) '4(',F12,' ',F13,' ',F23,' ',F22,')'
        WRITE(4,231) NCENTER
        WRITE(5,231) NUM(I1)
    ENDIF  
    IF (INDEX_S(I1).EQ.4) THEN
        WRITE(2,235) '4(',F11,' ',F21,' ',F24,' ',F14,')'
        WRITE(4,231) NCENTER
        WRITE(5,231) NUM(I1)
    ENDIF
    IF (INDEX_S(I1).EQ.5) THEN
        WRITE(2,235) '4(',F14,' ',F24,' ',F23,' ',F13,')'
        WRITE(4,231) NCENTER
        WRITE(5,231) NUM(I1)
    ENDIF
    IF (INDEX_S(I1).EQ.6) THEN
        WRITE(2,235) '4(',F11,' ',F12,' ',F22,' ',F21,')'
        WRITE(4,231) NCENTER
        WRITE(5,231) NUM(I1)
    ENDIF       
    ENDIF  
    ENDDO
        
    ELSE
    !CASE 3-2
    F11=(CNJ+1)*(CNI+1)*(KP-1)+(CNI+1)*(JP-1)+IP-1
    F12=F11+1
    F13=F12+CNI+1
    F14=F13-1
    F21=F11+(CNJ+1)*(CNI+1)
    IF (KP.EQ.NK) F21=(CNI+1)*(JP-1)+IP-1
    F22=F21+1
    F23=F22+CNI+1
    F24=F23-1
    CALL CAL_FIJ_CY(F11,F12,F13,F14,F21,F22,F23,F24)
    NUM(7)=CELL_ID_CY(IP,JP,KP)
    NUM(1)=CELL_ID_CY(IP,JP,KP-1) !FRONT FACE
    NUM(2)=CELL_ID_CY(IP,JP,KP+1) !BACK  FACE
    NUM(3)=-1 !LEFT  FACE
    NUM(4)=CELL_ID_CY(IP-1,JP,KP) !RIGHT FACE
    NUM(5)=CELL_ID_CY(IP,JP+1,KP) !UP    FACE
    NUM(6)=CELL_ID_CY(IP,JP-1,KP) !DOWN  FACE
    CALL BUBBLE_SORT(NUM,INDEX_S)
    DO I1=1,7
    IF (NUM(I1).GT.NCENTER) THEN   
    IF (INDEX_S(I1).EQ.1) THEN
        WRITE(2,235) '4(',F11,' ',F14,' ',F13,' ',F12,')'
        WRITE(4,231) NCENTER
        WRITE(5,231) NUM(I1)
    ENDIF
    IF (INDEX_S(I1).EQ.2) THEN
        WRITE(2,235) '4(',F21,' ',F22,' ',F23,' ',F24,')'
        WRITE(4,231) NCENTER
        WRITE(5,231) NUM(I1)
    ENDIF
    IF (INDEX_S(I1).EQ.3) THEN
        WRITE(2,235) '4(',F12,' ',F13,' ',F23,' ',F22,')'
        WRITE(4,231) NCENTER
        WRITE(5,231) NUM(I1)
    ENDIF  
    IF (INDEX_S(I1).EQ.4) THEN
        WRITE(2,235) '4(',F11,' ',F21,' ',F24,' ',F14,')'
        WRITE(4,231) NCENTER
        WRITE(5,231) NUM(I1)
    ENDIF
    IF (INDEX_S(I1).EQ.5) THEN
        WRITE(2,235) '4(',F14,' ',F24,' ',F23,' ',F13,')'
        WRITE(4,231) NCENTER
        WRITE(5,231) NUM(I1)
    ENDIF
    IF (INDEX_S(I1).EQ.6) THEN
        WRITE(2,235) '4(',F11,' ',F12,' ',F22,' ',F21,')'
        WRITE(4,231) NCENTER
        WRITE(5,231) NUM(I1)
    ENDIF       
    ENDIF  
    ENDDO
    
    ENDIF    
    ELSE
    IF (JP.EQ.1) THEN 
    !CASE 2-1 
    F11=(CNJ+1)*(CNI+1)*(KP-1)+(CNI+1)*(JP-1)+IP-1
    F12=F11+1
    F13=F12+CNI+1
    F14=F13-1
    F21=F11+(CNJ+1)*(CNI+1)
    IF (KP.EQ.NK) F21=(CNI+1)*(JP-1)+IP-1
    F22=F21+1
    F23=F22+CNI+1
    F24=F23-1
    CALL CAL_FIJ_CY(F11,F12,F13,F14,F21,F22,F23,F24)
    NUM(7)=CELL_ID_CY(IP,JP,KP)
    NUM(1)=CELL_ID_CY(IP,JP,KP-1) !FRONT FACE
    NUM(2)=CELL_ID_CY(IP,JP,KP+1) !BACK  FACE
    NUM(3)=CELL_ID_CY(IP+1,JP,KP) !LEFT  FACE
    NUM(4)=CELL_ID_CY(IP-1,JP,KP) !RIGHT FACE
    NUM(5)=CELL_ID_CY(IP,JP+1,KP) !UP    FACE
    NUM(6)=-1 !DOWN  FACE
    CALL BUBBLE_SORT(NUM,INDEX_S)
    DO I1=1,7
    IF (NUM(I1).GT.NCENTER) THEN   
    IF (INDEX_S(I1).EQ.1) THEN
        WRITE(2,235) '4(',F11,' ',F14,' ',F13,' ',F12,')'
        WRITE(4,231) NCENTER
        WRITE(5,231) NUM(I1)
    ENDIF
    IF (INDEX_S(I1).EQ.2) THEN
        WRITE(2,235) '4(',F21,' ',F22,' ',F23,' ',F24,')'
        WRITE(4,231) NCENTER
        WRITE(5,231) NUM(I1)
    ENDIF
    IF (INDEX_S(I1).EQ.3) THEN
        WRITE(2,235) '4(',F12,' ',F13,' ',F23,' ',F22,')'
        WRITE(4,231) NCENTER
        WRITE(5,231) NUM(I1)
    ENDIF  
    IF (INDEX_S(I1).EQ.4) THEN
        WRITE(2,235) '4(',F11,' ',F21,' ',F24,' ',F14,')'
        WRITE(4,231) NCENTER
        WRITE(5,231) NUM(I1)
    ENDIF
    IF (INDEX_S(I1).EQ.5) THEN
        WRITE(2,235) '4(',F14,' ',F24,' ',F23,' ',F13,')'
        WRITE(4,231) NCENTER
        WRITE(5,231) NUM(I1)
    ENDIF
    IF (INDEX_S(I1).EQ.6) THEN
        WRITE(2,235) '4(',F11,' ',F12,' ',F22,' ',F21,')'
        WRITE(4,231) NCENTER
        WRITE(5,231) NUM(I1)
    ENDIF       
    ENDIF  
    ENDDO
        
    ELSEIF (JP.EQ.CNJ) THEN
    !CASE 2-3
    F11=(CNJ+1)*(CNI+1)*(KP-1)+(CNI+1)*(JP-1)+IP-1
    F12=F11+1
    F13=F12+CNI+1
    F14=F13-1
    F21=F11+(CNJ+1)*(CNI+1)
    IF (KP.EQ.NK) F21=(CNI+1)*(JP-1)+IP-1
    F22=F21+1
    F23=F22+CNI+1
    F24=F23-1
    CALL CAL_FIJ_CY(F11,F12,F13,F14,F21,F22,F23,F24)
    NUM(7)=CELL_ID_CY(IP,JP,KP)
    NUM(1)=CELL_ID_CY(IP,JP,KP-1) !FRONT FACE
    NUM(2)=CELL_ID_CY(IP,JP,KP+1) !BACK  FACE
    NUM(3)=CELL_ID_CY(IP+1,JP,KP) !LEFT  FACE
    NUM(4)=CELL_ID_CY(IP-1,JP,KP) !RIGHT FACE
    NUM(5)=-1 !UP    FACE
    NUM(6)=CELL_ID_CY(IP,JP-1,KP) !DOWN  FACE
    CALL BUBBLE_SORT(NUM,INDEX_S)
    DO I1=1,7
    IF (NUM(I1).GT.NCENTER) THEN   
    IF (INDEX_S(I1).EQ.1) THEN
        WRITE(2,235) '4(',F11,' ',F14,' ',F13,' ',F12,')'
        WRITE(4,231) NCENTER
        WRITE(5,231) NUM(I1)
    ENDIF
    IF (INDEX_S(I1).EQ.2) THEN
        WRITE(2,235) '4(',F21,' ',F22,' ',F23,' ',F24,')'
        WRITE(4,231) NCENTER
        WRITE(5,231) NUM(I1)
    ENDIF
    IF (INDEX_S(I1).EQ.3) THEN
        WRITE(2,235) '4(',F12,' ',F13,' ',F23,' ',F22,')'
        WRITE(4,231) NCENTER
        WRITE(5,231) NUM(I1)
    ENDIF  
    IF (INDEX_S(I1).EQ.4) THEN
        WRITE(2,235) '4(',F11,' ',F21,' ',F24,' ',F14,')'
        WRITE(4,231) NCENTER
        WRITE(5,231) NUM(I1)
    ENDIF
    IF (INDEX_S(I1).EQ.5) THEN
        WRITE(2,235) '4(',F14,' ',F24,' ',F23,' ',F13,')'
        WRITE(4,231) NCENTER
        WRITE(5,231) NUM(I1)
    ENDIF
    IF (INDEX_S(I1).EQ.6) THEN
        WRITE(2,235) '4(',F11,' ',F12,' ',F22,' ',F21,')'
        WRITE(4,231) NCENTER
        WRITE(5,231) NUM(I1)
    ENDIF       
    ENDIF  
    ENDDO
        
    ELSE
    !CASE 2-2
    F11=(CNJ+1)*(CNI+1)*(KP-1)+(CNI+1)*(JP-1)+IP-1
    F12=F11+1
    F13=F12+CNI+1
    F14=F13-1
    F21=F11+(CNJ+1)*(CNI+1)
    IF (KP.EQ.NK) F21=(CNI+1)*(JP-1)+IP-1
    F22=F21+1
    F23=F22+CNI+1
    F24=F23-1
    CALL CAL_FIJ_CY(F11,F12,F13,F14,F21,F22,F23,F24)
    NUM(7)=CELL_ID_CY(IP,JP,KP)
    NUM(1)=CELL_ID_CY(IP,JP,KP-1) !FRONT FACE
    NUM(2)=CELL_ID_CY(IP,JP,KP+1) !BACK  FACE
    NUM(3)=CELL_ID_CY(IP+1,JP,KP) !LEFT  FACE
    NUM(4)=CELL_ID_CY(IP-1,JP,KP) !RIGHT FACE
    NUM(5)=CELL_ID_CY(IP,JP+1,KP) !UP    FACE
    NUM(6)=CELL_ID_CY(IP,JP-1,KP) !DOWN  FACE
    CALL BUBBLE_SORT(NUM,INDEX_S)
    DO I1=1,7
    IF (NUM(I1).GT.NCENTER) THEN   
    IF (INDEX_S(I1).EQ.1) THEN
        WRITE(2,235) '4(',F11,' ',F14,' ',F13,' ',F12,')'
        WRITE(4,231) NCENTER
        WRITE(5,231) NUM(I1)
    ENDIF
    IF (INDEX_S(I1).EQ.2) THEN
        WRITE(2,235) '4(',F21,' ',F22,' ',F23,' ',F24,')'
        WRITE(4,231) NCENTER
        WRITE(5,231) NUM(I1)
    ENDIF
    IF (INDEX_S(I1).EQ.3) THEN
        WRITE(2,235) '4(',F12,' ',F13,' ',F23,' ',F22,')'
        WRITE(4,231) NCENTER
        WRITE(5,231) NUM(I1)
    ENDIF  
    IF (INDEX_S(I1).EQ.4) THEN
        WRITE(2,235) '4(',F11,' ',F21,' ',F24,' ',F14,')'
        WRITE(4,231) NCENTER
        WRITE(5,231) NUM(I1)
    ENDIF
    IF (INDEX_S(I1).EQ.5) THEN
        WRITE(2,235) '4(',F14,' ',F24,' ',F23,' ',F13,')'
        WRITE(4,231) NCENTER
        WRITE(5,231) NUM(I1)
    ENDIF
    IF (INDEX_S(I1).EQ.6) THEN
        WRITE(2,235) '4(',F11,' ',F12,' ',F22,' ',F21,')'
        WRITE(4,231) NCENTER
        WRITE(5,231) NUM(I1)
    ENDIF       
    ENDIF  
    ENDDO
        
    ENDIF    
    ENDIF
    
        
    ENDIF
    ENDIF
    
    
    ELSEIF (IF_CY.EQ.0) THEN
!*******************RECTANGLE***************  
    IP=I-CNI*(MI-1)
    JP=J-CNL*(MJ1-1)
    KP=K-CNL*(MJ2-1)
    IF (IP.LE.CNI.AND.IP.GE.1) THEN
    IF (JP.LE.CNL.AND.JP.GE.1) THEN
    IF (KP.LE.CNL.AND.KP.GE.1) THEN
        NCENTER=CELL_ID_RECT(IP,JP,KP)  
    IF (IP.EQ.1) THEN
    IF (JP.EQ.1) THEN
    IF (KP.EQ.1) THEN
    !CASE 1-1-1 
    F11=(CNL+1)*(CNI+1)*(KP-1)+(CNI+1)*(JP-1)+IP-1
    F12=F11+1
    F13=F12+CNI+1
    F14=F13-1
    F21=F11+(CNL+1)*(CNI+1)
    F22=F21+1
    F23=F22+CNI+1
    F24=F23-1 
    CALL CAL_FIJ_RECT(F11,F12,F13,F14,F21,F22,F23,F24)
    NUM(7)=CELL_ID_RECT(IP,JP,KP)
    NUM(1)=-1 !FRONT FACE
    NUM(2)=CELL_ID_RECT(IP,JP,KP+1) !BACK  FACE
    NUM(3)=CELL_ID_RECT(IP+1,JP,KP) !LEFT  FACE
    NUM(4)=-1 !RIGHT FACE
    NUM(5)=CELL_ID_RECT(IP,JP+1,KP) !UP    FACE
    NUM(6)=-1 !DOWN  FACE
    CALL BUBBLE_SORT(NUM,INDEX_S)
    DO I1=1,7
    IF (NUM(I1).GT.NCENTER) THEN   
    IF (INDEX_S(I1).EQ.1) THEN
        WRITE(2,235) '4(',F11,' ',F14,' ',F13,' ',F12,')'
        WRITE(4,231) NCENTER
        WRITE(5,231) NUM(I1)
    ENDIF
    IF (INDEX_S(I1).EQ.2) THEN
        WRITE(2,235) '4(',F21,' ',F22,' ',F23,' ',F24,')'
        WRITE(4,231) NCENTER
        WRITE(5,231) NUM(I1)
    ENDIF
    IF (INDEX_S(I1).EQ.3) THEN
        WRITE(2,235) '4(',F12,' ',F13,' ',F23,' ',F22,')'
        WRITE(4,231) NCENTER
        WRITE(5,231) NUM(I1)
    ENDIF  
    IF (INDEX_S(I1).EQ.4) THEN
        WRITE(2,235) '4(',F11,' ',F21,' ',F24,' ',F14,')'
        WRITE(4,231) NCENTER
        WRITE(5,231) NUM(I1)
    ENDIF
    IF (INDEX_S(I1).EQ.5) THEN
        WRITE(2,235) '4(',F14,' ',F24,' ',F23,' ',F13,')'
        WRITE(4,231) NCENTER
        WRITE(5,231) NUM(I1)
    ENDIF
    IF (INDEX_S(I1).EQ.6) THEN
        WRITE(2,235) '4(',F11,' ',F12,' ',F22,' ',F21,')'
        WRITE(4,231) NCENTER
        WRITE(5,231) NUM(I1)
    ENDIF       
    ENDIF  
    ENDDO
    
    ELSEIF (KP.EQ.CNL) THEN
    !CASE 1-1-3 
    F11=(CNL+1)*(CNI+1)*(KP-1)+(CNI+1)*(JP-1)+IP-1
    F12=F11+1
    F13=F12+CNI+1
    F14=F13-1
    F21=F11+(CNL+1)*(CNI+1)
    F22=F21+1
    F23=F22+CNI+1
    F24=F23-1   
    CALL CAL_FIJ_RECT(F11,F12,F13,F14,F21,F22,F23,F24)
    NUM(7)=CELL_ID_RECT(IP,JP,KP)
    NUM(1)=CELL_ID_RECT(IP,JP,KP-1) !FRONT FACE
    NUM(2)=-1 !BACK  FACE
    NUM(3)=CELL_ID_RECT(IP+1,JP,KP) !LEFT  FACE
    NUM(4)=-1 !RIGHT FACE
    NUM(5)=CELL_ID_RECT(IP,JP+1,KP) !UP    FACE
    NUM(6)=-1 !DOWN  FACE
    CALL BUBBLE_SORT(NUM,INDEX_S)
    DO I1=1,7
    IF (NUM(I1).GT.NCENTER) THEN   
    IF (INDEX_S(I1).EQ.1) THEN
        WRITE(2,235) '4(',F11,' ',F14,' ',F13,' ',F12,')'
        WRITE(4,231) NCENTER
        WRITE(5,231) NUM(I1)
    ENDIF
    IF (INDEX_S(I1).EQ.2) THEN
        WRITE(2,235) '4(',F21,' ',F22,' ',F23,' ',F24,')'
        WRITE(4,231) NCENTER
        WRITE(5,231) NUM(I1)
    ENDIF
    IF (INDEX_S(I1).EQ.3) THEN
        WRITE(2,235) '4(',F12,' ',F13,' ',F23,' ',F22,')'
        WRITE(4,231) NCENTER
        WRITE(5,231) NUM(I1)
    ENDIF  
    IF (INDEX_S(I1).EQ.4) THEN
        WRITE(2,235) '4(',F11,' ',F21,' ',F24,' ',F14,')'
        WRITE(4,231) NCENTER
        WRITE(5,231) NUM(I1)
    ENDIF
    IF (INDEX_S(I1).EQ.5) THEN
        WRITE(2,235) '4(',F14,' ',F24,' ',F23,' ',F13,')'
        WRITE(4,231) NCENTER
        WRITE(5,231) NUM(I1)
    ENDIF
    IF (INDEX_S(I1).EQ.6) THEN
        WRITE(2,235) '4(',F11,' ',F12,' ',F22,' ',F21,')'
        WRITE(4,231) NCENTER
        WRITE(5,231) NUM(I1)
    ENDIF       
    ENDIF  
    ENDDO
    
    ELSE
    !CASE 1-1-2 
    F11=(CNL+1)*(CNI+1)*(KP-1)+(CNI+1)*(JP-1)+IP-1
    F12=F11+1
    F13=F12+CNI+1
    F14=F13-1
    F21=F11+(CNL+1)*(CNI+1)
    F22=F21+1
    F23=F22+CNI+1
    F24=F23-1  
    CALL CAL_FIJ_RECT(F11,F12,F13,F14,F21,F22,F23,F24)
    NUM(7)=CELL_ID_RECT(IP,JP,KP)
    NUM(1)=CELL_ID_RECT(IP,JP,KP-1) !FRONT FACE
    NUM(2)=CELL_ID_RECT(IP,JP,KP+1) !BACK  FACE
    NUM(3)=CELL_ID_RECT(IP+1,JP,KP) !LEFT  FACE
    NUM(4)=-1 !RIGHT FACE
    NUM(5)=CELL_ID_RECT(IP,JP+1,KP) !UP    FACE
    NUM(6)=-1 !DOWN  FACE
    CALL BUBBLE_SORT(NUM,INDEX_S)
    DO I1=1,7
    IF (NUM(I1).GT.NCENTER) THEN   
    IF (INDEX_S(I1).EQ.1) THEN
        WRITE(2,235) '4(',F11,' ',F14,' ',F13,' ',F12,')'
        WRITE(4,231) NCENTER
        WRITE(5,231) NUM(I1)
    ENDIF
    IF (INDEX_S(I1).EQ.2) THEN
        WRITE(2,235) '4(',F21,' ',F22,' ',F23,' ',F24,')'
        WRITE(4,231) NCENTER
        WRITE(5,231) NUM(I1)
    ENDIF
    IF (INDEX_S(I1).EQ.3) THEN
        WRITE(2,235) '4(',F12,' ',F13,' ',F23,' ',F22,')'
        WRITE(4,231) NCENTER
        WRITE(5,231) NUM(I1)
    ENDIF  
    IF (INDEX_S(I1).EQ.4) THEN
        WRITE(2,235) '4(',F11,' ',F21,' ',F24,' ',F14,')'
        WRITE(4,231) NCENTER
        WRITE(5,231) NUM(I1)
    ENDIF
    IF (INDEX_S(I1).EQ.5) THEN
        WRITE(2,235) '4(',F14,' ',F24,' ',F23,' ',F13,')'
        WRITE(4,231) NCENTER
        WRITE(5,231) NUM(I1)
    ENDIF
    IF (INDEX_S(I1).EQ.6) THEN
        WRITE(2,235) '4(',F11,' ',F12,' ',F22,' ',F21,')'
        WRITE(4,231) NCENTER
        WRITE(5,231) NUM(I1)
    ENDIF       
    ENDIF  
    ENDDO
    
    ENDIF    
    ELSEIF (JP.EQ.CNL) THEN
    IF (KP.EQ.1) THEN
    !CASE 1-3-1  
    F11=(CNL+1)*(CNI+1)*(KP-1)+(CNI+1)*(JP-1)+IP-1
    F12=F11+1
    F13=F12+CNI+1
    F14=F13-1
    F21=F11+(CNL+1)*(CNI+1)
    F22=F21+1
    F23=F22+CNI+1
    F24=F23-1  
    CALL CAL_FIJ_RECT(F11,F12,F13,F14,F21,F22,F23,F24)
    NUM(7)=CELL_ID_RECT(IP,JP,KP)
    NUM(1)=-1 !FRONT FACE
    NUM(2)=CELL_ID_RECT(IP,JP,KP+1) !BACK  FACE
    NUM(3)=CELL_ID_RECT(IP+1,JP,KP) !LEFT  FACE
    NUM(4)=-1 !RIGHT FACE
    NUM(5)=-1 !UP    FACE
    NUM(6)=CELL_ID_RECT(IP,JP-1,KP) !DOWN  FACE
    CALL BUBBLE_SORT(NUM,INDEX_S)
    DO I1=1,7
    IF (NUM(I1).GT.NCENTER) THEN   
    IF (INDEX_S(I1).EQ.1) THEN
        WRITE(2,235) '4(',F11,' ',F14,' ',F13,' ',F12,')'
        WRITE(4,231) NCENTER
        WRITE(5,231) NUM(I1)
    ENDIF
    IF (INDEX_S(I1).EQ.2) THEN
        WRITE(2,235) '4(',F21,' ',F22,' ',F23,' ',F24,')'
        WRITE(4,231) NCENTER
        WRITE(5,231) NUM(I1)
    ENDIF
    IF (INDEX_S(I1).EQ.3) THEN
        WRITE(2,235) '4(',F12,' ',F13,' ',F23,' ',F22,')'
        WRITE(4,231) NCENTER
        WRITE(5,231) NUM(I1)
    ENDIF  
    IF (INDEX_S(I1).EQ.4) THEN
        WRITE(2,235) '4(',F11,' ',F21,' ',F24,' ',F14,')'
        WRITE(4,231) NCENTER
        WRITE(5,231) NUM(I1)
    ENDIF
    IF (INDEX_S(I1).EQ.5) THEN
        WRITE(2,235) '4(',F14,' ',F24,' ',F23,' ',F13,')'
        WRITE(4,231) NCENTER
        WRITE(5,231) NUM(I1)
    ENDIF
    IF (INDEX_S(I1).EQ.6) THEN
        WRITE(2,235) '4(',F11,' ',F12,' ',F22,' ',F21,')'
        WRITE(4,231) NCENTER
        WRITE(5,231) NUM(I1)
    ENDIF       
    ENDIF  
    ENDDO
    
    ELSEIF (KP.EQ.CNL) THEN
    !CASE 1-3-3  
    F11=(CNL+1)*(CNI+1)*(KP-1)+(CNI+1)*(JP-1)+IP-1
    F12=F11+1
    F13=F12+CNI+1
    F14=F13-1
    F21=F11+(CNL+1)*(CNI+1)
    F22=F21+1
    F23=F22+CNI+1
    F24=F23-1 
    CALL CAL_FIJ_RECT(F11,F12,F13,F14,F21,F22,F23,F24)
    NUM(7)=CELL_ID_RECT(IP,JP,KP)
    NUM(1)=CELL_ID_RECT(IP,JP,KP-1) !FRONT FACE
    NUM(2)=-1 !BACK  FACE
    NUM(3)=CELL_ID_RECT(IP+1,JP,KP) !LEFT  FACE
    NUM(4)=-1 !RIGHT FACE
    NUM(5)=-1 !UP    FACE
    NUM(6)=CELL_ID_RECT(IP,JP-1,KP) !DOWN  FACE
    CALL BUBBLE_SORT(NUM,INDEX_S)
    DO I1=1,7
    IF (NUM(I1).GT.NCENTER) THEN   
    IF (INDEX_S(I1).EQ.1) THEN
        WRITE(2,235) '4(',F11,' ',F14,' ',F13,' ',F12,')'
        WRITE(4,231) NCENTER
        WRITE(5,231) NUM(I1)
    ENDIF
    IF (INDEX_S(I1).EQ.2) THEN
        WRITE(2,235) '4(',F21,' ',F22,' ',F23,' ',F24,')'
        WRITE(4,231) NCENTER
        WRITE(5,231) NUM(I1)
    ENDIF
    IF (INDEX_S(I1).EQ.3) THEN
        WRITE(2,235) '4(',F12,' ',F13,' ',F23,' ',F22,')'
        WRITE(4,231) NCENTER
        WRITE(5,231) NUM(I1)
    ENDIF  
    IF (INDEX_S(I1).EQ.4) THEN
        WRITE(2,235) '4(',F11,' ',F21,' ',F24,' ',F14,')'
        WRITE(4,231) NCENTER
        WRITE(5,231) NUM(I1)
    ENDIF
    IF (INDEX_S(I1).EQ.5) THEN
        WRITE(2,235) '4(',F14,' ',F24,' ',F23,' ',F13,')'
        WRITE(4,231) NCENTER
        WRITE(5,231) NUM(I1)
    ENDIF
    IF (INDEX_S(I1).EQ.6) THEN
        WRITE(2,235) '4(',F11,' ',F12,' ',F22,' ',F21,')'
        WRITE(4,231) NCENTER
        WRITE(5,231) NUM(I1)
    ENDIF       
    ENDIF  
    ENDDO
    
    ELSE
    !CASE 1-3-2 
    F11=(CNL+1)*(CNI+1)*(KP-1)+(CNI+1)*(JP-1)+IP-1
    F12=F11+1
    F13=F12+CNI+1
    F14=F13-1
    F21=F11+(CNL+1)*(CNI+1)
    F22=F21+1
    F23=F22+CNI+1
    F24=F23-1   
    CALL CAL_FIJ_RECT(F11,F12,F13,F14,F21,F22,F23,F24)
    NUM(7)=CELL_ID_RECT(IP,JP,KP)
    NUM(1)=CELL_ID_RECT(IP,JP,KP-1) !FRONT FACE
    NUM(2)=CELL_ID_RECT(IP,JP,KP+1) !BACK  FACE
    NUM(3)=CELL_ID_RECT(IP+1,JP,KP) !LEFT  FACE
    NUM(4)=-1 !RIGHT FACE
    NUM(5)=-1 !UP    FACE
    NUM(6)=CELL_ID_RECT(IP,JP-1,KP) !DOWN  FACE
    CALL BUBBLE_SORT(NUM,INDEX_S)
    DO I1=1,7
    IF (NUM(I1).GT.NCENTER) THEN   
    IF (INDEX_S(I1).EQ.1) THEN
        WRITE(2,235) '4(',F11,' ',F14,' ',F13,' ',F12,')'
        WRITE(4,231) NCENTER
        WRITE(5,231) NUM(I1)
    ENDIF
    IF (INDEX_S(I1).EQ.2) THEN
        WRITE(2,235) '4(',F21,' ',F22,' ',F23,' ',F24,')'
        WRITE(4,231) NCENTER
        WRITE(5,231) NUM(I1)
    ENDIF
    IF (INDEX_S(I1).EQ.3) THEN
        WRITE(2,235) '4(',F12,' ',F13,' ',F23,' ',F22,')'
        WRITE(4,231) NCENTER
        WRITE(5,231) NUM(I1)
    ENDIF  
    IF (INDEX_S(I1).EQ.4) THEN
        WRITE(2,235) '4(',F11,' ',F21,' ',F24,' ',F14,')'
        WRITE(4,231) NCENTER
        WRITE(5,231) NUM(I1)
    ENDIF
    IF (INDEX_S(I1).EQ.5) THEN
        WRITE(2,235) '4(',F14,' ',F24,' ',F23,' ',F13,')'
        WRITE(4,231) NCENTER
        WRITE(5,231) NUM(I1)
    ENDIF
    IF (INDEX_S(I1).EQ.6) THEN
        WRITE(2,235) '4(',F11,' ',F12,' ',F22,' ',F21,')'
        WRITE(4,231) NCENTER
        WRITE(5,231) NUM(I1)
    ENDIF       
    ENDIF  
    ENDDO
    
    ENDIF    
    ELSE
    IF (KP.EQ.1) THEN
    !CASE 1-2-1 
    F11=(CNL+1)*(CNI+1)*(KP-1)+(CNI+1)*(JP-1)+IP-1
    F12=F11+1
    F13=F12+CNI+1
    F14=F13-1
    F21=F11+(CNL+1)*(CNI+1)
    F22=F21+1
    F23=F22+CNI+1
    F24=F23-1   
    CALL CAL_FIJ_RECT(F11,F12,F13,F14,F21,F22,F23,F24)
    NUM(7)=CELL_ID_RECT(IP,JP,KP)
    NUM(1)=-1 !FRONT FACE
    NUM(2)=CELL_ID_RECT(IP,JP,KP+1) !BACK  FACE
    NUM(3)=CELL_ID_RECT(IP+1,JP,KP) !LEFT  FACE
    NUM(4)=-1 !RIGHT FACE
    NUM(5)=CELL_ID_RECT(IP,JP+1,KP) !UP    FACE
    NUM(6)=CELL_ID_RECT(IP,JP-1,KP) !DOWN  FACE
    CALL BUBBLE_SORT(NUM,INDEX_S)
    DO I1=1,7
    IF (NUM(I1).GT.NCENTER) THEN   
    IF (INDEX_S(I1).EQ.1) THEN
        WRITE(2,235) '4(',F11,' ',F14,' ',F13,' ',F12,')'
        WRITE(4,231) NCENTER
        WRITE(5,231) NUM(I1)
    ENDIF
    IF (INDEX_S(I1).EQ.2) THEN
        WRITE(2,235) '4(',F21,' ',F22,' ',F23,' ',F24,')'
        WRITE(4,231) NCENTER
        WRITE(5,231) NUM(I1)
    ENDIF
    IF (INDEX_S(I1).EQ.3) THEN
        WRITE(2,235) '4(',F12,' ',F13,' ',F23,' ',F22,')'
        WRITE(4,231) NCENTER
        WRITE(5,231) NUM(I1)
    ENDIF  
    IF (INDEX_S(I1).EQ.4) THEN
        WRITE(2,235) '4(',F11,' ',F21,' ',F24,' ',F14,')'
        WRITE(4,231) NCENTER
        WRITE(5,231) NUM(I1)
    ENDIF
    IF (INDEX_S(I1).EQ.5) THEN
        WRITE(2,235) '4(',F14,' ',F24,' ',F23,' ',F13,')'
        WRITE(4,231) NCENTER
        WRITE(5,231) NUM(I1)
    ENDIF
    IF (INDEX_S(I1).EQ.6) THEN
        WRITE(2,235) '4(',F11,' ',F12,' ',F22,' ',F21,')'
        WRITE(4,231) NCENTER
        WRITE(5,231) NUM(I1)
    ENDIF       
    ENDIF  
    ENDDO
    
    ELSEIF (KP.EQ.CNL) THEN
    !CASE 1-2-3  
    F11=(CNL+1)*(CNI+1)*(KP-1)+(CNI+1)*(JP-1)+IP-1
    F12=F11+1
    F13=F12+CNI+1
    F14=F13-1
    F21=F11+(CNL+1)*(CNI+1)
    F22=F21+1
    F23=F22+CNI+1
    F24=F23-1   
    CALL CAL_FIJ_RECT(F11,F12,F13,F14,F21,F22,F23,F24)
    NUM(7)=CELL_ID_RECT(IP,JP,KP)
    NUM(1)=CELL_ID_RECT(IP,JP,KP-1) !FRONT FACE
    NUM(2)=-1 !BACK  FACE
    NUM(3)=CELL_ID_RECT(IP+1,JP,KP) !LEFT  FACE
    NUM(4)=-1 !RIGHT FACE
    NUM(5)=CELL_ID_RECT(IP,JP+1,KP) !UP    FACE
    NUM(6)=CELL_ID_RECT(IP,JP-1,KP) !DOWN  FACE
    CALL BUBBLE_SORT(NUM,INDEX_S)
    DO I1=1,7
    IF (NUM(I1).GT.NCENTER) THEN   
    IF (INDEX_S(I1).EQ.1) THEN
        WRITE(2,235) '4(',F11,' ',F14,' ',F13,' ',F12,')'
        WRITE(4,231) NCENTER
        WRITE(5,231) NUM(I1)
    ENDIF
    IF (INDEX_S(I1).EQ.2) THEN
        WRITE(2,235) '4(',F21,' ',F22,' ',F23,' ',F24,')'
        WRITE(4,231) NCENTER
        WRITE(5,231) NUM(I1)
    ENDIF
    IF (INDEX_S(I1).EQ.3) THEN
        WRITE(2,235) '4(',F12,' ',F13,' ',F23,' ',F22,')'
        WRITE(4,231) NCENTER
        WRITE(5,231) NUM(I1)
    ENDIF  
    IF (INDEX_S(I1).EQ.4) THEN
        WRITE(2,235) '4(',F11,' ',F21,' ',F24,' ',F14,')'
        WRITE(4,231) NCENTER
        WRITE(5,231) NUM(I1)
    ENDIF
    IF (INDEX_S(I1).EQ.5) THEN
        WRITE(2,235) '4(',F14,' ',F24,' ',F23,' ',F13,')'
        WRITE(4,231) NCENTER
        WRITE(5,231) NUM(I1)
    ENDIF
    IF (INDEX_S(I1).EQ.6) THEN
        WRITE(2,235) '4(',F11,' ',F12,' ',F22,' ',F21,')'
        WRITE(4,231) NCENTER
        WRITE(5,231) NUM(I1)
    ENDIF       
    ENDIF  
    ENDDO
    
    ELSE
    !CASE 1-2-2 
    F11=(CNL+1)*(CNI+1)*(KP-1)+(CNI+1)*(JP-1)+IP-1
    F12=F11+1
    F13=F12+CNI+1
    F14=F13-1
    F21=F11+(CNL+1)*(CNI+1)
    F22=F21+1
    F23=F22+CNI+1
    F24=F23-1    
    CALL CAL_FIJ_RECT(F11,F12,F13,F14,F21,F22,F23,F24)
    NUM(7)=CELL_ID_RECT(IP,JP,KP)
    NUM(1)=CELL_ID_RECT(IP,JP,KP-1) !FRONT FACE
    NUM(2)=CELL_ID_RECT(IP,JP,KP+1) !BACK  FACE
    NUM(3)=CELL_ID_RECT(IP+1,JP,KP) !LEFT  FACE
    NUM(4)=-1 !RIGHT FACE
    NUM(5)=CELL_ID_RECT(IP,JP+1,KP) !UP    FACE
    NUM(6)=CELL_ID_RECT(IP,JP-1,KP) !DOWN  FACE
    CALL BUBBLE_SORT(NUM,INDEX_S)
    DO I1=1,7
    IF (NUM(I1).GT.NCENTER) THEN   
    IF (INDEX_S(I1).EQ.1) THEN
        WRITE(2,235) '4(',F11,' ',F14,' ',F13,' ',F12,')'
        WRITE(4,231) NCENTER
        WRITE(5,231) NUM(I1)
    ENDIF
    IF (INDEX_S(I1).EQ.2) THEN
        WRITE(2,235) '4(',F21,' ',F22,' ',F23,' ',F24,')'
        WRITE(4,231) NCENTER
        WRITE(5,231) NUM(I1)
    ENDIF
    IF (INDEX_S(I1).EQ.3) THEN
        WRITE(2,235) '4(',F12,' ',F13,' ',F23,' ',F22,')'
        WRITE(4,231) NCENTER
        WRITE(5,231) NUM(I1)
    ENDIF  
    IF (INDEX_S(I1).EQ.4) THEN
        WRITE(2,235) '4(',F11,' ',F21,' ',F24,' ',F14,')'
        WRITE(4,231) NCENTER
        WRITE(5,231) NUM(I1)
    ENDIF
    IF (INDEX_S(I1).EQ.5) THEN
        WRITE(2,235) '4(',F14,' ',F24,' ',F23,' ',F13,')'
        WRITE(4,231) NCENTER
        WRITE(5,231) NUM(I1)
    ENDIF
    IF (INDEX_S(I1).EQ.6) THEN
        WRITE(2,235) '4(',F11,' ',F12,' ',F22,' ',F21,')'
        WRITE(4,231) NCENTER
        WRITE(5,231) NUM(I1)
    ENDIF       
    ENDIF  
    ENDDO
    
    ENDIF    
    ENDIF    
    ELSEIF (IP.EQ.CNI) THEN
    IF (JP.EQ.1) THEN
    IF (KP.EQ.1) THEN
    !CASE 3-1-1  
    F11=(CNL+1)*(CNI+1)*(KP-1)+(CNI+1)*(JP-1)+IP-1
    F12=F11+1
    F13=F12+CNI+1
    F14=F13-1
    F21=F11+(CNL+1)*(CNI+1)
    F22=F21+1
    F23=F22+CNI+1
    F24=F23-1  
    CALL CAL_FIJ_RECT(F11,F12,F13,F14,F21,F22,F23,F24)
    NUM(7)=CELL_ID_RECT(IP,JP,KP)
    NUM(1)=-1 !FRONT FACE
    NUM(2)=CELL_ID_RECT(IP,JP,KP+1) !BACK  FACE
    NUM(3)=-1 !LEFT  FACE
    NUM(4)=CELL_ID_RECT(IP-1,JP,KP) !RIGHT FACE
    NUM(5)=CELL_ID_RECT(IP,JP+1,KP) !UP    FACE
    NUM(6)=-1 !DOWN  FACE
    CALL BUBBLE_SORT(NUM,INDEX_S)
    DO I1=1,7
    IF (NUM(I1).GT.NCENTER) THEN   
    IF (INDEX_S(I1).EQ.1) THEN
        WRITE(2,235) '4(',F11,' ',F14,' ',F13,' ',F12,')'
        WRITE(4,231) NCENTER
        WRITE(5,231) NUM(I1)
    ENDIF
    IF (INDEX_S(I1).EQ.2) THEN
        WRITE(2,235) '4(',F21,' ',F22,' ',F23,' ',F24,')'
        WRITE(4,231) NCENTER
        WRITE(5,231) NUM(I1)
    ENDIF
    IF (INDEX_S(I1).EQ.3) THEN
        WRITE(2,235) '4(',F12,' ',F13,' ',F23,' ',F22,')'
        WRITE(4,231) NCENTER
        WRITE(5,231) NUM(I1)
    ENDIF  
    IF (INDEX_S(I1).EQ.4) THEN
        WRITE(2,235) '4(',F11,' ',F21,' ',F24,' ',F14,')'
        WRITE(4,231) NCENTER
        WRITE(5,231) NUM(I1)
    ENDIF
    IF (INDEX_S(I1).EQ.5) THEN
        WRITE(2,235) '4(',F14,' ',F24,' ',F23,' ',F13,')'
        WRITE(4,231) NCENTER
        WRITE(5,231) NUM(I1)
    ENDIF
    IF (INDEX_S(I1).EQ.6) THEN
        WRITE(2,235) '4(',F11,' ',F12,' ',F22,' ',F21,')'
        WRITE(4,231) NCENTER
        WRITE(5,231) NUM(I1)
    ENDIF       
    ENDIF  
    ENDDO
    
    ELSEIF (KP.EQ.CNL) THEN
    !CASE 3-1-3 
    F11=(CNL+1)*(CNI+1)*(KP-1)+(CNI+1)*(JP-1)+IP-1
    F12=F11+1
    F13=F12+CNI+1
    F14=F13-1
    F21=F11+(CNL+1)*(CNI+1)
    F22=F21+1
    F23=F22+CNI+1
    F24=F23-1   
    CALL CAL_FIJ_RECT(F11,F12,F13,F14,F21,F22,F23,F24)
    NUM(7)=CELL_ID_RECT(IP,JP,KP)
    NUM(1)=CELL_ID_RECT(IP,JP,KP-1) !FRONT FACE
    NUM(2)=-1 !BACK  FACE
    NUM(3)=-1 !LEFT  FACE
    NUM(4)=CELL_ID_RECT(IP-1,JP,KP) !RIGHT FACE
    NUM(5)=CELL_ID_RECT(IP,JP+1,KP) !UP    FACE
    NUM(6)=-1 !DOWN  FACE
    CALL BUBBLE_SORT(NUM,INDEX_S)
    DO I1=1,7
    IF (NUM(I1).GT.NCENTER) THEN   
    IF (INDEX_S(I1).EQ.1) THEN
        WRITE(2,235) '4(',F11,' ',F14,' ',F13,' ',F12,')'
        WRITE(4,231) NCENTER
        WRITE(5,231) NUM(I1)
    ENDIF
    IF (INDEX_S(I1).EQ.2) THEN
        WRITE(2,235) '4(',F21,' ',F22,' ',F23,' ',F24,')'
        WRITE(4,231) NCENTER
        WRITE(5,231) NUM(I1)
    ENDIF
    IF (INDEX_S(I1).EQ.3) THEN
        WRITE(2,235) '4(',F12,' ',F13,' ',F23,' ',F22,')'
        WRITE(4,231) NCENTER
        WRITE(5,231) NUM(I1)
    ENDIF  
    IF (INDEX_S(I1).EQ.4) THEN
        WRITE(2,235) '4(',F11,' ',F21,' ',F24,' ',F14,')'
        WRITE(4,231) NCENTER
        WRITE(5,231) NUM(I1)
    ENDIF
    IF (INDEX_S(I1).EQ.5) THEN
        WRITE(2,235) '4(',F14,' ',F24,' ',F23,' ',F13,')'
        WRITE(4,231) NCENTER
        WRITE(5,231) NUM(I1)
    ENDIF
    IF (INDEX_S(I1).EQ.6) THEN
        WRITE(2,235) '4(',F11,' ',F12,' ',F22,' ',F21,')'
        WRITE(4,231) NCENTER
        WRITE(5,231) NUM(I1)
    ENDIF       
    ENDIF  
    ENDDO
    
    ELSE
    !CASE 3-1-2 
    F11=(CNL+1)*(CNI+1)*(KP-1)+(CNI+1)*(JP-1)+IP-1
    F12=F11+1
    F13=F12+CNI+1
    F14=F13-1
    F21=F11+(CNL+1)*(CNI+1)
    F22=F21+1
    F23=F22+CNI+1
    F24=F23-1    
    CALL CAL_FIJ_RECT(F11,F12,F13,F14,F21,F22,F23,F24)
    NUM(7)=CELL_ID_RECT(IP,JP,KP)
    NUM(1)=CELL_ID_RECT(IP,JP,KP-1) !FRONT FACE
    NUM(2)=CELL_ID_RECT(IP,JP,KP+1) !BACK  FACE
    NUM(3)=-1 !LEFT  FACE
    NUM(4)=CELL_ID_RECT(IP-1,JP,KP) !RIGHT FACE
    NUM(5)=CELL_ID_RECT(IP,JP+1,KP) !UP    FACE
    NUM(6)=-1 !DOWN  FACE
    CALL BUBBLE_SORT(NUM,INDEX_S)
    DO I1=1,7
    IF (NUM(I1).GT.NCENTER) THEN   
    IF (INDEX_S(I1).EQ.1) THEN
        WRITE(2,235) '4(',F11,' ',F14,' ',F13,' ',F12,')'
        WRITE(4,231) NCENTER
        WRITE(5,231) NUM(I1)
    ENDIF
    IF (INDEX_S(I1).EQ.2) THEN
        WRITE(2,235) '4(',F21,' ',F22,' ',F23,' ',F24,')'
        WRITE(4,231) NCENTER
        WRITE(5,231) NUM(I1)
    ENDIF
    IF (INDEX_S(I1).EQ.3) THEN
        WRITE(2,235) '4(',F12,' ',F13,' ',F23,' ',F22,')'
        WRITE(4,231) NCENTER
        WRITE(5,231) NUM(I1)
    ENDIF  
    IF (INDEX_S(I1).EQ.4) THEN
        WRITE(2,235) '4(',F11,' ',F21,' ',F24,' ',F14,')'
        WRITE(4,231) NCENTER
        WRITE(5,231) NUM(I1)
    ENDIF
    IF (INDEX_S(I1).EQ.5) THEN
        WRITE(2,235) '4(',F14,' ',F24,' ',F23,' ',F13,')'
        WRITE(4,231) NCENTER
        WRITE(5,231) NUM(I1)
    ENDIF
    IF (INDEX_S(I1).EQ.6) THEN
        WRITE(2,235) '4(',F11,' ',F12,' ',F22,' ',F21,')'
        WRITE(4,231) NCENTER
        WRITE(5,231) NUM(I1)
    ENDIF       
    ENDIF  
    ENDDO
    
    ENDIF    
    ELSEIF (JP.EQ.CNL) THEN
    IF (KP.EQ.1) THEN
    !CASE 3-3-1  
    F11=(CNL+1)*(CNI+1)*(KP-1)+(CNI+1)*(JP-1)+IP-1
    F12=F11+1
    F13=F12+CNI+1
    F14=F13-1
    F21=F11+(CNL+1)*(CNI+1)
    F22=F21+1
    F23=F22+CNI+1
    F24=F23-1   
    CALL CAL_FIJ_RECT(F11,F12,F13,F14,F21,F22,F23,F24)
    NUM(7)=CELL_ID_RECT(IP,JP,KP)
    NUM(1)=-1 !FRONT FACE
    NUM(2)=CELL_ID_RECT(IP,JP,KP+1) !BACK  FACE
    NUM(3)=-1 !LEFT  FACE
    NUM(4)=CELL_ID_RECT(IP-1,JP,KP) !RIGHT FACE
    NUM(5)=-1 !UP    FACE
    NUM(6)=CELL_ID_RECT(IP,JP-1,KP) !DOWN  FACE
    CALL BUBBLE_SORT(NUM,INDEX_S)
    DO I1=1,7
    IF (NUM(I1).GT.NCENTER) THEN   
    IF (INDEX_S(I1).EQ.1) THEN
        WRITE(2,235) '4(',F11,' ',F14,' ',F13,' ',F12,')'
        WRITE(4,231) NCENTER
        WRITE(5,231) NUM(I1)
    ENDIF
    IF (INDEX_S(I1).EQ.2) THEN
        WRITE(2,235) '4(',F21,' ',F22,' ',F23,' ',F24,')'
        WRITE(4,231) NCENTER
        WRITE(5,231) NUM(I1)
    ENDIF
    IF (INDEX_S(I1).EQ.3) THEN
        WRITE(2,235) '4(',F12,' ',F13,' ',F23,' ',F22,')'
        WRITE(4,231) NCENTER
        WRITE(5,231) NUM(I1)
    ENDIF  
    IF (INDEX_S(I1).EQ.4) THEN
        WRITE(2,235) '4(',F11,' ',F21,' ',F24,' ',F14,')'
        WRITE(4,231) NCENTER
        WRITE(5,231) NUM(I1)
    ENDIF
    IF (INDEX_S(I1).EQ.5) THEN
        WRITE(2,235) '4(',F14,' ',F24,' ',F23,' ',F13,')'
        WRITE(4,231) NCENTER
        WRITE(5,231) NUM(I1)
    ENDIF
    IF (INDEX_S(I1).EQ.6) THEN
        WRITE(2,235) '4(',F11,' ',F12,' ',F22,' ',F21,')'
        WRITE(4,231) NCENTER
        WRITE(5,231) NUM(I1)
    ENDIF       
    ENDIF  
    ENDDO
    
    ELSEIF (KP.EQ.CNL) THEN
    !CASE 3-3-3 
    F11=(CNL+1)*(CNI+1)*(KP-1)+(CNI+1)*(JP-1)+IP-1
    F12=F11+1
    F13=F12+CNI+1
    F14=F13-1
    F21=F11+(CNL+1)*(CNI+1)
    F22=F21+1
    F23=F22+CNI+1
    F24=F23-1 
    CALL CAL_FIJ_RECT(F11,F12,F13,F14,F21,F22,F23,F24)
    NUM(7)=CELL_ID_RECT(IP,JP,KP)
    NUM(1)=CELL_ID_RECT(IP,JP,KP-1) !FRONT FACE
    NUM(2)=-1 !BACK  FACE
    NUM(3)=-1 !LEFT  FACE
    NUM(4)=CELL_ID_RECT(IP-1,JP,KP) !RIGHT FACE
    NUM(5)=-1 !UP    FACE
    NUM(6)=CELL_ID_RECT(IP,JP-1,KP) !DOWN  FACE
    CALL BUBBLE_SORT(NUM,INDEX_S)
    DO I1=1,7
    IF (NUM(I1).GT.NCENTER) THEN   
    IF (INDEX_S(I1).EQ.1) THEN
        WRITE(2,235) '4(',F11,' ',F14,' ',F13,' ',F12,')'
        WRITE(4,231) NCENTER
        WRITE(5,231) NUM(I1)
    ENDIF
    IF (INDEX_S(I1).EQ.2) THEN
        WRITE(2,235) '4(',F21,' ',F22,' ',F23,' ',F24,')'
        WRITE(4,231) NCENTER
        WRITE(5,231) NUM(I1)
    ENDIF
    IF (INDEX_S(I1).EQ.3) THEN
        WRITE(2,235) '4(',F12,' ',F13,' ',F23,' ',F22,')'
        WRITE(4,231) NCENTER
        WRITE(5,231) NUM(I1)
    ENDIF  
    IF (INDEX_S(I1).EQ.4) THEN
        WRITE(2,235) '4(',F11,' ',F21,' ',F24,' ',F14,')'
        WRITE(4,231) NCENTER
        WRITE(5,231) NUM(I1)
    ENDIF
    IF (INDEX_S(I1).EQ.5) THEN
        WRITE(2,235) '4(',F14,' ',F24,' ',F23,' ',F13,')'
        WRITE(4,231) NCENTER
        WRITE(5,231) NUM(I1)
    ENDIF
    IF (INDEX_S(I1).EQ.6) THEN
        WRITE(2,235) '4(',F11,' ',F12,' ',F22,' ',F21,')'
        WRITE(4,231) NCENTER
        WRITE(5,231) NUM(I1)
    ENDIF       
    ENDIF  
    ENDDO
    
    ELSE
    !CASE 3-3-2 
    F11=(CNL+1)*(CNI+1)*(KP-1)+(CNI+1)*(JP-1)+IP-1
    F12=F11+1
    F13=F12+CNI+1
    F14=F13-1
    F21=F11+(CNL+1)*(CNI+1)
    F22=F21+1
    F23=F22+CNI+1
    F24=F23-1    
    CALL CAL_FIJ_RECT(F11,F12,F13,F14,F21,F22,F23,F24)
    NUM(7)=CELL_ID_RECT(IP,JP,KP)
    NUM(1)=CELL_ID_RECT(IP,JP,KP-1) !FRONT FACE
    NUM(2)=CELL_ID_RECT(IP,JP,KP+1) !BACK  FACE
    NUM(3)=-1 !LEFT  FACE
    NUM(4)=CELL_ID_RECT(IP-1,JP,KP) !RIGHT FACE
    NUM(5)=-1 !UP    FACE
    NUM(6)=CELL_ID_RECT(IP,JP-1,KP) !DOWN  FACE
    CALL BUBBLE_SORT(NUM,INDEX_S)
    DO I1=1,7
    IF (NUM(I1).GT.NCENTER) THEN   
    IF (INDEX_S(I1).EQ.1) THEN
        WRITE(2,235) '4(',F11,' ',F14,' ',F13,' ',F12,')'
        WRITE(4,231) NCENTER
        WRITE(5,231) NUM(I1)
    ENDIF
    IF (INDEX_S(I1).EQ.2) THEN
        WRITE(2,235) '4(',F21,' ',F22,' ',F23,' ',F24,')'
        WRITE(4,231) NCENTER
        WRITE(5,231) NUM(I1)
    ENDIF
    IF (INDEX_S(I1).EQ.3) THEN
        WRITE(2,235) '4(',F12,' ',F13,' ',F23,' ',F22,')'
        WRITE(4,231) NCENTER
        WRITE(5,231) NUM(I1)
    ENDIF  
    IF (INDEX_S(I1).EQ.4) THEN
        WRITE(2,235) '4(',F11,' ',F21,' ',F24,' ',F14,')'
        WRITE(4,231) NCENTER
        WRITE(5,231) NUM(I1)
    ENDIF
    IF (INDEX_S(I1).EQ.5) THEN
        WRITE(2,235) '4(',F14,' ',F24,' ',F23,' ',F13,')'
        WRITE(4,231) NCENTER
        WRITE(5,231) NUM(I1)
    ENDIF
    IF (INDEX_S(I1).EQ.6) THEN
        WRITE(2,235) '4(',F11,' ',F12,' ',F22,' ',F21,')'
        WRITE(4,231) NCENTER
        WRITE(5,231) NUM(I1)
    ENDIF       
    ENDIF  
    ENDDO
    
    ENDIF    
    ELSE
    IF (KP.EQ.1) THEN
    !CASE 3-2-1  
    F11=(CNL+1)*(CNI+1)*(KP-1)+(CNI+1)*(JP-1)+IP-1
    F12=F11+1
    F13=F12+CNI+1
    F14=F13-1
    F21=F11+(CNL+1)*(CNI+1)
    F22=F21+1
    F23=F22+CNI+1
    F24=F23-1   
    CALL CAL_FIJ_RECT(F11,F12,F13,F14,F21,F22,F23,F24)
    NUM(7)=CELL_ID_RECT(IP,JP,KP)
    NUM(1)=-1 !FRONT FACE
    NUM(2)=CELL_ID_RECT(IP,JP,KP+1) !BACK  FACE
    NUM(3)=-1 !LEFT  FACE
    NUM(4)=CELL_ID_RECT(IP-1,JP,KP) !RIGHT FACE
    NUM(5)=CELL_ID_RECT(IP,JP+1,KP) !UP    FACE
    NUM(6)=CELL_ID_RECT(IP,JP-1,KP) !DOWN  FACE
    CALL BUBBLE_SORT(NUM,INDEX_S)
    DO I1=1,7
    IF (NUM(I1).GT.NCENTER) THEN   
    IF (INDEX_S(I1).EQ.1) THEN
        WRITE(2,235) '4(',F11,' ',F14,' ',F13,' ',F12,')'
        WRITE(4,231) NCENTER
        WRITE(5,231) NUM(I1)
    ENDIF
    IF (INDEX_S(I1).EQ.2) THEN
        WRITE(2,235) '4(',F21,' ',F22,' ',F23,' ',F24,')'
        WRITE(4,231) NCENTER
        WRITE(5,231) NUM(I1)
    ENDIF
    IF (INDEX_S(I1).EQ.3) THEN
        WRITE(2,235) '4(',F12,' ',F13,' ',F23,' ',F22,')'
        WRITE(4,231) NCENTER
        WRITE(5,231) NUM(I1)
    ENDIF  
    IF (INDEX_S(I1).EQ.4) THEN
        WRITE(2,235) '4(',F11,' ',F21,' ',F24,' ',F14,')'
        WRITE(4,231) NCENTER
        WRITE(5,231) NUM(I1)
    ENDIF
    IF (INDEX_S(I1).EQ.5) THEN
        WRITE(2,235) '4(',F14,' ',F24,' ',F23,' ',F13,')'
        WRITE(4,231) NCENTER
        WRITE(5,231) NUM(I1)
    ENDIF
    IF (INDEX_S(I1).EQ.6) THEN
        WRITE(2,235) '4(',F11,' ',F12,' ',F22,' ',F21,')'
        WRITE(4,231) NCENTER
        WRITE(5,231) NUM(I1)
    ENDIF       
    ENDIF  
    ENDDO
    
    ELSEIF (KP.EQ.CNL) THEN
    !CASE 3-2-3  
    F11=(CNL+1)*(CNI+1)*(KP-1)+(CNI+1)*(JP-1)+IP-1
    F12=F11+1
    F13=F12+CNI+1
    F14=F13-1
    F21=F11+(CNL+1)*(CNI+1)
    F22=F21+1
    F23=F22+CNI+1
    F24=F23-1   
    CALL CAL_FIJ_RECT(F11,F12,F13,F14,F21,F22,F23,F24)
    NUM(7)=CELL_ID_RECT(IP,JP,KP)
    NUM(1)=CELL_ID_RECT(IP,JP,KP-1) !FRONT FACE
    NUM(2)=-1 !BACK  FACE
    NUM(3)=-1 !LEFT  FACE
    NUM(4)=CELL_ID_RECT(IP-1,JP,KP) !RIGHT FACE
    NUM(5)=CELL_ID_RECT(IP,JP+1,KP) !UP    FACE
    NUM(6)=CELL_ID_RECT(IP,JP-1,KP) !DOWN  FACE
    CALL BUBBLE_SORT(NUM,INDEX_S)
    DO I1=1,7
    IF (NUM(I1).GT.NCENTER) THEN   
    IF (INDEX_S(I1).EQ.1) THEN
        WRITE(2,235) '4(',F11,' ',F14,' ',F13,' ',F12,')'
        WRITE(4,231) NCENTER
        WRITE(5,231) NUM(I1)
    ENDIF
    IF (INDEX_S(I1).EQ.2) THEN
        WRITE(2,235) '4(',F21,' ',F22,' ',F23,' ',F24,')'
        WRITE(4,231) NCENTER
        WRITE(5,231) NUM(I1)
    ENDIF
    IF (INDEX_S(I1).EQ.3) THEN
        WRITE(2,235) '4(',F12,' ',F13,' ',F23,' ',F22,')'
        WRITE(4,231) NCENTER
        WRITE(5,231) NUM(I1)
    ENDIF  
    IF (INDEX_S(I1).EQ.4) THEN
        WRITE(2,235) '4(',F11,' ',F21,' ',F24,' ',F14,')'
        WRITE(4,231) NCENTER
        WRITE(5,231) NUM(I1)
    ENDIF
    IF (INDEX_S(I1).EQ.5) THEN
        WRITE(2,235) '4(',F14,' ',F24,' ',F23,' ',F13,')'
        WRITE(4,231) NCENTER
        WRITE(5,231) NUM(I1)
    ENDIF
    IF (INDEX_S(I1).EQ.6) THEN
        WRITE(2,235) '4(',F11,' ',F12,' ',F22,' ',F21,')'
        WRITE(4,231) NCENTER
        WRITE(5,231) NUM(I1)
    ENDIF       
    ENDIF  
    ENDDO
    
    ELSE
    !CASE 3-2-2  
    F11=(CNL+1)*(CNI+1)*(KP-1)+(CNI+1)*(JP-1)+IP-1
    F12=F11+1
    F13=F12+CNI+1
    F14=F13-1
    F21=F11+(CNL+1)*(CNI+1)
    F22=F21+1
    F23=F22+CNI+1
    F24=F23-1    
    CALL CAL_FIJ_RECT(F11,F12,F13,F14,F21,F22,F23,F24)
    NUM(7)=CELL_ID_RECT(IP,JP,KP)
    NUM(1)=CELL_ID_RECT(IP,JP,KP-1) !FRONT FACE
    NUM(2)=CELL_ID_RECT(IP,JP,KP+1) !BACK  FACE
    NUM(3)=-1 !LEFT  FACE
    NUM(4)=CELL_ID_RECT(IP-1,JP,KP) !RIGHT FACE
    NUM(5)=CELL_ID_RECT(IP,JP+1,KP) !UP    FACE
    NUM(6)=CELL_ID_RECT(IP,JP-1,KP) !DOWN  FACE
    CALL BUBBLE_SORT(NUM,INDEX_S)
    DO I1=1,7
    IF (NUM(I1).GT.NCENTER) THEN   
    IF (INDEX_S(I1).EQ.1) THEN
        WRITE(2,235) '4(',F11,' ',F14,' ',F13,' ',F12,')'
        WRITE(4,231) NCENTER
        WRITE(5,231) NUM(I1)
    ENDIF
    IF (INDEX_S(I1).EQ.2) THEN
        WRITE(2,235) '4(',F21,' ',F22,' ',F23,' ',F24,')'
        WRITE(4,231) NCENTER
        WRITE(5,231) NUM(I1)
    ENDIF
    IF (INDEX_S(I1).EQ.3) THEN
        WRITE(2,235) '4(',F12,' ',F13,' ',F23,' ',F22,')'
        WRITE(4,231) NCENTER
        WRITE(5,231) NUM(I1)
    ENDIF  
    IF (INDEX_S(I1).EQ.4) THEN
        WRITE(2,235) '4(',F11,' ',F21,' ',F24,' ',F14,')'
        WRITE(4,231) NCENTER
        WRITE(5,231) NUM(I1)
    ENDIF
    IF (INDEX_S(I1).EQ.5) THEN
        WRITE(2,235) '4(',F14,' ',F24,' ',F23,' ',F13,')'
        WRITE(4,231) NCENTER
        WRITE(5,231) NUM(I1)
    ENDIF
    IF (INDEX_S(I1).EQ.6) THEN
        WRITE(2,235) '4(',F11,' ',F12,' ',F22,' ',F21,')'
        WRITE(4,231) NCENTER
        WRITE(5,231) NUM(I1)
    ENDIF       
    ENDIF  
    ENDDO
    
    ENDIF    
    ENDIF    
    ELSE
    IF (JP.EQ.1) THEN
    IF (KP.EQ.1) THEN
    !CASE 2-1-1 
    F11=(CNL+1)*(CNI+1)*(KP-1)+(CNI+1)*(JP-1)+IP-1
    F12=F11+1
    F13=F12+CNI+1
    F14=F13-1
    F21=F11+(CNL+1)*(CNI+1)
    F22=F21+1
    F23=F22+CNI+1
    F24=F23-1 
    CALL CAL_FIJ_RECT(F11,F12,F13,F14,F21,F22,F23,F24)
    NUM(7)=CELL_ID_RECT(IP,JP,KP)
    NUM(1)=-1 !FRONT FACE
    NUM(2)=CELL_ID_RECT(IP,JP,KP+1) !BACK  FACE
    NUM(3)=CELL_ID_RECT(IP+1,JP,KP) !LEFT  FACE
    NUM(4)=CELL_ID_RECT(IP-1,JP,KP) !RIGHT FACE
    NUM(5)=CELL_ID_RECT(IP,JP+1,KP) !UP    FACE
    NUM(6)=-1 !DOWN  FACE
    CALL BUBBLE_SORT(NUM,INDEX_S)
    DO I1=1,7
    IF (NUM(I1).GT.NCENTER) THEN   
    IF (INDEX_S(I1).EQ.1) THEN
        WRITE(2,235) '4(',F11,' ',F14,' ',F13,' ',F12,')'
        WRITE(4,231) NCENTER
        WRITE(5,231) NUM(I1)
    ENDIF
    IF (INDEX_S(I1).EQ.2) THEN
        WRITE(2,235) '4(',F21,' ',F22,' ',F23,' ',F24,')'
        WRITE(4,231) NCENTER
        WRITE(5,231) NUM(I1)
    ENDIF
    IF (INDEX_S(I1).EQ.3) THEN
        WRITE(2,235) '4(',F12,' ',F13,' ',F23,' ',F22,')'
        WRITE(4,231) NCENTER
        WRITE(5,231) NUM(I1)
    ENDIF  
    IF (INDEX_S(I1).EQ.4) THEN
        WRITE(2,235) '4(',F11,' ',F21,' ',F24,' ',F14,')'
        WRITE(4,231) NCENTER
        WRITE(5,231) NUM(I1)
    ENDIF
    IF (INDEX_S(I1).EQ.5) THEN
        WRITE(2,235) '4(',F14,' ',F24,' ',F23,' ',F13,')'
        WRITE(4,231) NCENTER
        WRITE(5,231) NUM(I1)
    ENDIF
    IF (INDEX_S(I1).EQ.6) THEN
        WRITE(2,235) '4(',F11,' ',F12,' ',F22,' ',F21,')'
        WRITE(4,231) NCENTER
        WRITE(5,231) NUM(I1)
    ENDIF       
    ENDIF  
    ENDDO
    
    ELSEIF (KP.EQ.CNL) THEN
    !CASE 2-1-3 
    F11=(CNL+1)*(CNI+1)*(KP-1)+(CNI+1)*(JP-1)+IP-1
    F12=F11+1
    F13=F12+CNI+1
    F14=F13-1
    F21=F11+(CNL+1)*(CNI+1)
    F22=F21+1
    F23=F22+CNI+1
    F24=F23-1  
    CALL CAL_FIJ_RECT(F11,F12,F13,F14,F21,F22,F23,F24)
    NUM(7)=CELL_ID_RECT(IP,JP,KP)
    NUM(1)=CELL_ID_RECT(IP,JP,KP-1) !FRONT FACE
    NUM(2)=-1 !BACK  FACE
    NUM(3)=CELL_ID_RECT(IP+1,JP,KP) !LEFT  FACE
    NUM(4)=CELL_ID_RECT(IP-1,JP,KP) !RIGHT FACE
    NUM(5)=CELL_ID_RECT(IP,JP+1,KP) !UP    FACE
    NUM(6)=-1 !DOWN  FACE
    CALL BUBBLE_SORT(NUM,INDEX_S)
    DO I1=1,7
    IF (NUM(I1).GT.NCENTER) THEN   
    IF (INDEX_S(I1).EQ.1) THEN
        WRITE(2,235) '4(',F11,' ',F14,' ',F13,' ',F12,')'
        WRITE(4,231) NCENTER
        WRITE(5,231) NUM(I1)
    ENDIF
    IF (INDEX_S(I1).EQ.2) THEN
        WRITE(2,235) '4(',F21,' ',F22,' ',F23,' ',F24,')'
        WRITE(4,231) NCENTER
        WRITE(5,231) NUM(I1)
    ENDIF
    IF (INDEX_S(I1).EQ.3) THEN
        WRITE(2,235) '4(',F12,' ',F13,' ',F23,' ',F22,')'
        WRITE(4,231) NCENTER
        WRITE(5,231) NUM(I1)
    ENDIF  
    IF (INDEX_S(I1).EQ.4) THEN
        WRITE(2,235) '4(',F11,' ',F21,' ',F24,' ',F14,')'
        WRITE(4,231) NCENTER
        WRITE(5,231) NUM(I1)
    ENDIF
    IF (INDEX_S(I1).EQ.5) THEN
        WRITE(2,235) '4(',F14,' ',F24,' ',F23,' ',F13,')'
        WRITE(4,231) NCENTER
        WRITE(5,231) NUM(I1)
    ENDIF
    IF (INDEX_S(I1).EQ.6) THEN
        WRITE(2,235) '4(',F11,' ',F12,' ',F22,' ',F21,')'
        WRITE(4,231) NCENTER
        WRITE(5,231) NUM(I1)
    ENDIF       
    ENDIF  
    ENDDO
    
    ELSE
    !CASE 2-1-2 
    F11=(CNL+1)*(CNI+1)*(KP-1)+(CNI+1)*(JP-1)+IP-1
    F12=F11+1
    F13=F12+CNI+1
    F14=F13-1
    F21=F11+(CNL+1)*(CNI+1)
    F22=F21+1
    F23=F22+CNI+1
    F24=F23-1    
    CALL CAL_FIJ_RECT(F11,F12,F13,F14,F21,F22,F23,F24)
    NUM(7)=CELL_ID_RECT(IP,JP,KP)
    NUM(1)=CELL_ID_RECT(IP,JP,KP-1) !FRONT FACE
    NUM(2)=CELL_ID_RECT(IP,JP,KP+1) !BACK  FACE
    NUM(3)=CELL_ID_RECT(IP+1,JP,KP) !LEFT  FACE
    NUM(4)=CELL_ID_RECT(IP-1,JP,KP) !RIGHT FACE
    NUM(5)=CELL_ID_RECT(IP,JP+1,KP) !UP    FACE
    NUM(6)=-1 !DOWN  FACE
    CALL BUBBLE_SORT(NUM,INDEX_S)
    DO I1=1,7
    IF (NUM(I1).GT.NCENTER) THEN   
    IF (INDEX_S(I1).EQ.1) THEN
        WRITE(2,235) '4(',F11,' ',F14,' ',F13,' ',F12,')'
        WRITE(4,231) NCENTER
        WRITE(5,231) NUM(I1)
    ENDIF
    IF (INDEX_S(I1).EQ.2) THEN
        WRITE(2,235) '4(',F21,' ',F22,' ',F23,' ',F24,')'
        WRITE(4,231) NCENTER
        WRITE(5,231) NUM(I1)
    ENDIF
    IF (INDEX_S(I1).EQ.3) THEN
        WRITE(2,235) '4(',F12,' ',F13,' ',F23,' ',F22,')'
        WRITE(4,231) NCENTER
        WRITE(5,231) NUM(I1)
    ENDIF  
    IF (INDEX_S(I1).EQ.4) THEN
        WRITE(2,235) '4(',F11,' ',F21,' ',F24,' ',F14,')'
        WRITE(4,231) NCENTER
        WRITE(5,231) NUM(I1)
    ENDIF
    IF (INDEX_S(I1).EQ.5) THEN
        WRITE(2,235) '4(',F14,' ',F24,' ',F23,' ',F13,')'
        WRITE(4,231) NCENTER
        WRITE(5,231) NUM(I1)
    ENDIF
    IF (INDEX_S(I1).EQ.6) THEN
        WRITE(2,235) '4(',F11,' ',F12,' ',F22,' ',F21,')'
        WRITE(4,231) NCENTER
        WRITE(5,231) NUM(I1)
    ENDIF       
    ENDIF  
    ENDDO
    
    ENDIF    
    ELSEIF (JP.EQ.CNL) THEN
    IF (KP.EQ.1) THEN
    !CASE 2-3-1 
    F11=(CNL+1)*(CNI+1)*(KP-1)+(CNI+1)*(JP-1)+IP-1
    F12=F11+1
    F13=F12+CNI+1
    F14=F13-1
    F21=F11+(CNL+1)*(CNI+1)
    F22=F21+1
    F23=F22+CNI+1
    F24=F23-1    
    CALL CAL_FIJ_RECT(F11,F12,F13,F14,F21,F22,F23,F24)
    NUM(7)=CELL_ID_RECT(IP,JP,KP)
    NUM(1)=-1 !FRONT FACE
    NUM(2)=CELL_ID_RECT(IP,JP,KP+1) !BACK  FACE
    NUM(3)=CELL_ID_RECT(IP+1,JP,KP) !LEFT  FACE
    NUM(4)=CELL_ID_RECT(IP-1,JP,KP) !RIGHT FACE
    NUM(5)=-1 !UP    FACE
    NUM(6)=CELL_ID_RECT(IP,JP-1,KP) !DOWN  FACE
    CALL BUBBLE_SORT(NUM,INDEX_S)
    DO I1=1,7
    IF (NUM(I1).GT.NCENTER) THEN   
    IF (INDEX_S(I1).EQ.1) THEN
        WRITE(2,235) '4(',F11,' ',F14,' ',F13,' ',F12,')'
        WRITE(4,231) NCENTER
        WRITE(5,231) NUM(I1)
    ENDIF
    IF (INDEX_S(I1).EQ.2) THEN
        WRITE(2,235) '4(',F21,' ',F22,' ',F23,' ',F24,')'
        WRITE(4,231) NCENTER
        WRITE(5,231) NUM(I1)
    ENDIF
    IF (INDEX_S(I1).EQ.3) THEN
        WRITE(2,235) '4(',F12,' ',F13,' ',F23,' ',F22,')'
        WRITE(4,231) NCENTER
        WRITE(5,231) NUM(I1)
    ENDIF  
    IF (INDEX_S(I1).EQ.4) THEN
        WRITE(2,235) '4(',F11,' ',F21,' ',F24,' ',F14,')'
        WRITE(4,231) NCENTER
        WRITE(5,231) NUM(I1)
    ENDIF
    IF (INDEX_S(I1).EQ.5) THEN
        WRITE(2,235) '4(',F14,' ',F24,' ',F23,' ',F13,')'
        WRITE(4,231) NCENTER
        WRITE(5,231) NUM(I1)
    ENDIF
    IF (INDEX_S(I1).EQ.6) THEN
        WRITE(2,235) '4(',F11,' ',F12,' ',F22,' ',F21,')'
        WRITE(4,231) NCENTER
        WRITE(5,231) NUM(I1)
    ENDIF       
    ENDIF  
    ENDDO
    
    ELSEIF (KP.EQ.CNL) THEN
    !CASE 2-3-3 
    F11=(CNL+1)*(CNI+1)*(KP-1)+(CNI+1)*(JP-1)+IP-1
    F12=F11+1
    F13=F12+CNI+1
    F14=F13-1
    F21=F11+(CNL+1)*(CNI+1)
    F22=F21+1
    F23=F22+CNI+1
    F24=F23-1    
    CALL CAL_FIJ_RECT(F11,F12,F13,F14,F21,F22,F23,F24)
    NUM(7)=CELL_ID_RECT(IP,JP,KP)
    NUM(1)=CELL_ID_RECT(IP,JP,KP-1) !FRONT FACE
    NUM(2)=-1 !BACK  FACE
    NUM(3)=CELL_ID_RECT(IP+1,JP,KP) !LEFT  FACE
    NUM(4)=CELL_ID_RECT(IP-1,JP,KP) !RIGHT FACE
    NUM(5)=-1 !UP    FACE
    NUM(6)=CELL_ID_RECT(IP,JP-1,KP) !DOWN  FACE
    CALL BUBBLE_SORT(NUM,INDEX_S)
    DO I1=1,7
    IF (NUM(I1).GT.NCENTER) THEN   
    IF (INDEX_S(I1).EQ.1) THEN
        WRITE(2,235) '4(',F11,' ',F14,' ',F13,' ',F12,')'
        WRITE(4,231) NCENTER
        WRITE(5,231) NUM(I1)
    ENDIF
    IF (INDEX_S(I1).EQ.2) THEN
        WRITE(2,235) '4(',F21,' ',F22,' ',F23,' ',F24,')'
        WRITE(4,231) NCENTER
        WRITE(5,231) NUM(I1)
    ENDIF
    IF (INDEX_S(I1).EQ.3) THEN
        WRITE(2,235) '4(',F12,' ',F13,' ',F23,' ',F22,')'
        WRITE(4,231) NCENTER
        WRITE(5,231) NUM(I1)
    ENDIF  
    IF (INDEX_S(I1).EQ.4) THEN
        WRITE(2,235) '4(',F11,' ',F21,' ',F24,' ',F14,')'
        WRITE(4,231) NCENTER
        WRITE(5,231) NUM(I1)
    ENDIF
    IF (INDEX_S(I1).EQ.5) THEN
        WRITE(2,235) '4(',F14,' ',F24,' ',F23,' ',F13,')'
        WRITE(4,231) NCENTER
        WRITE(5,231) NUM(I1)
    ENDIF
    IF (INDEX_S(I1).EQ.6) THEN
        WRITE(2,235) '4(',F11,' ',F12,' ',F22,' ',F21,')'
        WRITE(4,231) NCENTER
        WRITE(5,231) NUM(I1)
    ENDIF       
    ENDIF  
    ENDDO
    
    ELSE
    !CASE 2-3-2 
    F11=(CNL+1)*(CNI+1)*(KP-1)+(CNI+1)*(JP-1)+IP-1
    F12=F11+1
    F13=F12+CNI+1
    F14=F13-1
    F21=F11+(CNL+1)*(CNI+1)
    F22=F21+1
    F23=F22+CNI+1
    F24=F23-1    
    CALL CAL_FIJ_RECT(F11,F12,F13,F14,F21,F22,F23,F24)
    NUM(7)=CELL_ID_RECT(IP,JP,KP)
    NUM(1)=CELL_ID_RECT(IP,JP,KP-1) !FRONT FACE
    NUM(2)=CELL_ID_RECT(IP,JP,KP+1) !BACK  FACE
    NUM(3)=CELL_ID_RECT(IP+1,JP,KP) !LEFT  FACE
    NUM(4)=CELL_ID_RECT(IP-1,JP,KP) !RIGHT FACE
    NUM(5)=-1 !UP    FACE
    NUM(6)=CELL_ID_RECT(IP,JP-1,KP) !DOWN  FACE
    CALL BUBBLE_SORT(NUM,INDEX_S)
    DO I1=1,7
    IF (NUM(I1).GT.NCENTER) THEN   
    IF (INDEX_S(I1).EQ.1) THEN
        WRITE(2,235) '4(',F11,' ',F14,' ',F13,' ',F12,')'
        WRITE(4,231) NCENTER
        WRITE(5,231) NUM(I1)
    ENDIF
    IF (INDEX_S(I1).EQ.2) THEN
        WRITE(2,235) '4(',F21,' ',F22,' ',F23,' ',F24,')'
        WRITE(4,231) NCENTER
        WRITE(5,231) NUM(I1)
    ENDIF
    IF (INDEX_S(I1).EQ.3) THEN
        WRITE(2,235) '4(',F12,' ',F13,' ',F23,' ',F22,')'
        WRITE(4,231) NCENTER
        WRITE(5,231) NUM(I1)
    ENDIF  
    IF (INDEX_S(I1).EQ.4) THEN
        WRITE(2,235) '4(',F11,' ',F21,' ',F24,' ',F14,')'
        WRITE(4,231) NCENTER
        WRITE(5,231) NUM(I1)
    ENDIF
    IF (INDEX_S(I1).EQ.5) THEN
        WRITE(2,235) '4(',F14,' ',F24,' ',F23,' ',F13,')'
        WRITE(4,231) NCENTER
        WRITE(5,231) NUM(I1)
    ENDIF
    IF (INDEX_S(I1).EQ.6) THEN
        WRITE(2,235) '4(',F11,' ',F12,' ',F22,' ',F21,')'
        WRITE(4,231) NCENTER
        WRITE(5,231) NUM(I1)
    ENDIF       
    ENDIF  
    ENDDO
    
    ENDIF    
    ELSE
    IF (KP.EQ.1) THEN
    !CASE 2-2-1 
    F11=(CNL+1)*(CNI+1)*(KP-1)+(CNI+1)*(JP-1)+IP-1
    F12=F11+1
    F13=F12+CNI+1
    F14=F13-1
    F21=F11+(CNL+1)*(CNI+1)
    F22=F21+1
    F23=F22+CNI+1
    F24=F23-1   
    CALL CAL_FIJ_RECT(F11,F12,F13,F14,F21,F22,F23,F24)
    NUM(7)=CELL_ID_RECT(IP,JP,KP)
    NUM(1)=-1 !FRONT FACE
    NUM(2)=CELL_ID_RECT(IP,JP,KP+1) !BACK  FACE
    NUM(3)=CELL_ID_RECT(IP+1,JP,KP) !LEFT  FACE
    NUM(4)=CELL_ID_RECT(IP-1,JP,KP) !RIGHT FACE
    NUM(5)=CELL_ID_RECT(IP,JP+1,KP) !UP    FACE
    NUM(6)=CELL_ID_RECT(IP,JP-1,KP) !DOWN  FACE
    CALL BUBBLE_SORT(NUM,INDEX_S)
    DO I1=1,7
    IF (NUM(I1).GT.NCENTER) THEN   
    IF (INDEX_S(I1).EQ.1) THEN
        WRITE(2,235) '4(',F11,' ',F14,' ',F13,' ',F12,')'
        WRITE(4,231) NCENTER
        WRITE(5,231) NUM(I1)
    ENDIF
    IF (INDEX_S(I1).EQ.2) THEN
        WRITE(2,235) '4(',F21,' ',F22,' ',F23,' ',F24,')'
        WRITE(4,231) NCENTER
        WRITE(5,231) NUM(I1)
    ENDIF
    IF (INDEX_S(I1).EQ.3) THEN
        WRITE(2,235) '4(',F12,' ',F13,' ',F23,' ',F22,')'
        WRITE(4,231) NCENTER
        WRITE(5,231) NUM(I1)
    ENDIF  
    IF (INDEX_S(I1).EQ.4) THEN
        WRITE(2,235) '4(',F11,' ',F21,' ',F24,' ',F14,')'
        WRITE(4,231) NCENTER
        WRITE(5,231) NUM(I1)
    ENDIF
    IF (INDEX_S(I1).EQ.5) THEN
        WRITE(2,235) '4(',F14,' ',F24,' ',F23,' ',F13,')'
        WRITE(4,231) NCENTER
        WRITE(5,231) NUM(I1)
    ENDIF
    IF (INDEX_S(I1).EQ.6) THEN
        WRITE(2,235) '4(',F11,' ',F12,' ',F22,' ',F21,')'
        WRITE(4,231) NCENTER
        WRITE(5,231) NUM(I1)
    ENDIF       
    ENDIF  
    ENDDO
    
    ELSEIF (KP.EQ.CNL) THEN
    !CASE 2-2-3  
    F11=(CNL+1)*(CNI+1)*(KP-1)+(CNI+1)*(JP-1)+IP-1
    F12=F11+1
    F13=F12+CNI+1
    F14=F13-1
    F21=F11+(CNL+1)*(CNI+1)
    F22=F21+1
    F23=F22+CNI+1
    F24=F23-1   
    CALL CAL_FIJ_RECT(F11,F12,F13,F14,F21,F22,F23,F24)
    NUM(7)=CELL_ID_RECT(IP,JP,KP)
    NUM(1)=CELL_ID_RECT(IP,JP,KP-1) !FRONT FACE
    NUM(2)=-1 !BACK  FACE
    NUM(3)=CELL_ID_RECT(IP+1,JP,KP) !LEFT  FACE
    NUM(4)=CELL_ID_RECT(IP-1,JP,KP) !RIGHT FACE
    NUM(5)=CELL_ID_RECT(IP,JP+1,KP) !UP    FACE
    NUM(6)=CELL_ID_RECT(IP,JP-1,KP) !DOWN  FACE
    CALL BUBBLE_SORT(NUM,INDEX_S)
    DO I1=1,7
    IF (NUM(I1).GT.NCENTER) THEN   
    IF (INDEX_S(I1).EQ.1) THEN
        WRITE(2,235) '4(',F11,' ',F14,' ',F13,' ',F12,')'
        WRITE(4,231) NCENTER
        WRITE(5,231) NUM(I1)
    ENDIF
    IF (INDEX_S(I1).EQ.2) THEN
        WRITE(2,235) '4(',F21,' ',F22,' ',F23,' ',F24,')'
        WRITE(4,231) NCENTER
        WRITE(5,231) NUM(I1)
    ENDIF
    IF (INDEX_S(I1).EQ.3) THEN
        WRITE(2,235) '4(',F12,' ',F13,' ',F23,' ',F22,')'
        WRITE(4,231) NCENTER
        WRITE(5,231) NUM(I1)
    ENDIF  
    IF (INDEX_S(I1).EQ.4) THEN
        WRITE(2,235) '4(',F11,' ',F21,' ',F24,' ',F14,')'
        WRITE(4,231) NCENTER
        WRITE(5,231) NUM(I1)
    ENDIF
    IF (INDEX_S(I1).EQ.5) THEN
        WRITE(2,235) '4(',F14,' ',F24,' ',F23,' ',F13,')'
        WRITE(4,231) NCENTER
        WRITE(5,231) NUM(I1)
    ENDIF
    IF (INDEX_S(I1).EQ.6) THEN
        WRITE(2,235) '4(',F11,' ',F12,' ',F22,' ',F21,')'
        WRITE(4,231) NCENTER
        WRITE(5,231) NUM(I1)
    ENDIF       
    ENDIF  
    ENDDO
    
    ELSE
    !CASE 2-2-2 
    F11=(CNL+1)*(CNI+1)*(KP-1)+(CNI+1)*(JP-1)+IP-1
    F12=F11+1
    F13=F12+CNI+1
    F14=F13-1
    F21=F11+(CNL+1)*(CNI+1)
    F22=F21+1
    F23=F22+CNI+1
    F24=F23-1 
    CALL CAL_FIJ_RECT(F11,F12,F13,F14,F21,F22,F23,F24)
    NUM(7)=CELL_ID_RECT(IP,JP,KP)
    NUM(1)=CELL_ID_RECT(IP,JP,KP-1) !FRONT FACE
    NUM(2)=CELL_ID_RECT(IP,JP,KP+1) !BACK  FACE
    NUM(3)=CELL_ID_RECT(IP+1,JP,KP) !LEFT  FACE
    NUM(4)=CELL_ID_RECT(IP-1,JP,KP) !RIGHT FACE
    NUM(5)=CELL_ID_RECT(IP,JP+1,KP) !UP    FACE
    NUM(6)=CELL_ID_RECT(IP,JP-1,KP) !DOWN  FACE
    CALL BUBBLE_SORT(NUM,INDEX_S)
    DO I1=1,7
    IF (NUM(I1).GT.NCENTER) THEN   
    IF (INDEX_S(I1).EQ.1) THEN
        WRITE(2,235) '4(',F11,' ',F14,' ',F13,' ',F12,')'
        WRITE(4,231) NCENTER
        WRITE(5,231) NUM(I1)
    ENDIF
    IF (INDEX_S(I1).EQ.2) THEN
        WRITE(2,235) '4(',F21,' ',F22,' ',F23,' ',F24,')'
        WRITE(4,231) NCENTER
        WRITE(5,231) NUM(I1)
    ENDIF
    IF (INDEX_S(I1).EQ.3) THEN
        WRITE(2,235) '4(',F12,' ',F13,' ',F23,' ',F22,')'
        WRITE(4,231) NCENTER
        WRITE(5,231) NUM(I1)
    ENDIF  
    IF (INDEX_S(I1).EQ.4) THEN
        WRITE(2,235) '4(',F11,' ',F21,' ',F24,' ',F14,')'
        WRITE(4,231) NCENTER
        WRITE(5,231) NUM(I1)
    ENDIF
    IF (INDEX_S(I1).EQ.5) THEN
        WRITE(2,235) '4(',F14,' ',F24,' ',F23,' ',F13,')'
        WRITE(4,231) NCENTER
        WRITE(5,231) NUM(I1)
    ENDIF
    IF (INDEX_S(I1).EQ.6) THEN
        WRITE(2,235) '4(',F11,' ',F12,' ',F22,' ',F21,')'
        WRITE(4,231) NCENTER
        WRITE(5,231) NUM(I1)
    ENDIF       
    ENDIF  
    ENDDO
    
    ENDIF    
    ENDIF    
    ENDIF
        
    ENDIF
    ENDIF
    ENDIF
    
    ENDIF
    ENDDO
    
    If ( iMPI_MyID.EQ.0 ) THEN
    write(*,*) 'test 03'
    ENDIF
    !**********************************************************
    
    !�ж����λ��
    !1.�ý��̵��ڲ��棺own�����У���neighbour������
    !2.�ý������������̵Ľ����棺own�����У���neighbour�������У���ʱ����Ϊ+��������own�������У���neighbour�����У���ʱ����Ϊ-��
    !3.�����棺own�������У���neighbour��������
    !�ı䣺CALL GET_INDEX(RENUM2(NUM(I1)),I_NEIGH,J_NEIGH,K_NEIGH,IF_CY)
    !��ɣ�CALL GET_INDEX(NUM(I1),I_NEIGH,J_NEIGH,K_NEIGH,IF_CY)
    FACE_ID=0
    DO NCENTER=0,C_NUM_ALL-1
    IF_OWN=0
    !CALL GET_INDEX(RENUM2(NCENTER),I,J,K,IF_CY)
    CALL GET_INDEX(NCENTER,I,J,K,IF_CY)
    IF (IF_CY.EQ.1) THEN
!*******************CYLINDER***************        
    IP=I-CNI*(MI-1)
    JP=J-CNJ*(MJ-1)
    KP=K
    IF (IP.LE.CNI.AND.IP.GE.1) THEN
    IF (JP.LE.CNJ.AND.JP.GE.1) THEN
        IF_OWN=1
    ENDIF
    ENDIF
    
    IF (I.EQ.1) THEN
    IF (J.EQ.1) THEN
    IF (K.LE.NL) THEN
    !CASE 1-1-1 
    NUM(7)=NCENTER
    CALL GET_CELLID(I,J,K-1,NUM(1),1)
    CALL GET_CELLID(I,J,K+1,NUM(2),1)
    CALL GET_CELLID(I+1,J,K,NUM(3),1)
    NUM(4)=-1
    CALL GET_CELLID(I,J+1,K,NUM(5),1)
    CALL GET_CELLID(I,K,1,NUM(6),0)
    CALL BUBBLE_SORT(NUM,INDEX_S)
    DO I1=1,7
    IF (NUM(I1).GT.NCENTER) THEN  
    IF_NEIGH=0
    CALL GET_INDEX(NUM(I1),I_NEIGH,J_NEIGH,K_NEIGH,IF_CY)
    IF (IF_CY.EQ.1) THEN
    IP_NEIGH=I_NEIGH-CNI*(MI-1)
    JP_NEIGH=J_NEIGH-CNJ*(MJ-1)
    KP_NEIGH=K_NEIGH
    IF (IP_NEIGH.LE.CNI.AND.IP_NEIGH.GE.1) THEN
    IF (JP_NEIGH.LE.CNJ.AND.JP_NEIGH.GE.1) THEN
        IF_NEIGH=1
    ENDIF
    ENDIF 
    ELSE
    IP_NEIGH=I_NEIGH-CNI*(MI-1)
    JP_NEIGH=J_NEIGH-CNL*(MJ1-1)
    KP_NEIGH=K_NEIGH-CNL*(MJ2-1)
    IF (IP_NEIGH.LE.CNI.AND.IP_NEIGH.GE.1) THEN
    IF (JP_NEIGH.LE.CNL.AND.JP_NEIGH.GE.1) THEN
    IF (KP_NEIGH.LE.CNL.AND.KP_NEIGH.GE.1) THEN
        IF_NEIGH=1 
    ENDIF
    ENDIF
    ENDIF    
    ENDIF
    IF (IF_CY.EQ.1) THEN
    IF (INDEX_S(I1).EQ.1) THEN
        FACE_ID=FACE_ID+1
        IF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.1) THEN
        WRITE(3,231) FACE_ID
        ENDIF
    ENDIF
    IF (INDEX_S(I1).EQ.2) THEN
        FACE_ID=FACE_ID+1
        IF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.1) THEN
        WRITE(3,231) FACE_ID
        ENDIF
    ENDIF
    IF (INDEX_S(I1).EQ.3) THEN
        FACE_ID=FACE_ID+1
        IF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.1) THEN
        WRITE(3,231) FACE_ID
        ELSEIF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.0) THEN
        F_LEFT_CY(JP,KP)=FACE_ID
        ELSEIF (IF_OWN.EQ.0.AND.IF_NEIGH.EQ.1) THEN
        F_RIGHT_CY(JP_NEIGH,KP_NEIGH)=-FACE_ID
        ENDIF
    ENDIF  
    IF (INDEX_S(I1).EQ.4) THEN
        FACE_ID=FACE_ID+1
        IF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.1) THEN
        WRITE(3,231) FACE_ID
        ELSEIF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.0) THEN
        F_RIGHT_CY(JP,KP)=FACE_ID
        ELSEIF (IF_OWN.EQ.0.AND.IF_NEIGH.EQ.1) THEN
        F_LEFT_CY(JP_NEIGH,KP_NEIGH)=-FACE_ID
        ENDIF
    ENDIF
    IF (INDEX_S(I1).EQ.5) THEN
        FACE_ID=FACE_ID+1
        IF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.1) THEN
        WRITE(3,231) FACE_ID
        ELSEIF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.0) THEN
        F_UP_CY(IP,KP)=FACE_ID
        ELSEIF (IF_OWN.EQ.0.AND.IF_NEIGH.EQ.1) THEN
        F_DOWN_CY(IP_NEIGH,KP_NEIGH)=-FACE_ID
        ENDIF
    ENDIF
    IF (INDEX_S(I1).EQ.6) THEN
        FACE_ID=FACE_ID+1
        IF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.1) THEN
        WRITE(3,231) FACE_ID
        ELSEIF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.0) THEN
        F_DOWN_CY(IP,KP)=FACE_ID
        ELSEIF (IF_OWN.EQ.0.AND.IF_NEIGH.EQ.1) THEN
        F_UP_CY(IP_NEIGH,KP_NEIGH)=-FACE_ID
        ENDIF
    ENDIF       
     
    ELSE
    IF (INDEX_S(I1).EQ.6) THEN
        FACE_ID=FACE_ID+1
        IF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.0) THEN
        F_DOWN_CY(IP,KP)=FACE_ID
        ELSEIF (IF_OWN.EQ.0.AND.IF_NEIGH.EQ.1) THEN
        F_FRONT_RECT(IP_NEIGH,JP_NEIGH)=-FACE_ID
        ENDIF
    ENDIF
        
    ENDIF
    
    ENDIF 
    ENDDO
    
    
    
    
    ELSEIF (K.GT.NL .AND. K.LE.(2*NL)) THEN
    !CASE 1-1-2
    NUM(7)=NCENTER
    CALL GET_CELLID(I,J,K-1,NUM(1),1)
    CALL GET_CELLID(I,J,K+1,NUM(2),1)
    CALL GET_CELLID(I+1,J,K,NUM(3),1)
    NUM(4)=-1
    CALL GET_CELLID(I,J+1,K,NUM(5),1)
    CALL GET_CELLID(I,NL,K-NL,NUM(6),0)
    CALL BUBBLE_SORT(NUM,INDEX_S)
    DO I1=1,7
    IF (NUM(I1).GT.NCENTER) THEN  
    IF_NEIGH=0
    CALL GET_INDEX(NUM(I1),I_NEIGH,J_NEIGH,K_NEIGH,IF_CY)
    IF (IF_CY.EQ.1) THEN
    IP_NEIGH=I_NEIGH-CNI*(MI-1)
    JP_NEIGH=J_NEIGH-CNJ*(MJ-1)
    KP_NEIGH=K_NEIGH
    IF (IP_NEIGH.LE.CNI.AND.IP_NEIGH.GE.1) THEN
    IF (JP_NEIGH.LE.CNJ.AND.JP_NEIGH.GE.1) THEN
        IF_NEIGH=1
    ENDIF
    ENDIF 
    ELSE
    IP_NEIGH=I_NEIGH-CNI*(MI-1)
    JP_NEIGH=J_NEIGH-CNL*(MJ1-1)
    KP_NEIGH=K_NEIGH-CNL*(MJ2-1)
    IF (IP_NEIGH.LE.CNI.AND.IP_NEIGH.GE.1) THEN
    IF (JP_NEIGH.LE.CNL.AND.JP_NEIGH.GE.1) THEN
    IF (KP_NEIGH.LE.CNL.AND.KP_NEIGH.GE.1) THEN
        IF_NEIGH=1 
    ENDIF
    ENDIF
    ENDIF   
    ENDIF
    
    IF (IF_CY.EQ.1) THEN
    IF (INDEX_S(I1).EQ.1) THEN
        FACE_ID=FACE_ID+1
        IF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.1) THEN
        WRITE(3,231) FACE_ID
        ENDIF
    ENDIF
    IF (INDEX_S(I1).EQ.2) THEN
        FACE_ID=FACE_ID+1
        IF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.1) THEN
        WRITE(3,231) FACE_ID
        ENDIF
    ENDIF
    IF (INDEX_S(I1).EQ.3) THEN
        FACE_ID=FACE_ID+1
        IF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.1) THEN
        WRITE(3,231) FACE_ID
        ELSEIF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.0) THEN
        F_LEFT_CY(JP,KP)=FACE_ID
        ELSEIF (IF_OWN.EQ.0.AND.IF_NEIGH.EQ.1) THEN
        F_RIGHT_CY(JP_NEIGH,KP_NEIGH)=-FACE_ID
        ENDIF
    ENDIF  
    IF (INDEX_S(I1).EQ.4) THEN
        FACE_ID=FACE_ID+1
        IF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.1) THEN
        WRITE(3,231) FACE_ID
        ELSEIF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.0) THEN
        F_RIGHT_CY(JP,KP)=FACE_ID
        ELSEIF (IF_OWN.EQ.0.AND.IF_NEIGH.EQ.1) THEN
        F_LEFT_CY(JP_NEIGH,KP_NEIGH)=-FACE_ID
        ENDIF
    ENDIF
    IF (INDEX_S(I1).EQ.5) THEN
        FACE_ID=FACE_ID+1
        IF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.1) THEN
        WRITE(3,231) FACE_ID
        ELSEIF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.0) THEN
        F_UP_CY(IP,KP)=FACE_ID
        ELSEIF (IF_OWN.EQ.0.AND.IF_NEIGH.EQ.1) THEN
        F_DOWN_CY(IP_NEIGH,KP_NEIGH)=-FACE_ID
        ENDIF
    ENDIF
    IF (INDEX_S(I1).EQ.6) THEN
        FACE_ID=FACE_ID+1
        IF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.1) THEN
        WRITE(3,231) FACE_ID
        ELSEIF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.0) THEN
        F_DOWN_CY(IP,KP)=FACE_ID
        ELSEIF (IF_OWN.EQ.0.AND.IF_NEIGH.EQ.1) THEN
        F_UP_CY(IP_NEIGH,KP_NEIGH)=-FACE_ID
        ENDIF
    ENDIF       
     
    ELSE
    IF (INDEX_S(I1).EQ.6) THEN
        FACE_ID=FACE_ID+1
        IF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.0) THEN
        F_DOWN_CY(IP,KP)=FACE_ID
        ELSEIF (IF_OWN.EQ.0.AND.IF_NEIGH.EQ.1) THEN
        F_UP_RECT(IP_NEIGH,KP_NEIGH)=-FACE_ID
        ENDIF
    ENDIF
        
    ENDIF
    
    ENDIF 
    ENDDO
    ELSEIF (K.GT.(2*NL) .AND. K.LE.(3*NL)) THEN
    !CASE 1-1-3
    NUM(7)=NCENTER
    CALL GET_CELLID(I,J,K-1,NUM(1),1)
    CALL GET_CELLID(I,J,K+1,NUM(2),1)
    CALL GET_CELLID(I+1,J,K,NUM(3),1)
    NUM(4)=-1
    CALL GET_CELLID(I,J+1,K,NUM(5),1)
    CALL GET_CELLID(I,3*NL+1-K,NL,NUM(6),0)
    CALL BUBBLE_SORT(NUM,INDEX_S)
    DO I1=1,7
    IF (NUM(I1).GT.NCENTER) THEN  
    IF_NEIGH=0
    CALL GET_INDEX(NUM(I1),I_NEIGH,J_NEIGH,K_NEIGH,IF_CY)
    IF (IF_CY.EQ.1) THEN
    IP_NEIGH=I_NEIGH-CNI*(MI-1)
    JP_NEIGH=J_NEIGH-CNJ*(MJ-1)
    KP_NEIGH=K_NEIGH
    IF (IP_NEIGH.LE.CNI.AND.IP_NEIGH.GE.1) THEN
    IF (JP_NEIGH.LE.CNJ.AND.JP_NEIGH.GE.1) THEN
        IF_NEIGH=1
    ENDIF
    ENDIF 
    ELSE
    IP_NEIGH=I_NEIGH-CNI*(MI-1)
    JP_NEIGH=J_NEIGH-CNL*(MJ1-1)
    KP_NEIGH=K_NEIGH-CNL*(MJ2-1)
    IF (IP_NEIGH.LE.CNI.AND.IP_NEIGH.GE.1) THEN
    IF (JP_NEIGH.LE.CNL.AND.JP_NEIGH.GE.1) THEN
    IF (KP_NEIGH.LE.CNL.AND.KP_NEIGH.GE.1) THEN
        IF_NEIGH=1 
    ENDIF
    ENDIF
    ENDIF   
    ENDIF
    
    IF (IF_CY.EQ.1) THEN
    IF (INDEX_S(I1).EQ.1) THEN
        FACE_ID=FACE_ID+1
        IF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.1) THEN
        WRITE(3,231) FACE_ID
        ENDIF
    ENDIF
    IF (INDEX_S(I1).EQ.2) THEN
        FACE_ID=FACE_ID+1
        IF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.1) THEN
        WRITE(3,231) FACE_ID
        ENDIF
    ENDIF
    IF (INDEX_S(I1).EQ.3) THEN
        FACE_ID=FACE_ID+1
        IF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.1) THEN
        WRITE(3,231) FACE_ID
        ELSEIF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.0) THEN
        F_LEFT_CY(JP,KP)=FACE_ID
        ELSEIF (IF_OWN.EQ.0.AND.IF_NEIGH.EQ.1) THEN
        F_RIGHT_CY(JP_NEIGH,KP_NEIGH)=-FACE_ID
        ENDIF
    ENDIF  
    IF (INDEX_S(I1).EQ.4) THEN
        FACE_ID=FACE_ID+1
        IF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.1) THEN
        WRITE(3,231) FACE_ID
        ELSEIF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.0) THEN
        F_RIGHT_CY(JP,KP)=FACE_ID
        ELSEIF (IF_OWN.EQ.0.AND.IF_NEIGH.EQ.1) THEN
        F_LEFT_CY(JP_NEIGH,KP_NEIGH)=-FACE_ID
        ENDIF
    ENDIF
    IF (INDEX_S(I1).EQ.5) THEN
        FACE_ID=FACE_ID+1
        IF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.1) THEN
        WRITE(3,231) FACE_ID
        ELSEIF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.0) THEN
        F_UP_CY(IP,KP)=FACE_ID
        ELSEIF (IF_OWN.EQ.0.AND.IF_NEIGH.EQ.1) THEN
        F_DOWN_CY(IP_NEIGH,KP_NEIGH)=-FACE_ID
        ENDIF
    ENDIF
    IF (INDEX_S(I1).EQ.6) THEN
        FACE_ID=FACE_ID+1
        IF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.1) THEN
        WRITE(3,231) FACE_ID
        ELSEIF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.0) THEN
        F_DOWN_CY(IP,KP)=FACE_ID
        ELSEIF (IF_OWN.EQ.0.AND.IF_NEIGH.EQ.1) THEN
        F_UP_CY(IP_NEIGH,KP_NEIGH)=-FACE_ID
        ENDIF
    ENDIF       
     
    ELSE
    IF (INDEX_S(I1).EQ.6) THEN
        FACE_ID=FACE_ID+1
        IF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.0) THEN
        F_DOWN_CY(IP,KP)=FACE_ID
        ELSEIF (IF_OWN.EQ.0.AND.IF_NEIGH.EQ.1) THEN
        F_BACK_RECT(IP_NEIGH,JP_NEIGH)=-FACE_ID
        ENDIF
    ENDIF
        
    ENDIF
    
    ENDIF 
    ENDDO
    ELSE
    !CASE 1-1-4
    NUM(7)=NCENTER
    CALL GET_CELLID(I,J,K-1,NUM(1),1)
    CALL GET_CELLID(I,J,K+1,NUM(2),1)
    CALL GET_CELLID(I+1,J,K,NUM(3),1)
    NUM(4)=-1
    CALL GET_CELLID(I,J+1,K,NUM(5),1)
    CALL GET_CELLID(I,1,4*NL+1-K,NUM(6),0)
    CALL BUBBLE_SORT(NUM,INDEX_S)
    DO I1=1,7
    IF (NUM(I1).GT.NCENTER) THEN  
    IF_NEIGH=0
    CALL GET_INDEX(NUM(I1),I_NEIGH,J_NEIGH,K_NEIGH,IF_CY)
    IF (IF_CY.EQ.1) THEN
    IP_NEIGH=I_NEIGH-CNI*(MI-1)
    JP_NEIGH=J_NEIGH-CNJ*(MJ-1)
    KP_NEIGH=K_NEIGH
    IF (IP_NEIGH.LE.CNI.AND.IP_NEIGH.GE.1) THEN
    IF (JP_NEIGH.LE.CNJ.AND.JP_NEIGH.GE.1) THEN
        IF_NEIGH=1
    ENDIF
    ENDIF 
    ELSE
    IP_NEIGH=I_NEIGH-CNI*(MI-1)
    JP_NEIGH=J_NEIGH-CNL*(MJ1-1)
    KP_NEIGH=K_NEIGH-CNL*(MJ2-1)
    IF (IP_NEIGH.LE.CNI.AND.IP_NEIGH.GE.1) THEN
    IF (JP_NEIGH.LE.CNL.AND.JP_NEIGH.GE.1) THEN
    IF (KP_NEIGH.LE.CNL.AND.KP_NEIGH.GE.1) THEN
        IF_NEIGH=1 
    ENDIF
    ENDIF
    ENDIF    
    ENDIF
    
    IF (IF_CY.EQ.1) THEN
    IF (INDEX_S(I1).EQ.1) THEN
        FACE_ID=FACE_ID+1
        IF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.1) THEN
        WRITE(3,231) FACE_ID
        ENDIF
    ENDIF
    IF (INDEX_S(I1).EQ.2) THEN
        FACE_ID=FACE_ID+1
        IF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.1) THEN
        WRITE(3,231) FACE_ID
        ENDIF
    ENDIF
    IF (INDEX_S(I1).EQ.3) THEN
        FACE_ID=FACE_ID+1
        IF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.1) THEN
        WRITE(3,231) FACE_ID
        ELSEIF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.0) THEN
        F_LEFT_CY(JP,KP)=FACE_ID
        ELSEIF (IF_OWN.EQ.0.AND.IF_NEIGH.EQ.1) THEN
        F_RIGHT_CY(JP_NEIGH,KP_NEIGH)=-FACE_ID
        ENDIF
    ENDIF  
    IF (INDEX_S(I1).EQ.4) THEN
        FACE_ID=FACE_ID+1
        IF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.1) THEN
        WRITE(3,231) FACE_ID
        ELSEIF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.0) THEN
        F_RIGHT_CY(JP,KP)=FACE_ID
        ELSEIF (IF_OWN.EQ.0.AND.IF_NEIGH.EQ.1) THEN
        F_LEFT_CY(JP_NEIGH,KP_NEIGH)=-FACE_ID
        ENDIF
    ENDIF
    IF (INDEX_S(I1).EQ.5) THEN
        FACE_ID=FACE_ID+1
        IF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.1) THEN
        WRITE(3,231) FACE_ID
        ELSEIF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.0) THEN
        F_UP_CY(IP,KP)=FACE_ID
        ELSEIF (IF_OWN.EQ.0.AND.IF_NEIGH.EQ.1) THEN
        F_DOWN_CY(IP_NEIGH,KP_NEIGH)=-FACE_ID
        ENDIF
    ENDIF
    IF (INDEX_S(I1).EQ.6) THEN
        FACE_ID=FACE_ID+1
        IF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.1) THEN
        WRITE(3,231) FACE_ID
        ELSEIF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.0) THEN
        F_DOWN_CY(IP,KP)=FACE_ID
        ELSEIF (IF_OWN.EQ.0.AND.IF_NEIGH.EQ.1) THEN
        F_UP_CY(IP_NEIGH,KP_NEIGH)=-FACE_ID
        ENDIF
    ENDIF       
     
    ELSE
    IF (INDEX_S(I1).EQ.6) THEN
        FACE_ID=FACE_ID+1
        IF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.0) THEN
        F_DOWN_CY(IP,KP)=FACE_ID
        ELSEIF (IF_OWN.EQ.0.AND.IF_NEIGH.EQ.1) THEN
        F_DOWN_RECT(IP_NEIGH,KP_NEIGH)=-FACE_ID
        ENDIF
    ENDIF
        
    ENDIF
    
    ENDIF 
    ENDDO    
    ENDIF     
    
    ELSEIF (J.EQ.NJ) THEN
    !CASE 1-3  
    NUM(7)=NCENTER
    CALL GET_CELLID(I,J,K-1,NUM(1),1)
    CALL GET_CELLID(I,J,K+1,NUM(2),1)
    CALL GET_CELLID(I+1,J,K,NUM(3),1)
    NUM(4)=-1
    NUM(5)=-1
    CALL GET_CELLID(I,J-1,K,NUM(6),1)
    CALL BUBBLE_SORT(NUM,INDEX_S)
    DO I1=1,7
    IF (NUM(I1).GT.NCENTER) THEN  
    IF_NEIGH=0
    CALL GET_INDEX(NUM(I1),I_NEIGH,J_NEIGH,K_NEIGH,IF_CY)
    IP_NEIGH=I_NEIGH-CNI*(MI-1)
    JP_NEIGH=J_NEIGH-CNJ*(MJ-1)
    KP_NEIGH=K_NEIGH
    IF (IP_NEIGH.LE.CNI.AND.IP_NEIGH.GE.1) THEN
    IF (JP_NEIGH.LE.CNJ.AND.JP_NEIGH.GE.1) THEN
        IF_NEIGH=1
    ENDIF
    ENDIF    
    
    IF (INDEX_S(I1).EQ.1) THEN
        FACE_ID=FACE_ID+1
        IF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.1) THEN
        WRITE(3,231) FACE_ID
        ENDIF
    ENDIF
    IF (INDEX_S(I1).EQ.2) THEN
        FACE_ID=FACE_ID+1
        IF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.1) THEN
        WRITE(3,231) FACE_ID
        ENDIF
    ENDIF
    IF (INDEX_S(I1).EQ.3) THEN
        FACE_ID=FACE_ID+1
        IF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.1) THEN
        WRITE(3,231) FACE_ID
        ELSEIF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.0) THEN
        F_LEFT_CY(JP,KP)=FACE_ID
        ELSEIF (IF_OWN.EQ.0.AND.IF_NEIGH.EQ.1) THEN
        F_RIGHT_CY(JP_NEIGH,KP_NEIGH)=-FACE_ID
        ENDIF
    ENDIF  
    IF (INDEX_S(I1).EQ.4) THEN
        FACE_ID=FACE_ID+1
        IF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.1) THEN
        WRITE(3,231) FACE_ID
        ELSEIF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.0) THEN
        F_RIGHT_CY(JP,KP)=FACE_ID
        ELSEIF (IF_OWN.EQ.0.AND.IF_NEIGH.EQ.1) THEN
        F_LEFT_CY(JP_NEIGH,KP_NEIGH)=-FACE_ID
        ENDIF
    ENDIF
    IF (INDEX_S(I1).EQ.5) THEN
        FACE_ID=FACE_ID+1
        IF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.1) THEN
        WRITE(3,231) FACE_ID
        ELSEIF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.0) THEN
        F_UP_CY(IP,KP)=FACE_ID
        ELSEIF (IF_OWN.EQ.0.AND.IF_NEIGH.EQ.1) THEN
        F_DOWN_CY(IP_NEIGH,KP_NEIGH)=-FACE_ID
        ENDIF
    ENDIF
    IF (INDEX_S(I1).EQ.6) THEN
        FACE_ID=FACE_ID+1
        IF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.1) THEN
        WRITE(3,231) FACE_ID
        ELSEIF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.0) THEN
        F_DOWN_CY(IP,KP)=FACE_ID
        ELSEIF (IF_OWN.EQ.0.AND.IF_NEIGH.EQ.1) THEN
        F_UP_CY(IP_NEIGH,KP_NEIGH)=-FACE_ID
        ENDIF
    ENDIF       
    
    ENDIF 
    ENDDO
    
    ELSE
    !CASE 1-2 
    NUM(7)=NCENTER
    CALL GET_CELLID(I,J,K-1,NUM(1),1)
    CALL GET_CELLID(I,J,K+1,NUM(2),1)
    CALL GET_CELLID(I+1,J,K,NUM(3),1)
    NUM(4)=-1
    CALL GET_CELLID(I,J+1,K,NUM(5),1)
    CALL GET_CELLID(I,J-1,K,NUM(6),1)
    CALL BUBBLE_SORT(NUM,INDEX_S)
    DO I1=1,7
    IF (NUM(I1).GT.NCENTER) THEN  
    IF_NEIGH=0
    CALL GET_INDEX(NUM(I1),I_NEIGH,J_NEIGH,K_NEIGH,IF_CY)
    IP_NEIGH=I_NEIGH-CNI*(MI-1)
    JP_NEIGH=J_NEIGH-CNJ*(MJ-1)
    KP_NEIGH=K_NEIGH
    IF (IP_NEIGH.LE.CNI.AND.IP_NEIGH.GE.1) THEN
    IF (JP_NEIGH.LE.CNJ.AND.JP_NEIGH.GE.1) THEN
        IF_NEIGH=1
    ENDIF
    ENDIF    
    
    IF (INDEX_S(I1).EQ.1) THEN
        FACE_ID=FACE_ID+1
        IF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.1) THEN
        WRITE(3,231) FACE_ID
        ENDIF
    ENDIF
    IF (INDEX_S(I1).EQ.2) THEN
        FACE_ID=FACE_ID+1
        IF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.1) THEN
        WRITE(3,231) FACE_ID
        ENDIF
    ENDIF
    IF (INDEX_S(I1).EQ.3) THEN
        FACE_ID=FACE_ID+1
        IF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.1) THEN
        WRITE(3,231) FACE_ID
        ELSEIF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.0) THEN
        F_LEFT_CY(JP,KP)=FACE_ID
        ELSEIF (IF_OWN.EQ.0.AND.IF_NEIGH.EQ.1) THEN
        F_RIGHT_CY(JP_NEIGH,KP_NEIGH)=-FACE_ID
        ENDIF
    ENDIF  
    IF (INDEX_S(I1).EQ.4) THEN
        FACE_ID=FACE_ID+1
        IF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.1) THEN
        WRITE(3,231) FACE_ID
        ELSEIF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.0) THEN
        F_RIGHT_CY(JP,KP)=FACE_ID
        ELSEIF (IF_OWN.EQ.0.AND.IF_NEIGH.EQ.1) THEN
        F_LEFT_CY(JP_NEIGH,KP_NEIGH)=-FACE_ID
        ENDIF
    ENDIF
    IF (INDEX_S(I1).EQ.5) THEN
        FACE_ID=FACE_ID+1
        IF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.1) THEN
        WRITE(3,231) FACE_ID
        ELSEIF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.0) THEN
        F_UP_CY(IP,KP)=FACE_ID
        ELSEIF (IF_OWN.EQ.0.AND.IF_NEIGH.EQ.1) THEN
        F_DOWN_CY(IP_NEIGH,KP_NEIGH)=-FACE_ID
        ENDIF
    ENDIF
    IF (INDEX_S(I1).EQ.6) THEN
        FACE_ID=FACE_ID+1
        IF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.1) THEN
        WRITE(3,231) FACE_ID
        ELSEIF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.0) THEN
        F_DOWN_CY(IP,KP)=FACE_ID
        ELSEIF (IF_OWN.EQ.0.AND.IF_NEIGH.EQ.1) THEN
        F_UP_CY(IP_NEIGH,KP_NEIGH)=-FACE_ID
        ENDIF
    ENDIF       
    
    ENDIF 
    ENDDO
    
    ENDIF    
    ELSEIF (I.EQ.NI) THEN
    IF (J.EQ.1) THEN
    IF (K.LE.NL) THEN
    !CASE 3-1-1 
    NUM(7)=NCENTER
    CALL GET_CELLID(I,J,K-1,NUM(1),1)
    CALL GET_CELLID(I,J,K+1,NUM(2),1)
    NUM(3)=-1
    CALL GET_CELLID(I-1,J,K,NUM(4),1)
    CALL GET_CELLID(I,J+1,K,NUM(5),1)
    CALL GET_CELLID(I,K,1,NUM(6),0)
    CALL BUBBLE_SORT(NUM,INDEX_S)
    DO I1=1,7
    IF (NUM(I1).GT.NCENTER) THEN  
    IF_NEIGH=0
    CALL GET_INDEX(NUM(I1),I_NEIGH,J_NEIGH,K_NEIGH,IF_CY)
    IF (IF_CY.EQ.1) THEN
    IP_NEIGH=I_NEIGH-CNI*(MI-1)
    JP_NEIGH=J_NEIGH-CNJ*(MJ-1)
    KP_NEIGH=K_NEIGH
    IF (IP_NEIGH.LE.CNI.AND.IP_NEIGH.GE.1) THEN
    IF (JP_NEIGH.LE.CNJ.AND.JP_NEIGH.GE.1) THEN
        IF_NEIGH=1
    ENDIF
    ENDIF 
    ELSE
    IP_NEIGH=I_NEIGH-CNI*(MI-1)
    JP_NEIGH=J_NEIGH-CNL*(MJ1-1)
    KP_NEIGH=K_NEIGH-CNL*(MJ2-1)
    IF (IP_NEIGH.LE.CNI.AND.IP_NEIGH.GE.1) THEN
    IF (JP_NEIGH.LE.CNL.AND.JP_NEIGH.GE.1) THEN
    IF (KP_NEIGH.LE.CNL.AND.KP_NEIGH.GE.1) THEN
        IF_NEIGH=1 
    ENDIF
    ENDIF
    ENDIF    
    ENDIF
    IF (IF_CY.EQ.1) THEN
    IF (INDEX_S(I1).EQ.1) THEN
        FACE_ID=FACE_ID+1
        IF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.1) THEN
        WRITE(3,231) FACE_ID
        ENDIF
    ENDIF
    IF (INDEX_S(I1).EQ.2) THEN
        FACE_ID=FACE_ID+1
        IF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.1) THEN
        WRITE(3,231) FACE_ID
        ENDIF
    ENDIF
    IF (INDEX_S(I1).EQ.3) THEN
        FACE_ID=FACE_ID+1
        IF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.1) THEN
        WRITE(3,231) FACE_ID
        ELSEIF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.0) THEN
        F_LEFT_CY(JP,KP)=FACE_ID
        ELSEIF (IF_OWN.EQ.0.AND.IF_NEIGH.EQ.1) THEN
        F_RIGHT_CY(JP_NEIGH,KP_NEIGH)=-FACE_ID
        ENDIF
    ENDIF  
    IF (INDEX_S(I1).EQ.4) THEN
        FACE_ID=FACE_ID+1
        IF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.1) THEN
        WRITE(3,231) FACE_ID
        ELSEIF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.0) THEN
        F_RIGHT_CY(JP,KP)=FACE_ID
        ELSEIF (IF_OWN.EQ.0.AND.IF_NEIGH.EQ.1) THEN
        F_LEFT_CY(JP_NEIGH,KP_NEIGH)=-FACE_ID
        ENDIF
    ENDIF
    IF (INDEX_S(I1).EQ.5) THEN
        FACE_ID=FACE_ID+1
        IF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.1) THEN
        WRITE(3,231) FACE_ID
        ELSEIF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.0) THEN
        F_UP_CY(IP,KP)=FACE_ID
        ELSEIF (IF_OWN.EQ.0.AND.IF_NEIGH.EQ.1) THEN
        F_DOWN_CY(IP_NEIGH,KP_NEIGH)=-FACE_ID
        ENDIF
    ENDIF
    IF (INDEX_S(I1).EQ.6) THEN
        FACE_ID=FACE_ID+1
        IF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.1) THEN
        WRITE(3,231) FACE_ID
        ELSEIF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.0) THEN
        F_DOWN_CY(IP,KP)=FACE_ID
        ELSEIF (IF_OWN.EQ.0.AND.IF_NEIGH.EQ.1) THEN
        F_UP_CY(IP_NEIGH,KP_NEIGH)=-FACE_ID
        ENDIF
    ENDIF       
     
    ELSE
    IF (INDEX_S(I1).EQ.6) THEN
        FACE_ID=FACE_ID+1
        IF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.0) THEN
        F_DOWN_CY(IP,KP)=FACE_ID
        ELSEIF (IF_OWN.EQ.0.AND.IF_NEIGH.EQ.1) THEN
        F_FRONT_RECT(IP_NEIGH,JP_NEIGH)=-FACE_ID
        ENDIF
    ENDIF
        
    ENDIF
    
    ENDIF 
    ENDDO
    
    
    
    
    ELSEIF (K.GT.NL .AND. K.LE.(2*NL)) THEN
    !CASE 3-1-2
    NUM(7)=NCENTER
    CALL GET_CELLID(I,J,K-1,NUM(1),1)
    CALL GET_CELLID(I,J,K+1,NUM(2),1)
    NUM(3)=-1
    CALL GET_CELLID(I-1,J,K,NUM(4),1)
    CALL GET_CELLID(I,J+1,K,NUM(5),1)
    CALL GET_CELLID(I,NL,K-NL,NUM(6),0)
    CALL BUBBLE_SORT(NUM,INDEX_S)
    DO I1=1,7
    IF (NUM(I1).GT.NCENTER) THEN  
    IF_NEIGH=0
    CALL GET_INDEX(NUM(I1),I_NEIGH,J_NEIGH,K_NEIGH,IF_CY)
    IF (IF_CY.EQ.1) THEN
    IP_NEIGH=I_NEIGH-CNI*(MI-1)
    JP_NEIGH=J_NEIGH-CNJ*(MJ-1)
    KP_NEIGH=K_NEIGH
    IF (IP_NEIGH.LE.CNI.AND.IP_NEIGH.GE.1) THEN
    IF (JP_NEIGH.LE.CNJ.AND.JP_NEIGH.GE.1) THEN
        IF_NEIGH=1
    ENDIF
    ENDIF 
    ELSE
    IP_NEIGH=I_NEIGH-CNI*(MI-1)
    JP_NEIGH=J_NEIGH-CNL*(MJ1-1)
    KP_NEIGH=K_NEIGH-CNL*(MJ2-1)
    IF (IP_NEIGH.LE.CNI.AND.IP_NEIGH.GE.1) THEN
    IF (JP_NEIGH.LE.CNL.AND.JP_NEIGH.GE.1) THEN
    IF (KP_NEIGH.LE.CNL.AND.KP_NEIGH.GE.1) THEN
        IF_NEIGH=1 
    ENDIF
    ENDIF
    ENDIF   
    ENDIF
    
    IF (IF_CY.EQ.1) THEN
    IF (INDEX_S(I1).EQ.1) THEN
        FACE_ID=FACE_ID+1
        IF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.1) THEN
        WRITE(3,231) FACE_ID
        ENDIF
    ENDIF
    IF (INDEX_S(I1).EQ.2) THEN
        FACE_ID=FACE_ID+1
        IF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.1) THEN
        WRITE(3,231) FACE_ID
        ENDIF
    ENDIF
    IF (INDEX_S(I1).EQ.3) THEN
        FACE_ID=FACE_ID+1
        IF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.1) THEN
        WRITE(3,231) FACE_ID
        ELSEIF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.0) THEN
        F_LEFT_CY(JP,KP)=FACE_ID
        ELSEIF (IF_OWN.EQ.0.AND.IF_NEIGH.EQ.1) THEN
        F_RIGHT_CY(JP_NEIGH,KP_NEIGH)=-FACE_ID
        ENDIF
    ENDIF  
    IF (INDEX_S(I1).EQ.4) THEN
        FACE_ID=FACE_ID+1
        IF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.1) THEN
        WRITE(3,231) FACE_ID
        ELSEIF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.0) THEN
        F_RIGHT_CY(JP,KP)=FACE_ID
        ELSEIF (IF_OWN.EQ.0.AND.IF_NEIGH.EQ.1) THEN
        F_LEFT_CY(JP_NEIGH,KP_NEIGH)=-FACE_ID
        ENDIF
    ENDIF
    IF (INDEX_S(I1).EQ.5) THEN
        FACE_ID=FACE_ID+1
        IF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.1) THEN
        WRITE(3,231) FACE_ID
        ELSEIF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.0) THEN
        F_UP_CY(IP,KP)=FACE_ID
        ELSEIF (IF_OWN.EQ.0.AND.IF_NEIGH.EQ.1) THEN
        F_DOWN_CY(IP_NEIGH,KP_NEIGH)=-FACE_ID
        ENDIF
    ENDIF
    IF (INDEX_S(I1).EQ.6) THEN
        FACE_ID=FACE_ID+1
        IF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.1) THEN
        WRITE(3,231) FACE_ID
        ELSEIF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.0) THEN
        F_DOWN_CY(IP,KP)=FACE_ID
        ELSEIF (IF_OWN.EQ.0.AND.IF_NEIGH.EQ.1) THEN
        F_UP_CY(IP_NEIGH,KP_NEIGH)=-FACE_ID
        ENDIF
    ENDIF       
     
    ELSE
    IF (INDEX_S(I1).EQ.6) THEN
        FACE_ID=FACE_ID+1
        IF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.0) THEN
        F_DOWN_CY(IP,KP)=FACE_ID
        ELSEIF (IF_OWN.EQ.0.AND.IF_NEIGH.EQ.1) THEN
        F_UP_RECT(IP_NEIGH,KP_NEIGH)=-FACE_ID
        ENDIF
    ENDIF
        
    ENDIF
    
    ENDIF 
    ENDDO
    ELSEIF (K.GT.(2*NL) .AND. K.LE.(3*NL)) THEN
    !CASE 3-1-3
    NUM(7)=NCENTER
    CALL GET_CELLID(I,J,K-1,NUM(1),1)
    CALL GET_CELLID(I,J,K+1,NUM(2),1)
    NUM(3)=-1
    CALL GET_CELLID(I-1,J,K,NUM(4),1)
    CALL GET_CELLID(I,J+1,K,NUM(5),1)
    CALL GET_CELLID(I,3*NL+1-K,NL,NUM(6),0)
    CALL BUBBLE_SORT(NUM,INDEX_S)
    DO I1=1,7
    IF (NUM(I1).GT.NCENTER) THEN  
    IF_NEIGH=0
    CALL GET_INDEX(NUM(I1),I_NEIGH,J_NEIGH,K_NEIGH,IF_CY)
    IF (IF_CY.EQ.1) THEN
    IP_NEIGH=I_NEIGH-CNI*(MI-1)
    JP_NEIGH=J_NEIGH-CNJ*(MJ-1)
    KP_NEIGH=K_NEIGH
    IF (IP_NEIGH.LE.CNI.AND.IP_NEIGH.GE.1) THEN
    IF (JP_NEIGH.LE.CNJ.AND.JP_NEIGH.GE.1) THEN
        IF_NEIGH=1
    ENDIF
    ENDIF 
    ELSE
    IP_NEIGH=I_NEIGH-CNI*(MI-1)
    JP_NEIGH=J_NEIGH-CNL*(MJ1-1)
    KP_NEIGH=K_NEIGH-CNL*(MJ2-1)
    IF (IP_NEIGH.LE.CNI.AND.IP_NEIGH.GE.1) THEN
    IF (JP_NEIGH.LE.CNL.AND.JP_NEIGH.GE.1) THEN
    IF (KP_NEIGH.LE.CNL.AND.KP_NEIGH.GE.1) THEN
        IF_NEIGH=1 
    ENDIF
    ENDIF
    ENDIF    
    ENDIF
    
    IF (IF_CY.EQ.1) THEN
    IF (INDEX_S(I1).EQ.1) THEN
        FACE_ID=FACE_ID+1
        IF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.1) THEN
        WRITE(3,231) FACE_ID
        ENDIF
    ENDIF
    IF (INDEX_S(I1).EQ.2) THEN
        FACE_ID=FACE_ID+1
        IF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.1) THEN
        WRITE(3,231) FACE_ID
        ENDIF
    ENDIF
    IF (INDEX_S(I1).EQ.3) THEN
        FACE_ID=FACE_ID+1
        IF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.1) THEN
        WRITE(3,231) FACE_ID
        ELSEIF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.0) THEN
        F_LEFT_CY(JP,KP)=FACE_ID
        ELSEIF (IF_OWN.EQ.0.AND.IF_NEIGH.EQ.1) THEN
        F_RIGHT_CY(JP_NEIGH,KP_NEIGH)=-FACE_ID
        ENDIF
    ENDIF  
    IF (INDEX_S(I1).EQ.4) THEN
        FACE_ID=FACE_ID+1
        IF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.1) THEN
        WRITE(3,231) FACE_ID
        ELSEIF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.0) THEN
        F_RIGHT_CY(JP,KP)=FACE_ID
        ELSEIF (IF_OWN.EQ.0.AND.IF_NEIGH.EQ.1) THEN
        F_LEFT_CY(JP_NEIGH,KP_NEIGH)=-FACE_ID
        ENDIF
    ENDIF
    IF (INDEX_S(I1).EQ.5) THEN
        FACE_ID=FACE_ID+1
        IF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.1) THEN
        WRITE(3,231) FACE_ID
        ELSEIF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.0) THEN
        F_UP_CY(IP,KP)=FACE_ID
        ELSEIF (IF_OWN.EQ.0.AND.IF_NEIGH.EQ.1) THEN
        F_DOWN_CY(IP_NEIGH,KP_NEIGH)=-FACE_ID
        ENDIF
    ENDIF
    IF (INDEX_S(I1).EQ.6) THEN
        FACE_ID=FACE_ID+1
        IF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.1) THEN
        WRITE(3,231) FACE_ID
        ELSEIF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.0) THEN
        F_DOWN_CY(IP,KP)=FACE_ID
        ELSEIF (IF_OWN.EQ.0.AND.IF_NEIGH.EQ.1) THEN
        F_UP_CY(IP_NEIGH,KP_NEIGH)=-FACE_ID
        ENDIF
    ENDIF       
     
    ELSE
    IF (INDEX_S(I1).EQ.6) THEN
        FACE_ID=FACE_ID+1
        IF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.0) THEN
        F_DOWN_CY(IP,KP)=FACE_ID
        ELSEIF (IF_OWN.EQ.0.AND.IF_NEIGH.EQ.1) THEN
        F_BACK_RECT(IP_NEIGH,JP_NEIGH)=-FACE_ID
        ENDIF
    ENDIF
        
    ENDIF
    
    ENDIF 
    ENDDO
    ELSE
    !CASE 3-1-4
    NUM(7)=NCENTER
    CALL GET_CELLID(I,J,K-1,NUM(1),1)
    CALL GET_CELLID(I,J,K+1,NUM(2),1)
    NUM(3)=-1
    CALL GET_CELLID(I-1,J,K,NUM(4),1)
    CALL GET_CELLID(I,J+1,K,NUM(5),1)
    CALL GET_CELLID(I,1,4*NL+1-K,NUM(6),0)
    CALL BUBBLE_SORT(NUM,INDEX_S)
    DO I1=1,7
    IF (NUM(I1).GT.NCENTER) THEN  
    IF_NEIGH=0
    CALL GET_INDEX(NUM(I1),I_NEIGH,J_NEIGH,K_NEIGH,IF_CY)
    IF (IF_CY.EQ.1) THEN
    IP_NEIGH=I_NEIGH-CNI*(MI-1)
    JP_NEIGH=J_NEIGH-CNJ*(MJ-1)
    KP_NEIGH=K_NEIGH
    IF (IP_NEIGH.LE.CNI.AND.IP_NEIGH.GE.1) THEN
    IF (JP_NEIGH.LE.CNJ.AND.JP_NEIGH.GE.1) THEN
        IF_NEIGH=1
    ENDIF
    ENDIF 
    ELSE
    IP_NEIGH=I_NEIGH-CNI*(MI-1)
    JP_NEIGH=J_NEIGH-CNL*(MJ1-1)
    KP_NEIGH=K_NEIGH-CNL*(MJ2-1)
    IF (IP_NEIGH.LE.CNI.AND.IP_NEIGH.GE.1) THEN
    IF (JP_NEIGH.LE.CNL.AND.JP_NEIGH.GE.1) THEN
    IF (KP_NEIGH.LE.CNL.AND.KP_NEIGH.GE.1) THEN
        IF_NEIGH=1 
    ENDIF
    ENDIF
    ENDIF   
    ENDIF
    
    IF (IF_CY.EQ.1) THEN
    IF (INDEX_S(I1).EQ.1) THEN
        FACE_ID=FACE_ID+1
        IF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.1) THEN
        WRITE(3,231) FACE_ID
        ENDIF
    ENDIF
    IF (INDEX_S(I1).EQ.2) THEN
        FACE_ID=FACE_ID+1
        IF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.1) THEN
        WRITE(3,231) FACE_ID
        ENDIF
    ENDIF
    IF (INDEX_S(I1).EQ.3) THEN
        FACE_ID=FACE_ID+1
        IF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.1) THEN
        WRITE(3,231) FACE_ID
        ELSEIF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.0) THEN
        F_LEFT_CY(JP,KP)=FACE_ID
        ELSEIF (IF_OWN.EQ.0.AND.IF_NEIGH.EQ.1) THEN
        F_RIGHT_CY(JP_NEIGH,KP_NEIGH)=-FACE_ID
        ENDIF
    ENDIF  
    IF (INDEX_S(I1).EQ.4) THEN
        FACE_ID=FACE_ID+1
        IF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.1) THEN
        WRITE(3,231) FACE_ID
        ELSEIF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.0) THEN
        F_RIGHT_CY(JP,KP)=FACE_ID
        ELSEIF (IF_OWN.EQ.0.AND.IF_NEIGH.EQ.1) THEN
        F_LEFT_CY(JP_NEIGH,KP_NEIGH)=-FACE_ID
        ENDIF
    ENDIF
    IF (INDEX_S(I1).EQ.5) THEN
        FACE_ID=FACE_ID+1
        IF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.1) THEN
        WRITE(3,231) FACE_ID
        ELSEIF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.0) THEN
        F_UP_CY(IP,KP)=FACE_ID
        ELSEIF (IF_OWN.EQ.0.AND.IF_NEIGH.EQ.1) THEN
        F_DOWN_CY(IP_NEIGH,KP_NEIGH)=-FACE_ID
        ENDIF
    ENDIF
    IF (INDEX_S(I1).EQ.6) THEN
        FACE_ID=FACE_ID+1
        IF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.1) THEN
        WRITE(3,231) FACE_ID
        ELSEIF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.0) THEN
        F_DOWN_CY(IP,KP)=FACE_ID
        ELSEIF (IF_OWN.EQ.0.AND.IF_NEIGH.EQ.1) THEN
        F_UP_CY(IP_NEIGH,KP_NEIGH)=-FACE_ID
        ENDIF
    ENDIF       
     
    ELSE
    IF (INDEX_S(I1).EQ.6) THEN
        FACE_ID=FACE_ID+1
        IF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.0) THEN
        F_DOWN_CY(IP,KP)=FACE_ID
        ELSEIF (IF_OWN.EQ.0.AND.IF_NEIGH.EQ.1) THEN
        F_DOWN_RECT(IP_NEIGH,KP_NEIGH)=-FACE_ID
        ENDIF
    ENDIF
        
    ENDIF
    
    ENDIF 
    ENDDO    
    ENDIF
    
    ELSEIF (J.EQ.NJ) THEN
    !CASE 3-3 
    NUM(7)=NCENTER
    CALL GET_CELLID(I,J,K-1,NUM(1),1)
    CALL GET_CELLID(I,J,K+1,NUM(2),1)
    NUM(3)=-1
    CALL GET_CELLID(I-1,J,K,NUM(4),1)
    NUM(5)=-1
    CALL GET_CELLID(I,J-1,K,NUM(6),1)
    CALL BUBBLE_SORT(NUM,INDEX_S)
    DO I1=1,7
    IF (NUM(I1).GT.NCENTER) THEN  
    IF_NEIGH=0
    CALL GET_INDEX(NUM(I1),I_NEIGH,J_NEIGH,K_NEIGH,IF_CY)
    IP_NEIGH=I_NEIGH-CNI*(MI-1)
    JP_NEIGH=J_NEIGH-CNJ*(MJ-1)
    KP_NEIGH=K_NEIGH
    IF (IP_NEIGH.LE.CNI.AND.IP_NEIGH.GE.1) THEN
    IF (JP_NEIGH.LE.CNJ.AND.JP_NEIGH.GE.1) THEN
        IF_NEIGH=1
    ENDIF
    ENDIF    
    
    IF (INDEX_S(I1).EQ.1) THEN
        FACE_ID=FACE_ID+1
        IF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.1) THEN
        WRITE(3,231) FACE_ID
        ENDIF
    ENDIF
    IF (INDEX_S(I1).EQ.2) THEN
        FACE_ID=FACE_ID+1
        IF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.1) THEN
        WRITE(3,231) FACE_ID
        ENDIF
    ENDIF
    IF (INDEX_S(I1).EQ.3) THEN
        FACE_ID=FACE_ID+1
        IF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.1) THEN
        WRITE(3,231) FACE_ID
        ELSEIF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.0) THEN
        F_LEFT_CY(JP,KP)=FACE_ID
        ELSEIF (IF_OWN.EQ.0.AND.IF_NEIGH.EQ.1) THEN
        F_RIGHT_CY(JP_NEIGH,KP_NEIGH)=-FACE_ID
        ENDIF
    ENDIF  
    IF (INDEX_S(I1).EQ.4) THEN
        FACE_ID=FACE_ID+1
        IF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.1) THEN
        WRITE(3,231) FACE_ID
        ELSEIF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.0) THEN
        F_RIGHT_CY(JP,KP)=FACE_ID
        ELSEIF (IF_OWN.EQ.0.AND.IF_NEIGH.EQ.1) THEN
        F_LEFT_CY(JP_NEIGH,KP_NEIGH)=-FACE_ID
        ENDIF
    ENDIF
    IF (INDEX_S(I1).EQ.5) THEN
        FACE_ID=FACE_ID+1
        IF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.1) THEN
        WRITE(3,231) FACE_ID
        ELSEIF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.0) THEN
        F_UP_CY(IP,KP)=FACE_ID
        ELSEIF (IF_OWN.EQ.0.AND.IF_NEIGH.EQ.1) THEN
        F_DOWN_CY(IP_NEIGH,KP_NEIGH)=-FACE_ID
        ENDIF
    ENDIF
    IF (INDEX_S(I1).EQ.6) THEN
        FACE_ID=FACE_ID+1
        IF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.1) THEN
        WRITE(3,231) FACE_ID
        ELSEIF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.0) THEN
        F_DOWN_CY(IP,KP)=FACE_ID
        ELSEIF (IF_OWN.EQ.0.AND.IF_NEIGH.EQ.1) THEN
        F_UP_CY(IP_NEIGH,KP_NEIGH)=-FACE_ID
        ENDIF
    ENDIF       
    
    ENDIF 
    ENDDO
    ELSE
    !CASE 3-2  
    NUM(7)=NCENTER
    CALL GET_CELLID(I,J,K-1,NUM(1),1)
    CALL GET_CELLID(I,J,K+1,NUM(2),1)
    NUM(3)=-1
    CALL GET_CELLID(I-1,J,K,NUM(4),1)
    CALL GET_CELLID(I,J+1,K,NUM(5),1)
    CALL GET_CELLID(I,J-1,K,NUM(6),1)
    CALL BUBBLE_SORT(NUM,INDEX_S)
    DO I1=1,7
    IF (NUM(I1).GT.NCENTER) THEN  
    IF_NEIGH=0
    CALL GET_INDEX(NUM(I1),I_NEIGH,J_NEIGH,K_NEIGH,IF_CY)
    IP_NEIGH=I_NEIGH-CNI*(MI-1)
    JP_NEIGH=J_NEIGH-CNJ*(MJ-1)
    KP_NEIGH=K_NEIGH
    IF (IP_NEIGH.LE.CNI.AND.IP_NEIGH.GE.1) THEN
    IF (JP_NEIGH.LE.CNJ.AND.JP_NEIGH.GE.1) THEN
        IF_NEIGH=1
    ENDIF
    ENDIF    
    
    IF (INDEX_S(I1).EQ.1) THEN
        FACE_ID=FACE_ID+1
        IF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.1) THEN
        WRITE(3,231) FACE_ID
        ENDIF
    ENDIF
    IF (INDEX_S(I1).EQ.2) THEN
        FACE_ID=FACE_ID+1
        IF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.1) THEN
        WRITE(3,231) FACE_ID
        ENDIF
    ENDIF
    IF (INDEX_S(I1).EQ.3) THEN
        FACE_ID=FACE_ID+1
        IF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.1) THEN
        WRITE(3,231) FACE_ID
        ELSEIF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.0) THEN
        F_LEFT_CY(JP,KP)=FACE_ID
        ELSEIF (IF_OWN.EQ.0.AND.IF_NEIGH.EQ.1) THEN
        F_RIGHT_CY(JP_NEIGH,KP_NEIGH)=-FACE_ID
        ENDIF
    ENDIF  
    IF (INDEX_S(I1).EQ.4) THEN
        FACE_ID=FACE_ID+1
        IF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.1) THEN
        WRITE(3,231) FACE_ID
        ELSEIF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.0) THEN
        F_RIGHT_CY(JP,KP)=FACE_ID
        ELSEIF (IF_OWN.EQ.0.AND.IF_NEIGH.EQ.1) THEN
        F_LEFT_CY(JP_NEIGH,KP_NEIGH)=-FACE_ID
        ENDIF
    ENDIF
    IF (INDEX_S(I1).EQ.5) THEN
        FACE_ID=FACE_ID+1
        IF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.1) THEN
        WRITE(3,231) FACE_ID
        ELSEIF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.0) THEN
        F_UP_CY(IP,KP)=FACE_ID
        ELSEIF (IF_OWN.EQ.0.AND.IF_NEIGH.EQ.1) THEN
        F_DOWN_CY(IP_NEIGH,KP_NEIGH)=-FACE_ID
        ENDIF
    ENDIF
    IF (INDEX_S(I1).EQ.6) THEN
        FACE_ID=FACE_ID+1
        IF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.1) THEN
        WRITE(3,231) FACE_ID
        ELSEIF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.0) THEN
        F_DOWN_CY(IP,KP)=FACE_ID
        ELSEIF (IF_OWN.EQ.0.AND.IF_NEIGH.EQ.1) THEN
        F_UP_CY(IP_NEIGH,KP_NEIGH)=-FACE_ID
        ENDIF
    ENDIF       
    
    ENDIF 
    ENDDO
    ENDIF    
    ELSE
    IF (J.EQ.1) THEN
    IF (K.LE.NL) THEN
    !CASE 2-1-1 
    NUM(7)=NCENTER
    CALL GET_CELLID(I,J,K-1,NUM(1),1)
    CALL GET_CELLID(I,J,K+1,NUM(2),1)
    CALL GET_CELLID(I+1,J,K,NUM(3),1)
    CALL GET_CELLID(I-1,J,K,NUM(4),1)
    CALL GET_CELLID(I,J+1,K,NUM(5),1)
    CALL GET_CELLID(I,K,1,NUM(6),0)
    CALL BUBBLE_SORT(NUM,INDEX_S)
    DO I1=1,7
    IF (NUM(I1).GT.NCENTER) THEN  
    IF_NEIGH=0
    CALL GET_INDEX(NUM(I1),I_NEIGH,J_NEIGH,K_NEIGH,IF_CY)
    IF (IF_CY.EQ.1) THEN
    IP_NEIGH=I_NEIGH-CNI*(MI-1)
    JP_NEIGH=J_NEIGH-CNJ*(MJ-1)
    KP_NEIGH=K_NEIGH
    IF (IP_NEIGH.LE.CNI.AND.IP_NEIGH.GE.1) THEN
    IF (JP_NEIGH.LE.CNJ.AND.JP_NEIGH.GE.1) THEN
        IF_NEIGH=1
    ENDIF
    ENDIF 
    ELSE
    IP_NEIGH=I_NEIGH-CNI*(MI-1)
    JP_NEIGH=J_NEIGH-CNL*(MJ1-1)
    KP_NEIGH=K_NEIGH-CNL*(MJ2-1)
    IF (IP_NEIGH.LE.CNI.AND.IP_NEIGH.GE.1) THEN
    IF (JP_NEIGH.LE.CNL.AND.JP_NEIGH.GE.1) THEN
    IF (KP_NEIGH.LE.CNL.AND.KP_NEIGH.GE.1) THEN
        IF_NEIGH=1 
    ENDIF
    ENDIF
    ENDIF   
    ENDIF
    
    IF (IF_CY.EQ.1) THEN
    IF (INDEX_S(I1).EQ.1) THEN
        FACE_ID=FACE_ID+1
        IF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.1) THEN
        WRITE(3,231) FACE_ID
        ENDIF
    ENDIF
    IF (INDEX_S(I1).EQ.2) THEN
        FACE_ID=FACE_ID+1
        IF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.1) THEN
        WRITE(3,231) FACE_ID
        ENDIF
    ENDIF
    IF (INDEX_S(I1).EQ.3) THEN
        FACE_ID=FACE_ID+1
        IF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.1) THEN
        WRITE(3,231) FACE_ID
        ELSEIF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.0) THEN
        F_LEFT_CY(JP,KP)=FACE_ID
        ELSEIF (IF_OWN.EQ.0.AND.IF_NEIGH.EQ.1) THEN
        F_RIGHT_CY(JP_NEIGH,KP_NEIGH)=-FACE_ID
        ENDIF
    ENDIF  
    IF (INDEX_S(I1).EQ.4) THEN
        FACE_ID=FACE_ID+1
        IF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.1) THEN
        WRITE(3,231) FACE_ID
        ELSEIF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.0) THEN
        F_RIGHT_CY(JP,KP)=FACE_ID
        ELSEIF (IF_OWN.EQ.0.AND.IF_NEIGH.EQ.1) THEN
        F_LEFT_CY(JP_NEIGH,KP_NEIGH)=-FACE_ID
        ENDIF
    ENDIF
    IF (INDEX_S(I1).EQ.5) THEN
        FACE_ID=FACE_ID+1
        IF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.1) THEN
        WRITE(3,231) FACE_ID
        ELSEIF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.0) THEN
        F_UP_CY(IP,KP)=FACE_ID
        ELSEIF (IF_OWN.EQ.0.AND.IF_NEIGH.EQ.1) THEN
        F_DOWN_CY(IP_NEIGH,KP_NEIGH)=-FACE_ID
        ENDIF
    ENDIF
    IF (INDEX_S(I1).EQ.6) THEN
        FACE_ID=FACE_ID+1
        IF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.1) THEN
        WRITE(3,231) FACE_ID
        ELSEIF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.0) THEN
        F_DOWN_CY(IP,KP)=FACE_ID
        ELSEIF (IF_OWN.EQ.0.AND.IF_NEIGH.EQ.1) THEN
        F_UP_CY(IP_NEIGH,KP_NEIGH)=-FACE_ID
        ENDIF
    ENDIF       
     
    ELSE
    IF (INDEX_S(I1).EQ.6) THEN
        FACE_ID=FACE_ID+1
        IF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.0) THEN
        F_DOWN_CY(IP,KP)=FACE_ID
        ELSEIF (IF_OWN.EQ.0.AND.IF_NEIGH.EQ.1) THEN
        F_FRONT_RECT(IP_NEIGH,JP_NEIGH)=-FACE_ID
        ENDIF
    ENDIF
        
    ENDIF
    
    ENDIF 
    ENDDO
    
    
    
    
    ELSEIF (K.GT.NL .AND. K.LE.(2*NL)) THEN
    !CASE 2-1-2
    NUM(7)=NCENTER
    CALL GET_CELLID(I,J,K-1,NUM(1),1)
    CALL GET_CELLID(I,J,K+1,NUM(2),1)
    CALL GET_CELLID(I+1,J,K,NUM(3),1)
    CALL GET_CELLID(I-1,J,K,NUM(4),1)
    CALL GET_CELLID(I,J+1,K,NUM(5),1)
    CALL GET_CELLID(I,NL,K-NL,NUM(6),0)
    CALL BUBBLE_SORT(NUM,INDEX_S)
    DO I1=1,7
    IF (NUM(I1).GT.NCENTER) THEN  
    IF_NEIGH=0
    CALL GET_INDEX(NUM(I1),I_NEIGH,J_NEIGH,K_NEIGH,IF_CY)
    IF (IF_CY.EQ.1) THEN
    IP_NEIGH=I_NEIGH-CNI*(MI-1)
    JP_NEIGH=J_NEIGH-CNJ*(MJ-1)
    KP_NEIGH=K_NEIGH
    IF (IP_NEIGH.LE.CNI.AND.IP_NEIGH.GE.1) THEN
    IF (JP_NEIGH.LE.CNJ.AND.JP_NEIGH.GE.1) THEN
        IF_NEIGH=1
    ENDIF
    ENDIF 
    ELSE
    IP_NEIGH=I_NEIGH-CNI*(MI-1)
    JP_NEIGH=J_NEIGH-CNL*(MJ1-1)
    KP_NEIGH=K_NEIGH-CNL*(MJ2-1)
    IF (IP_NEIGH.LE.CNI.AND.IP_NEIGH.GE.1) THEN
    IF (JP_NEIGH.LE.CNL.AND.JP_NEIGH.GE.1) THEN
    IF (KP_NEIGH.LE.CNL.AND.KP_NEIGH.GE.1) THEN
        IF_NEIGH=1 
    ENDIF
    ENDIF
    ENDIF   
    ENDIF
    
    IF (IF_CY.EQ.1) THEN
    IF (INDEX_S(I1).EQ.1) THEN
        FACE_ID=FACE_ID+1
        IF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.1) THEN
        WRITE(3,231) FACE_ID
        ENDIF
    ENDIF
    IF (INDEX_S(I1).EQ.2) THEN
        FACE_ID=FACE_ID+1
        IF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.1) THEN
        WRITE(3,231) FACE_ID
        ENDIF
    ENDIF
    IF (INDEX_S(I1).EQ.3) THEN
        FACE_ID=FACE_ID+1
        IF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.1) THEN
        WRITE(3,231) FACE_ID
        ELSEIF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.0) THEN
        F_LEFT_CY(JP,KP)=FACE_ID
        ELSEIF (IF_OWN.EQ.0.AND.IF_NEIGH.EQ.1) THEN
        F_RIGHT_CY(JP_NEIGH,KP_NEIGH)=-FACE_ID
        ENDIF
    ENDIF  
    IF (INDEX_S(I1).EQ.4) THEN
        FACE_ID=FACE_ID+1
        IF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.1) THEN
        WRITE(3,231) FACE_ID
        ELSEIF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.0) THEN
        F_RIGHT_CY(JP,KP)=FACE_ID
        ELSEIF (IF_OWN.EQ.0.AND.IF_NEIGH.EQ.1) THEN
        F_LEFT_CY(JP_NEIGH,KP_NEIGH)=-FACE_ID
        ENDIF
    ENDIF
    IF (INDEX_S(I1).EQ.5) THEN
        FACE_ID=FACE_ID+1
        IF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.1) THEN
        WRITE(3,231) FACE_ID
        ELSEIF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.0) THEN
        F_UP_CY(IP,KP)=FACE_ID
        ELSEIF (IF_OWN.EQ.0.AND.IF_NEIGH.EQ.1) THEN
        F_DOWN_CY(IP_NEIGH,KP_NEIGH)=-FACE_ID
        ENDIF
    ENDIF
    IF (INDEX_S(I1).EQ.6) THEN
        FACE_ID=FACE_ID+1
        IF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.1) THEN
        WRITE(3,231) FACE_ID
        ELSEIF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.0) THEN
        F_DOWN_CY(IP,KP)=FACE_ID
        ELSEIF (IF_OWN.EQ.0.AND.IF_NEIGH.EQ.1) THEN
        F_UP_CY(IP_NEIGH,KP_NEIGH)=-FACE_ID
        ENDIF
    ENDIF       
     
    ELSE
    IF (INDEX_S(I1).EQ.6) THEN
        FACE_ID=FACE_ID+1
        IF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.0) THEN
        F_DOWN_CY(IP,KP)=FACE_ID
        ELSEIF (IF_OWN.EQ.0.AND.IF_NEIGH.EQ.1) THEN
        F_UP_RECT(IP_NEIGH,KP_NEIGH)=-FACE_ID
        ENDIF
    ENDIF
        
    ENDIF
    
    ENDIF 
    ENDDO
    ELSEIF (K.GT.(2*NL) .AND. K.LE.(3*NL)) THEN
    !CASE 2-1-3
    NUM(7)=NCENTER
    CALL GET_CELLID(I,J,K-1,NUM(1),1)
    CALL GET_CELLID(I,J,K+1,NUM(2),1)
    CALL GET_CELLID(I+1,J,K,NUM(3),1)
    CALL GET_CELLID(I-1,J,K,NUM(4),1)
    CALL GET_CELLID(I,J+1,K,NUM(5),1)
    CALL GET_CELLID(I,3*NL+1-K,NL,NUM(6),0)
    CALL BUBBLE_SORT(NUM,INDEX_S)
    DO I1=1,7
    IF (NUM(I1).GT.NCENTER) THEN  
    IF_NEIGH=0
    CALL GET_INDEX(NUM(I1),I_NEIGH,J_NEIGH,K_NEIGH,IF_CY)
    IF (IF_CY.EQ.1) THEN
    IP_NEIGH=I_NEIGH-CNI*(MI-1)
    JP_NEIGH=J_NEIGH-CNJ*(MJ-1)
    KP_NEIGH=K_NEIGH
    IF (IP_NEIGH.LE.CNI.AND.IP_NEIGH.GE.1) THEN
    IF (JP_NEIGH.LE.CNJ.AND.JP_NEIGH.GE.1) THEN
        IF_NEIGH=1
    ENDIF
    ENDIF 
    ELSE
    IP_NEIGH=I_NEIGH-CNI*(MI-1)
    JP_NEIGH=J_NEIGH-CNL*(MJ1-1)
    KP_NEIGH=K_NEIGH-CNL*(MJ2-1)
    IF (IP_NEIGH.LE.CNI.AND.IP_NEIGH.GE.1) THEN
    IF (JP_NEIGH.LE.CNL.AND.JP_NEIGH.GE.1) THEN
    IF (KP_NEIGH.LE.CNL.AND.KP_NEIGH.GE.1) THEN
        IF_NEIGH=1 
    ENDIF
    ENDIF
    ENDIF    
    ENDIF
    
    IF (IF_CY.EQ.1) THEN
    IF (INDEX_S(I1).EQ.1) THEN
        FACE_ID=FACE_ID+1
        IF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.1) THEN
        WRITE(3,231) FACE_ID
        ENDIF
    ENDIF
    IF (INDEX_S(I1).EQ.2) THEN
        FACE_ID=FACE_ID+1
        IF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.1) THEN
        WRITE(3,231) FACE_ID
        ENDIF
    ENDIF
    IF (INDEX_S(I1).EQ.3) THEN
        FACE_ID=FACE_ID+1
        IF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.1) THEN
        WRITE(3,231) FACE_ID
        ELSEIF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.0) THEN
        F_LEFT_CY(JP,KP)=FACE_ID
        ELSEIF (IF_OWN.EQ.0.AND.IF_NEIGH.EQ.1) THEN
        F_RIGHT_CY(JP_NEIGH,KP_NEIGH)=-FACE_ID
        ENDIF
    ENDIF  
    IF (INDEX_S(I1).EQ.4) THEN
        FACE_ID=FACE_ID+1
        IF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.1) THEN
        WRITE(3,231) FACE_ID
        ELSEIF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.0) THEN
        F_RIGHT_CY(JP,KP)=FACE_ID
        ELSEIF (IF_OWN.EQ.0.AND.IF_NEIGH.EQ.1) THEN
        F_LEFT_CY(JP_NEIGH,KP_NEIGH)=-FACE_ID
        ENDIF
    ENDIF
    IF (INDEX_S(I1).EQ.5) THEN
        FACE_ID=FACE_ID+1
        IF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.1) THEN
        WRITE(3,231) FACE_ID
        ELSEIF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.0) THEN
        F_UP_CY(IP,KP)=FACE_ID
        ELSEIF (IF_OWN.EQ.0.AND.IF_NEIGH.EQ.1) THEN
        F_DOWN_CY(IP_NEIGH,KP_NEIGH)=-FACE_ID
        ENDIF
    ENDIF
    IF (INDEX_S(I1).EQ.6) THEN
        FACE_ID=FACE_ID+1
        IF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.1) THEN
        WRITE(3,231) FACE_ID
        ELSEIF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.0) THEN
        F_DOWN_CY(IP,KP)=FACE_ID
        ELSEIF (IF_OWN.EQ.0.AND.IF_NEIGH.EQ.1) THEN
        F_UP_CY(IP_NEIGH,KP_NEIGH)=-FACE_ID
        ENDIF
    ENDIF       
     
    ELSE
    IF (INDEX_S(I1).EQ.6) THEN
        FACE_ID=FACE_ID+1
        IF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.0) THEN
        F_DOWN_CY(IP,KP)=FACE_ID
        ELSEIF (IF_OWN.EQ.0.AND.IF_NEIGH.EQ.1) THEN
        F_BACK_RECT(IP_NEIGH,JP_NEIGH)=-FACE_ID
        ENDIF
    ENDIF
        
    ENDIF
    
    ENDIF 
    ENDDO
    ELSE
    !CASE 2-1-4
    NUM(7)=NCENTER
    CALL GET_CELLID(I,J,K-1,NUM(1),1)
    CALL GET_CELLID(I,J,K+1,NUM(2),1)
    CALL GET_CELLID(I+1,J,K,NUM(3),1)
    CALL GET_CELLID(I-1,J,K,NUM(4),1)
    CALL GET_CELLID(I,J+1,K,NUM(5),1)
    CALL GET_CELLID(I,1,4*NL+1-K,NUM(6),0)
    CALL BUBBLE_SORT(NUM,INDEX_S)
    DO I1=1,7
    IF (NUM(I1).GT.NCENTER) THEN  
    IF_NEIGH=0
    CALL GET_INDEX(NUM(I1),I_NEIGH,J_NEIGH,K_NEIGH,IF_CY)
    IF (IF_CY.EQ.1) THEN
    IP_NEIGH=I_NEIGH-CNI*(MI-1)
    JP_NEIGH=J_NEIGH-CNJ*(MJ-1)
    KP_NEIGH=K_NEIGH
    IF (IP_NEIGH.LE.CNI.AND.IP_NEIGH.GE.1) THEN
    IF (JP_NEIGH.LE.CNJ.AND.JP_NEIGH.GE.1) THEN
        IF_NEIGH=1
    ENDIF
    ENDIF 
    ELSE
    IP_NEIGH=I_NEIGH-CNI*(MI-1)
    JP_NEIGH=J_NEIGH-CNL*(MJ1-1)
    KP_NEIGH=K_NEIGH-CNL*(MJ2-1)
    IF (IP_NEIGH.LE.CNI.AND.IP_NEIGH.GE.1) THEN
    IF (JP_NEIGH.LE.CNL.AND.JP_NEIGH.GE.1) THEN
    IF (KP_NEIGH.LE.CNL.AND.KP_NEIGH.GE.1) THEN
        IF_NEIGH=1 
    ENDIF
    ENDIF
    ENDIF    
    ENDIF
    
    IF (IF_CY.EQ.1) THEN
    IF (INDEX_S(I1).EQ.1) THEN
        FACE_ID=FACE_ID+1
        IF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.1) THEN
        WRITE(3,231) FACE_ID
        ENDIF
    ENDIF
    IF (INDEX_S(I1).EQ.2) THEN
        FACE_ID=FACE_ID+1
        IF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.1) THEN
        WRITE(3,231) FACE_ID
        ENDIF
    ENDIF
    IF (INDEX_S(I1).EQ.3) THEN
        FACE_ID=FACE_ID+1
        IF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.1) THEN
        WRITE(3,231) FACE_ID
        ELSEIF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.0) THEN
        F_LEFT_CY(JP,KP)=FACE_ID
        ELSEIF (IF_OWN.EQ.0.AND.IF_NEIGH.EQ.1) THEN
        F_RIGHT_CY(JP_NEIGH,KP_NEIGH)=-FACE_ID
        ENDIF
    ENDIF  
    IF (INDEX_S(I1).EQ.4) THEN
        FACE_ID=FACE_ID+1
        IF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.1) THEN
        WRITE(3,231) FACE_ID
        ELSEIF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.0) THEN
        F_RIGHT_CY(JP,KP)=FACE_ID
        ELSEIF (IF_OWN.EQ.0.AND.IF_NEIGH.EQ.1) THEN
        F_LEFT_CY(JP_NEIGH,KP_NEIGH)=-FACE_ID
        ENDIF
    ENDIF
    IF (INDEX_S(I1).EQ.5) THEN
        FACE_ID=FACE_ID+1
        IF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.1) THEN
        WRITE(3,231) FACE_ID
        ELSEIF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.0) THEN
        F_UP_CY(IP,KP)=FACE_ID
        ELSEIF (IF_OWN.EQ.0.AND.IF_NEIGH.EQ.1) THEN
        F_DOWN_CY(IP_NEIGH,KP_NEIGH)=-FACE_ID
        ENDIF
    ENDIF
    IF (INDEX_S(I1).EQ.6) THEN
        FACE_ID=FACE_ID+1
        IF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.1) THEN
        WRITE(3,231) FACE_ID
        ELSEIF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.0) THEN
        F_DOWN_CY(IP,KP)=FACE_ID
        ELSEIF (IF_OWN.EQ.0.AND.IF_NEIGH.EQ.1) THEN
        F_UP_CY(IP_NEIGH,KP_NEIGH)=-FACE_ID
        ENDIF
    ENDIF       
     
    ELSE
    IF (INDEX_S(I1).EQ.6) THEN
        FACE_ID=FACE_ID+1
        IF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.0) THEN
        F_DOWN_CY(IP,KP)=FACE_ID
        ELSEIF (IF_OWN.EQ.0.AND.IF_NEIGH.EQ.1) THEN
        F_DOWN_RECT(IP_NEIGH,KP_NEIGH)=-FACE_ID
        ENDIF
    ENDIF
        
    ENDIF
    
    ENDIF 
    ENDDO    
    ENDIF
    ELSEIF (J.EQ.NJ) THEN
    !CASE 2-3 
    NUM(7)=NCENTER
    CALL GET_CELLID(I,J,K-1,NUM(1),1)
    CALL GET_CELLID(I,J,K+1,NUM(2),1)
    CALL GET_CELLID(I+1,J,K,NUM(3),1)
    CALL GET_CELLID(I-1,J,K,NUM(4),1)
    NUM(5)=-1
    CALL GET_CELLID(I,J-1,K,NUM(6),1)
    CALL BUBBLE_SORT(NUM,INDEX_S)
    DO I1=1,7
    IF (NUM(I1).GT.NCENTER) THEN  
    IF_NEIGH=0
    CALL GET_INDEX(NUM(I1),I_NEIGH,J_NEIGH,K_NEIGH,IF_CY)
    IP_NEIGH=I_NEIGH-CNI*(MI-1)
    JP_NEIGH=J_NEIGH-CNJ*(MJ-1)
    KP_NEIGH=K_NEIGH
    IF (IP_NEIGH.LE.CNI.AND.IP_NEIGH.GE.1) THEN
    IF (JP_NEIGH.LE.CNJ.AND.JP_NEIGH.GE.1) THEN
        IF_NEIGH=1
    ENDIF
    ENDIF    
    
    IF (INDEX_S(I1).EQ.1) THEN
        FACE_ID=FACE_ID+1
        IF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.1) THEN
        WRITE(3,231) FACE_ID
        ENDIF
    ENDIF
    IF (INDEX_S(I1).EQ.2) THEN
        FACE_ID=FACE_ID+1
        IF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.1) THEN
        WRITE(3,231) FACE_ID
        ENDIF
    ENDIF
    IF (INDEX_S(I1).EQ.3) THEN
        FACE_ID=FACE_ID+1
        IF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.1) THEN
        WRITE(3,231) FACE_ID
        ELSEIF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.0) THEN
        F_LEFT_CY(JP,KP)=FACE_ID
        ELSEIF (IF_OWN.EQ.0.AND.IF_NEIGH.EQ.1) THEN
        F_RIGHT_CY(JP_NEIGH,KP_NEIGH)=-FACE_ID
        ENDIF
    ENDIF  
    IF (INDEX_S(I1).EQ.4) THEN
        FACE_ID=FACE_ID+1
        IF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.1) THEN
        WRITE(3,231) FACE_ID
        ELSEIF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.0) THEN
        F_RIGHT_CY(JP,KP)=FACE_ID
        ELSEIF (IF_OWN.EQ.0.AND.IF_NEIGH.EQ.1) THEN
        F_LEFT_CY(JP_NEIGH,KP_NEIGH)=-FACE_ID
        ENDIF
    ENDIF
    IF (INDEX_S(I1).EQ.5) THEN
        FACE_ID=FACE_ID+1
        IF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.1) THEN
        WRITE(3,231) FACE_ID
        ELSEIF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.0) THEN
        F_UP_CY(IP,KP)=FACE_ID
        ELSEIF (IF_OWN.EQ.0.AND.IF_NEIGH.EQ.1) THEN
        F_DOWN_CY(IP_NEIGH,KP_NEIGH)=-FACE_ID
        ENDIF
    ENDIF
    IF (INDEX_S(I1).EQ.6) THEN
        FACE_ID=FACE_ID+1
        IF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.1) THEN
        WRITE(3,231) FACE_ID
        ELSEIF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.0) THEN
        F_DOWN_CY(IP,KP)=FACE_ID
        ELSEIF (IF_OWN.EQ.0.AND.IF_NEIGH.EQ.1) THEN
        F_UP_CY(IP_NEIGH,KP_NEIGH)=-FACE_ID
        ENDIF
    ENDIF       
    
    ENDIF 
    ENDDO
    ELSE
    !CASE 2-2 
    NUM(7)=NCENTER
    CALL GET_CELLID(I,J,K-1,NUM(1),1)
    CALL GET_CELLID(I,J,K+1,NUM(2),1)
    CALL GET_CELLID(I+1,J,K,NUM(3),1)
    CALL GET_CELLID(I-1,J,K,NUM(4),1)
    CALL GET_CELLID(I,J+1,K,NUM(5),1)
    CALL GET_CELLID(I,J-1,K,NUM(6),1)
    CALL BUBBLE_SORT(NUM,INDEX_S)
    DO I1=1,7
    IF (NUM(I1).GT.NCENTER) THEN  
    IF_NEIGH=0
    CALL GET_INDEX(NUM(I1),I_NEIGH,J_NEIGH,K_NEIGH,IF_CY)
    IP_NEIGH=I_NEIGH-CNI*(MI-1)
    JP_NEIGH=J_NEIGH-CNJ*(MJ-1)
    KP_NEIGH=K_NEIGH
    IF (IP_NEIGH.LE.CNI.AND.IP_NEIGH.GE.1) THEN
    IF (JP_NEIGH.LE.CNJ.AND.JP_NEIGH.GE.1) THEN
        IF_NEIGH=1
    ENDIF
    ENDIF    
    
    IF (INDEX_S(I1).EQ.1) THEN
        FACE_ID=FACE_ID+1
        IF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.1) THEN
        WRITE(3,231) FACE_ID
        ENDIF
    ENDIF
    IF (INDEX_S(I1).EQ.2) THEN
        FACE_ID=FACE_ID+1
        IF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.1) THEN
        WRITE(3,231) FACE_ID
        ENDIF
    ENDIF
    IF (INDEX_S(I1).EQ.3) THEN
        FACE_ID=FACE_ID+1
        IF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.1) THEN
        WRITE(3,231) FACE_ID
        ELSEIF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.0) THEN
        F_LEFT_CY(JP,KP)=FACE_ID
        ELSEIF (IF_OWN.EQ.0.AND.IF_NEIGH.EQ.1) THEN
        F_RIGHT_CY(JP_NEIGH,KP_NEIGH)=-FACE_ID
        ENDIF
    ENDIF  
    IF (INDEX_S(I1).EQ.4) THEN
        FACE_ID=FACE_ID+1
        IF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.1) THEN
        WRITE(3,231) FACE_ID
        ELSEIF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.0) THEN
        F_RIGHT_CY(JP,KP)=FACE_ID
        ELSEIF (IF_OWN.EQ.0.AND.IF_NEIGH.EQ.1) THEN
        F_LEFT_CY(JP_NEIGH,KP_NEIGH)=-FACE_ID
        ENDIF
    ENDIF
    IF (INDEX_S(I1).EQ.5) THEN
        FACE_ID=FACE_ID+1
        IF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.1) THEN
        WRITE(3,231) FACE_ID
        ELSEIF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.0) THEN
        F_UP_CY(IP,KP)=FACE_ID
        ELSEIF (IF_OWN.EQ.0.AND.IF_NEIGH.EQ.1) THEN
        F_DOWN_CY(IP_NEIGH,KP_NEIGH)=-FACE_ID
        ENDIF
    ENDIF
    IF (INDEX_S(I1).EQ.6) THEN
        FACE_ID=FACE_ID+1
        IF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.1) THEN
        WRITE(3,231) FACE_ID
        ELSEIF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.0) THEN
        F_DOWN_CY(IP,KP)=FACE_ID
        ELSEIF (IF_OWN.EQ.0.AND.IF_NEIGH.EQ.1) THEN
        F_UP_CY(IP_NEIGH,KP_NEIGH)=-FACE_ID
        ENDIF
    ENDIF       
    
    ENDIF 
    ENDDO
    ENDIF    
    ENDIF
    
    
    
    
    
    
    
    ELSEIF (IF_CY.EQ.0) THEN
!*******************RECTANGLE*************** 
    IP=I-CNI*(MI-1)
    JP=J-CNL*(MJ1-1)
    KP=K-CNL*(MJ2-1)
    IF (IP.LE.CNI.AND.IP.GE.1) THEN
    IF (JP.LE.CNL.AND.JP.GE.1) THEN
    IF (KP.LE.CNL.AND.KP.GE.1) THEN
        IF_OWN=1  
    ENDIF
    ENDIF
    ENDIF
    
    IF (I.EQ.1) THEN
    IF (J.EQ.1) THEN
    IF (K.EQ.1) THEN
    !CASE 1-1-1 
    NUM(7)=NCENTER
    CALL GET_CELLID(I,1,1,NUM(1),1)
    CALL GET_CELLID(I,J,K+1,NUM(2),0)
    CALL GET_CELLID(I+1,J,K,NUM(3),0)
    NUM(4)=-1
    CALL GET_CELLID(I,J+1,K,NUM(5),0)
    CALL GET_CELLID(I,1,NK,NUM(6),1)
    CALL BUBBLE_SORT(NUM,INDEX_S)
    DO I1=1,7
    IF (NUM(I1).GT.NCENTER) THEN  
    IF_NEIGH=0
    CALL GET_INDEX(NUM(I1),I_NEIGH,J_NEIGH,K_NEIGH,IF_CY)
    IF (IF_CY.EQ.1) THEN
    IP_NEIGH=I_NEIGH-CNI*(MI-1)
    JP_NEIGH=J_NEIGH-CNJ*(MJ-1)
    KP_NEIGH=K_NEIGH
    IF (IP_NEIGH.LE.CNI.AND.IP_NEIGH.GE.1) THEN
    IF (JP_NEIGH.LE.CNJ.AND.JP_NEIGH.GE.1) THEN
        IF_NEIGH=1
    ENDIF
    ENDIF 
    ELSE
    IP_NEIGH=I_NEIGH-CNI*(MI-1)
    JP_NEIGH=J_NEIGH-CNL*(MJ1-1)
    KP_NEIGH=K_NEIGH-CNL*(MJ2-1)
    IF (IP_NEIGH.LE.CNI.AND.IP_NEIGH.GE.1) THEN
    IF (JP_NEIGH.LE.CNL.AND.JP_NEIGH.GE.1) THEN
    IF (KP_NEIGH.LE.CNL.AND.KP_NEIGH.GE.1) THEN
        IF_NEIGH=1 
    ENDIF
    ENDIF
    ENDIF    
    ENDIF
    
    IF (IF_CY.EQ.0) THEN
    IF (INDEX_S(I1).EQ.1) THEN
        FACE_ID=FACE_ID+1
        IF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.1) THEN
        WRITE(3,231) FACE_ID
        ELSEIF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.0) THEN
        F_FRONT_RECT(IP,JP)=FACE_ID
        ELSEIF (IF_OWN.EQ.0.AND.IF_NEIGH.EQ.1) THEN
        F_BACK_RECT(IP_NEIGH,JP_NEIGH)=-FACE_ID
        ENDIF
    ENDIF
    IF (INDEX_S(I1).EQ.2) THEN
        FACE_ID=FACE_ID+1
        IF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.1) THEN
        WRITE(3,231) FACE_ID
        ELSEIF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.0) THEN
        F_BACK_RECT(IP,JP)=FACE_ID
        ELSEIF (IF_OWN.EQ.0.AND.IF_NEIGH.EQ.1) THEN
        F_FRONT_RECT(IP_NEIGH,JP_NEIGH)=-FACE_ID
        ENDIF
    ENDIF
    IF (INDEX_S(I1).EQ.3) THEN
        FACE_ID=FACE_ID+1
        IF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.1) THEN
        WRITE(3,231) FACE_ID
        ELSEIF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.0) THEN
        F_LEFT_RECT(JP,KP)=FACE_ID
        ELSEIF (IF_OWN.EQ.0.AND.IF_NEIGH.EQ.1) THEN
        F_RIGHT_RECT(JP_NEIGH,KP_NEIGH)=-FACE_ID
        ENDIF
    ENDIF  
    IF (INDEX_S(I1).EQ.4) THEN
        FACE_ID=FACE_ID+1
        IF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.1) THEN
        WRITE(3,231) FACE_ID
        ELSEIF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.0) THEN
        F_RIGHT_RECT(JP,KP)=FACE_ID
        ELSEIF (IF_OWN.EQ.0.AND.IF_NEIGH.EQ.1) THEN
        F_LEFT_RECT(JP_NEIGH,KP_NEIGH)=-FACE_ID
        ENDIF
    ENDIF
    IF (INDEX_S(I1).EQ.5) THEN
        FACE_ID=FACE_ID+1
        IF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.1) THEN
        WRITE(3,231) FACE_ID
        ELSEIF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.0) THEN
        F_UP_RECT(IP,KP)=FACE_ID
        ELSEIF (IF_OWN.EQ.0.AND.IF_NEIGH.EQ.1) THEN
        F_DOWN_RECT(IP_NEIGH,KP_NEIGH)=-FACE_ID
        ENDIF
    ENDIF
    IF (INDEX_S(I1).EQ.6) THEN
        FACE_ID=FACE_ID+1
        IF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.1) THEN
        WRITE(3,231) FACE_ID
        ELSEIF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.0) THEN
        F_DOWN_RECT(IP,KP)=FACE_ID
        ELSEIF (IF_OWN.EQ.0.AND.IF_NEIGH.EQ.1) THEN
        F_UP_RECT(IP_NEIGH,KP_NEIGH)=-FACE_ID
        ENDIF
    ENDIF       
     
    ELSE
    IF (INDEX_S(I1).EQ.1) THEN
        FACE_ID=FACE_ID+1
        IF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.0) THEN
        F_FRONT_RECT(IP,JP)=FACE_ID
        ELSEIF (IF_OWN.EQ.0.AND.IF_NEIGH.EQ.1) THEN
        F_DOWN_CY(IP_NEIGH,KP_NEIGH)=-FACE_ID
        ENDIF
    ENDIF
    IF (INDEX_S(I1).EQ.6) THEN
        FACE_ID=FACE_ID+1
        IF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.0) THEN
        F_DOWN_RECT(IP,KP)=FACE_ID
        ELSEIF (IF_OWN.EQ.0.AND.IF_NEIGH.EQ.1) THEN
        F_DOWN_CY(IP_NEIGH,KP_NEIGH)=-FACE_ID
        ENDIF
    ENDIF
        
    ENDIF
    
    ENDIF
    ENDDO    
        
    ELSEIF (K.EQ.NL) THEN
    !CASE 1-1-3
    NUM(7)=NCENTER
    CALL GET_CELLID(I,J,K-1,NUM(1),0)
    CALL GET_CELLID(I,1,3*NL,NUM(2),1)
    CALL GET_CELLID(I+1,J,K,NUM(3),0)
    NUM(4)=-1
    CALL GET_CELLID(I,J+1,K,NUM(5),0)
    CALL GET_CELLID(I,1,3*NL+1,NUM(6),1)
    CALL BUBBLE_SORT(NUM,INDEX_S)
    DO I1=1,7
    IF (NUM(I1).GT.NCENTER) THEN  
    IF_NEIGH=0
    CALL GET_INDEX(NUM(I1),I_NEIGH,J_NEIGH,K_NEIGH,IF_CY)
    IF (IF_CY.EQ.1) THEN
    IP_NEIGH=I_NEIGH-CNI*(MI-1)
    JP_NEIGH=J_NEIGH-CNJ*(MJ-1)
    KP_NEIGH=K_NEIGH
    IF (IP_NEIGH.LE.CNI.AND.IP_NEIGH.GE.1) THEN
    IF (JP_NEIGH.LE.CNJ.AND.JP_NEIGH.GE.1) THEN
        IF_NEIGH=1
    ENDIF
    ENDIF 
    ELSE
    IP_NEIGH=I_NEIGH-CNI*(MI-1)
    JP_NEIGH=J_NEIGH-CNL*(MJ1-1)
    KP_NEIGH=K_NEIGH-CNL*(MJ2-1)
    IF (IP_NEIGH.LE.CNI.AND.IP_NEIGH.GE.1) THEN
    IF (JP_NEIGH.LE.CNL.AND.JP_NEIGH.GE.1) THEN
    IF (KP_NEIGH.LE.CNL.AND.KP_NEIGH.GE.1) THEN
        IF_NEIGH=1 
    ENDIF
    ENDIF
    ENDIF    
    ENDIF
    
    IF (IF_CY.EQ.0) THEN
    IF (INDEX_S(I1).EQ.1) THEN
        FACE_ID=FACE_ID+1
        IF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.1) THEN
        WRITE(3,231) FACE_ID
        ELSEIF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.0) THEN
        F_FRONT_RECT(IP,JP)=FACE_ID
        ELSEIF (IF_OWN.EQ.0.AND.IF_NEIGH.EQ.1) THEN
        F_BACK_RECT(IP_NEIGH,JP_NEIGH)=-FACE_ID
        ENDIF
    ENDIF
    IF (INDEX_S(I1).EQ.2) THEN
        FACE_ID=FACE_ID+1
        IF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.1) THEN
        WRITE(3,231) FACE_ID
        ELSEIF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.0) THEN
        F_BACK_RECT(IP,JP)=FACE_ID
        ELSEIF (IF_OWN.EQ.0.AND.IF_NEIGH.EQ.1) THEN
        F_FRONT_RECT(IP_NEIGH,JP_NEIGH)=-FACE_ID
        ENDIF
    ENDIF
    IF (INDEX_S(I1).EQ.3) THEN
        FACE_ID=FACE_ID+1
        IF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.1) THEN
        WRITE(3,231) FACE_ID
        ELSEIF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.0) THEN
        F_LEFT_RECT(JP,KP)=FACE_ID
        ELSEIF (IF_OWN.EQ.0.AND.IF_NEIGH.EQ.1) THEN
        F_RIGHT_RECT(JP_NEIGH,KP_NEIGH)=-FACE_ID
        ENDIF
    ENDIF  
    IF (INDEX_S(I1).EQ.4) THEN
        FACE_ID=FACE_ID+1
        IF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.1) THEN
        WRITE(3,231) FACE_ID
        ELSEIF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.0) THEN
        F_RIGHT_RECT(JP,KP)=FACE_ID
        ELSEIF (IF_OWN.EQ.0.AND.IF_NEIGH.EQ.1) THEN
        F_LEFT_RECT(JP_NEIGH,KP_NEIGH)=-FACE_ID
        ENDIF
    ENDIF
    IF (INDEX_S(I1).EQ.5) THEN
        FACE_ID=FACE_ID+1
        IF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.1) THEN
        WRITE(3,231) FACE_ID
        ELSEIF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.0) THEN
        F_UP_RECT(IP,KP)=FACE_ID
        ELSEIF (IF_OWN.EQ.0.AND.IF_NEIGH.EQ.1) THEN
        F_DOWN_RECT(IP_NEIGH,KP_NEIGH)=-FACE_ID
        ENDIF
    ENDIF
    IF (INDEX_S(I1).EQ.6) THEN
        FACE_ID=FACE_ID+1
        IF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.1) THEN
        WRITE(3,231) FACE_ID
        ELSEIF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.0) THEN
        F_DOWN_RECT(IP,KP)=FACE_ID
        ELSEIF (IF_OWN.EQ.0.AND.IF_NEIGH.EQ.1) THEN
        F_UP_RECT(IP_NEIGH,KP_NEIGH)=-FACE_ID
        ENDIF
    ENDIF       
     
    ELSE
    IF (INDEX_S(I1).EQ.2) THEN
        FACE_ID=FACE_ID+1
        IF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.0) THEN
        F_BACK_RECT(IP,JP)=FACE_ID
        ELSEIF (IF_OWN.EQ.0.AND.IF_NEIGH.EQ.1) THEN
        F_DOWN_CY(IP_NEIGH,KP_NEIGH)=-FACE_ID
        ENDIF
    ENDIF
    IF (INDEX_S(I1).EQ.6) THEN
        FACE_ID=FACE_ID+1
        IF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.0) THEN
        F_DOWN_RECT(IP,KP)=FACE_ID
        ELSEIF (IF_OWN.EQ.0.AND.IF_NEIGH.EQ.1) THEN
        F_DOWN_CY(IP_NEIGH,KP_NEIGH)=-FACE_ID
        ENDIF
    ENDIF
        
    ENDIF
    
    ENDIF
    ENDDO    
    
    ELSE
    !CASE 1-1-2 
    NUM(7)=NCENTER
    CALL GET_CELLID(I,J,K-1,NUM(1),0)
    CALL GET_CELLID(I,J,K+1,NUM(2),0)
    CALL GET_CELLID(I+1,J,K,NUM(3),0)
    NUM(4)=-1
    CALL GET_CELLID(I,J+1,K,NUM(5),0)
    CALL GET_CELLID(I,1,NK+1-K,NUM(6),1)
    CALL BUBBLE_SORT(NUM,INDEX_S)
    DO I1=1,7
    IF (NUM(I1).GT.NCENTER) THEN  
    IF_NEIGH=0
    CALL GET_INDEX(NUM(I1),I_NEIGH,J_NEIGH,K_NEIGH,IF_CY)
    IF (IF_CY.EQ.1) THEN
    IP_NEIGH=I_NEIGH-CNI*(MI-1)
    JP_NEIGH=J_NEIGH-CNJ*(MJ-1)
    KP_NEIGH=K_NEIGH
    IF (IP_NEIGH.LE.CNI.AND.IP_NEIGH.GE.1) THEN
    IF (JP_NEIGH.LE.CNJ.AND.JP_NEIGH.GE.1) THEN
        IF_NEIGH=1
    ENDIF
    ENDIF 
    ELSE
    IP_NEIGH=I_NEIGH-CNI*(MI-1)
    JP_NEIGH=J_NEIGH-CNL*(MJ1-1)
    KP_NEIGH=K_NEIGH-CNL*(MJ2-1)
    IF (IP_NEIGH.LE.CNI.AND.IP_NEIGH.GE.1) THEN
    IF (JP_NEIGH.LE.CNL.AND.JP_NEIGH.GE.1) THEN
    IF (KP_NEIGH.LE.CNL.AND.KP_NEIGH.GE.1) THEN
        IF_NEIGH=1 
    ENDIF
    ENDIF
    ENDIF    
    ENDIF
    
    IF (IF_CY.EQ.0) THEN
    IF (INDEX_S(I1).EQ.1) THEN
        FACE_ID=FACE_ID+1
        IF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.1) THEN
        WRITE(3,231) FACE_ID
        ELSEIF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.0) THEN
        F_FRONT_RECT(IP,JP)=FACE_ID
        ELSEIF (IF_OWN.EQ.0.AND.IF_NEIGH.EQ.1) THEN
        F_BACK_RECT(IP_NEIGH,JP_NEIGH)=-FACE_ID
        ENDIF
    ENDIF
    IF (INDEX_S(I1).EQ.2) THEN
        FACE_ID=FACE_ID+1
        IF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.1) THEN
        WRITE(3,231) FACE_ID
        ELSEIF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.0) THEN
        F_BACK_RECT(IP,JP)=FACE_ID
        ELSEIF (IF_OWN.EQ.0.AND.IF_NEIGH.EQ.1) THEN
        F_FRONT_RECT(IP_NEIGH,JP_NEIGH)=-FACE_ID
        ENDIF
    ENDIF
    IF (INDEX_S(I1).EQ.3) THEN
        FACE_ID=FACE_ID+1
        IF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.1) THEN
        WRITE(3,231) FACE_ID
        ELSEIF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.0) THEN
        F_LEFT_RECT(JP,KP)=FACE_ID
        ELSEIF (IF_OWN.EQ.0.AND.IF_NEIGH.EQ.1) THEN
        F_RIGHT_RECT(JP_NEIGH,KP_NEIGH)=-FACE_ID
        ENDIF
    ENDIF  
    IF (INDEX_S(I1).EQ.4) THEN
        FACE_ID=FACE_ID+1
        IF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.1) THEN
        WRITE(3,231) FACE_ID
        ELSEIF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.0) THEN
        F_RIGHT_RECT(JP,KP)=FACE_ID
        ELSEIF (IF_OWN.EQ.0.AND.IF_NEIGH.EQ.1) THEN
        F_LEFT_RECT(JP_NEIGH,KP_NEIGH)=-FACE_ID
        ENDIF
    ENDIF
    IF (INDEX_S(I1).EQ.5) THEN
        FACE_ID=FACE_ID+1
        IF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.1) THEN
        WRITE(3,231) FACE_ID
        ELSEIF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.0) THEN
        F_UP_RECT(IP,KP)=FACE_ID
        ELSEIF (IF_OWN.EQ.0.AND.IF_NEIGH.EQ.1) THEN
        F_DOWN_RECT(IP_NEIGH,KP_NEIGH)=-FACE_ID
        ENDIF
    ENDIF
    IF (INDEX_S(I1).EQ.6) THEN
        FACE_ID=FACE_ID+1
        IF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.1) THEN
        WRITE(3,231) FACE_ID
        ELSEIF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.0) THEN
        F_DOWN_RECT(IP,KP)=FACE_ID
        ELSEIF (IF_OWN.EQ.0.AND.IF_NEIGH.EQ.1) THEN
        F_UP_RECT(IP_NEIGH,KP_NEIGH)=-FACE_ID
        ENDIF
    ENDIF       
     
    ELSE
    IF (INDEX_S(I1).EQ.6) THEN
        FACE_ID=FACE_ID+1
        IF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.0) THEN
        F_DOWN_RECT(IP,KP)=FACE_ID
        ELSEIF (IF_OWN.EQ.0.AND.IF_NEIGH.EQ.1) THEN
        F_DOWN_CY(IP_NEIGH,KP_NEIGH)=-FACE_ID
        ENDIF
    ENDIF
        
    ENDIF
    
    ENDIF
    ENDDO
    
    ENDIF    
    ELSEIF (J.EQ.NL) THEN
    IF (K.EQ.1) THEN
    !CASE 1-3-1  
    NUM(7)=NCENTER
    CALL GET_CELLID(I,1,NL,NUM(1),1)
    CALL GET_CELLID(I,J,K+1,NUM(2),0)
    CALL GET_CELLID(I+1,J,K,NUM(3),0)
    NUM(4)=-1
    CALL GET_CELLID(I,1,NL+1,NUM(5),1)
    CALL GET_CELLID(I,J-1,K,NUM(6),0)
    CALL BUBBLE_SORT(NUM,INDEX_S)
    DO I1=1,7
    IF (NUM(I1).GT.NCENTER) THEN  
    IF_NEIGH=0
    CALL GET_INDEX(NUM(I1),I_NEIGH,J_NEIGH,K_NEIGH,IF_CY)
    IF (IF_CY.EQ.1) THEN
    IP_NEIGH=I_NEIGH-CNI*(MI-1)
    JP_NEIGH=J_NEIGH-CNJ*(MJ-1)
    KP_NEIGH=K_NEIGH
    IF (IP_NEIGH.LE.CNI.AND.IP_NEIGH.GE.1) THEN
    IF (JP_NEIGH.LE.CNJ.AND.JP_NEIGH.GE.1) THEN
        IF_NEIGH=1
    ENDIF
    ENDIF 
    ELSE
    IP_NEIGH=I_NEIGH-CNI*(MI-1)
    JP_NEIGH=J_NEIGH-CNL*(MJ1-1)
    KP_NEIGH=K_NEIGH-CNL*(MJ2-1)
    IF (IP_NEIGH.LE.CNI.AND.IP_NEIGH.GE.1) THEN
    IF (JP_NEIGH.LE.CNL.AND.JP_NEIGH.GE.1) THEN
    IF (KP_NEIGH.LE.CNL.AND.KP_NEIGH.GE.1) THEN
        IF_NEIGH=1 
    ENDIF
    ENDIF
    ENDIF    
    ENDIF
    
    IF (IF_CY.EQ.0) THEN
    IF (INDEX_S(I1).EQ.1) THEN
        FACE_ID=FACE_ID+1
        IF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.1) THEN
        WRITE(3,231) FACE_ID
        ELSEIF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.0) THEN
        F_FRONT_RECT(IP,JP)=FACE_ID
        ELSEIF (IF_OWN.EQ.0.AND.IF_NEIGH.EQ.1) THEN
        F_BACK_RECT(IP_NEIGH,JP_NEIGH)=-FACE_ID
        ENDIF
    ENDIF
    IF (INDEX_S(I1).EQ.2) THEN
        FACE_ID=FACE_ID+1
        IF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.1) THEN
        WRITE(3,231) FACE_ID
        ELSEIF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.0) THEN
        F_BACK_RECT(IP,JP)=FACE_ID
        ELSEIF (IF_OWN.EQ.0.AND.IF_NEIGH.EQ.1) THEN
        F_FRONT_RECT(IP_NEIGH,JP_NEIGH)=-FACE_ID
        ENDIF
    ENDIF
    IF (INDEX_S(I1).EQ.3) THEN
        FACE_ID=FACE_ID+1
        IF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.1) THEN
        WRITE(3,231) FACE_ID
        ELSEIF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.0) THEN
        F_LEFT_RECT(JP,KP)=FACE_ID
        ELSEIF (IF_OWN.EQ.0.AND.IF_NEIGH.EQ.1) THEN
        F_RIGHT_RECT(JP_NEIGH,KP_NEIGH)=-FACE_ID
        ENDIF
    ENDIF  
    IF (INDEX_S(I1).EQ.4) THEN
        FACE_ID=FACE_ID+1
        IF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.1) THEN
        WRITE(3,231) FACE_ID
        ELSEIF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.0) THEN
        F_RIGHT_RECT(JP,KP)=FACE_ID
        ELSEIF (IF_OWN.EQ.0.AND.IF_NEIGH.EQ.1) THEN
        F_LEFT_RECT(JP_NEIGH,KP_NEIGH)=-FACE_ID
        ENDIF
    ENDIF
    IF (INDEX_S(I1).EQ.5) THEN
        FACE_ID=FACE_ID+1
        IF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.1) THEN
        WRITE(3,231) FACE_ID
        ELSEIF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.0) THEN
        F_UP_RECT(IP,KP)=FACE_ID
        ELSEIF (IF_OWN.EQ.0.AND.IF_NEIGH.EQ.1) THEN
        F_DOWN_RECT(IP_NEIGH,KP_NEIGH)=-FACE_ID
        ENDIF
    ENDIF
    IF (INDEX_S(I1).EQ.6) THEN
        FACE_ID=FACE_ID+1
        IF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.1) THEN
        WRITE(3,231) FACE_ID
        ELSEIF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.0) THEN
        F_DOWN_RECT(IP,KP)=FACE_ID
        ELSEIF (IF_OWN.EQ.0.AND.IF_NEIGH.EQ.1) THEN
        F_UP_RECT(IP_NEIGH,KP_NEIGH)=-FACE_ID
        ENDIF
    ENDIF       
     
    ELSE
    IF (INDEX_S(I1).EQ.1) THEN
        FACE_ID=FACE_ID+1
        IF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.0) THEN
        F_FRONT_RECT(IP,JP)=FACE_ID
        ELSEIF (IF_OWN.EQ.0.AND.IF_NEIGH.EQ.1) THEN
        F_DOWN_CY(IP_NEIGH,KP_NEIGH)=-FACE_ID
        ENDIF
    ENDIF
    IF (INDEX_S(I1).EQ.5) THEN
        FACE_ID=FACE_ID+1
        IF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.0) THEN
        F_UP_RECT(IP,KP)=FACE_ID
        ELSEIF (IF_OWN.EQ.0.AND.IF_NEIGH.EQ.1) THEN
        F_DOWN_CY(IP_NEIGH,KP_NEIGH)=-FACE_ID
        ENDIF
    ENDIF
        
    ENDIF
    
    ENDIF
    ENDDO    
     
    ELSEIF (K.EQ.NL) THEN
    !CASE 1-3-3 
    NUM(7)=NCENTER
    CALL GET_CELLID(I,J,K-1,NUM(1),0)
    CALL GET_CELLID(I,1,2*NL+1,NUM(2),1)
    CALL GET_CELLID(I+1,J,K,NUM(3),0)
    NUM(4)=-1
    CALL GET_CELLID(I,1,2*NL,NUM(5),1)
    CALL GET_CELLID(I,J-1,K,NUM(6),0)
    CALL BUBBLE_SORT(NUM,INDEX_S)
    DO I1=1,7
    IF (NUM(I1).GT.NCENTER) THEN  
    IF_NEIGH=0
    CALL GET_INDEX(NUM(I1),I_NEIGH,J_NEIGH,K_NEIGH,IF_CY)
    IF (IF_CY.EQ.1) THEN
    IP_NEIGH=I_NEIGH-CNI*(MI-1)
    JP_NEIGH=J_NEIGH-CNJ*(MJ-1)
    KP_NEIGH=K_NEIGH
    IF (IP_NEIGH.LE.CNI.AND.IP_NEIGH.GE.1) THEN
    IF (JP_NEIGH.LE.CNJ.AND.JP_NEIGH.GE.1) THEN
        IF_NEIGH=1
    ENDIF
    ENDIF 
    ELSE
    IP_NEIGH=I_NEIGH-CNI*(MI-1)
    JP_NEIGH=J_NEIGH-CNL*(MJ1-1)
    KP_NEIGH=K_NEIGH-CNL*(MJ2-1)
    IF (IP_NEIGH.LE.CNI.AND.IP_NEIGH.GE.1) THEN
    IF (JP_NEIGH.LE.CNL.AND.JP_NEIGH.GE.1) THEN
    IF (KP_NEIGH.LE.CNL.AND.KP_NEIGH.GE.1) THEN
        IF_NEIGH=1 
    ENDIF
    ENDIF
    ENDIF    
    ENDIF
    
    IF (IF_CY.EQ.0) THEN
    IF (INDEX_S(I1).EQ.1) THEN
        FACE_ID=FACE_ID+1
        IF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.1) THEN
        WRITE(3,231) FACE_ID
        ELSEIF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.0) THEN
        F_FRONT_RECT(IP,JP)=FACE_ID
        ELSEIF (IF_OWN.EQ.0.AND.IF_NEIGH.EQ.1) THEN
        F_BACK_RECT(IP_NEIGH,JP_NEIGH)=-FACE_ID
        ENDIF
    ENDIF
    IF (INDEX_S(I1).EQ.2) THEN
        FACE_ID=FACE_ID+1
        IF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.1) THEN
        WRITE(3,231) FACE_ID
        ELSEIF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.0) THEN
        F_BACK_RECT(IP,JP)=FACE_ID
        ELSEIF (IF_OWN.EQ.0.AND.IF_NEIGH.EQ.1) THEN
        F_FRONT_RECT(IP_NEIGH,JP_NEIGH)=-FACE_ID
        ENDIF
    ENDIF
    IF (INDEX_S(I1).EQ.3) THEN
        FACE_ID=FACE_ID+1
        IF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.1) THEN
        WRITE(3,231) FACE_ID
        ELSEIF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.0) THEN
        F_LEFT_RECT(JP,KP)=FACE_ID
        ELSEIF (IF_OWN.EQ.0.AND.IF_NEIGH.EQ.1) THEN
        F_RIGHT_RECT(JP_NEIGH,KP_NEIGH)=-FACE_ID
        ENDIF
    ENDIF  
    IF (INDEX_S(I1).EQ.4) THEN
        FACE_ID=FACE_ID+1
        IF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.1) THEN
        WRITE(3,231) FACE_ID
        ELSEIF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.0) THEN
        F_RIGHT_RECT(JP,KP)=FACE_ID
        ELSEIF (IF_OWN.EQ.0.AND.IF_NEIGH.EQ.1) THEN
        F_LEFT_RECT(JP_NEIGH,KP_NEIGH)=-FACE_ID
        ENDIF
    ENDIF
    IF (INDEX_S(I1).EQ.5) THEN
        FACE_ID=FACE_ID+1
        IF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.1) THEN
        WRITE(3,231) FACE_ID
        ELSEIF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.0) THEN
        F_UP_RECT(IP,KP)=FACE_ID
        ELSEIF (IF_OWN.EQ.0.AND.IF_NEIGH.EQ.1) THEN
        F_DOWN_RECT(IP_NEIGH,KP_NEIGH)=-FACE_ID
        ENDIF
    ENDIF
    IF (INDEX_S(I1).EQ.6) THEN
        FACE_ID=FACE_ID+1
        IF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.1) THEN
        WRITE(3,231) FACE_ID
        ELSEIF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.0) THEN
        F_DOWN_RECT(IP,KP)=FACE_ID
        ELSEIF (IF_OWN.EQ.0.AND.IF_NEIGH.EQ.1) THEN
        F_UP_RECT(IP_NEIGH,KP_NEIGH)=-FACE_ID
        ENDIF
    ENDIF       
     
    ELSE
    IF (INDEX_S(I1).EQ.2) THEN
        FACE_ID=FACE_ID+1
        IF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.0) THEN
        F_BACK_RECT(IP,JP)=FACE_ID
        ELSEIF (IF_OWN.EQ.0.AND.IF_NEIGH.EQ.1) THEN
        F_DOWN_CY(IP_NEIGH,KP_NEIGH)=-FACE_ID
        ENDIF
    ENDIF
    IF (INDEX_S(I1).EQ.5) THEN
        FACE_ID=FACE_ID+1
        IF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.0) THEN
        F_UP_RECT(IP,KP)=FACE_ID
        ELSEIF (IF_OWN.EQ.0.AND.IF_NEIGH.EQ.1) THEN
        F_DOWN_CY(IP_NEIGH,KP_NEIGH)=-FACE_ID
        ENDIF
    ENDIF
        
    ENDIF
    
    ENDIF
    ENDDO    
    
    ELSE
    !CASE 1-3-2  
    NUM(7)=NCENTER
    CALL GET_CELLID(I,J,K-1,NUM(1),0)
    CALL GET_CELLID(I,J,K+1,NUM(2),0)
    CALL GET_CELLID(I+1,J,K,NUM(3),0)
    NUM(4)=-1
    CALL GET_CELLID(I,1,NL+K,NUM(5),1)
    CALL GET_CELLID(I,J-1,K,NUM(6),0)
    CALL BUBBLE_SORT(NUM,INDEX_S)
    DO I1=1,7
    IF (NUM(I1).GT.NCENTER) THEN  
    IF_NEIGH=0
    CALL GET_INDEX(NUM(I1),I_NEIGH,J_NEIGH,K_NEIGH,IF_CY)
    IF (IF_CY.EQ.1) THEN
    IP_NEIGH=I_NEIGH-CNI*(MI-1)
    JP_NEIGH=J_NEIGH-CNJ*(MJ-1)
    KP_NEIGH=K_NEIGH
    IF (IP_NEIGH.LE.CNI.AND.IP_NEIGH.GE.1) THEN
    IF (JP_NEIGH.LE.CNJ.AND.JP_NEIGH.GE.1) THEN
        IF_NEIGH=1
    ENDIF
    ENDIF 
    ELSE
    IP_NEIGH=I_NEIGH-CNI*(MI-1)
    JP_NEIGH=J_NEIGH-CNL*(MJ1-1)
    KP_NEIGH=K_NEIGH-CNL*(MJ2-1)
    IF (IP_NEIGH.LE.CNI.AND.IP_NEIGH.GE.1) THEN
    IF (JP_NEIGH.LE.CNL.AND.JP_NEIGH.GE.1) THEN
    IF (KP_NEIGH.LE.CNL.AND.KP_NEIGH.GE.1) THEN
        IF_NEIGH=1 
    ENDIF
    ENDIF
    ENDIF    
    ENDIF
    
    IF (IF_CY.EQ.0) THEN
    IF (INDEX_S(I1).EQ.1) THEN
        FACE_ID=FACE_ID+1
        IF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.1) THEN
        WRITE(3,231) FACE_ID
        ELSEIF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.0) THEN
        F_FRONT_RECT(IP,JP)=FACE_ID
        ELSEIF (IF_OWN.EQ.0.AND.IF_NEIGH.EQ.1) THEN
        F_BACK_RECT(IP_NEIGH,JP_NEIGH)=-FACE_ID
        ENDIF
    ENDIF
    IF (INDEX_S(I1).EQ.2) THEN
        FACE_ID=FACE_ID+1
        IF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.1) THEN
        WRITE(3,231) FACE_ID
        ELSEIF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.0) THEN
        F_BACK_RECT(IP,JP)=FACE_ID
        ELSEIF (IF_OWN.EQ.0.AND.IF_NEIGH.EQ.1) THEN
        F_FRONT_RECT(IP_NEIGH,JP_NEIGH)=-FACE_ID
        ENDIF
    ENDIF
    IF (INDEX_S(I1).EQ.3) THEN
        FACE_ID=FACE_ID+1
        IF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.1) THEN
        WRITE(3,231) FACE_ID
        ELSEIF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.0) THEN
        F_LEFT_RECT(JP,KP)=FACE_ID
        ELSEIF (IF_OWN.EQ.0.AND.IF_NEIGH.EQ.1) THEN
        F_RIGHT_RECT(JP_NEIGH,KP_NEIGH)=-FACE_ID
        ENDIF
    ENDIF  
    IF (INDEX_S(I1).EQ.4) THEN
        FACE_ID=FACE_ID+1
        IF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.1) THEN
        WRITE(3,231) FACE_ID
        ELSEIF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.0) THEN
        F_RIGHT_RECT(JP,KP)=FACE_ID
        ELSEIF (IF_OWN.EQ.0.AND.IF_NEIGH.EQ.1) THEN
        F_LEFT_RECT(JP_NEIGH,KP_NEIGH)=-FACE_ID
        ENDIF
    ENDIF
    IF (INDEX_S(I1).EQ.5) THEN
        FACE_ID=FACE_ID+1
        IF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.1) THEN
        WRITE(3,231) FACE_ID
        ELSEIF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.0) THEN
        F_UP_RECT(IP,KP)=FACE_ID
        ELSEIF (IF_OWN.EQ.0.AND.IF_NEIGH.EQ.1) THEN
        F_DOWN_RECT(IP_NEIGH,KP_NEIGH)=-FACE_ID
        ENDIF
    ENDIF
    IF (INDEX_S(I1).EQ.6) THEN
        FACE_ID=FACE_ID+1
        IF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.1) THEN
        WRITE(3,231) FACE_ID
        ELSEIF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.0) THEN
        F_DOWN_RECT(IP,KP)=FACE_ID
        ELSEIF (IF_OWN.EQ.0.AND.IF_NEIGH.EQ.1) THEN
        F_UP_RECT(IP_NEIGH,KP_NEIGH)=-FACE_ID
        ENDIF
    ENDIF       
     
    ELSE
    IF (INDEX_S(I1).EQ.5) THEN
        FACE_ID=FACE_ID+1
        IF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.0) THEN
        F_UP_RECT(IP,KP)=FACE_ID
        ELSEIF (IF_OWN.EQ.0.AND.IF_NEIGH.EQ.1) THEN
        F_DOWN_CY(IP_NEIGH,KP_NEIGH)=-FACE_ID
        ENDIF
    ENDIF
        
    ENDIF
    
    ENDIF
    ENDDO    
    
    ENDIF    
    ELSE
    IF (K.EQ.1) THEN
    !CASE 1-2-1 
    NUM(7)=NCENTER
    CALL GET_CELLID(I,1,J,NUM(1),1)
    CALL GET_CELLID(I,J,K+1,NUM(2),0)
    CALL GET_CELLID(I+1,J,K,NUM(3),0)
    NUM(4)=-1
    CALL GET_CELLID(I,J+1,K,NUM(5),0)
    CALL GET_CELLID(I,J-1,K,NUM(6),0)
    CALL BUBBLE_SORT(NUM,INDEX_S)
    DO I1=1,7
    IF (NUM(I1).GT.NCENTER) THEN  
    IF_NEIGH=0
    CALL GET_INDEX(NUM(I1),I_NEIGH,J_NEIGH,K_NEIGH,IF_CY)
    IF (IF_CY.EQ.1) THEN
    IP_NEIGH=I_NEIGH-CNI*(MI-1)
    JP_NEIGH=J_NEIGH-CNJ*(MJ-1)
    KP_NEIGH=K_NEIGH
    IF (IP_NEIGH.LE.CNI.AND.IP_NEIGH.GE.1) THEN
    IF (JP_NEIGH.LE.CNJ.AND.JP_NEIGH.GE.1) THEN
        IF_NEIGH=1
    ENDIF
    ENDIF 
    ELSE
    IP_NEIGH=I_NEIGH-CNI*(MI-1)
    JP_NEIGH=J_NEIGH-CNL*(MJ1-1)
    KP_NEIGH=K_NEIGH-CNL*(MJ2-1)
    IF (IP_NEIGH.LE.CNI.AND.IP_NEIGH.GE.1) THEN
    IF (JP_NEIGH.LE.CNL.AND.JP_NEIGH.GE.1) THEN
    IF (KP_NEIGH.LE.CNL.AND.KP_NEIGH.GE.1) THEN
        IF_NEIGH=1 
    ENDIF
    ENDIF
    ENDIF    
    ENDIF
    
    IF (IF_CY.EQ.0) THEN
    IF (INDEX_S(I1).EQ.1) THEN
        FACE_ID=FACE_ID+1
        IF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.1) THEN
        WRITE(3,231) FACE_ID
        ELSEIF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.0) THEN
        F_FRONT_RECT(IP,JP)=FACE_ID
        ELSEIF (IF_OWN.EQ.0.AND.IF_NEIGH.EQ.1) THEN
        F_BACK_RECT(IP_NEIGH,JP_NEIGH)=-FACE_ID
        ENDIF
    ENDIF
    IF (INDEX_S(I1).EQ.2) THEN
        FACE_ID=FACE_ID+1
        IF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.1) THEN
        WRITE(3,231) FACE_ID
        ELSEIF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.0) THEN
        F_BACK_RECT(IP,JP)=FACE_ID
        ELSEIF (IF_OWN.EQ.0.AND.IF_NEIGH.EQ.1) THEN
        F_FRONT_RECT(IP_NEIGH,JP_NEIGH)=-FACE_ID
        ENDIF
    ENDIF
    IF (INDEX_S(I1).EQ.3) THEN
        FACE_ID=FACE_ID+1
        IF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.1) THEN
        WRITE(3,231) FACE_ID
        ELSEIF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.0) THEN
        F_LEFT_RECT(JP,KP)=FACE_ID
        ELSEIF (IF_OWN.EQ.0.AND.IF_NEIGH.EQ.1) THEN
        F_RIGHT_RECT(JP_NEIGH,KP_NEIGH)=-FACE_ID
        ENDIF
    ENDIF  
    IF (INDEX_S(I1).EQ.4) THEN
        FACE_ID=FACE_ID+1
        IF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.1) THEN
        WRITE(3,231) FACE_ID
        ELSEIF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.0) THEN
        F_RIGHT_RECT(JP,KP)=FACE_ID
        ELSEIF (IF_OWN.EQ.0.AND.IF_NEIGH.EQ.1) THEN
        F_LEFT_RECT(JP_NEIGH,KP_NEIGH)=-FACE_ID
        ENDIF
    ENDIF
    IF (INDEX_S(I1).EQ.5) THEN
        FACE_ID=FACE_ID+1
        IF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.1) THEN
        WRITE(3,231) FACE_ID
        ELSEIF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.0) THEN
        F_UP_RECT(IP,KP)=FACE_ID
        ELSEIF (IF_OWN.EQ.0.AND.IF_NEIGH.EQ.1) THEN
        F_DOWN_RECT(IP_NEIGH,KP_NEIGH)=-FACE_ID
        ENDIF
    ENDIF
    IF (INDEX_S(I1).EQ.6) THEN
        FACE_ID=FACE_ID+1
        IF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.1) THEN
        WRITE(3,231) FACE_ID
        ELSEIF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.0) THEN
        F_DOWN_RECT(IP,KP)=FACE_ID
        ELSEIF (IF_OWN.EQ.0.AND.IF_NEIGH.EQ.1) THEN
        F_UP_RECT(IP_NEIGH,KP_NEIGH)=-FACE_ID
        ENDIF
    ENDIF       
     
    ELSE
    IF (INDEX_S(I1).EQ.1) THEN
        FACE_ID=FACE_ID+1
        IF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.0) THEN
        F_FRONT_RECT(IP,JP)=FACE_ID
        ELSEIF (IF_OWN.EQ.0.AND.IF_NEIGH.EQ.1) THEN
        F_DOWN_CY(IP_NEIGH,KP_NEIGH)=-FACE_ID
        ENDIF
    ENDIF
        
    ENDIF
    
    ENDIF
    ENDDO    
    
    ELSEIF (K.EQ.NL) THEN
    !CASE 1-2-3 
    NUM(7)=NCENTER
    CALL GET_CELLID(I,J,K-1,NUM(1),0)
    CALL GET_CELLID(I,1,3*NL+1-J,NUM(2),1)
    CALL GET_CELLID(I+1,J,K,NUM(3),0)
    NUM(4)=-1
    CALL GET_CELLID(I,J+1,K,NUM(5),0)
    CALL GET_CELLID(I,J-1,K,NUM(6),0)
    CALL BUBBLE_SORT(NUM,INDEX_S)
    DO I1=1,7
    IF (NUM(I1).GT.NCENTER) THEN  
    IF_NEIGH=0
    CALL GET_INDEX(NUM(I1),I_NEIGH,J_NEIGH,K_NEIGH,IF_CY)
    IF (IF_CY.EQ.1) THEN
    IP_NEIGH=I_NEIGH-CNI*(MI-1)
    JP_NEIGH=J_NEIGH-CNJ*(MJ-1)
    KP_NEIGH=K_NEIGH
    IF (IP_NEIGH.LE.CNI.AND.IP_NEIGH.GE.1) THEN
    IF (JP_NEIGH.LE.CNJ.AND.JP_NEIGH.GE.1) THEN
        IF_NEIGH=1
    ENDIF
    ENDIF 
    ELSE
    IP_NEIGH=I_NEIGH-CNI*(MI-1)
    JP_NEIGH=J_NEIGH-CNL*(MJ1-1)
    KP_NEIGH=K_NEIGH-CNL*(MJ2-1)
    IF (IP_NEIGH.LE.CNI.AND.IP_NEIGH.GE.1) THEN
    IF (JP_NEIGH.LE.CNL.AND.JP_NEIGH.GE.1) THEN
    IF (KP_NEIGH.LE.CNL.AND.KP_NEIGH.GE.1) THEN
        IF_NEIGH=1 
    ENDIF
    ENDIF
    ENDIF    
    ENDIF
    
    IF (IF_CY.EQ.0) THEN
    IF (INDEX_S(I1).EQ.1) THEN
        FACE_ID=FACE_ID+1
        IF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.1) THEN
        WRITE(3,231) FACE_ID
        ELSEIF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.0) THEN
        F_FRONT_RECT(IP,JP)=FACE_ID
        ELSEIF (IF_OWN.EQ.0.AND.IF_NEIGH.EQ.1) THEN
        F_BACK_RECT(IP_NEIGH,JP_NEIGH)=-FACE_ID
        ENDIF
    ENDIF
    IF (INDEX_S(I1).EQ.2) THEN
        FACE_ID=FACE_ID+1
        IF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.1) THEN
        WRITE(3,231) FACE_ID
        ELSEIF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.0) THEN
        F_BACK_RECT(IP,JP)=FACE_ID
        ELSEIF (IF_OWN.EQ.0.AND.IF_NEIGH.EQ.1) THEN
        F_FRONT_RECT(IP_NEIGH,JP_NEIGH)=-FACE_ID
        ENDIF
    ENDIF
    IF (INDEX_S(I1).EQ.3) THEN
        FACE_ID=FACE_ID+1
        IF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.1) THEN
        WRITE(3,231) FACE_ID
        ELSEIF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.0) THEN
        F_LEFT_RECT(JP,KP)=FACE_ID
        ELSEIF (IF_OWN.EQ.0.AND.IF_NEIGH.EQ.1) THEN
        F_RIGHT_RECT(JP_NEIGH,KP_NEIGH)=-FACE_ID
        ENDIF
    ENDIF  
    IF (INDEX_S(I1).EQ.4) THEN
        FACE_ID=FACE_ID+1
        IF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.1) THEN
        WRITE(3,231) FACE_ID
        ELSEIF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.0) THEN
        F_RIGHT_RECT(JP,KP)=FACE_ID
        ELSEIF (IF_OWN.EQ.0.AND.IF_NEIGH.EQ.1) THEN
        F_LEFT_RECT(JP_NEIGH,KP_NEIGH)=-FACE_ID
        ENDIF
    ENDIF
    IF (INDEX_S(I1).EQ.5) THEN
        FACE_ID=FACE_ID+1
        IF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.1) THEN
        WRITE(3,231) FACE_ID
        ELSEIF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.0) THEN
        F_UP_RECT(IP,KP)=FACE_ID
        ELSEIF (IF_OWN.EQ.0.AND.IF_NEIGH.EQ.1) THEN
        F_DOWN_RECT(IP_NEIGH,KP_NEIGH)=-FACE_ID
        ENDIF
    ENDIF
    IF (INDEX_S(I1).EQ.6) THEN
        FACE_ID=FACE_ID+1
        IF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.1) THEN
        WRITE(3,231) FACE_ID
        ELSEIF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.0) THEN
        F_DOWN_RECT(IP,KP)=FACE_ID
        ELSEIF (IF_OWN.EQ.0.AND.IF_NEIGH.EQ.1) THEN
        F_UP_RECT(IP_NEIGH,KP_NEIGH)=-FACE_ID
        ENDIF
    ENDIF       
     
    ELSE
    IF (INDEX_S(I1).EQ.2) THEN
        FACE_ID=FACE_ID+1
        IF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.0) THEN
        F_BACK_RECT(IP,JP)=FACE_ID
        ELSEIF (IF_OWN.EQ.0.AND.IF_NEIGH.EQ.1) THEN
        F_DOWN_CY(IP_NEIGH,KP_NEIGH)=-FACE_ID
        ENDIF
    ENDIF
        
    ENDIF
    
    ENDIF
    ENDDO
    
    ELSE
    !CASE 1-2-2
    NUM(7)=NCENTER
    CALL GET_CELLID(I,J,K-1,NUM(1),0)
    CALL GET_CELLID(I,J,K+1,NUM(2),0)
    CALL GET_CELLID(I+1,J,K,NUM(3),0)
    NUM(4)=-1
    CALL GET_CELLID(I,J+1,K,NUM(5),0)
    CALL GET_CELLID(I,J-1,K,NUM(6),0)
    CALL BUBBLE_SORT(NUM,INDEX_S)
    DO I1=1,7
    IF (NUM(I1).GT.NCENTER) THEN  
    IF_NEIGH=0
    CALL GET_INDEX(NUM(I1),I_NEIGH,J_NEIGH,K_NEIGH,IF_CY)
    IF (IF_CY.EQ.1) THEN
    IP_NEIGH=I_NEIGH-CNI*(MI-1)
    JP_NEIGH=J_NEIGH-CNJ*(MJ-1)
    KP_NEIGH=K_NEIGH
    IF (IP_NEIGH.LE.CNI.AND.IP_NEIGH.GE.1) THEN
    IF (JP_NEIGH.LE.CNJ.AND.JP_NEIGH.GE.1) THEN
        IF_NEIGH=1
    ENDIF
    ENDIF 
    ELSE
    IP_NEIGH=I_NEIGH-CNI*(MI-1)
    JP_NEIGH=J_NEIGH-CNL*(MJ1-1)
    KP_NEIGH=K_NEIGH-CNL*(MJ2-1)
    IF (IP_NEIGH.LE.CNI.AND.IP_NEIGH.GE.1) THEN
    IF (JP_NEIGH.LE.CNL.AND.JP_NEIGH.GE.1) THEN
    IF (KP_NEIGH.LE.CNL.AND.KP_NEIGH.GE.1) THEN
        IF_NEIGH=1 
    ENDIF
    ENDIF
    ENDIF    
    ENDIF
    
    IF (INDEX_S(I1).EQ.1) THEN
        FACE_ID=FACE_ID+1
        IF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.1) THEN
        WRITE(3,231) FACE_ID
        ELSEIF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.0) THEN
        F_FRONT_RECT(IP,JP)=FACE_ID
        ELSEIF (IF_OWN.EQ.0.AND.IF_NEIGH.EQ.1) THEN
        F_BACK_RECT(IP_NEIGH,JP_NEIGH)=-FACE_ID
        ENDIF
    ENDIF
    IF (INDEX_S(I1).EQ.2) THEN
        FACE_ID=FACE_ID+1
        IF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.1) THEN
        WRITE(3,231) FACE_ID
        ELSEIF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.0) THEN
        F_BACK_RECT(IP,JP)=FACE_ID
        ELSEIF (IF_OWN.EQ.0.AND.IF_NEIGH.EQ.1) THEN
        F_FRONT_RECT(IP_NEIGH,JP_NEIGH)=-FACE_ID
        ENDIF
    ENDIF
    IF (INDEX_S(I1).EQ.3) THEN
        FACE_ID=FACE_ID+1
        IF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.1) THEN
        WRITE(3,231) FACE_ID
        ELSEIF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.0) THEN
        F_LEFT_RECT(JP,KP)=FACE_ID
        ELSEIF (IF_OWN.EQ.0.AND.IF_NEIGH.EQ.1) THEN
        F_RIGHT_RECT(JP_NEIGH,KP_NEIGH)=-FACE_ID
        ENDIF
    ENDIF  
    IF (INDEX_S(I1).EQ.4) THEN
        FACE_ID=FACE_ID+1
        IF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.1) THEN
        WRITE(3,231) FACE_ID
        ELSEIF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.0) THEN
        F_RIGHT_RECT(JP,KP)=FACE_ID
        ELSEIF (IF_OWN.EQ.0.AND.IF_NEIGH.EQ.1) THEN
        F_LEFT_RECT(JP_NEIGH,KP_NEIGH)=-FACE_ID
        ENDIF
    ENDIF
    IF (INDEX_S(I1).EQ.5) THEN
        FACE_ID=FACE_ID+1
        IF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.1) THEN
        WRITE(3,231) FACE_ID
        ELSEIF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.0) THEN
        F_UP_RECT(IP,KP)=FACE_ID
        ELSEIF (IF_OWN.EQ.0.AND.IF_NEIGH.EQ.1) THEN
        F_DOWN_RECT(IP_NEIGH,KP_NEIGH)=-FACE_ID
        ENDIF
    ENDIF
    IF (INDEX_S(I1).EQ.6) THEN
        FACE_ID=FACE_ID+1
        IF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.1) THEN
        WRITE(3,231) FACE_ID
        ELSEIF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.0) THEN
        F_DOWN_RECT(IP,KP)=FACE_ID
        ELSEIF (IF_OWN.EQ.0.AND.IF_NEIGH.EQ.1) THEN
        F_UP_RECT(IP_NEIGH,KP_NEIGH)=-FACE_ID
        ENDIF
    ENDIF       
    
    ENDIF
    ENDDO
    
    ENDIF    
    ENDIF    
    ELSEIF (I.EQ.NI) THEN
    IF (J.EQ.1) THEN
    IF (K.EQ.1) THEN
    !CASE 3-1-1
    NUM(7)=NCENTER
    CALL GET_CELLID(I,1,1,NUM(1),1)
    CALL GET_CELLID(I,J,K+1,NUM(2),0)
    NUM(3)=-1
    CALL GET_CELLID(I-1,J,K,NUM(4),0)
    CALL GET_CELLID(I,J+1,K,NUM(5),0)
    CALL GET_CELLID(I,1,NK,NUM(6),1)
    CALL BUBBLE_SORT(NUM,INDEX_S)
    DO I1=1,7
    IF (NUM(I1).GT.NCENTER) THEN  
    IF_NEIGH=0
    CALL GET_INDEX(NUM(I1),I_NEIGH,J_NEIGH,K_NEIGH,IF_CY)
    IF (IF_CY.EQ.1) THEN
    IP_NEIGH=I_NEIGH-CNI*(MI-1)
    JP_NEIGH=J_NEIGH-CNJ*(MJ-1)
    KP_NEIGH=K_NEIGH
    IF (IP_NEIGH.LE.CNI.AND.IP_NEIGH.GE.1) THEN
    IF (JP_NEIGH.LE.CNJ.AND.JP_NEIGH.GE.1) THEN
        IF_NEIGH=1
    ENDIF
    ENDIF 
    ELSE
    IP_NEIGH=I_NEIGH-CNI*(MI-1)
    JP_NEIGH=J_NEIGH-CNL*(MJ1-1)
    KP_NEIGH=K_NEIGH-CNL*(MJ2-1)
    IF (IP_NEIGH.LE.CNI.AND.IP_NEIGH.GE.1) THEN
    IF (JP_NEIGH.LE.CNL.AND.JP_NEIGH.GE.1) THEN
    IF (KP_NEIGH.LE.CNL.AND.KP_NEIGH.GE.1) THEN
        IF_NEIGH=1 
    ENDIF
    ENDIF
    ENDIF    
    ENDIF
    
    IF (IF_CY.EQ.0) THEN
    IF (INDEX_S(I1).EQ.1) THEN
        FACE_ID=FACE_ID+1
        IF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.1) THEN
        WRITE(3,231) FACE_ID
        ELSEIF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.0) THEN
        F_FRONT_RECT(IP,JP)=FACE_ID
        ELSEIF (IF_OWN.EQ.0.AND.IF_NEIGH.EQ.1) THEN
        F_BACK_RECT(IP_NEIGH,JP_NEIGH)=-FACE_ID
        ENDIF
    ENDIF
    IF (INDEX_S(I1).EQ.2) THEN
        FACE_ID=FACE_ID+1
        IF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.1) THEN
        WRITE(3,231) FACE_ID
        ELSEIF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.0) THEN
        F_BACK_RECT(IP,JP)=FACE_ID
        ELSEIF (IF_OWN.EQ.0.AND.IF_NEIGH.EQ.1) THEN
        F_FRONT_RECT(IP_NEIGH,JP_NEIGH)=-FACE_ID
        ENDIF
    ENDIF
    IF (INDEX_S(I1).EQ.3) THEN
        FACE_ID=FACE_ID+1
        IF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.1) THEN
        WRITE(3,231) FACE_ID
        ELSEIF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.0) THEN
        F_LEFT_RECT(JP,KP)=FACE_ID
        ELSEIF (IF_OWN.EQ.0.AND.IF_NEIGH.EQ.1) THEN
        F_RIGHT_RECT(JP_NEIGH,KP_NEIGH)=-FACE_ID
        ENDIF
    ENDIF  
    IF (INDEX_S(I1).EQ.4) THEN
        FACE_ID=FACE_ID+1
        IF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.1) THEN
        WRITE(3,231) FACE_ID
        ELSEIF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.0) THEN
        F_RIGHT_RECT(JP,KP)=FACE_ID
        ELSEIF (IF_OWN.EQ.0.AND.IF_NEIGH.EQ.1) THEN
        F_LEFT_RECT(JP_NEIGH,KP_NEIGH)=-FACE_ID
        ENDIF
    ENDIF
    IF (INDEX_S(I1).EQ.5) THEN
        FACE_ID=FACE_ID+1
        IF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.1) THEN
        WRITE(3,231) FACE_ID
        ELSEIF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.0) THEN
        F_UP_RECT(IP,KP)=FACE_ID
        ELSEIF (IF_OWN.EQ.0.AND.IF_NEIGH.EQ.1) THEN
        F_DOWN_RECT(IP_NEIGH,KP_NEIGH)=-FACE_ID
        ENDIF
    ENDIF
    IF (INDEX_S(I1).EQ.6) THEN
        FACE_ID=FACE_ID+1
        IF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.1) THEN
        WRITE(3,231) FACE_ID
        ELSEIF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.0) THEN
        F_DOWN_RECT(IP,KP)=FACE_ID
        ELSEIF (IF_OWN.EQ.0.AND.IF_NEIGH.EQ.1) THEN
        F_UP_RECT(IP_NEIGH,KP_NEIGH)=-FACE_ID
        ENDIF
    ENDIF       
     
    ELSE
    IF (INDEX_S(I1).EQ.1) THEN
        FACE_ID=FACE_ID+1
        IF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.0) THEN
        F_FRONT_RECT(IP,JP)=FACE_ID
        ELSEIF (IF_OWN.EQ.0.AND.IF_NEIGH.EQ.1) THEN
        F_DOWN_CY(IP_NEIGH,KP_NEIGH)=-FACE_ID
        ENDIF
    ENDIF
    IF (INDEX_S(I1).EQ.6) THEN
        FACE_ID=FACE_ID+1
        IF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.0) THEN
        F_DOWN_RECT(IP,KP)=FACE_ID
        ELSEIF (IF_OWN.EQ.0.AND.IF_NEIGH.EQ.1) THEN
        F_DOWN_CY(IP_NEIGH,KP_NEIGH)=-FACE_ID
        ENDIF
    ENDIF
        
    ENDIF
    
    ENDIF
    ENDDO    
    
    ELSEIF (K.EQ.NL) THEN
    !CASE 3-1-3 
    NUM(7)=NCENTER
    CALL GET_CELLID(I,J,K-1,NUM(1),0)
    CALL GET_CELLID(I,1,3*NL,NUM(2),1)
    NUM(3)=-1
    CALL GET_CELLID(I-1,J,K,NUM(4),0)
    CALL GET_CELLID(I,J+1,K,NUM(5),0)
    CALL GET_CELLID(I,1,3*NL+1,NUM(6),1)
    CALL BUBBLE_SORT(NUM,INDEX_S)
    DO I1=1,7
    IF (NUM(I1).GT.NCENTER) THEN  
    IF_NEIGH=0
    CALL GET_INDEX(NUM(I1),I_NEIGH,J_NEIGH,K_NEIGH,IF_CY)
    IF (IF_CY.EQ.1) THEN
    IP_NEIGH=I_NEIGH-CNI*(MI-1)
    JP_NEIGH=J_NEIGH-CNJ*(MJ-1)
    KP_NEIGH=K_NEIGH
    IF (IP_NEIGH.LE.CNI.AND.IP_NEIGH.GE.1) THEN
    IF (JP_NEIGH.LE.CNJ.AND.JP_NEIGH.GE.1) THEN
        IF_NEIGH=1
    ENDIF
    ENDIF 
    ELSE
    IP_NEIGH=I_NEIGH-CNI*(MI-1)
    JP_NEIGH=J_NEIGH-CNL*(MJ1-1)
    KP_NEIGH=K_NEIGH-CNL*(MJ2-1)
    IF (IP_NEIGH.LE.CNI.AND.IP_NEIGH.GE.1) THEN
    IF (JP_NEIGH.LE.CNL.AND.JP_NEIGH.GE.1) THEN
    IF (KP_NEIGH.LE.CNL.AND.KP_NEIGH.GE.1) THEN
        IF_NEIGH=1 
    ENDIF
    ENDIF
    ENDIF    
    ENDIF
    
    IF (IF_CY.EQ.0) THEN
    IF (INDEX_S(I1).EQ.1) THEN
        FACE_ID=FACE_ID+1
        IF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.1) THEN
        WRITE(3,231) FACE_ID
        ELSEIF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.0) THEN
        F_FRONT_RECT(IP,JP)=FACE_ID
        ELSEIF (IF_OWN.EQ.0.AND.IF_NEIGH.EQ.1) THEN
        F_BACK_RECT(IP_NEIGH,JP_NEIGH)=-FACE_ID
        ENDIF
    ENDIF
    IF (INDEX_S(I1).EQ.2) THEN
        FACE_ID=FACE_ID+1
        IF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.1) THEN
        WRITE(3,231) FACE_ID
        ELSEIF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.0) THEN
        F_BACK_RECT(IP,JP)=FACE_ID
        ELSEIF (IF_OWN.EQ.0.AND.IF_NEIGH.EQ.1) THEN
        F_FRONT_RECT(IP_NEIGH,JP_NEIGH)=-FACE_ID
        ENDIF
    ENDIF
    IF (INDEX_S(I1).EQ.3) THEN
        FACE_ID=FACE_ID+1
        IF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.1) THEN
        WRITE(3,231) FACE_ID
        ELSEIF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.0) THEN
        F_LEFT_RECT(JP,KP)=FACE_ID
        ELSEIF (IF_OWN.EQ.0.AND.IF_NEIGH.EQ.1) THEN
        F_RIGHT_RECT(JP_NEIGH,KP_NEIGH)=-FACE_ID
        ENDIF
    ENDIF  
    IF (INDEX_S(I1).EQ.4) THEN
        FACE_ID=FACE_ID+1
        IF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.1) THEN
        WRITE(3,231) FACE_ID
        ELSEIF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.0) THEN
        F_RIGHT_RECT(JP,KP)=FACE_ID
        ELSEIF (IF_OWN.EQ.0.AND.IF_NEIGH.EQ.1) THEN
        F_LEFT_RECT(JP_NEIGH,KP_NEIGH)=-FACE_ID
        ENDIF
    ENDIF
    IF (INDEX_S(I1).EQ.5) THEN
        FACE_ID=FACE_ID+1
        IF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.1) THEN
        WRITE(3,231) FACE_ID
        ELSEIF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.0) THEN
        F_UP_RECT(IP,KP)=FACE_ID
        ELSEIF (IF_OWN.EQ.0.AND.IF_NEIGH.EQ.1) THEN
        F_DOWN_RECT(IP_NEIGH,KP_NEIGH)=-FACE_ID
        ENDIF
    ENDIF
    IF (INDEX_S(I1).EQ.6) THEN
        FACE_ID=FACE_ID+1
        IF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.1) THEN
        WRITE(3,231) FACE_ID
        ELSEIF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.0) THEN
        F_DOWN_RECT(IP,KP)=FACE_ID
        ELSEIF (IF_OWN.EQ.0.AND.IF_NEIGH.EQ.1) THEN
        F_UP_RECT(IP_NEIGH,KP_NEIGH)=-FACE_ID
        ENDIF
    ENDIF       
     
    ELSE
    IF (INDEX_S(I1).EQ.2) THEN
        FACE_ID=FACE_ID+1
        IF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.0) THEN
        F_BACK_RECT(IP,JP)=FACE_ID
        ELSEIF (IF_OWN.EQ.0.AND.IF_NEIGH.EQ.1) THEN
        F_DOWN_CY(IP_NEIGH,KP_NEIGH)=-FACE_ID
        ENDIF
    ENDIF
    IF (INDEX_S(I1).EQ.6) THEN
        FACE_ID=FACE_ID+1
        IF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.0) THEN
        F_DOWN_RECT(IP,KP)=FACE_ID
        ELSEIF (IF_OWN.EQ.0.AND.IF_NEIGH.EQ.1) THEN
        F_DOWN_CY(IP_NEIGH,KP_NEIGH)=-FACE_ID
        ENDIF
    ENDIF
        
    ENDIF
    
    ENDIF
    ENDDO    
    
    ELSE
    !CASE 3-1-2 
    NUM(7)=NCENTER
    CALL GET_CELLID(I,J,K-1,NUM(1),0)
    CALL GET_CELLID(I,J,K+1,NUM(2),0)
    NUM(3)=-1
    CALL GET_CELLID(I-1,J,K,NUM(4),0)
    CALL GET_CELLID(I,J+1,K,NUM(5),0)
    CALL GET_CELLID(I,1,NK+1-K,NUM(6),1)
    CALL BUBBLE_SORT(NUM,INDEX_S)
    DO I1=1,7
    IF (NUM(I1).GT.NCENTER) THEN  
    IF_NEIGH=0
    CALL GET_INDEX(NUM(I1),I_NEIGH,J_NEIGH,K_NEIGH,IF_CY)
    IF (IF_CY.EQ.1) THEN
    IP_NEIGH=I_NEIGH-CNI*(MI-1)
    JP_NEIGH=J_NEIGH-CNJ*(MJ-1)
    KP_NEIGH=K_NEIGH
    IF (IP_NEIGH.LE.CNI.AND.IP_NEIGH.GE.1) THEN
    IF (JP_NEIGH.LE.CNJ.AND.JP_NEIGH.GE.1) THEN
        IF_NEIGH=1
    ENDIF
    ENDIF 
    ELSE
    IP_NEIGH=I_NEIGH-CNI*(MI-1)
    JP_NEIGH=J_NEIGH-CNL*(MJ1-1)
    KP_NEIGH=K_NEIGH-CNL*(MJ2-1)
    IF (IP_NEIGH.LE.CNI.AND.IP_NEIGH.GE.1) THEN
    IF (JP_NEIGH.LE.CNL.AND.JP_NEIGH.GE.1) THEN
    IF (KP_NEIGH.LE.CNL.AND.KP_NEIGH.GE.1) THEN
        IF_NEIGH=1 
    ENDIF
    ENDIF
    ENDIF    
    ENDIF
    
    IF (IF_CY.EQ.0) THEN
    IF (INDEX_S(I1).EQ.1) THEN
        FACE_ID=FACE_ID+1
        IF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.1) THEN
        WRITE(3,231) FACE_ID
        ELSEIF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.0) THEN
        F_FRONT_RECT(IP,JP)=FACE_ID
        ELSEIF (IF_OWN.EQ.0.AND.IF_NEIGH.EQ.1) THEN
        F_BACK_RECT(IP_NEIGH,JP_NEIGH)=-FACE_ID
        ENDIF
    ENDIF
    IF (INDEX_S(I1).EQ.2) THEN
        FACE_ID=FACE_ID+1
        IF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.1) THEN
        WRITE(3,231) FACE_ID
        ELSEIF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.0) THEN
        F_BACK_RECT(IP,JP)=FACE_ID
        ELSEIF (IF_OWN.EQ.0.AND.IF_NEIGH.EQ.1) THEN
        F_FRONT_RECT(IP_NEIGH,JP_NEIGH)=-FACE_ID
        ENDIF
    ENDIF
    IF (INDEX_S(I1).EQ.3) THEN
        FACE_ID=FACE_ID+1
        IF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.1) THEN
        WRITE(3,231) FACE_ID
        ELSEIF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.0) THEN
        F_LEFT_RECT(JP,KP)=FACE_ID
        ELSEIF (IF_OWN.EQ.0.AND.IF_NEIGH.EQ.1) THEN
        F_RIGHT_RECT(JP_NEIGH,KP_NEIGH)=-FACE_ID
        ENDIF
    ENDIF  
    IF (INDEX_S(I1).EQ.4) THEN
        FACE_ID=FACE_ID+1
        IF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.1) THEN
        WRITE(3,231) FACE_ID
        ELSEIF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.0) THEN
        F_RIGHT_RECT(JP,KP)=FACE_ID
        ELSEIF (IF_OWN.EQ.0.AND.IF_NEIGH.EQ.1) THEN
        F_LEFT_RECT(JP_NEIGH,KP_NEIGH)=-FACE_ID
        ENDIF
    ENDIF
    IF (INDEX_S(I1).EQ.5) THEN
        FACE_ID=FACE_ID+1
        IF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.1) THEN
        WRITE(3,231) FACE_ID
        ELSEIF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.0) THEN
        F_UP_RECT(IP,KP)=FACE_ID
        ELSEIF (IF_OWN.EQ.0.AND.IF_NEIGH.EQ.1) THEN
        F_DOWN_RECT(IP_NEIGH,KP_NEIGH)=-FACE_ID
        ENDIF
    ENDIF
    IF (INDEX_S(I1).EQ.6) THEN
        FACE_ID=FACE_ID+1
        IF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.1) THEN
        WRITE(3,231) FACE_ID
        ELSEIF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.0) THEN
        F_DOWN_RECT(IP,KP)=FACE_ID
        ELSEIF (IF_OWN.EQ.0.AND.IF_NEIGH.EQ.1) THEN
        F_UP_RECT(IP_NEIGH,KP_NEIGH)=-FACE_ID
        ENDIF
    ENDIF       
     
    ELSE
    IF (INDEX_S(I1).EQ.6) THEN
        FACE_ID=FACE_ID+1
        IF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.0) THEN
        F_DOWN_RECT(IP,KP)=FACE_ID
        ELSEIF (IF_OWN.EQ.0.AND.IF_NEIGH.EQ.1) THEN
        F_DOWN_CY(IP_NEIGH,KP_NEIGH)=-FACE_ID
        ENDIF
    ENDIF
        
    ENDIF
    
    ENDIF
    ENDDO    
    
    ENDIF    
    ELSEIF (J.EQ.NL) THEN
    IF (K.EQ.1) THEN
    !CASE 3-3-1 
    NUM(7)=NCENTER
    CALL GET_CELLID(I,1,NL,NUM(1),1)
    CALL GET_CELLID(I,J,K+1,NUM(2),0)
    NUM(3)=-1
    CALL GET_CELLID(I-1,J,K,NUM(4),0)
    CALL GET_CELLID(I,1,NL+1,NUM(5),1)
    CALL GET_CELLID(I,J-1,K,NUM(6),0)
    CALL BUBBLE_SORT(NUM,INDEX_S)
    DO I1=1,7
    IF (NUM(I1).GT.NCENTER) THEN  
    IF_NEIGH=0
    CALL GET_INDEX(NUM(I1),I_NEIGH,J_NEIGH,K_NEIGH,IF_CY)
    IF (IF_CY.EQ.1) THEN
    IP_NEIGH=I_NEIGH-CNI*(MI-1)
    JP_NEIGH=J_NEIGH-CNJ*(MJ-1)
    KP_NEIGH=K_NEIGH
    IF (IP_NEIGH.LE.CNI.AND.IP_NEIGH.GE.1) THEN
    IF (JP_NEIGH.LE.CNJ.AND.JP_NEIGH.GE.1) THEN
        IF_NEIGH=1
    ENDIF
    ENDIF 
    ELSE
    IP_NEIGH=I_NEIGH-CNI*(MI-1)
    JP_NEIGH=J_NEIGH-CNL*(MJ1-1)
    KP_NEIGH=K_NEIGH-CNL*(MJ2-1)
    IF (IP_NEIGH.LE.CNI.AND.IP_NEIGH.GE.1) THEN
    IF (JP_NEIGH.LE.CNL.AND.JP_NEIGH.GE.1) THEN
    IF (KP_NEIGH.LE.CNL.AND.KP_NEIGH.GE.1) THEN
        IF_NEIGH=1 
    ENDIF
    ENDIF
    ENDIF    
    ENDIF
    
    IF (IF_CY.EQ.0) THEN
    IF (INDEX_S(I1).EQ.1) THEN
        FACE_ID=FACE_ID+1
        IF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.1) THEN
        WRITE(3,231) FACE_ID
        ELSEIF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.0) THEN
        F_FRONT_RECT(IP,JP)=FACE_ID
        ELSEIF (IF_OWN.EQ.0.AND.IF_NEIGH.EQ.1) THEN
        F_BACK_RECT(IP_NEIGH,JP_NEIGH)=-FACE_ID
        ENDIF
    ENDIF
    IF (INDEX_S(I1).EQ.2) THEN
        FACE_ID=FACE_ID+1
        IF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.1) THEN
        WRITE(3,231) FACE_ID
        ELSEIF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.0) THEN
        F_BACK_RECT(IP,JP)=FACE_ID
        ELSEIF (IF_OWN.EQ.0.AND.IF_NEIGH.EQ.1) THEN
        F_FRONT_RECT(IP_NEIGH,JP_NEIGH)=-FACE_ID
        ENDIF
    ENDIF
    IF (INDEX_S(I1).EQ.3) THEN
        FACE_ID=FACE_ID+1
        IF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.1) THEN
        WRITE(3,231) FACE_ID
        ELSEIF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.0) THEN
        F_LEFT_RECT(JP,KP)=FACE_ID
        ELSEIF (IF_OWN.EQ.0.AND.IF_NEIGH.EQ.1) THEN
        F_RIGHT_RECT(JP_NEIGH,KP_NEIGH)=-FACE_ID
        ENDIF
    ENDIF  
    IF (INDEX_S(I1).EQ.4) THEN
        FACE_ID=FACE_ID+1
        IF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.1) THEN
        WRITE(3,231) FACE_ID
        ELSEIF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.0) THEN
        F_RIGHT_RECT(JP,KP)=FACE_ID
        ELSEIF (IF_OWN.EQ.0.AND.IF_NEIGH.EQ.1) THEN
        F_LEFT_RECT(JP_NEIGH,KP_NEIGH)=-FACE_ID
        ENDIF
    ENDIF
    IF (INDEX_S(I1).EQ.5) THEN
        FACE_ID=FACE_ID+1
        IF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.1) THEN
        WRITE(3,231) FACE_ID
        ELSEIF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.0) THEN
        F_UP_RECT(IP,KP)=FACE_ID
        ELSEIF (IF_OWN.EQ.0.AND.IF_NEIGH.EQ.1) THEN
        F_DOWN_RECT(IP_NEIGH,KP_NEIGH)=-FACE_ID
        ENDIF
    ENDIF
    IF (INDEX_S(I1).EQ.6) THEN
        FACE_ID=FACE_ID+1
        IF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.1) THEN
        WRITE(3,231) FACE_ID
        ELSEIF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.0) THEN
        F_DOWN_RECT(IP,KP)=FACE_ID
        ELSEIF (IF_OWN.EQ.0.AND.IF_NEIGH.EQ.1) THEN
        F_UP_RECT(IP_NEIGH,KP_NEIGH)=-FACE_ID
        ENDIF
    ENDIF       
     
    ELSE
    IF (INDEX_S(I1).EQ.1) THEN
        FACE_ID=FACE_ID+1
        IF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.0) THEN
        F_FRONT_RECT(IP,JP)=FACE_ID
        ELSEIF (IF_OWN.EQ.0.AND.IF_NEIGH.EQ.1) THEN
        F_DOWN_CY(IP_NEIGH,KP_NEIGH)=-FACE_ID
        ENDIF
    ENDIF
    IF (INDEX_S(I1).EQ.5) THEN
        FACE_ID=FACE_ID+1
        IF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.0) THEN
        F_UP_RECT(IP,KP)=FACE_ID
        ELSEIF (IF_OWN.EQ.0.AND.IF_NEIGH.EQ.1) THEN
        F_DOWN_CY(IP_NEIGH,KP_NEIGH)=-FACE_ID
        ENDIF
    ENDIF
        
    ENDIF
    
    ENDIF
    ENDDO    
    
    ELSEIF (K.EQ.NL) THEN
    !CASE 3-3-3 
    NUM(7)=NCENTER
    CALL GET_CELLID(I,J,K-1,NUM(1),0)
    CALL GET_CELLID(I,1,2*NL+1,NUM(2),1)
    NUM(3)=-1
    CALL GET_CELLID(I-1,J,K,NUM(4),0)
    CALL GET_CELLID(I,1,2*NL,NUM(5),1)
    CALL GET_CELLID(I,J-1,K,NUM(6),0)
    CALL BUBBLE_SORT(NUM,INDEX_S)
    DO I1=1,7
    IF (NUM(I1).GT.NCENTER) THEN  
    IF_NEIGH=0
    CALL GET_INDEX(NUM(I1),I_NEIGH,J_NEIGH,K_NEIGH,IF_CY)
    IF (IF_CY.EQ.1) THEN
    IP_NEIGH=I_NEIGH-CNI*(MI-1)
    JP_NEIGH=J_NEIGH-CNJ*(MJ-1)
    KP_NEIGH=K_NEIGH
    IF (IP_NEIGH.LE.CNI.AND.IP_NEIGH.GE.1) THEN
    IF (JP_NEIGH.LE.CNJ.AND.JP_NEIGH.GE.1) THEN
        IF_NEIGH=1
    ENDIF
    ENDIF 
    ELSE
    IP_NEIGH=I_NEIGH-CNI*(MI-1)
    JP_NEIGH=J_NEIGH-CNL*(MJ1-1)
    KP_NEIGH=K_NEIGH-CNL*(MJ2-1)
    IF (IP_NEIGH.LE.CNI.AND.IP_NEIGH.GE.1) THEN
    IF (JP_NEIGH.LE.CNL.AND.JP_NEIGH.GE.1) THEN
    IF (KP_NEIGH.LE.CNL.AND.KP_NEIGH.GE.1) THEN
        IF_NEIGH=1 
    ENDIF
    ENDIF
    ENDIF    
    ENDIF
    
    IF (IF_CY.EQ.0) THEN
    IF (INDEX_S(I1).EQ.1) THEN
        FACE_ID=FACE_ID+1
        IF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.1) THEN
        WRITE(3,231) FACE_ID
        ELSEIF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.0) THEN
        F_FRONT_RECT(IP,JP)=FACE_ID
        ELSEIF (IF_OWN.EQ.0.AND.IF_NEIGH.EQ.1) THEN
        F_BACK_RECT(IP_NEIGH,JP_NEIGH)=-FACE_ID
        ENDIF
    ENDIF
    IF (INDEX_S(I1).EQ.2) THEN
        FACE_ID=FACE_ID+1
        IF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.1) THEN
        WRITE(3,231) FACE_ID
        ELSEIF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.0) THEN
        F_BACK_RECT(IP,JP)=FACE_ID
        ELSEIF (IF_OWN.EQ.0.AND.IF_NEIGH.EQ.1) THEN
        F_FRONT_RECT(IP_NEIGH,JP_NEIGH)=-FACE_ID
        ENDIF
    ENDIF
    IF (INDEX_S(I1).EQ.3) THEN
        FACE_ID=FACE_ID+1
        IF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.1) THEN
        WRITE(3,231) FACE_ID
        ELSEIF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.0) THEN
        F_LEFT_RECT(JP,KP)=FACE_ID
        ELSEIF (IF_OWN.EQ.0.AND.IF_NEIGH.EQ.1) THEN
        F_RIGHT_RECT(JP_NEIGH,KP_NEIGH)=-FACE_ID
        ENDIF
    ENDIF  
    IF (INDEX_S(I1).EQ.4) THEN
        FACE_ID=FACE_ID+1
        IF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.1) THEN
        WRITE(3,231) FACE_ID
        ELSEIF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.0) THEN
        F_RIGHT_RECT(JP,KP)=FACE_ID
        ELSEIF (IF_OWN.EQ.0.AND.IF_NEIGH.EQ.1) THEN
        F_LEFT_RECT(JP_NEIGH,KP_NEIGH)=-FACE_ID
        ENDIF
    ENDIF
    IF (INDEX_S(I1).EQ.5) THEN
        FACE_ID=FACE_ID+1
        IF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.1) THEN
        WRITE(3,231) FACE_ID
        ELSEIF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.0) THEN
        F_UP_RECT(IP,KP)=FACE_ID
        ELSEIF (IF_OWN.EQ.0.AND.IF_NEIGH.EQ.1) THEN
        F_DOWN_RECT(IP_NEIGH,KP_NEIGH)=-FACE_ID
        ENDIF
    ENDIF
    IF (INDEX_S(I1).EQ.6) THEN
        FACE_ID=FACE_ID+1
        IF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.1) THEN
        WRITE(3,231) FACE_ID
        ELSEIF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.0) THEN
        F_DOWN_RECT(IP,KP)=FACE_ID
        ELSEIF (IF_OWN.EQ.0.AND.IF_NEIGH.EQ.1) THEN
        F_UP_RECT(IP_NEIGH,KP_NEIGH)=-FACE_ID
        ENDIF
    ENDIF       
     
    ELSE
    IF (INDEX_S(I1).EQ.2) THEN
        FACE_ID=FACE_ID+1
        IF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.0) THEN
        F_BACK_RECT(IP,JP)=FACE_ID
        ELSEIF (IF_OWN.EQ.0.AND.IF_NEIGH.EQ.1) THEN
        F_DOWN_CY(IP_NEIGH,KP_NEIGH)=-FACE_ID
        ENDIF
    ENDIF
    IF (INDEX_S(I1).EQ.5) THEN
        FACE_ID=FACE_ID+1
        IF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.0) THEN
        F_UP_RECT(IP,KP)=FACE_ID
        ELSEIF (IF_OWN.EQ.0.AND.IF_NEIGH.EQ.1) THEN
        F_DOWN_CY(IP_NEIGH,KP_NEIGH)=-FACE_ID
        ENDIF
    ENDIF
        
    ENDIF
    
    ENDIF
    ENDDO    
    
    ELSE
    !CASE 3-3-2
    NUM(7)=NCENTER
    CALL GET_CELLID(I,J,K-1,NUM(1),0)
    CALL GET_CELLID(I,J,K+1,NUM(2),0)
    NUM(3)=-1
    CALL GET_CELLID(I-1,J,K,NUM(4),0)
    CALL GET_CELLID(I,1,NL+K,NUM(5),1)
    CALL GET_CELLID(I,J-1,K,NUM(6),0)
    CALL BUBBLE_SORT(NUM,INDEX_S)
    DO I1=1,7
    IF (NUM(I1).GT.NCENTER) THEN  
    IF_NEIGH=0
    CALL GET_INDEX(NUM(I1),I_NEIGH,J_NEIGH,K_NEIGH,IF_CY)
    IF (IF_CY.EQ.1) THEN
    IP_NEIGH=I_NEIGH-CNI*(MI-1)
    JP_NEIGH=J_NEIGH-CNJ*(MJ-1)
    KP_NEIGH=K_NEIGH
    IF (IP_NEIGH.LE.CNI.AND.IP_NEIGH.GE.1) THEN
    IF (JP_NEIGH.LE.CNJ.AND.JP_NEIGH.GE.1) THEN
        IF_NEIGH=1
    ENDIF
    ENDIF 
    ELSE
    IP_NEIGH=I_NEIGH-CNI*(MI-1)
    JP_NEIGH=J_NEIGH-CNL*(MJ1-1)
    KP_NEIGH=K_NEIGH-CNL*(MJ2-1)
    IF (IP_NEIGH.LE.CNI.AND.IP_NEIGH.GE.1) THEN
    IF (JP_NEIGH.LE.CNL.AND.JP_NEIGH.GE.1) THEN
    IF (KP_NEIGH.LE.CNL.AND.KP_NEIGH.GE.1) THEN
        IF_NEIGH=1 
    ENDIF
    ENDIF
    ENDIF    
    ENDIF
    
    IF (IF_CY.EQ.0) THEN
    IF (INDEX_S(I1).EQ.1) THEN
        FACE_ID=FACE_ID+1
        IF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.1) THEN
        WRITE(3,231) FACE_ID
        ELSEIF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.0) THEN
        F_FRONT_RECT(IP,JP)=FACE_ID
        ELSEIF (IF_OWN.EQ.0.AND.IF_NEIGH.EQ.1) THEN
        F_BACK_RECT(IP_NEIGH,JP_NEIGH)=-FACE_ID
        ENDIF
    ENDIF
    IF (INDEX_S(I1).EQ.2) THEN
        FACE_ID=FACE_ID+1
        IF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.1) THEN
        WRITE(3,231) FACE_ID
        ELSEIF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.0) THEN
        F_BACK_RECT(IP,JP)=FACE_ID
        ELSEIF (IF_OWN.EQ.0.AND.IF_NEIGH.EQ.1) THEN
        F_FRONT_RECT(IP_NEIGH,JP_NEIGH)=-FACE_ID
        ENDIF
    ENDIF
    IF (INDEX_S(I1).EQ.3) THEN
        FACE_ID=FACE_ID+1
        IF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.1) THEN
        WRITE(3,231) FACE_ID
        ELSEIF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.0) THEN
        F_LEFT_RECT(JP,KP)=FACE_ID
        ELSEIF (IF_OWN.EQ.0.AND.IF_NEIGH.EQ.1) THEN
        F_RIGHT_RECT(JP_NEIGH,KP_NEIGH)=-FACE_ID
        ENDIF
    ENDIF  
    IF (INDEX_S(I1).EQ.4) THEN
        FACE_ID=FACE_ID+1
        IF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.1) THEN
        WRITE(3,231) FACE_ID
        ELSEIF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.0) THEN
        F_RIGHT_RECT(JP,KP)=FACE_ID
        ELSEIF (IF_OWN.EQ.0.AND.IF_NEIGH.EQ.1) THEN
        F_LEFT_RECT(JP_NEIGH,KP_NEIGH)=-FACE_ID
        ENDIF
    ENDIF
    IF (INDEX_S(I1).EQ.5) THEN
        FACE_ID=FACE_ID+1
        IF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.1) THEN
        WRITE(3,231) FACE_ID
        ELSEIF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.0) THEN
        F_UP_RECT(IP,KP)=FACE_ID
        ELSEIF (IF_OWN.EQ.0.AND.IF_NEIGH.EQ.1) THEN
        F_DOWN_RECT(IP_NEIGH,KP_NEIGH)=-FACE_ID
        ENDIF
    ENDIF
    IF (INDEX_S(I1).EQ.6) THEN
        FACE_ID=FACE_ID+1
        IF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.1) THEN
        WRITE(3,231) FACE_ID
        ELSEIF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.0) THEN
        F_DOWN_RECT(IP,KP)=FACE_ID
        ELSEIF (IF_OWN.EQ.0.AND.IF_NEIGH.EQ.1) THEN
        F_UP_RECT(IP_NEIGH,KP_NEIGH)=-FACE_ID
        ENDIF
    ENDIF       
     
    ELSE
    IF (INDEX_S(I1).EQ.5) THEN
        FACE_ID=FACE_ID+1
        IF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.0) THEN
        F_UP_RECT(IP,KP)=FACE_ID
        ELSEIF (IF_OWN.EQ.0.AND.IF_NEIGH.EQ.1) THEN
        F_DOWN_CY(IP_NEIGH,KP_NEIGH)=-FACE_ID
        ENDIF
    ENDIF
        
    ENDIF
    
    ENDIF
    ENDDO    
    
    ENDIF    
    ELSE
    IF (K.EQ.1) THEN
    !CASE 3-2-1
    NUM(7)=NCENTER
    CALL GET_CELLID(I,1,J,NUM(1),1)
    CALL GET_CELLID(I,J,K+1,NUM(2),0)
    NUM(3)=-1
    CALL GET_CELLID(I-1,J,K,NUM(4),0)
    CALL GET_CELLID(I,J+1,K,NUM(5),0)
    CALL GET_CELLID(I,J-1,K,NUM(6),0)
    CALL BUBBLE_SORT(NUM,INDEX_S)
    DO I1=1,7
    IF (NUM(I1).GT.NCENTER) THEN  
    IF_NEIGH=0
    CALL GET_INDEX(NUM(I1),I_NEIGH,J_NEIGH,K_NEIGH,IF_CY)
    IF (IF_CY.EQ.1) THEN
    IP_NEIGH=I_NEIGH-CNI*(MI-1)
    JP_NEIGH=J_NEIGH-CNJ*(MJ-1)
    KP_NEIGH=K_NEIGH
    IF (IP_NEIGH.LE.CNI.AND.IP_NEIGH.GE.1) THEN
    IF (JP_NEIGH.LE.CNJ.AND.JP_NEIGH.GE.1) THEN
        IF_NEIGH=1
    ENDIF
    ENDIF 
    ELSE
    IP_NEIGH=I_NEIGH-CNI*(MI-1)
    JP_NEIGH=J_NEIGH-CNL*(MJ1-1)
    KP_NEIGH=K_NEIGH-CNL*(MJ2-1)
    IF (IP_NEIGH.LE.CNI.AND.IP_NEIGH.GE.1) THEN
    IF (JP_NEIGH.LE.CNL.AND.JP_NEIGH.GE.1) THEN
    IF (KP_NEIGH.LE.CNL.AND.KP_NEIGH.GE.1) THEN
        IF_NEIGH=1 
    ENDIF
    ENDIF
    ENDIF    
    ENDIF
    
    IF (IF_CY.EQ.0) THEN
    IF (INDEX_S(I1).EQ.1) THEN
        FACE_ID=FACE_ID+1
        IF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.1) THEN
        WRITE(3,231) FACE_ID
        ELSEIF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.0) THEN
        F_FRONT_RECT(IP,JP)=FACE_ID
        ELSEIF (IF_OWN.EQ.0.AND.IF_NEIGH.EQ.1) THEN
        F_BACK_RECT(IP_NEIGH,JP_NEIGH)=-FACE_ID
        ENDIF
    ENDIF
    IF (INDEX_S(I1).EQ.2) THEN
        FACE_ID=FACE_ID+1
        IF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.1) THEN
        WRITE(3,231) FACE_ID
        ELSEIF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.0) THEN
        F_BACK_RECT(IP,JP)=FACE_ID
        ELSEIF (IF_OWN.EQ.0.AND.IF_NEIGH.EQ.1) THEN
        F_FRONT_RECT(IP_NEIGH,JP_NEIGH)=-FACE_ID
        ENDIF
    ENDIF
    IF (INDEX_S(I1).EQ.3) THEN
        FACE_ID=FACE_ID+1
        IF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.1) THEN
        WRITE(3,231) FACE_ID
        ELSEIF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.0) THEN
        F_LEFT_RECT(JP,KP)=FACE_ID
        ELSEIF (IF_OWN.EQ.0.AND.IF_NEIGH.EQ.1) THEN
        F_RIGHT_RECT(JP_NEIGH,KP_NEIGH)=-FACE_ID
        ENDIF
    ENDIF  
    IF (INDEX_S(I1).EQ.4) THEN
        FACE_ID=FACE_ID+1
        IF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.1) THEN
        WRITE(3,231) FACE_ID
        ELSEIF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.0) THEN
        F_RIGHT_RECT(JP,KP)=FACE_ID
        ELSEIF (IF_OWN.EQ.0.AND.IF_NEIGH.EQ.1) THEN
        F_LEFT_RECT(JP_NEIGH,KP_NEIGH)=-FACE_ID
        ENDIF
    ENDIF
    IF (INDEX_S(I1).EQ.5) THEN
        FACE_ID=FACE_ID+1
        IF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.1) THEN
        WRITE(3,231) FACE_ID
        ELSEIF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.0) THEN
        F_UP_RECT(IP,KP)=FACE_ID
        ELSEIF (IF_OWN.EQ.0.AND.IF_NEIGH.EQ.1) THEN
        F_DOWN_RECT(IP_NEIGH,KP_NEIGH)=-FACE_ID
        ENDIF
    ENDIF
    IF (INDEX_S(I1).EQ.6) THEN
        FACE_ID=FACE_ID+1
        IF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.1) THEN
        WRITE(3,231) FACE_ID
        ELSEIF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.0) THEN
        F_DOWN_RECT(IP,KP)=FACE_ID
        ELSEIF (IF_OWN.EQ.0.AND.IF_NEIGH.EQ.1) THEN
        F_UP_RECT(IP_NEIGH,KP_NEIGH)=-FACE_ID
        ENDIF
    ENDIF       
     
    ELSE
    IF (INDEX_S(I1).EQ.1) THEN
        FACE_ID=FACE_ID+1
        IF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.0) THEN
        F_FRONT_RECT(IP,JP)=FACE_ID
        ELSEIF (IF_OWN.EQ.0.AND.IF_NEIGH.EQ.1) THEN
        F_DOWN_CY(IP_NEIGH,KP_NEIGH)=-FACE_ID
        ENDIF
    ENDIF
        
    ENDIF
    
    ENDIF
    ENDDO    
    
    ELSEIF (K.EQ.NL) THEN
    !CASE 3-2-3
    NUM(7)=NCENTER
    CALL GET_CELLID(I,J,K-1,NUM(1),0)
    CALL GET_CELLID(I,1,3*NL+1-J,NUM(2),1)
    NUM(3)=-1
    CALL GET_CELLID(I-1,J,K,NUM(4),0)
    CALL GET_CELLID(I,J+1,K,NUM(5),0)
    CALL GET_CELLID(I,J-1,K,NUM(6),0)
    CALL BUBBLE_SORT(NUM,INDEX_S)
    DO I1=1,7
    IF (NUM(I1).GT.NCENTER) THEN  
    IF_NEIGH=0
    CALL GET_INDEX(NUM(I1),I_NEIGH,J_NEIGH,K_NEIGH,IF_CY)
    IF (IF_CY.EQ.1) THEN
    IP_NEIGH=I_NEIGH-CNI*(MI-1)
    JP_NEIGH=J_NEIGH-CNJ*(MJ-1)
    KP_NEIGH=K_NEIGH
    IF (IP_NEIGH.LE.CNI.AND.IP_NEIGH.GE.1) THEN
    IF (JP_NEIGH.LE.CNJ.AND.JP_NEIGH.GE.1) THEN
        IF_NEIGH=1
    ENDIF
    ENDIF 
    ELSE
    IP_NEIGH=I_NEIGH-CNI*(MI-1)
    JP_NEIGH=J_NEIGH-CNL*(MJ1-1)
    KP_NEIGH=K_NEIGH-CNL*(MJ2-1)
    IF (IP_NEIGH.LE.CNI.AND.IP_NEIGH.GE.1) THEN
    IF (JP_NEIGH.LE.CNL.AND.JP_NEIGH.GE.1) THEN
    IF (KP_NEIGH.LE.CNL.AND.KP_NEIGH.GE.1) THEN
        IF_NEIGH=1 
    ENDIF
    ENDIF
    ENDIF    
    ENDIF
    
    IF (IF_CY.EQ.0) THEN
    IF (INDEX_S(I1).EQ.1) THEN
        FACE_ID=FACE_ID+1
        IF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.1) THEN
        WRITE(3,231) FACE_ID
        ELSEIF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.0) THEN
        F_FRONT_RECT(IP,JP)=FACE_ID
        ELSEIF (IF_OWN.EQ.0.AND.IF_NEIGH.EQ.1) THEN
        F_BACK_RECT(IP_NEIGH,JP_NEIGH)=-FACE_ID
        ENDIF
    ENDIF
    IF (INDEX_S(I1).EQ.2) THEN
        FACE_ID=FACE_ID+1
        IF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.1) THEN
        WRITE(3,231) FACE_ID
        ELSEIF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.0) THEN
        F_BACK_RECT(IP,JP)=FACE_ID
        ELSEIF (IF_OWN.EQ.0.AND.IF_NEIGH.EQ.1) THEN
        F_FRONT_RECT(IP_NEIGH,JP_NEIGH)=-FACE_ID
        ENDIF
    ENDIF
    IF (INDEX_S(I1).EQ.3) THEN
        FACE_ID=FACE_ID+1
        IF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.1) THEN
        WRITE(3,231) FACE_ID
        ELSEIF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.0) THEN
        F_LEFT_RECT(JP,KP)=FACE_ID
        ELSEIF (IF_OWN.EQ.0.AND.IF_NEIGH.EQ.1) THEN
        F_RIGHT_RECT(JP_NEIGH,KP_NEIGH)=-FACE_ID
        ENDIF
    ENDIF  
    IF (INDEX_S(I1).EQ.4) THEN
        FACE_ID=FACE_ID+1
        IF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.1) THEN
        WRITE(3,231) FACE_ID
        ELSEIF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.0) THEN
        F_RIGHT_RECT(JP,KP)=FACE_ID
        ELSEIF (IF_OWN.EQ.0.AND.IF_NEIGH.EQ.1) THEN
        F_LEFT_RECT(JP_NEIGH,KP_NEIGH)=-FACE_ID
        ENDIF
    ENDIF
    IF (INDEX_S(I1).EQ.5) THEN
        FACE_ID=FACE_ID+1
        IF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.1) THEN
        WRITE(3,231) FACE_ID
        ELSEIF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.0) THEN
        F_UP_RECT(IP,KP)=FACE_ID
        ELSEIF (IF_OWN.EQ.0.AND.IF_NEIGH.EQ.1) THEN
        F_DOWN_RECT(IP_NEIGH,KP_NEIGH)=-FACE_ID
        ENDIF
    ENDIF
    IF (INDEX_S(I1).EQ.6) THEN
        FACE_ID=FACE_ID+1
        IF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.1) THEN
        WRITE(3,231) FACE_ID
        ELSEIF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.0) THEN
        F_DOWN_RECT(IP,KP)=FACE_ID
        ELSEIF (IF_OWN.EQ.0.AND.IF_NEIGH.EQ.1) THEN
        F_UP_RECT(IP_NEIGH,KP_NEIGH)=-FACE_ID
        ENDIF
    ENDIF       
     
    ELSE
    IF (INDEX_S(I1).EQ.2) THEN
        FACE_ID=FACE_ID+1
        IF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.0) THEN
        F_BACK_RECT(IP,JP)=FACE_ID
        ELSEIF (IF_OWN.EQ.0.AND.IF_NEIGH.EQ.1) THEN
        F_DOWN_CY(IP_NEIGH,KP_NEIGH)=-FACE_ID
        ENDIF
    ENDIF
        
    ENDIF
    
    ENDIF
    ENDDO    
    
    ELSE
    !CASE 3-2-2 
    NUM(7)=NCENTER
    CALL GET_CELLID(I,J,K-1,NUM(1),0)
    CALL GET_CELLID(I,J,K+1,NUM(2),0)
    NUM(3)=-1
    CALL GET_CELLID(I-1,J,K,NUM(4),0)
    CALL GET_CELLID(I,J+1,K,NUM(5),0)
    CALL GET_CELLID(I,J-1,K,NUM(6),0)
    CALL BUBBLE_SORT(NUM,INDEX_S)
    DO I1=1,7
    IF (NUM(I1).GT.NCENTER) THEN  
    IF_NEIGH=0
    CALL GET_INDEX(NUM(I1),I_NEIGH,J_NEIGH,K_NEIGH,IF_CY)
    IF (IF_CY.EQ.1) THEN
    IP_NEIGH=I_NEIGH-CNI*(MI-1)
    JP_NEIGH=J_NEIGH-CNJ*(MJ-1)
    KP_NEIGH=K_NEIGH
    IF (IP_NEIGH.LE.CNI.AND.IP_NEIGH.GE.1) THEN
    IF (JP_NEIGH.LE.CNJ.AND.JP_NEIGH.GE.1) THEN
        IF_NEIGH=1
    ENDIF
    ENDIF 
    ELSE
    IP_NEIGH=I_NEIGH-CNI*(MI-1)
    JP_NEIGH=J_NEIGH-CNL*(MJ1-1)
    KP_NEIGH=K_NEIGH-CNL*(MJ2-1)
    IF (IP_NEIGH.LE.CNI.AND.IP_NEIGH.GE.1) THEN
    IF (JP_NEIGH.LE.CNL.AND.JP_NEIGH.GE.1) THEN
    IF (KP_NEIGH.LE.CNL.AND.KP_NEIGH.GE.1) THEN
        IF_NEIGH=1 
    ENDIF
    ENDIF
    ENDIF    
    ENDIF
    
    IF (INDEX_S(I1).EQ.1) THEN
        FACE_ID=FACE_ID+1
        IF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.1) THEN
        WRITE(3,231) FACE_ID
        ELSEIF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.0) THEN
        F_FRONT_RECT(IP,JP)=FACE_ID
        ELSEIF (IF_OWN.EQ.0.AND.IF_NEIGH.EQ.1) THEN
        F_BACK_RECT(IP_NEIGH,JP_NEIGH)=-FACE_ID
        ENDIF
    ENDIF
    IF (INDEX_S(I1).EQ.2) THEN
        FACE_ID=FACE_ID+1
        IF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.1) THEN
        WRITE(3,231) FACE_ID
        ELSEIF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.0) THEN
        F_BACK_RECT(IP,JP)=FACE_ID
        ELSEIF (IF_OWN.EQ.0.AND.IF_NEIGH.EQ.1) THEN
        F_FRONT_RECT(IP_NEIGH,JP_NEIGH)=-FACE_ID
        ENDIF
    ENDIF
    IF (INDEX_S(I1).EQ.3) THEN
        FACE_ID=FACE_ID+1
        IF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.1) THEN
        WRITE(3,231) FACE_ID
        ELSEIF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.0) THEN
        F_LEFT_RECT(JP,KP)=FACE_ID
        ELSEIF (IF_OWN.EQ.0.AND.IF_NEIGH.EQ.1) THEN
        F_RIGHT_RECT(JP_NEIGH,KP_NEIGH)=-FACE_ID
        ENDIF
    ENDIF  
    IF (INDEX_S(I1).EQ.4) THEN
        FACE_ID=FACE_ID+1
        IF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.1) THEN
        WRITE(3,231) FACE_ID
        ELSEIF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.0) THEN
        F_RIGHT_RECT(JP,KP)=FACE_ID
        ELSEIF (IF_OWN.EQ.0.AND.IF_NEIGH.EQ.1) THEN
        F_LEFT_RECT(JP_NEIGH,KP_NEIGH)=-FACE_ID
        ENDIF
    ENDIF
    IF (INDEX_S(I1).EQ.5) THEN
        FACE_ID=FACE_ID+1
        IF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.1) THEN
        WRITE(3,231) FACE_ID
        ELSEIF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.0) THEN
        F_UP_RECT(IP,KP)=FACE_ID
        ELSEIF (IF_OWN.EQ.0.AND.IF_NEIGH.EQ.1) THEN
        F_DOWN_RECT(IP_NEIGH,KP_NEIGH)=-FACE_ID
        ENDIF
    ENDIF
    IF (INDEX_S(I1).EQ.6) THEN
        FACE_ID=FACE_ID+1
        IF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.1) THEN
        WRITE(3,231) FACE_ID
        ELSEIF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.0) THEN
        F_DOWN_RECT(IP,KP)=FACE_ID
        ELSEIF (IF_OWN.EQ.0.AND.IF_NEIGH.EQ.1) THEN
        F_UP_RECT(IP_NEIGH,KP_NEIGH)=-FACE_ID
        ENDIF
    ENDIF       
    
    ENDIF
    ENDDO    
    
    ENDIF    
    ENDIF    
    ELSE
    IF (J.EQ.1) THEN
    IF (K.EQ.1) THEN
    !CASE 2-1-1  
    NUM(7)=NCENTER
    CALL GET_CELLID(I,1,1,NUM(1),1)
    CALL GET_CELLID(I,J,K+1,NUM(2),0)
    CALL GET_CELLID(I+1,J,K,NUM(3),0)
    CALL GET_CELLID(I-1,J,K,NUM(4),0)
    CALL GET_CELLID(I,J+1,K,NUM(5),0)
    CALL GET_CELLID(I,1,NK,NUM(6),1)
    CALL BUBBLE_SORT(NUM,INDEX_S)
    DO I1=1,7
    IF (NUM(I1).GT.NCENTER) THEN  
    IF_NEIGH=0
    CALL GET_INDEX(NUM(I1),I_NEIGH,J_NEIGH,K_NEIGH,IF_CY)
    IF (IF_CY.EQ.1) THEN
    IP_NEIGH=I_NEIGH-CNI*(MI-1)
    JP_NEIGH=J_NEIGH-CNJ*(MJ-1)
    KP_NEIGH=K_NEIGH
    IF (IP_NEIGH.LE.CNI.AND.IP_NEIGH.GE.1) THEN
    IF (JP_NEIGH.LE.CNJ.AND.JP_NEIGH.GE.1) THEN
        IF_NEIGH=1
    ENDIF
    ENDIF 
    ELSE
    IP_NEIGH=I_NEIGH-CNI*(MI-1)
    JP_NEIGH=J_NEIGH-CNL*(MJ1-1)
    KP_NEIGH=K_NEIGH-CNL*(MJ2-1)
    IF (IP_NEIGH.LE.CNI.AND.IP_NEIGH.GE.1) THEN
    IF (JP_NEIGH.LE.CNL.AND.JP_NEIGH.GE.1) THEN
    IF (KP_NEIGH.LE.CNL.AND.KP_NEIGH.GE.1) THEN
        IF_NEIGH=1 
    ENDIF
    ENDIF
    ENDIF    
    ENDIF
    
    IF (IF_CY.EQ.0) THEN
    IF (INDEX_S(I1).EQ.1) THEN
        FACE_ID=FACE_ID+1
        IF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.1) THEN
        WRITE(3,231) FACE_ID
        ELSEIF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.0) THEN
        F_FRONT_RECT(IP,JP)=FACE_ID
        ELSEIF (IF_OWN.EQ.0.AND.IF_NEIGH.EQ.1) THEN
        F_BACK_RECT(IP_NEIGH,JP_NEIGH)=-FACE_ID
        ENDIF
    ENDIF
    IF (INDEX_S(I1).EQ.2) THEN
        FACE_ID=FACE_ID+1
        IF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.1) THEN
        WRITE(3,231) FACE_ID
        ELSEIF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.0) THEN
        F_BACK_RECT(IP,JP)=FACE_ID
        ELSEIF (IF_OWN.EQ.0.AND.IF_NEIGH.EQ.1) THEN
        F_FRONT_RECT(IP_NEIGH,JP_NEIGH)=-FACE_ID
        ENDIF
    ENDIF
    IF (INDEX_S(I1).EQ.3) THEN
        FACE_ID=FACE_ID+1
        IF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.1) THEN
        WRITE(3,231) FACE_ID
        ELSEIF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.0) THEN
        F_LEFT_RECT(JP,KP)=FACE_ID
        ELSEIF (IF_OWN.EQ.0.AND.IF_NEIGH.EQ.1) THEN
        F_RIGHT_RECT(JP_NEIGH,KP_NEIGH)=-FACE_ID
        ENDIF
    ENDIF  
    IF (INDEX_S(I1).EQ.4) THEN
        FACE_ID=FACE_ID+1
        IF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.1) THEN
        WRITE(3,231) FACE_ID
        ELSEIF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.0) THEN
        F_RIGHT_RECT(JP,KP)=FACE_ID
        ELSEIF (IF_OWN.EQ.0.AND.IF_NEIGH.EQ.1) THEN
        F_LEFT_RECT(JP_NEIGH,KP_NEIGH)=-FACE_ID
        ENDIF
    ENDIF
    IF (INDEX_S(I1).EQ.5) THEN
        FACE_ID=FACE_ID+1
        IF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.1) THEN
        WRITE(3,231) FACE_ID
        ELSEIF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.0) THEN
        F_UP_RECT(IP,KP)=FACE_ID
        ELSEIF (IF_OWN.EQ.0.AND.IF_NEIGH.EQ.1) THEN
        F_DOWN_RECT(IP_NEIGH,KP_NEIGH)=-FACE_ID
        ENDIF
    ENDIF
    IF (INDEX_S(I1).EQ.6) THEN
        FACE_ID=FACE_ID+1
        IF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.1) THEN
        WRITE(3,231) FACE_ID
        ELSEIF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.0) THEN
        F_DOWN_RECT(IP,KP)=FACE_ID
        ELSEIF (IF_OWN.EQ.0.AND.IF_NEIGH.EQ.1) THEN
        F_UP_RECT(IP_NEIGH,KP_NEIGH)=-FACE_ID
        ENDIF
    ENDIF       
     
    ELSE
    IF (INDEX_S(I1).EQ.1) THEN
        FACE_ID=FACE_ID+1
        IF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.0) THEN
        F_FRONT_RECT(IP,JP)=FACE_ID
        ELSEIF (IF_OWN.EQ.0.AND.IF_NEIGH.EQ.1) THEN
        F_DOWN_CY(IP_NEIGH,KP_NEIGH)=-FACE_ID
        ENDIF
    ENDIF
    IF (INDEX_S(I1).EQ.6) THEN
        FACE_ID=FACE_ID+1
        IF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.0) THEN
        F_DOWN_RECT(IP,KP)=FACE_ID
        ELSEIF (IF_OWN.EQ.0.AND.IF_NEIGH.EQ.1) THEN
        F_DOWN_CY(IP_NEIGH,KP_NEIGH)=-FACE_ID
        ENDIF
    ENDIF
        
    ENDIF
    
    ENDIF
    ENDDO    
    
    ELSEIF (K.EQ.NL) THEN
    !CASE 2-1-3 
    NUM(7)=NCENTER
    CALL GET_CELLID(I,J,K-1,NUM(1),0)
    CALL GET_CELLID(I,1,3*NL,NUM(2),1)
    CALL GET_CELLID(I+1,J,K,NUM(3),0)
    CALL GET_CELLID(I-1,J,K,NUM(4),0)
    CALL GET_CELLID(I,J+1,K,NUM(5),0)
    CALL GET_CELLID(I,1,3*NL+1,NUM(6),1)
    CALL BUBBLE_SORT(NUM,INDEX_S)
    DO I1=1,7
    IF (NUM(I1).GT.NCENTER) THEN  
    IF_NEIGH=0
    CALL GET_INDEX(NUM(I1),I_NEIGH,J_NEIGH,K_NEIGH,IF_CY)
    IF (IF_CY.EQ.1) THEN
    IP_NEIGH=I_NEIGH-CNI*(MI-1)
    JP_NEIGH=J_NEIGH-CNJ*(MJ-1)
    KP_NEIGH=K_NEIGH
    IF (IP_NEIGH.LE.CNI.AND.IP_NEIGH.GE.1) THEN
    IF (JP_NEIGH.LE.CNJ.AND.JP_NEIGH.GE.1) THEN
        IF_NEIGH=1
    ENDIF
    ENDIF 
    ELSE
    IP_NEIGH=I_NEIGH-CNI*(MI-1)
    JP_NEIGH=J_NEIGH-CNL*(MJ1-1)
    KP_NEIGH=K_NEIGH-CNL*(MJ2-1)
    IF (IP_NEIGH.LE.CNI.AND.IP_NEIGH.GE.1) THEN
    IF (JP_NEIGH.LE.CNL.AND.JP_NEIGH.GE.1) THEN
    IF (KP_NEIGH.LE.CNL.AND.KP_NEIGH.GE.1) THEN
        IF_NEIGH=1 
    ENDIF
    ENDIF
    ENDIF    
    ENDIF
    
    IF (IF_CY.EQ.0) THEN
    IF (INDEX_S(I1).EQ.1) THEN
        FACE_ID=FACE_ID+1
        IF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.1) THEN
        WRITE(3,231) FACE_ID
        ELSEIF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.0) THEN
        F_FRONT_RECT(IP,JP)=FACE_ID
        ELSEIF (IF_OWN.EQ.0.AND.IF_NEIGH.EQ.1) THEN
        F_BACK_RECT(IP_NEIGH,JP_NEIGH)=-FACE_ID
        ENDIF
    ENDIF
    IF (INDEX_S(I1).EQ.2) THEN
        FACE_ID=FACE_ID+1
        IF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.1) THEN
        WRITE(3,231) FACE_ID
        ELSEIF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.0) THEN
        F_BACK_RECT(IP,JP)=FACE_ID
        ELSEIF (IF_OWN.EQ.0.AND.IF_NEIGH.EQ.1) THEN
        F_FRONT_RECT(IP_NEIGH,JP_NEIGH)=-FACE_ID
        ENDIF
    ENDIF
    IF (INDEX_S(I1).EQ.3) THEN
        FACE_ID=FACE_ID+1
        IF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.1) THEN
        WRITE(3,231) FACE_ID
        ELSEIF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.0) THEN
        F_LEFT_RECT(JP,KP)=FACE_ID
        ELSEIF (IF_OWN.EQ.0.AND.IF_NEIGH.EQ.1) THEN
        F_RIGHT_RECT(JP_NEIGH,KP_NEIGH)=-FACE_ID
        ENDIF
    ENDIF  
    IF (INDEX_S(I1).EQ.4) THEN
        FACE_ID=FACE_ID+1
        IF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.1) THEN
        WRITE(3,231) FACE_ID
        ELSEIF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.0) THEN
        F_RIGHT_RECT(JP,KP)=FACE_ID
        ELSEIF (IF_OWN.EQ.0.AND.IF_NEIGH.EQ.1) THEN
        F_LEFT_RECT(JP_NEIGH,KP_NEIGH)=-FACE_ID
        ENDIF
    ENDIF
    IF (INDEX_S(I1).EQ.5) THEN
        FACE_ID=FACE_ID+1
        IF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.1) THEN
        WRITE(3,231) FACE_ID
        ELSEIF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.0) THEN
        F_UP_RECT(IP,KP)=FACE_ID
        ELSEIF (IF_OWN.EQ.0.AND.IF_NEIGH.EQ.1) THEN
        F_DOWN_RECT(IP_NEIGH,KP_NEIGH)=-FACE_ID
        ENDIF
    ENDIF
    IF (INDEX_S(I1).EQ.6) THEN
        FACE_ID=FACE_ID+1
        IF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.1) THEN
        WRITE(3,231) FACE_ID
        ELSEIF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.0) THEN
        F_DOWN_RECT(IP,KP)=FACE_ID
        ELSEIF (IF_OWN.EQ.0.AND.IF_NEIGH.EQ.1) THEN
        F_UP_RECT(IP_NEIGH,KP_NEIGH)=-FACE_ID
        ENDIF
    ENDIF       
     
    ELSE
    IF (INDEX_S(I1).EQ.2) THEN
        FACE_ID=FACE_ID+1
        IF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.0) THEN
        F_BACK_RECT(IP,JP)=FACE_ID
        ELSEIF (IF_OWN.EQ.0.AND.IF_NEIGH.EQ.1) THEN
        F_DOWN_CY(IP_NEIGH,KP_NEIGH)=-FACE_ID
        ENDIF
    ENDIF
    IF (INDEX_S(I1).EQ.6) THEN
        FACE_ID=FACE_ID+1
        IF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.0) THEN
        F_DOWN_RECT(IP,KP)=FACE_ID
        ELSEIF (IF_OWN.EQ.0.AND.IF_NEIGH.EQ.1) THEN
        F_DOWN_CY(IP_NEIGH,KP_NEIGH)=-FACE_ID
        ENDIF
    ENDIF
        
    ENDIF
    
    ENDIF
    ENDDO    
    
    ELSE
    !CASE 2-1-2 
    NUM(7)=NCENTER
    CALL GET_CELLID(I,J,K-1,NUM(1),0)
    CALL GET_CELLID(I,J,K+1,NUM(2),0)
    CALL GET_CELLID(I+1,J,K,NUM(3),0)
    CALL GET_CELLID(I-1,J,K,NUM(4),0)
    CALL GET_CELLID(I,J+1,K,NUM(5),0)
    CALL GET_CELLID(I,1,NK+1-J,NUM(6),1)
    CALL BUBBLE_SORT(NUM,INDEX_S)
    DO I1=1,7
    IF (NUM(I1).GT.NCENTER) THEN  
    IF_NEIGH=0
    CALL GET_INDEX(NUM(I1),I_NEIGH,J_NEIGH,K_NEIGH,IF_CY)
    IF (IF_CY.EQ.1) THEN
    IP_NEIGH=I_NEIGH-CNI*(MI-1)
    JP_NEIGH=J_NEIGH-CNJ*(MJ-1)
    KP_NEIGH=K_NEIGH
    IF (IP_NEIGH.LE.CNI.AND.IP_NEIGH.GE.1) THEN
    IF (JP_NEIGH.LE.CNJ.AND.JP_NEIGH.GE.1) THEN
        IF_NEIGH=1
    ENDIF
    ENDIF 
    ELSE
    IP_NEIGH=I_NEIGH-CNI*(MI-1)
    JP_NEIGH=J_NEIGH-CNL*(MJ1-1)
    KP_NEIGH=K_NEIGH-CNL*(MJ2-1)
    IF (IP_NEIGH.LE.CNI.AND.IP_NEIGH.GE.1) THEN
    IF (JP_NEIGH.LE.CNL.AND.JP_NEIGH.GE.1) THEN
    IF (KP_NEIGH.LE.CNL.AND.KP_NEIGH.GE.1) THEN
        IF_NEIGH=1 
    ENDIF
    ENDIF
    ENDIF    
    ENDIF
    
    IF (IF_CY.EQ.0) THEN
    IF (INDEX_S(I1).EQ.1) THEN
        FACE_ID=FACE_ID+1
        IF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.1) THEN
        WRITE(3,231) FACE_ID
        ELSEIF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.0) THEN
        F_FRONT_RECT(IP,JP)=FACE_ID
        ELSEIF (IF_OWN.EQ.0.AND.IF_NEIGH.EQ.1) THEN
        F_BACK_RECT(IP_NEIGH,JP_NEIGH)=-FACE_ID
        ENDIF
    ENDIF
    IF (INDEX_S(I1).EQ.2) THEN
        FACE_ID=FACE_ID+1
        IF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.1) THEN
        WRITE(3,231) FACE_ID
        ELSEIF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.0) THEN
        F_BACK_RECT(IP,JP)=FACE_ID
        ELSEIF (IF_OWN.EQ.0.AND.IF_NEIGH.EQ.1) THEN
        F_FRONT_RECT(IP_NEIGH,JP_NEIGH)=-FACE_ID
        ENDIF
    ENDIF
    IF (INDEX_S(I1).EQ.3) THEN
        FACE_ID=FACE_ID+1
        IF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.1) THEN
        WRITE(3,231) FACE_ID
        ELSEIF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.0) THEN
        F_LEFT_RECT(JP,KP)=FACE_ID
        ELSEIF (IF_OWN.EQ.0.AND.IF_NEIGH.EQ.1) THEN
        F_RIGHT_RECT(JP_NEIGH,KP_NEIGH)=-FACE_ID
        ENDIF
    ENDIF  
    IF (INDEX_S(I1).EQ.4) THEN
        FACE_ID=FACE_ID+1
        IF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.1) THEN
        WRITE(3,231) FACE_ID
        ELSEIF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.0) THEN
        F_RIGHT_RECT(JP,KP)=FACE_ID
        ELSEIF (IF_OWN.EQ.0.AND.IF_NEIGH.EQ.1) THEN
        F_LEFT_RECT(JP_NEIGH,KP_NEIGH)=-FACE_ID
        ENDIF
    ENDIF
    IF (INDEX_S(I1).EQ.5) THEN
        FACE_ID=FACE_ID+1
        IF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.1) THEN
        WRITE(3,231) FACE_ID
        ELSEIF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.0) THEN
        F_UP_RECT(IP,KP)=FACE_ID
        ELSEIF (IF_OWN.EQ.0.AND.IF_NEIGH.EQ.1) THEN
        F_DOWN_RECT(IP_NEIGH,KP_NEIGH)=-FACE_ID
        ENDIF
    ENDIF
    IF (INDEX_S(I1).EQ.6) THEN
        FACE_ID=FACE_ID+1
        IF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.1) THEN
        WRITE(3,231) FACE_ID
        ELSEIF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.0) THEN
        F_DOWN_RECT(IP,KP)=FACE_ID
        ELSEIF (IF_OWN.EQ.0.AND.IF_NEIGH.EQ.1) THEN
        F_UP_RECT(IP_NEIGH,KP_NEIGH)=-FACE_ID
        ENDIF
    ENDIF       
     
    ELSE
    IF (INDEX_S(I1).EQ.6) THEN
        FACE_ID=FACE_ID+1
        IF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.0) THEN
        F_DOWN_RECT(IP,KP)=FACE_ID
        ELSEIF (IF_OWN.EQ.0.AND.IF_NEIGH.EQ.1) THEN
        F_DOWN_CY(IP_NEIGH,KP_NEIGH)=-FACE_ID
        ENDIF
    ENDIF
        
    ENDIF
    
    ENDIF
    ENDDO    
    
    ENDIF    
    ELSEIF (J.EQ.NL) THEN
    IF (K.EQ.1) THEN
    !CASE 2-3-1  
    NUM(7)=NCENTER
    CALL GET_CELLID(I,1,1,NUM(1),1)
    CALL GET_CELLID(I,J,K+1,NUM(2),0)
    CALL GET_CELLID(I+1,J,K,NUM(3),0)
    CALL GET_CELLID(I-1,J,K,NUM(4),0)
    CALL GET_CELLID(I,1,NL+1,NUM(5),1)
    CALL GET_CELLID(I,J-1,K,NUM(6),0)
    CALL BUBBLE_SORT(NUM,INDEX_S)
    DO I1=1,7
    IF (NUM(I1).GT.NCENTER) THEN  
    IF_NEIGH=0
    CALL GET_INDEX(NUM(I1),I_NEIGH,J_NEIGH,K_NEIGH,IF_CY)
    IF (IF_CY.EQ.1) THEN
    IP_NEIGH=I_NEIGH-CNI*(MI-1)
    JP_NEIGH=J_NEIGH-CNJ*(MJ-1)
    KP_NEIGH=K_NEIGH
    IF (IP_NEIGH.LE.CNI.AND.IP_NEIGH.GE.1) THEN
    IF (JP_NEIGH.LE.CNJ.AND.JP_NEIGH.GE.1) THEN
        IF_NEIGH=1
    ENDIF
    ENDIF 
    ELSE
    IP_NEIGH=I_NEIGH-CNI*(MI-1)
    JP_NEIGH=J_NEIGH-CNL*(MJ1-1)
    KP_NEIGH=K_NEIGH-CNL*(MJ2-1)
    IF (IP_NEIGH.LE.CNI.AND.IP_NEIGH.GE.1) THEN
    IF (JP_NEIGH.LE.CNL.AND.JP_NEIGH.GE.1) THEN
    IF (KP_NEIGH.LE.CNL.AND.KP_NEIGH.GE.1) THEN
        IF_NEIGH=1 
    ENDIF
    ENDIF
    ENDIF    
    ENDIF
    
    IF (IF_CY.EQ.0) THEN
    IF (INDEX_S(I1).EQ.1) THEN
        FACE_ID=FACE_ID+1
        IF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.1) THEN
        WRITE(3,231) FACE_ID
        ELSEIF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.0) THEN
        F_FRONT_RECT(IP,JP)=FACE_ID
        ELSEIF (IF_OWN.EQ.0.AND.IF_NEIGH.EQ.1) THEN
        F_BACK_RECT(IP_NEIGH,JP_NEIGH)=-FACE_ID
        ENDIF
    ENDIF
    IF (INDEX_S(I1).EQ.2) THEN
        FACE_ID=FACE_ID+1
        IF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.1) THEN
        WRITE(3,231) FACE_ID
        ELSEIF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.0) THEN
        F_BACK_RECT(IP,JP)=FACE_ID
        ELSEIF (IF_OWN.EQ.0.AND.IF_NEIGH.EQ.1) THEN
        F_FRONT_RECT(IP_NEIGH,JP_NEIGH)=-FACE_ID
        ENDIF
    ENDIF
    IF (INDEX_S(I1).EQ.3) THEN
        FACE_ID=FACE_ID+1
        IF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.1) THEN
        WRITE(3,231) FACE_ID
        ELSEIF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.0) THEN
        F_LEFT_RECT(JP,KP)=FACE_ID
        ELSEIF (IF_OWN.EQ.0.AND.IF_NEIGH.EQ.1) THEN
        F_RIGHT_RECT(JP_NEIGH,KP_NEIGH)=-FACE_ID
        ENDIF
    ENDIF  
    IF (INDEX_S(I1).EQ.4) THEN
        FACE_ID=FACE_ID+1
        IF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.1) THEN
        WRITE(3,231) FACE_ID
        ELSEIF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.0) THEN
        F_RIGHT_RECT(JP,KP)=FACE_ID
        ELSEIF (IF_OWN.EQ.0.AND.IF_NEIGH.EQ.1) THEN
        F_LEFT_RECT(JP_NEIGH,KP_NEIGH)=-FACE_ID
        ENDIF
    ENDIF
    IF (INDEX_S(I1).EQ.5) THEN
        FACE_ID=FACE_ID+1
        IF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.1) THEN
        WRITE(3,231) FACE_ID
        ELSEIF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.0) THEN
        F_UP_RECT(IP,KP)=FACE_ID
        ELSEIF (IF_OWN.EQ.0.AND.IF_NEIGH.EQ.1) THEN
        F_DOWN_RECT(IP_NEIGH,KP_NEIGH)=-FACE_ID
        ENDIF
    ENDIF
    IF (INDEX_S(I1).EQ.6) THEN
        FACE_ID=FACE_ID+1
        IF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.1) THEN
        WRITE(3,231) FACE_ID
        ELSEIF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.0) THEN
        F_DOWN_RECT(IP,KP)=FACE_ID
        ELSEIF (IF_OWN.EQ.0.AND.IF_NEIGH.EQ.1) THEN
        F_UP_RECT(IP_NEIGH,KP_NEIGH)=-FACE_ID
        ENDIF
    ENDIF       
     
    ELSE
    IF (INDEX_S(I1).EQ.1) THEN
        FACE_ID=FACE_ID+1
        IF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.0) THEN
        F_FRONT_RECT(IP,JP)=FACE_ID
        ELSEIF (IF_OWN.EQ.0.AND.IF_NEIGH.EQ.1) THEN
        F_DOWN_CY(IP_NEIGH,KP_NEIGH)=-FACE_ID
        ENDIF
    ENDIF
    IF (INDEX_S(I1).EQ.5) THEN
        FACE_ID=FACE_ID+1
        IF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.0) THEN
        F_UP_RECT(IP,KP)=FACE_ID
        ELSEIF (IF_OWN.EQ.0.AND.IF_NEIGH.EQ.1) THEN
        F_DOWN_CY(IP_NEIGH,KP_NEIGH)=-FACE_ID
        ENDIF
    ENDIF
        
    ENDIF
    
    ENDIF
    ENDDO    
    
    ELSEIF (K.EQ.NL) THEN
    !CASE 2-3-3 
    NUM(7)=NCENTER
    CALL GET_CELLID(I,J,K-1,NUM(1),0)
    CALL GET_CELLID(I,1,2*NL+1,NUM(2),1)
    CALL GET_CELLID(I+1,J,K,NUM(3),0)
    CALL GET_CELLID(I-1,J,K,NUM(4),0)
    CALL GET_CELLID(I,1,2*NL,NUM(5),1)
    CALL GET_CELLID(I,J-1,K,NUM(6),0)
    CALL BUBBLE_SORT(NUM,INDEX_S)
    DO I1=1,7
    IF (NUM(I1).GT.NCENTER) THEN  
    IF_NEIGH=0
    CALL GET_INDEX(NUM(I1),I_NEIGH,J_NEIGH,K_NEIGH,IF_CY)
    IF (IF_CY.EQ.1) THEN
    IP_NEIGH=I_NEIGH-CNI*(MI-1)
    JP_NEIGH=J_NEIGH-CNJ*(MJ-1)
    KP_NEIGH=K_NEIGH
    IF (IP_NEIGH.LE.CNI.AND.IP_NEIGH.GE.1) THEN
    IF (JP_NEIGH.LE.CNJ.AND.JP_NEIGH.GE.1) THEN
        IF_NEIGH=1
    ENDIF
    ENDIF 
    ELSE
    IP_NEIGH=I_NEIGH-CNI*(MI-1)
    JP_NEIGH=J_NEIGH-CNL*(MJ1-1)
    KP_NEIGH=K_NEIGH-CNL*(MJ2-1)
    IF (IP_NEIGH.LE.CNI.AND.IP_NEIGH.GE.1) THEN
    IF (JP_NEIGH.LE.CNL.AND.JP_NEIGH.GE.1) THEN
    IF (KP_NEIGH.LE.CNL.AND.KP_NEIGH.GE.1) THEN
        IF_NEIGH=1 
    ENDIF
    ENDIF
    ENDIF    
    ENDIF
    
    IF (IF_CY.EQ.0) THEN
    IF (INDEX_S(I1).EQ.1) THEN
        FACE_ID=FACE_ID+1
        IF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.1) THEN
        WRITE(3,231) FACE_ID
        ELSEIF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.0) THEN
        F_FRONT_RECT(IP,JP)=FACE_ID
        ELSEIF (IF_OWN.EQ.0.AND.IF_NEIGH.EQ.1) THEN
        F_BACK_RECT(IP_NEIGH,JP_NEIGH)=-FACE_ID
        ENDIF
    ENDIF
    IF (INDEX_S(I1).EQ.2) THEN
        FACE_ID=FACE_ID+1
        IF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.1) THEN
        WRITE(3,231) FACE_ID
        ELSEIF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.0) THEN
        F_BACK_RECT(IP,JP)=FACE_ID
        ELSEIF (IF_OWN.EQ.0.AND.IF_NEIGH.EQ.1) THEN
        F_FRONT_RECT(IP_NEIGH,JP_NEIGH)=-FACE_ID
        ENDIF
    ENDIF
    IF (INDEX_S(I1).EQ.3) THEN
        FACE_ID=FACE_ID+1
        IF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.1) THEN
        WRITE(3,231) FACE_ID
        ELSEIF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.0) THEN
        F_LEFT_RECT(JP,KP)=FACE_ID
        ELSEIF (IF_OWN.EQ.0.AND.IF_NEIGH.EQ.1) THEN
        F_RIGHT_RECT(JP_NEIGH,KP_NEIGH)=-FACE_ID
        ENDIF
    ENDIF  
    IF (INDEX_S(I1).EQ.4) THEN
        FACE_ID=FACE_ID+1
        IF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.1) THEN
        WRITE(3,231) FACE_ID
        ELSEIF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.0) THEN
        F_RIGHT_RECT(JP,KP)=FACE_ID
        ELSEIF (IF_OWN.EQ.0.AND.IF_NEIGH.EQ.1) THEN
        F_LEFT_RECT(JP_NEIGH,KP_NEIGH)=-FACE_ID
        ENDIF
    ENDIF
    IF (INDEX_S(I1).EQ.5) THEN
        FACE_ID=FACE_ID+1
        IF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.1) THEN
        WRITE(3,231) FACE_ID
        ELSEIF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.0) THEN
        F_UP_RECT(IP,KP)=FACE_ID
        ELSEIF (IF_OWN.EQ.0.AND.IF_NEIGH.EQ.1) THEN
        F_DOWN_RECT(IP_NEIGH,KP_NEIGH)=-FACE_ID
        ENDIF
    ENDIF
    IF (INDEX_S(I1).EQ.6) THEN
        FACE_ID=FACE_ID+1
        IF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.1) THEN
        WRITE(3,231) FACE_ID
        ELSEIF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.0) THEN
        F_DOWN_RECT(IP,KP)=FACE_ID
        ELSEIF (IF_OWN.EQ.0.AND.IF_NEIGH.EQ.1) THEN
        F_UP_RECT(IP_NEIGH,KP_NEIGH)=-FACE_ID
        ENDIF
    ENDIF       
     
    ELSE
    IF (INDEX_S(I1).EQ.2) THEN
        FACE_ID=FACE_ID+1
        IF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.0) THEN
        F_BACK_RECT(IP,JP)=FACE_ID
        ELSEIF (IF_OWN.EQ.0.AND.IF_NEIGH.EQ.1) THEN
        F_DOWN_CY(IP_NEIGH,KP_NEIGH)=-FACE_ID
        ENDIF
    ENDIF
    IF (INDEX_S(I1).EQ.5) THEN
        FACE_ID=FACE_ID+1
        IF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.0) THEN
        F_UP_RECT(IP,KP)=FACE_ID
        ELSEIF (IF_OWN.EQ.0.AND.IF_NEIGH.EQ.1) THEN
        F_DOWN_CY(IP_NEIGH,KP_NEIGH)=-FACE_ID
        ENDIF
    ENDIF
        
    ENDIF
    
    ENDIF
    ENDDO    
    
    ELSE
    !CASE 2-3-2  
    NUM(7)=NCENTER
    CALL GET_CELLID(I,J,K-1,NUM(1),0)
    CALL GET_CELLID(I,J,K+1,NUM(2),0)
    CALL GET_CELLID(I+1,J,K,NUM(3),0)
    CALL GET_CELLID(I-1,J,K,NUM(4),0)
    CALL GET_CELLID(I,1,NL+K,NUM(5),1)
    CALL GET_CELLID(I,J-1,K,NUM(6),0)
    CALL BUBBLE_SORT(NUM,INDEX_S)
    DO I1=1,7
    IF (NUM(I1).GT.NCENTER) THEN  
    IF_NEIGH=0
    CALL GET_INDEX(NUM(I1),I_NEIGH,J_NEIGH,K_NEIGH,IF_CY)
    IF (IF_CY.EQ.1) THEN
    IP_NEIGH=I_NEIGH-CNI*(MI-1)
    JP_NEIGH=J_NEIGH-CNJ*(MJ-1)
    KP_NEIGH=K_NEIGH
    IF (IP_NEIGH.LE.CNI.AND.IP_NEIGH.GE.1) THEN
    IF (JP_NEIGH.LE.CNJ.AND.JP_NEIGH.GE.1) THEN
        IF_NEIGH=1
    ENDIF
    ENDIF 
    ELSE
    IP_NEIGH=I_NEIGH-CNI*(MI-1)
    JP_NEIGH=J_NEIGH-CNL*(MJ1-1)
    KP_NEIGH=K_NEIGH-CNL*(MJ2-1)
    IF (IP_NEIGH.LE.CNI.AND.IP_NEIGH.GE.1) THEN
    IF (JP_NEIGH.LE.CNL.AND.JP_NEIGH.GE.1) THEN
    IF (KP_NEIGH.LE.CNL.AND.KP_NEIGH.GE.1) THEN
        IF_NEIGH=1 
    ENDIF
    ENDIF
    ENDIF    
    ENDIF
    
    IF (IF_CY.EQ.0) THEN
    IF (INDEX_S(I1).EQ.1) THEN
        FACE_ID=FACE_ID+1
        IF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.1) THEN
        WRITE(3,231) FACE_ID
        ELSEIF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.0) THEN
        F_FRONT_RECT(IP,JP)=FACE_ID
        ELSEIF (IF_OWN.EQ.0.AND.IF_NEIGH.EQ.1) THEN
        F_BACK_RECT(IP_NEIGH,JP_NEIGH)=-FACE_ID
        ENDIF
    ENDIF
    IF (INDEX_S(I1).EQ.2) THEN
        FACE_ID=FACE_ID+1
        IF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.1) THEN
        WRITE(3,231) FACE_ID
        ELSEIF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.0) THEN
        F_BACK_RECT(IP,JP)=FACE_ID
        ELSEIF (IF_OWN.EQ.0.AND.IF_NEIGH.EQ.1) THEN
        F_FRONT_RECT(IP_NEIGH,JP_NEIGH)=-FACE_ID
        ENDIF
    ENDIF
    IF (INDEX_S(I1).EQ.3) THEN
        FACE_ID=FACE_ID+1
        IF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.1) THEN
        WRITE(3,231) FACE_ID
        ELSEIF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.0) THEN
        F_LEFT_RECT(JP,KP)=FACE_ID
        ELSEIF (IF_OWN.EQ.0.AND.IF_NEIGH.EQ.1) THEN
        F_RIGHT_RECT(JP_NEIGH,KP_NEIGH)=-FACE_ID
        ENDIF
    ENDIF  
    IF (INDEX_S(I1).EQ.4) THEN
        FACE_ID=FACE_ID+1
        IF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.1) THEN
        WRITE(3,231) FACE_ID
        ELSEIF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.0) THEN
        F_RIGHT_RECT(JP,KP)=FACE_ID
        ELSEIF (IF_OWN.EQ.0.AND.IF_NEIGH.EQ.1) THEN
        F_LEFT_RECT(JP_NEIGH,KP_NEIGH)=-FACE_ID
        ENDIF
    ENDIF
    IF (INDEX_S(I1).EQ.5) THEN
        FACE_ID=FACE_ID+1
        IF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.1) THEN
        WRITE(3,231) FACE_ID
        ELSEIF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.0) THEN
        F_UP_RECT(IP,KP)=FACE_ID
        ELSEIF (IF_OWN.EQ.0.AND.IF_NEIGH.EQ.1) THEN
        F_DOWN_RECT(IP_NEIGH,KP_NEIGH)=-FACE_ID
        ENDIF
    ENDIF
    IF (INDEX_S(I1).EQ.6) THEN
        FACE_ID=FACE_ID+1
        IF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.1) THEN
        WRITE(3,231) FACE_ID
        ELSEIF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.0) THEN
        F_DOWN_RECT(IP,KP)=FACE_ID
        ELSEIF (IF_OWN.EQ.0.AND.IF_NEIGH.EQ.1) THEN
        F_UP_RECT(IP_NEIGH,KP_NEIGH)=-FACE_ID
        ENDIF
    ENDIF       
     
    ELSE
    IF (INDEX_S(I1).EQ.5) THEN
        FACE_ID=FACE_ID+1
        IF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.0) THEN
        F_UP_RECT(IP,KP)=FACE_ID
        ELSEIF (IF_OWN.EQ.0.AND.IF_NEIGH.EQ.1) THEN
        F_DOWN_CY(IP_NEIGH,KP_NEIGH)=-FACE_ID
        ENDIF
    ENDIF
        
    ENDIF
    
    ENDIF
    ENDDO    
    
    ENDIF    
    ELSE
    IF (K.EQ.1) THEN
    !CASE 2-2-1  
    NUM(7)=NCENTER
    CALL GET_CELLID(I,1,J,NUM(1),1)
    CALL GET_CELLID(I,J,K+1,NUM(2),0)
    CALL GET_CELLID(I+1,J,K,NUM(3),0)
    CALL GET_CELLID(I-1,J,K,NUM(4),0)
    CALL GET_CELLID(I,J+1,K,NUM(5),0)
    CALL GET_CELLID(I,J-1,K,NUM(6),0)
    CALL BUBBLE_SORT(NUM,INDEX_S)
    DO I1=1,7
    IF (NUM(I1).GT.NCENTER) THEN  
    IF_NEIGH=0
    CALL GET_INDEX(NUM(I1),I_NEIGH,J_NEIGH,K_NEIGH,IF_CY)
    IF (IF_CY.EQ.1) THEN
    IP_NEIGH=I_NEIGH-CNI*(MI-1)
    JP_NEIGH=J_NEIGH-CNJ*(MJ-1)
    KP_NEIGH=K_NEIGH
    IF (IP_NEIGH.LE.CNI.AND.IP_NEIGH.GE.1) THEN
    IF (JP_NEIGH.LE.CNJ.AND.JP_NEIGH.GE.1) THEN
        IF_NEIGH=1
    ENDIF
    ENDIF 
    ELSE
    IP_NEIGH=I_NEIGH-CNI*(MI-1)
    JP_NEIGH=J_NEIGH-CNL*(MJ1-1)
    KP_NEIGH=K_NEIGH-CNL*(MJ2-1)
    IF (IP_NEIGH.LE.CNI.AND.IP_NEIGH.GE.1) THEN
    IF (JP_NEIGH.LE.CNL.AND.JP_NEIGH.GE.1) THEN
    IF (KP_NEIGH.LE.CNL.AND.KP_NEIGH.GE.1) THEN
        IF_NEIGH=1 
    ENDIF
    ENDIF
    ENDIF    
    ENDIF
    
    IF (IF_CY.EQ.0) THEN
    IF (INDEX_S(I1).EQ.1) THEN
        FACE_ID=FACE_ID+1
        IF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.1) THEN
        WRITE(3,231) FACE_ID
        ELSEIF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.0) THEN
        F_FRONT_RECT(IP,JP)=FACE_ID
        ELSEIF (IF_OWN.EQ.0.AND.IF_NEIGH.EQ.1) THEN
        F_BACK_RECT(IP_NEIGH,JP_NEIGH)=-FACE_ID
        ENDIF
    ENDIF
    IF (INDEX_S(I1).EQ.2) THEN
        FACE_ID=FACE_ID+1
        IF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.1) THEN
        WRITE(3,231) FACE_ID
        ELSEIF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.0) THEN
        F_BACK_RECT(IP,JP)=FACE_ID
        ELSEIF (IF_OWN.EQ.0.AND.IF_NEIGH.EQ.1) THEN
        F_FRONT_RECT(IP_NEIGH,JP_NEIGH)=-FACE_ID
        ENDIF
    ENDIF
    IF (INDEX_S(I1).EQ.3) THEN
        FACE_ID=FACE_ID+1
        IF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.1) THEN
        WRITE(3,231) FACE_ID
        ELSEIF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.0) THEN
        F_LEFT_RECT(JP,KP)=FACE_ID
        ELSEIF (IF_OWN.EQ.0.AND.IF_NEIGH.EQ.1) THEN
        F_RIGHT_RECT(JP_NEIGH,KP_NEIGH)=-FACE_ID
        ENDIF
    ENDIF  
    IF (INDEX_S(I1).EQ.4) THEN
        FACE_ID=FACE_ID+1
        IF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.1) THEN
        WRITE(3,231) FACE_ID
        ELSEIF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.0) THEN
        F_RIGHT_RECT(JP,KP)=FACE_ID
        ELSEIF (IF_OWN.EQ.0.AND.IF_NEIGH.EQ.1) THEN
        F_LEFT_RECT(JP_NEIGH,KP_NEIGH)=-FACE_ID
        ENDIF
    ENDIF
    IF (INDEX_S(I1).EQ.5) THEN
        FACE_ID=FACE_ID+1
        IF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.1) THEN
        WRITE(3,231) FACE_ID
        ELSEIF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.0) THEN
        F_UP_RECT(IP,KP)=FACE_ID
        ELSEIF (IF_OWN.EQ.0.AND.IF_NEIGH.EQ.1) THEN
        F_DOWN_RECT(IP_NEIGH,KP_NEIGH)=-FACE_ID
        ENDIF
    ENDIF
    IF (INDEX_S(I1).EQ.6) THEN
        FACE_ID=FACE_ID+1
        IF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.1) THEN
        WRITE(3,231) FACE_ID
        ELSEIF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.0) THEN
        F_DOWN_RECT(IP,KP)=FACE_ID
        ELSEIF (IF_OWN.EQ.0.AND.IF_NEIGH.EQ.1) THEN
        F_UP_RECT(IP_NEIGH,KP_NEIGH)=-FACE_ID
        ENDIF
    ENDIF       
     
    ELSE
    IF (INDEX_S(I1).EQ.1) THEN
        FACE_ID=FACE_ID+1
        IF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.0) THEN
        F_FRONT_RECT(IP,JP)=FACE_ID
        ELSEIF (IF_OWN.EQ.0.AND.IF_NEIGH.EQ.1) THEN
        F_DOWN_CY(IP_NEIGH,KP_NEIGH)=-FACE_ID
        ENDIF
    ENDIF
        
    ENDIF
    
    ENDIF
    ENDDO    
    
    ELSEIF (K.EQ.NL) THEN
    !CASE 2-2-3 
    NUM(7)=NCENTER
    CALL GET_CELLID(I,J,K-1,NUM(1),0)
    CALL GET_CELLID(I,1,3*NL+1-J,NUM(2),1)
    CALL GET_CELLID(I+1,J,K,NUM(3),0)
    CALL GET_CELLID(I-1,J,K,NUM(4),0)
    CALL GET_CELLID(I,J+1,K,NUM(5),0)
    CALL GET_CELLID(I,J-1,K,NUM(6),0)
    CALL BUBBLE_SORT(NUM,INDEX_S)
    DO I1=1,7
    IF (NUM(I1).GT.NCENTER) THEN  
    IF_NEIGH=0
    CALL GET_INDEX(NUM(I1),I_NEIGH,J_NEIGH,K_NEIGH,IF_CY)
    IF (IF_CY.EQ.1) THEN
    IP_NEIGH=I_NEIGH-CNI*(MI-1)
    JP_NEIGH=J_NEIGH-CNJ*(MJ-1)
    KP_NEIGH=K_NEIGH
    IF (IP_NEIGH.LE.CNI.AND.IP_NEIGH.GE.1) THEN
    IF (JP_NEIGH.LE.CNJ.AND.JP_NEIGH.GE.1) THEN
        IF_NEIGH=1
    ENDIF
    ENDIF 
    ELSE
    IP_NEIGH=I_NEIGH-CNI*(MI-1)
    JP_NEIGH=J_NEIGH-CNL*(MJ1-1)
    KP_NEIGH=K_NEIGH-CNL*(MJ2-1)
    IF (IP_NEIGH.LE.CNI.AND.IP_NEIGH.GE.1) THEN
    IF (JP_NEIGH.LE.CNL.AND.JP_NEIGH.GE.1) THEN
    IF (KP_NEIGH.LE.CNL.AND.KP_NEIGH.GE.1) THEN
        IF_NEIGH=1 
    ENDIF
    ENDIF
    ENDIF    
    ENDIF
    
    IF (IF_CY.EQ.0) THEN
    IF (INDEX_S(I1).EQ.1) THEN
        FACE_ID=FACE_ID+1
        IF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.1) THEN
        WRITE(3,231) FACE_ID
        ELSEIF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.0) THEN
        F_FRONT_RECT(IP,JP)=FACE_ID
        ELSEIF (IF_OWN.EQ.0.AND.IF_NEIGH.EQ.1) THEN
        F_BACK_RECT(IP_NEIGH,JP_NEIGH)=-FACE_ID
        ENDIF
    ENDIF
    IF (INDEX_S(I1).EQ.2) THEN
        FACE_ID=FACE_ID+1
        IF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.1) THEN
        WRITE(3,231) FACE_ID
        ELSEIF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.0) THEN
        F_BACK_RECT(IP,JP)=FACE_ID
        ELSEIF (IF_OWN.EQ.0.AND.IF_NEIGH.EQ.1) THEN
        F_FRONT_RECT(IP_NEIGH,JP_NEIGH)=-FACE_ID
        ENDIF
    ENDIF
    IF (INDEX_S(I1).EQ.3) THEN
        FACE_ID=FACE_ID+1
        IF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.1) THEN
        WRITE(3,231) FACE_ID
        ELSEIF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.0) THEN
        F_LEFT_RECT(JP,KP)=FACE_ID
        ELSEIF (IF_OWN.EQ.0.AND.IF_NEIGH.EQ.1) THEN
        F_RIGHT_RECT(JP_NEIGH,KP_NEIGH)=-FACE_ID
        ENDIF
    ENDIF  
    IF (INDEX_S(I1).EQ.4) THEN
        FACE_ID=FACE_ID+1
        IF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.1) THEN
        WRITE(3,231) FACE_ID
        ELSEIF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.0) THEN
        F_RIGHT_RECT(JP,KP)=FACE_ID
        ELSEIF (IF_OWN.EQ.0.AND.IF_NEIGH.EQ.1) THEN
        F_LEFT_RECT(JP_NEIGH,KP_NEIGH)=-FACE_ID
        ENDIF
    ENDIF
    IF (INDEX_S(I1).EQ.5) THEN
        FACE_ID=FACE_ID+1
        IF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.1) THEN
        WRITE(3,231) FACE_ID
        ELSEIF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.0) THEN
        F_UP_RECT(IP,KP)=FACE_ID
        ELSEIF (IF_OWN.EQ.0.AND.IF_NEIGH.EQ.1) THEN
        F_DOWN_RECT(IP_NEIGH,KP_NEIGH)=-FACE_ID
        ENDIF
    ENDIF
    IF (INDEX_S(I1).EQ.6) THEN
        FACE_ID=FACE_ID+1
        IF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.1) THEN
        WRITE(3,231) FACE_ID
        ELSEIF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.0) THEN
        F_DOWN_RECT(IP,KP)=FACE_ID
        ELSEIF (IF_OWN.EQ.0.AND.IF_NEIGH.EQ.1) THEN
        F_UP_RECT(IP_NEIGH,KP_NEIGH)=-FACE_ID
        ENDIF
    ENDIF       
     
    ELSE
    IF (INDEX_S(I1).EQ.2) THEN
        FACE_ID=FACE_ID+1
        IF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.0) THEN
        F_BACK_RECT(IP,JP)=FACE_ID
        ELSEIF (IF_OWN.EQ.0.AND.IF_NEIGH.EQ.1) THEN
        F_DOWN_CY(IP_NEIGH,KP_NEIGH)=-FACE_ID
        ENDIF
    ENDIF
        
    ENDIF
    
    ENDIF
    ENDDO    
    
    ELSE
    !CASE 2-2-2
    NUM(7)=NCENTER
    CALL GET_CELLID(I,J,K-1,NUM(1),0)
    CALL GET_CELLID(I,J,K+1,NUM(2),0)
    CALL GET_CELLID(I+1,J,K,NUM(3),0)
    CALL GET_CELLID(I-1,J,K,NUM(4),0)
    CALL GET_CELLID(I,J+1,K,NUM(5),0)
    CALL GET_CELLID(I,J-1,K,NUM(6),0)
    CALL BUBBLE_SORT(NUM,INDEX_S)
    DO I1=1,7
    IF (NUM(I1).GT.NCENTER) THEN  
    IF_NEIGH=0
    CALL GET_INDEX(NUM(I1),I_NEIGH,J_NEIGH,K_NEIGH,IF_CY)
    IF (IF_CY.EQ.1) THEN
    IP_NEIGH=I_NEIGH-CNI*(MI-1)
    JP_NEIGH=J_NEIGH-CNJ*(MJ-1)
    KP_NEIGH=K_NEIGH
    IF (IP_NEIGH.LE.CNI.AND.IP_NEIGH.GE.1) THEN
    IF (JP_NEIGH.LE.CNJ.AND.JP_NEIGH.GE.1) THEN
        IF_NEIGH=1
    ENDIF
    ENDIF 
    ELSE
    IP_NEIGH=I_NEIGH-CNI*(MI-1)
    JP_NEIGH=J_NEIGH-CNL*(MJ1-1)
    KP_NEIGH=K_NEIGH-CNL*(MJ2-1)
    IF (IP_NEIGH.LE.CNI.AND.IP_NEIGH.GE.1) THEN
    IF (JP_NEIGH.LE.CNL.AND.JP_NEIGH.GE.1) THEN
    IF (KP_NEIGH.LE.CNL.AND.KP_NEIGH.GE.1) THEN
        IF_NEIGH=1 
    ENDIF
    ENDIF
    ENDIF    
    ENDIF
    
    IF (INDEX_S(I1).EQ.1) THEN
        FACE_ID=FACE_ID+1
        IF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.1) THEN
        WRITE(3,231) FACE_ID
        ELSEIF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.0) THEN
        F_FRONT_RECT(IP,JP)=FACE_ID
        ELSEIF (IF_OWN.EQ.0.AND.IF_NEIGH.EQ.1) THEN
        F_BACK_RECT(IP_NEIGH,JP_NEIGH)=-FACE_ID
        ENDIF
    ENDIF
    IF (INDEX_S(I1).EQ.2) THEN
        FACE_ID=FACE_ID+1
        IF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.1) THEN
        WRITE(3,231) FACE_ID
        ELSEIF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.0) THEN
        F_BACK_RECT(IP,JP)=FACE_ID
        ELSEIF (IF_OWN.EQ.0.AND.IF_NEIGH.EQ.1) THEN
        F_FRONT_RECT(IP_NEIGH,JP_NEIGH)=-FACE_ID
        ENDIF
    ENDIF
    IF (INDEX_S(I1).EQ.3) THEN
        FACE_ID=FACE_ID+1
        IF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.1) THEN
        WRITE(3,231) FACE_ID
        ELSEIF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.0) THEN
        F_LEFT_RECT(JP,KP)=FACE_ID
        ELSEIF (IF_OWN.EQ.0.AND.IF_NEIGH.EQ.1) THEN
        F_RIGHT_RECT(JP_NEIGH,KP_NEIGH)=-FACE_ID
        ENDIF
    ENDIF  
    IF (INDEX_S(I1).EQ.4) THEN
        FACE_ID=FACE_ID+1
        IF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.1) THEN
        WRITE(3,231) FACE_ID
        ELSEIF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.0) THEN
        F_RIGHT_RECT(JP,KP)=FACE_ID
        ELSEIF (IF_OWN.EQ.0.AND.IF_NEIGH.EQ.1) THEN
        F_LEFT_RECT(JP_NEIGH,KP_NEIGH)=-FACE_ID
        ENDIF
    ENDIF
    IF (INDEX_S(I1).EQ.5) THEN
        FACE_ID=FACE_ID+1
        IF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.1) THEN
        WRITE(3,231) FACE_ID
        ELSEIF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.0) THEN
        F_UP_RECT(IP,KP)=FACE_ID
        ELSEIF (IF_OWN.EQ.0.AND.IF_NEIGH.EQ.1) THEN
        F_DOWN_RECT(IP_NEIGH,KP_NEIGH)=-FACE_ID
        ENDIF
    ENDIF
    IF (INDEX_S(I1).EQ.6) THEN
        FACE_ID=FACE_ID+1
        IF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.1) THEN
        WRITE(3,231) FACE_ID
        ELSEIF (IF_OWN.EQ.1.AND.IF_NEIGH.EQ.0) THEN
        F_DOWN_RECT(IP,KP)=FACE_ID
        ELSEIF (IF_OWN.EQ.0.AND.IF_NEIGH.EQ.1) THEN
        F_UP_RECT(IP_NEIGH,KP_NEIGH)=-FACE_ID
        ENDIF
    ENDIF       
    
    ENDIF
    ENDDO    
    
    ENDIF    
    ENDIF    
    ENDIF      
        
    ENDIF
    
    
    ENDDO

    
    
    If ( iMPI_MyID.EQ.0 ) THEN
    write(*,*) 'test 04'
    ENDIF
    
    
    
    
!****************************************BOUNDARY FILES*****************************************    
    
    IF (MI.EQ.1) THEN 
!****************************************NI-START***************************************** 
    IF (MJ.EQ.1) THEN
    !CASE CYLINDER 1-1
        
        !RIGHT FACE--INLET2  !CASE CYLINDER 1-1
    IP=1
    I=IP+CNI*(MI-1)
    DO KP=1,NK
    K=KP
    DO JP=1,CNJ
    J=JP+CNJ*(MJ-1)
    F11=(CNJ+1)*(CNI+1)*(KP-1)+(CNI+1)*(JP-1)+IP-1
    F12=F11+1
    F13=F12+CNI+1
    F14=F13-1
    F21=F11+(CNJ+1)*(CNI+1)
    IF (KP.EQ.NK) F21=(CNI+1)*(JP-1)+IP-1
    F22=F21+1
    F23=F22+CNI+1
    F24=F23-1
    CALL CAL_FIJ_CY(F11,F12,F13,F14,F21,F22,F23,F24)
    WRITE(2,235) '4(',F11,' ',F21,' ',F24,' ',F14,')'
    FAD1=F_INTERNAL_NUM_ALL+NL*NL+NJ*(K-1)+J
    WRITE(3,231) FAD1
    OWN=CELL_ID_CY(IP,JP,KP)
    WRITE(4,231) OWN
    ENDDO   
    ENDDO
    
        !LEFT FACE--PROCESSOR_MI1_TO_MI2  !CASE CYLINDER 1-1
    IP=CNI
    I=IP+CNI*(MI-1)
    DO KP=1,NK
    K=KP
    DO JP=1,CNJ
    J=JP+CNJ*(MJ-1)
    F11=(CNJ+1)*(CNI+1)*(KP-1)+(CNI+1)*(JP-1)+IP-1
    F12=F11+1
    F13=F12+CNI+1
    F14=F13-1
    F21=F11+(CNJ+1)*(CNI+1)
    IF (KP.EQ.NK) F21=(CNI+1)*(JP-1)+IP-1
    F22=F21+1
    F23=F22+CNI+1
    F24=F23-1
    CALL CAL_FIJ_CY(F11,F12,F13,F14,F21,F22,F23,F24)
    WRITE(2,235) '4(',F12,' ',F13,' ',F23,' ',F22,')'
    FAD1=F_LEFT_CY(JP,KP)
    WRITE(3,231) FAD1
    OWN=CELL_ID_CY(IP,JP,KP)
    WRITE(4,231) OWN
    ENDDO   
    ENDDO
    
        !UP FACE--PROCESSOR_MI1_TO_MI2  !CASE CYLINDER 1-1
    JP=CNJ
    J=JP+CNJ*(MJ-1)
    DO KP=1,NK
    K=KP
    DO IP=1,CNI
    I=IP+CNI*(MI-1)
    F11=(CNJ+1)*(CNI+1)*(KP-1)+(CNI+1)*(JP-1)+IP-1
    F12=F11+1
    F13=F12+CNI+1
    F14=F13-1
    F21=F11+(CNJ+1)*(CNI+1)
    IF (KP.EQ.NK) F21=(CNI+1)*(JP-1)+IP-1
    F22=F21+1
    F23=F22+CNI+1
    F24=F23-1
    CALL CAL_FIJ_CY(F11,F12,F13,F14,F21,F22,F23,F24)
    WRITE(2,235) '4(',F14,' ',F24,' ',F23,' ',F13,')'
    FAD1=F_UP_CY(IP,KP)
    WRITE(3,231) FAD1
    OWN=CELL_ID_CY(IP,JP,KP)
    WRITE(4,231) OWN
    ENDDO   
    ENDDO
    
        !DOWN FACE--PROCESSOR_MI1_TO_MI2(RECTANGLE 1-1)  !CASE CYLINDER 1-1
    JP=1
    J=JP+CNJ*(MJ-1)
    DO KP=1,CNL
    K=KP
    DO IP=1,CNI
    I=IP+CNI*(MI-1)
    F11=(CNJ+1)*(CNI+1)*(KP-1)+(CNI+1)*(JP-1)+IP-1
    F12=F11+1
    F13=F12+CNI+1
    F14=F13-1
    F21=F11+(CNJ+1)*(CNI+1)
    IF (KP.EQ.NK) F21=(CNI+1)*(JP-1)+IP-1
    F22=F21+1
    F23=F22+CNI+1
    F24=F23-1
    CALL CAL_FIJ_CY(F11,F12,F13,F14,F21,F22,F23,F24)
    WRITE(2,235) '4(',F11,' ',F12,' ',F22,' ',F21,')'
    FAD1=F_DOWN_CY(IP,KP)
    WRITE(3,231) FAD1
    OWN=CELL_ID_CY(IP,JP,KP)
    WRITE(4,231) OWN
    ENDDO   
    ENDDO
    DO KP=NK-CNL+1,NK
    K=KP
    DO IP=1,CNI
    I=IP+CNI*(MI-1)
    F11=(CNJ+1)*(CNI+1)*(KP-1)+(CNI+1)*(JP-1)+IP-1
    F12=F11+1
    F13=F12+CNI+1
    F14=F13-1
    F21=F11+(CNJ+1)*(CNI+1)
    IF (KP.EQ.NK) F21=(CNI+1)*(JP-1)+IP-1
    F22=F21+1
    F23=F22+CNI+1
    F24=F23-1
    CALL CAL_FIJ_CY(F11,F12,F13,F14,F21,F22,F23,F24)
    WRITE(2,235) '4(',F11,' ',F12,' ',F22,' ',F21,')'
    FAD1=F_DOWN_CY(IP,KP)
    WRITE(3,231) FAD1
    OWN=CELL_ID_CY(IP,JP,KP)
    WRITE(4,231) OWN
    ENDDO   
    ENDDO
    
        !DOWN FACE--PROCESSOR_MI1_TO_MI2(RECTANGLE 2-1)  !CASE CYLINDER 1-1
    JP=1
    J=JP+CNJ*(MJ-1)
    DO NN=2,PNL-1
    DO KP=CNL*(NN-1)+1,CNL*NN
    K=KP
    DO IP=1,CNI
    I=IP+CNI*(MI-1)
    F11=(CNJ+1)*(CNI+1)*(KP-1)+(CNI+1)*(JP-1)+IP-1
    F12=F11+1
    F13=F12+CNI+1
    F14=F13-1
    F21=F11+(CNJ+1)*(CNI+1)
    IF (KP.EQ.NK) F21=(CNI+1)*(JP-1)+IP-1
    F22=F21+1
    F23=F22+CNI+1
    F24=F23-1
    CALL CAL_FIJ_CY(F11,F12,F13,F14,F21,F22,F23,F24)
    WRITE(2,235) '4(',F11,' ',F12,' ',F22,' ',F21,')'
    FAD1=F_DOWN_CY(IP,KP)
    WRITE(3,231) FAD1
    OWN=CELL_ID_CY(IP,JP,KP)
    WRITE(4,231) OWN
    ENDDO   
    ENDDO
    ENDDO
    
        !DOWN FACE--PROCESSOR_MI1_TO_MI2(RECTANGLE 3-1)  !CASE CYLINDER 1-1 
    JP=1
    J=JP+CNJ*(MJ-1)
    DO KP=NL-CNL+1,NL
    K=KP
    DO IP=1,CNI
    I=IP+CNI*(MI-1)
    F11=(CNJ+1)*(CNI+1)*(KP-1)+(CNI+1)*(JP-1)+IP-1
    F12=F11+1
    F13=F12+CNI+1
    F14=F13-1
    F21=F11+(CNJ+1)*(CNI+1)
    IF (KP.EQ.NK) F21=(CNI+1)*(JP-1)+IP-1
    F22=F21+1
    F23=F22+CNI+1
    F24=F23-1
    CALL CAL_FIJ_CY(F11,F12,F13,F14,F21,F22,F23,F24)
    WRITE(2,235) '4(',F11,' ',F12,' ',F22,' ',F21,')'
    FAD1=F_DOWN_CY(IP,KP)
    WRITE(3,231) FAD1
    OWN=CELL_ID_CY(IP,JP,KP)
    WRITE(4,231) OWN
    ENDDO   
    ENDDO
    DO KP=NL+1,NL+CNL
    K=KP
    DO IP=1,CNI
    I=IP+CNI*(MI-1)
    F11=(CNJ+1)*(CNI+1)*(KP-1)+(CNI+1)*(JP-1)+IP-1
    F12=F11+1
    F13=F12+CNI+1
    F14=F13-1
    F21=F11+(CNJ+1)*(CNI+1)
    IF (KP.EQ.NK) F21=(CNI+1)*(JP-1)+IP-1
    F22=F21+1
    F23=F22+CNI+1
    F24=F23-1
    CALL CAL_FIJ_CY(F11,F12,F13,F14,F21,F22,F23,F24)
    WRITE(2,235) '4(',F11,' ',F12,' ',F22,' ',F21,')'
    FAD1=F_DOWN_CY(IP,KP)
    WRITE(3,231) FAD1
    OWN=CELL_ID_CY(IP,JP,KP)
    WRITE(4,231) OWN
    ENDDO   
    ENDDO
        
        !DOWN FACE--PROCESSOR_MI1_TO_MI2(RECTANGLE 1-2 / 3-2)  !CASE CYLINDER 1-1  
    JP=1
    J=JP+CNJ*(MJ-1)
    DO NN=2,PNL-1
        
    DO KP=3*NL+CNL*(PNL-NN)+1,3*NL+CNL*(PNL-NN+1)
    K=KP
    DO IP=1,CNI
    I=IP+CNI*(MI-1)
    F11=(CNJ+1)*(CNI+1)*(KP-1)+(CNI+1)*(JP-1)+IP-1
    F12=F11+1
    F13=F12+CNI+1
    F14=F13-1
    F21=F11+(CNJ+1)*(CNI+1)
    IF (KP.EQ.NK) F21=(CNI+1)*(JP-1)+IP-1
    F22=F21+1
    F23=F22+CNI+1
    F24=F23-1
    CALL CAL_FIJ_CY(F11,F12,F13,F14,F21,F22,F23,F24)
    WRITE(2,235) '4(',F11,' ',F12,' ',F22,' ',F21,')'
    FAD1=F_DOWN_CY(IP,KP)
    WRITE(3,231) FAD1
    OWN=CELL_ID_CY(IP,JP,KP)
    WRITE(4,231) OWN
    ENDDO   
    ENDDO
    
    DO KP=NL+CNL*(NN-1)+1,NL+CNL*NN
    K=KP
    DO IP=1,CNI
    I=IP+CNI*(MI-1)
    F11=(CNJ+1)*(CNI+1)*(KP-1)+(CNI+1)*(JP-1)+IP-1
    F12=F11+1
    F13=F12+CNI+1
    F14=F13-1
    F21=F11+(CNJ+1)*(CNI+1)
    IF (KP.EQ.NK) F21=(CNI+1)*(JP-1)+IP-1
    F22=F21+1
    F23=F22+CNI+1
    F24=F23-1
    CALL CAL_FIJ_CY(F11,F12,F13,F14,F21,F22,F23,F24)
    WRITE(2,235) '4(',F11,' ',F12,' ',F22,' ',F21,')'
    FAD1=F_DOWN_CY(IP,KP)
    WRITE(3,231) FAD1
    OWN=CELL_ID_CY(IP,JP,KP)
    WRITE(4,231) OWN
    ENDDO   
    ENDDO
    
    ENDDO    
        
        !DOWN FACE--PROCESSOR_MI1_TO_MI2(RECTANGLE 1-3)  !CASE CYLINDER 1-1 
    JP=1
    J=JP+CNJ*(MJ-1)
    DO KP=3*NL-CNL+1,3*NL
    K=KP
    DO IP=1,CNI
    I=IP+CNI*(MI-1)
    F11=(CNJ+1)*(CNI+1)*(KP-1)+(CNI+1)*(JP-1)+IP-1
    F12=F11+1
    F13=F12+CNI+1
    F14=F13-1
    F21=F11+(CNJ+1)*(CNI+1)
    IF (KP.EQ.NK) F21=(CNI+1)*(JP-1)+IP-1
    F22=F21+1
    F23=F22+CNI+1
    F24=F23-1
    CALL CAL_FIJ_CY(F11,F12,F13,F14,F21,F22,F23,F24)
    WRITE(2,235) '4(',F11,' ',F12,' ',F22,' ',F21,')'
    FAD1=F_DOWN_CY(IP,KP)
    WRITE(3,231) FAD1
    OWN=CELL_ID_CY(IP,JP,KP)
    WRITE(4,231) OWN
    ENDDO   
    ENDDO
    DO KP=3*NL+1,3*NL+CNL
    K=KP
    DO IP=1,CNI
    I=IP+CNI*(MI-1)
    F11=(CNJ+1)*(CNI+1)*(KP-1)+(CNI+1)*(JP-1)+IP-1
    F12=F11+1
    F13=F12+CNI+1
    F14=F13-1
    F21=F11+(CNJ+1)*(CNI+1)
    IF (KP.EQ.NK) F21=(CNI+1)*(JP-1)+IP-1
    F22=F21+1
    F23=F22+CNI+1
    F24=F23-1
    CALL CAL_FIJ_CY(F11,F12,F13,F14,F21,F22,F23,F24)
    WRITE(2,235) '4(',F11,' ',F12,' ',F22,' ',F21,')'
    FAD1=F_DOWN_CY(IP,KP)
    WRITE(3,231) FAD1
    OWN=CELL_ID_CY(IP,JP,KP)
    WRITE(4,231) OWN
    ENDDO   
    ENDDO
    
        !DOWN FACE--PROCESSOR_MI1_TO_MI2(RECTANGLE 2-3)  !CASE CYLINDER 1-1  
    JP=1
    J=JP+CNJ*(MJ-1)
    DO NN=2,PNL-1
    DO KP=2*NL+CNL*(PNL-NN)+1,2*NL+CNL*(PNL-NN+1)
    K=KP
    DO IP=1,CNI
    I=IP+CNI*(MI-1)
    F11=(CNJ+1)*(CNI+1)*(KP-1)+(CNI+1)*(JP-1)+IP-1
    F12=F11+1
    F13=F12+CNI+1
    F14=F13-1
    F21=F11+(CNJ+1)*(CNI+1)
    IF (KP.EQ.NK) F21=(CNI+1)*(JP-1)+IP-1
    F22=F21+1
    F23=F22+CNI+1
    F24=F23-1
    CALL CAL_FIJ_CY(F11,F12,F13,F14,F21,F22,F23,F24)
    WRITE(2,235) '4(',F11,' ',F12,' ',F22,' ',F21,')'
    FAD1=F_DOWN_CY(IP,KP)
    WRITE(3,231) FAD1
    OWN=CELL_ID_CY(IP,JP,KP)
    WRITE(4,231) OWN
    ENDDO   
    ENDDO
    ENDDO
    
        !DOWN FACE--PROCESSOR_MI1_TO_MI2(RECTANGLE 3-3)  !CASE CYLINDER 1-1 
    JP=1
    J=JP+CNJ*(MJ-1)
    DO KP=2*NL-CNL+1,2*NL
    K=KP
    DO IP=1,CNI
    I=IP+CNI*(MI-1)
    F11=(CNJ+1)*(CNI+1)*(KP-1)+(CNI+1)*(JP-1)+IP-1
    F12=F11+1
    F13=F12+CNI+1
    F14=F13-1
    F21=F11+(CNJ+1)*(CNI+1)
    IF (KP.EQ.NK) F21=(CNI+1)*(JP-1)+IP-1
    F22=F21+1
    F23=F22+CNI+1
    F24=F23-1
    CALL CAL_FIJ_CY(F11,F12,F13,F14,F21,F22,F23,F24)
    WRITE(2,235) '4(',F11,' ',F12,' ',F22,' ',F21,')'
    FAD1=F_DOWN_CY(IP,KP)
    WRITE(3,231) FAD1
    OWN=CELL_ID_CY(IP,JP,KP)
    WRITE(4,231) OWN
    ENDDO   
    ENDDO
    DO KP=2*NL+1,2*NL+CNL
    K=KP
    DO IP=1,CNI
    I=IP+CNI*(MI-1)
    F11=(CNJ+1)*(CNI+1)*(KP-1)+(CNI+1)*(JP-1)+IP-1
    F12=F11+1
    F13=F12+CNI+1
    F14=F13-1
    F21=F11+(CNJ+1)*(CNI+1)
    IF (KP.EQ.NK) F21=(CNI+1)*(JP-1)+IP-1
    F22=F21+1
    F23=F22+CNI+1
    F24=F23-1
    CALL CAL_FIJ_CY(F11,F12,F13,F14,F21,F22,F23,F24)
    WRITE(2,235) '4(',F11,' ',F12,' ',F22,' ',F21,')'
    FAD1=F_DOWN_CY(IP,KP)
    WRITE(3,231) FAD1
    OWN=CELL_ID_CY(IP,JP,KP)
    WRITE(4,231) OWN
    ENDDO   
    ENDDO
    
    N=2+4*(PNL-1)
    WRITE(31,231) N+4
    WRITE(31,*) '('  
    WRITE(31,*) '0'
    WRITE(31,*) '1'
    WRITE(31,*) '2' 
    WRITE(31,*) '3' 
    DO I=1,N
    WRITE(31,*) '-1' 
    ENDDO
    WRITE(31,*) ')'
    
    WRITE(21,231) N+4
    WRITE(21,*) "("    
    WRITE(21,*) "    INLET1"
    WRITE(21,*) "    {"    
    WRITE(21,*) "        type            patch;"
    N=0
    WRITE(21,232) "        nFaces          ",N,";"
    N=F_INTERNAL_NUM
    WRITE(21,232) "        startFace       ",N,";"
    WRITE(21,*) "    }"
    WRITE(21,*) "    INLET2"
    WRITE(21,*) "    {"    
    WRITE(21,*) "        type            patch;"
    N=CNJ*NK
    WRITE(21,232) "        nFaces          ",N,";"
    N=F_INTERNAL_NUM
    WRITE(21,232) "        startFace       ",N,";"
    WRITE(21,*) "    }"
    WRITE(21,*) "    OUTLET"
    WRITE(21,*) "    {"
    WRITE(21,*) "        type            patch;"
    N=0
    WRITE(21,232) "        nFaces          ",N,";"
    N=F_INTERNAL_NUM+CNJ*NK
    WRITE(21,232) "        startFace       ",N,";"
    WRITE(21,*) "    }"
    WRITE(21,*) "    HULL"
    WRITE(21,*) "    {"
    WRITE(21,*) "        type            wall;"
    WRITE(21,*) "        inGroups        1(wall);"
    N=0
    WRITE(21,232) "        nFaces          ",N,";"
    N=F_INTERNAL_NUM+CNJ*NK
    WRITE(21,232) "        startFace       ",N,";"
    WRITE(21,*) "    }"
        !LEFT FACE--PROCESSOR_MI1_TO_MI2  !CASE CYLINDER 1-1 
    PRO1=PNI*(MJ-1)+MI-1
    PRO2=PNI*(MJ-1)+MI
    WRITE(21,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(21,*) '    {'
    WRITE(21,*) '        type            processor;'
    WRITE(21,*) '        inGroups        List<word> 1(processor);'
    N=CNJ*NK
    WRITE(21,232) "        nFaces          ",N,";"
    N=F_INTERNAL_NUM+CNJ*NK
    WRITE(21,232) "        startFace       ",N,";"
    WRITE(21,*) '        matchTolerance  0.0001;'
    WRITE(21,232) "        myProcNo        ",PRO1,";"
    WRITE(21,232) "        neighbProcNo    ",PRO2,";"
    WRITE(21,*) '    }'
    WRITE(41,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(41,*) '    {'
    WRITE(41,*) '        type            processor;'
    WRITE(41,*) '        value           nonuniform List<scalar> '
    N=CNJ*NK
    WRITE(41,'(I0)') N
    WRITE(41,*) '('
    DO KP=1,NK
    DO JP=1,CNJ
    WRITE(41,'(E21.14)') P_LEFT_CY(JP,KP)    
    ENDDO
    ENDDO
    WRITE(41,*) ')'
    WRITE(41,*) ';'
    WRITE(41,*) '    }'
    WRITE(51,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(51,*) '    {'
    WRITE(51,*) '        type            processor;'
    WRITE(51,*) '        value           nonuniform List<vector> '
    N=CNJ*NK
    WRITE(51,'(I0)') N
    WRITE(51,*) '('
    DO KP=1,NK
    DO JP=1,CNJ
    UU1=U_LEFT_CY(1,JP,KP)
    UU2=U_LEFT_CY(2,JP,KP)
    UU3=U_LEFT_CY(3,JP,KP)
    WRITE(51,230) '( ',UU1,' ',UU2,' ',UU3,' )'  
    ENDDO
    ENDDO
    WRITE(51,*) ')'
    WRITE(51,*) ';'
    WRITE(51,*) '    }'
    WRITE(61,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(61,*) '    {'
    WRITE(61,*) '        type            processor;'
    WRITE(61,*) "        value           uniform 0;"
    WRITE(61,*) '    }'
        !UP FACE--PROCESSOR_MI1_TO_MI2  !CASE CYLINDER 1-1
    PRO2=PNI*MJ+MI-1
    WRITE(21,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(21,*) '    {'
    WRITE(21,*) '        type            processor;'
    WRITE(21,*) '        inGroups        List<word> 1(processor);'
    N=CNI*NK
    WRITE(21,232) "        nFaces          ",N,";"
    N=F_INTERNAL_NUM+CNJ*NK+CNJ*NK
    WRITE(21,232) "        startFace       ",N,";"
    WRITE(21,*) '        matchTolerance  0.0001;'
    WRITE(21,232) "        myProcNo        ",PRO1,";"
    WRITE(21,232) "        neighbProcNo    ",PRO2,";"
    WRITE(21,*) '    }'
    WRITE(41,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(41,*) '    {'
    WRITE(41,*) '        type            processor;'
    WRITE(41,*) '        value           nonuniform List<scalar> '
    N=CNI*NK
    WRITE(41,'(I0)') N
    WRITE(41,*) '('
    DO KP=1,NK
    DO IP=1,CNI
    WRITE(41,'(E21.14)') P_LEFT_CY(IP,KP)    
    ENDDO
    ENDDO
    WRITE(41,*) ')'
    WRITE(41,*) ';'
    WRITE(41,*) '    }'
    WRITE(51,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(51,*) '    {'
    WRITE(51,*) '        type            processor;'
    WRITE(51,*) '        value           nonuniform List<vector> '
    N=CNI*NK
    WRITE(51,'(I0)') N
    WRITE(51,*) '('
    DO KP=1,NK
    DO IP=1,CNI
    UU1=U_UP_CY(1,IP,KP)
    UU2=U_UP_CY(2,IP,KP)
    UU3=U_UP_CY(3,IP,KP)
    WRITE(51,230) '( ',UU1,' ',UU2,' ',UU3,' )'  
    ENDDO
    ENDDO
    WRITE(51,*) ')'
    WRITE(51,*) ';'
    WRITE(51,*) '    }'
    WRITE(61,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(61,*) '    {'
    WRITE(61,*) '        type            processor;'
    WRITE(61,*) "        value           uniform 0;"
    WRITE(61,*) '    }'
        !DOWN FACE--PROCESSOR_MI1_TO_MI2(RECTANGLE 1-1)  !CASE CYLINDER 1-1
    MJ1=1
    MJ2=1
    PRO2=PNI*PNJ+PNL*PNI*(MJ2-1)+PNI*(MJ1-1)+MI-1
    WRITE(21,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(21,*) '    {'
    WRITE(21,*) '        type            processor;'
    WRITE(21,*) '        inGroups        List<word> 1(processor);'
    N=CNI*CNL*2
    WRITE(21,232) "        nFaces          ",N,";"
    N=F_INTERNAL_NUM+CNJ*NK+CNJ*NK+CNI*NK
    WRITE(21,232) "        startFace       ",N,";"
    WRITE(21,*) '        matchTolerance  0.0001;'
    WRITE(21,232) "        myProcNo        ",PRO1,";"
    WRITE(21,232) "        neighbProcNo    ",PRO2,";"
    WRITE(21,*) '    }'
    WRITE(41,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(41,*) '    {'
    WRITE(41,*) '        type            processor;'
    WRITE(41,*) '        value           nonuniform List<scalar> '
    N=CNI*CNL*2
    WRITE(41,'(I0)') N
    WRITE(41,*) '('
    DO KP=1,CNL
    DO IP=1,CNI
    WRITE(41,'(E21.14)') P_DOWN_CY(IP,KP)    
    ENDDO
    ENDDO
    DO KP=NK-CNL+1,NK
    DO IP=1,CNI
    WRITE(41,'(E21.14)') P_DOWN_CY(IP,KP)    
    ENDDO
    ENDDO
    WRITE(41,*) ')'
    WRITE(41,*) ';'
    WRITE(41,*) '    }'
    WRITE(51,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(51,*) '    {'
    WRITE(51,*) '        type            processor;'
    WRITE(51,*) '        value           nonuniform List<vector> '
    N=CNI*CNL*2
    WRITE(51,'(I0)') N
    WRITE(51,*) '('
    DO KP=1,CNL
    DO IP=1,CNI
    UU1=U_DOWN_CY(1,IP,KP)
    UU2=U_DOWN_CY(2,IP,KP)
    UU3=U_DOWN_CY(3,IP,KP)
    WRITE(51,230) '( ',UU1,' ',UU2,' ',UU3,' )'  
    ENDDO
    ENDDO
    DO KP=NK-CNL+1,NK
    DO IP=1,CNI
    UU1=U_DOWN_CY(1,IP,KP)
    UU2=U_DOWN_CY(2,IP,KP)
    UU3=U_DOWN_CY(3,IP,KP)
    WRITE(51,230) '( ',UU1,' ',UU2,' ',UU3,' )'  
    ENDDO
    ENDDO
    WRITE(51,*) ')'
    WRITE(51,*) ';'
    WRITE(51,*) '    }'
    WRITE(61,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(61,*) '    {'
    WRITE(61,*) '        type            processor;'
    WRITE(61,*) "        value           uniform 0;"
    WRITE(61,*) '    }'
    
        !DOWN FACE--PROCESSOR_MI1_TO_MI2(RECTANGLE 2-1)  !CASE CYLINDER 1-1
    DO MJ1=2,PNL-1
    MJ2=1
    PRO2=PNI*PNJ+PNL*PNI*(MJ2-1)+PNI*(MJ1-1)+MI-1
    WRITE(21,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(21,*) '    {'
    WRITE(21,*) '        type            processor;'
    WRITE(21,*) '        inGroups        List<word> 1(processor);'
    N=CNI*CNL
    WRITE(21,232) "        nFaces          ",N,";"
    N=F_INTERNAL_NUM+CNJ*NK+CNJ*NK+CNI*NK+CNI*CNL*2+(MJ1-2)*CNI*CNL
    WRITE(21,232) "        startFace       ",N,";"
    WRITE(21,*) '        matchTolerance  0.0001;'
    WRITE(21,232) "        myProcNo        ",PRO1,";"
    WRITE(21,232) "        neighbProcNo    ",PRO2,";"
    WRITE(21,*) '    }'
    WRITE(41,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(41,*) '    {'
    WRITE(41,*) '        type            processor;'
    WRITE(41,*) '        value           nonuniform List<scalar> '
    N=CNI*CNL
    WRITE(41,'(I0)') N
    WRITE(41,*) '('
    DO KP=CNL*(MJ1-1)+1,CNL*MJ1
    DO IP=1,CNI
    WRITE(41,'(E21.14)') P_DOWN_CY(IP,KP)    
    ENDDO
    ENDDO
    WRITE(41,*) ')'
    WRITE(41,*) ';'
    WRITE(41,*) '    }'
    WRITE(51,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(51,*) '    {'
    WRITE(51,*) '        type            processor;'
    WRITE(51,*) '        value           nonuniform List<vector> '
    N=CNI*CNL
    WRITE(51,'(I0)') N
    WRITE(51,*) '('
    DO KP=CNL*(MJ1-1)+1,CNL*MJ1
    DO IP=1,CNI
    UU1=U_DOWN_CY(1,IP,KP)
    UU2=U_DOWN_CY(2,IP,KP)
    UU3=U_DOWN_CY(3,IP,KP)
    WRITE(51,230) '( ',UU1,' ',UU2,' ',UU3,' )'  
    ENDDO
    ENDDO
    WRITE(51,*) ')'
    WRITE(51,*) ';'
    WRITE(51,*) '    }'
    WRITE(61,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(61,*) '    {'
    WRITE(61,*) '        type            processor;'
    WRITE(61,*) "        value           uniform 0;"
    WRITE(61,*) '    }'    
    ENDDO
    
        !DOWN FACE--PROCESSOR_MI1_TO_MI2(RECTANGLE 3-1)  !CASE CYLINDER 1-1
    MJ1=PNL
    MJ2=1
    PRO2=PNI*PNJ+PNL*PNI*(MJ2-1)+PNI*(MJ1-1)+MI-1
    WRITE(21,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(21,*) '    {'
    WRITE(21,*) '        type            processor;'
    WRITE(21,*) '        inGroups        List<word> 1(processor);'
    N=CNI*CNL*2
    WRITE(21,232) "        nFaces          ",N,";"
    N=F_INTERNAL_NUM+CNJ*NK+CNJ*NK+CNI*NK+CNI*CNL*2+(PNL-2)*CNI*CNL
    WRITE(21,232) "        startFace       ",N,";"
    WRITE(21,*) '        matchTolerance  0.0001;'
    WRITE(21,232) "        myProcNo        ",PRO1,";"
    WRITE(21,232) "        neighbProcNo    ",PRO2,";"
    WRITE(21,*) '    }'
    WRITE(41,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(41,*) '    {'
    WRITE(41,*) '        type            processor;'
    WRITE(41,*) '        value           nonuniform List<scalar> '
    N=CNI*CNL*2
    WRITE(41,'(I0)') N
    WRITE(41,*) '('
    DO KP=NL-CNL+1,NL
    DO IP=1,CNI
    WRITE(41,'(E21.14)') P_DOWN_CY(IP,KP)    
    ENDDO
    ENDDO
    DO KP=NL+1,NL+CNL
    DO IP=1,CNI
    WRITE(41,'(E21.14)') P_DOWN_CY(IP,KP)    
    ENDDO
    ENDDO
    WRITE(41,*) ')'
    WRITE(41,*) ';'
    WRITE(41,*) '    }'
    WRITE(51,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(51,*) '    {'
    WRITE(51,*) '        type            processor;'
    WRITE(51,*) '        value           nonuniform List<vector> '
    N=CNI*CNL*2
    WRITE(51,'(I0)') N
    WRITE(51,*) '('
    DO KP=NL-CNL+1,NL
    DO IP=1,CNI
    UU1=U_DOWN_CY(1,IP,KP)
    UU2=U_DOWN_CY(2,IP,KP)
    UU3=U_DOWN_CY(3,IP,KP)
    WRITE(51,230) '( ',UU1,' ',UU2,' ',UU3,' )'  
    ENDDO
    ENDDO
    DO KP=NL+1,NL+CNL
    DO IP=1,CNI
    UU1=U_DOWN_CY(1,IP,KP)
    UU2=U_DOWN_CY(2,IP,KP)
    UU3=U_DOWN_CY(3,IP,KP)
    WRITE(51,230) '( ',UU1,' ',UU2,' ',UU3,' )'  
    ENDDO
    ENDDO
    WRITE(51,*) ')'
    WRITE(51,*) ';'
    WRITE(51,*) '    }'
    WRITE(61,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(61,*) '    {'
    WRITE(61,*) '        type            processor;'
    WRITE(61,*) "        value           uniform 0;"
    WRITE(61,*) '    }'
    
        !DOWN FACE--PROCESSOR_MI1_TO_MI2(RECTANGLE 1-2 / 3-2)  !CASE CYLINDER 1-1
    DO MJ2=2,PNL-1
    MJ1=1
    PRO2=PNI*PNJ+PNL*PNI*(MJ2-1)+PNI*(MJ1-1)+MI-1
    WRITE(21,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(21,*) '    {'
    WRITE(21,*) '        type            processor;'
    WRITE(21,*) '        inGroups        List<word> 1(processor);'
    N=CNI*CNL
    WRITE(21,232) "        nFaces          ",N,";"
    N=F_INTERNAL_NUM+CNJ*NK+CNJ*NK+CNI*NK+CNI*CNL*2+(PNL-2)*CNI*CNL+CNI*CNL*2+(MJ2-2)*CNI*CNL*2
    WRITE(21,232) "        startFace       ",N,";"
    WRITE(21,*) '        matchTolerance  0.0001;'
    WRITE(21,232) "        myProcNo        ",PRO1,";"
    WRITE(21,232) "        neighbProcNo    ",PRO2,";"
    WRITE(21,*) '    }'
    WRITE(41,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(41,*) '    {'
    WRITE(41,*) '        type            processor;'
    WRITE(41,*) '        value           nonuniform List<scalar> '
    N=CNI*CNL
    WRITE(41,'(I0)') N
    WRITE(41,*) '('
    DO KP=3*NL+CNL*(PNL-MJ2)+1,3*NL+CNL*(PNL-MJ2+1)
    DO IP=1,CNI
    WRITE(41,'(E21.14)') P_DOWN_CY(IP,KP)    
    ENDDO
    ENDDO
    WRITE(41,*) ')'
    WRITE(41,*) ';'
    WRITE(41,*) '    }'
    WRITE(51,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(51,*) '    {'
    WRITE(51,*) '        type            processor;'
    WRITE(51,*) '        value           nonuniform List<vector> '
    N=CNI*CNL
    WRITE(51,'(I0)') N
    WRITE(51,*) '('
    DO KP=3*NL+CNL*(PNL-MJ2)+1,3*NL+CNL*(PNL-MJ2+1)
    DO IP=1,CNI
    UU1=U_DOWN_CY(1,IP,KP)
    UU2=U_DOWN_CY(2,IP,KP)
    UU3=U_DOWN_CY(3,IP,KP)
    WRITE(51,230) '( ',UU1,' ',UU2,' ',UU3,' )'  
    ENDDO
    ENDDO
    WRITE(51,*) ')'
    WRITE(51,*) ';'
    WRITE(51,*) '    }'
    WRITE(61,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(61,*) '    {'
    WRITE(61,*) '        type            processor;'
    WRITE(61,*) "        value           uniform 0;"
    WRITE(61,*) '    }' 
    
    MJ1=PNL
    PRO2=PNI*PNJ+PNL*PNI*(MJ2-1)+PNI*(MJ1-1)+MI-1
    WRITE(21,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(21,*) '    {'
    WRITE(21,*) '        type            processor;'
    WRITE(21,*) '        inGroups        List<word> 1(processor);'
    N=CNI*CNL
    WRITE(21,232) "        nFaces          ",N,";"
    N=F_INTERNAL_NUM+CNJ*NK+CNJ*NK+CNI*NK+CNI*CNL*2+(PNL-2)*CNI*CNL+CNI*CNL*2+(MJ2-2)*CNI*CNL*2+CNI*CNL
    WRITE(21,232) "        startFace       ",N,";"
    WRITE(21,*) '        matchTolerance  0.0001;'
    WRITE(21,232) "        myProcNo        ",PRO1,";"
    WRITE(21,232) "        neighbProcNo    ",PRO2,";"
    WRITE(21,*) '    }'
    WRITE(41,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(41,*) '    {'
    WRITE(41,*) '        type            processor;'
    WRITE(41,*) '        value           nonuniform List<scalar> '
    N=CNI*CNL
    WRITE(41,'(I0)') N
    WRITE(41,*) '('
    DO KP=NL+CNL*(MJ2-1)+1,NL+CNL*MJ2
    DO IP=1,CNI
    WRITE(41,'(E21.14)') P_DOWN_CY(IP,KP)    
    ENDDO
    ENDDO
    WRITE(41,*) ')'
    WRITE(41,*) ';'
    WRITE(41,*) '    }'
    WRITE(51,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(51,*) '    {'
    WRITE(51,*) '        type            processor;'
    WRITE(51,*) '        value           nonuniform List<vector> '
    N=CNI*CNL
    WRITE(51,'(I0)') N
    WRITE(51,*) '('
    DO KP=NL+CNL*(MJ2-1)+1,NL+CNL*MJ2
    DO IP=1,CNI
    UU1=U_DOWN_CY(1,IP,KP)
    UU2=U_DOWN_CY(2,IP,KP)
    UU3=U_DOWN_CY(3,IP,KP)
    WRITE(51,230) '( ',UU1,' ',UU2,' ',UU3,' )'  
    ENDDO
    ENDDO
    WRITE(51,*) ')'
    WRITE(51,*) ';'
    WRITE(51,*) '    }'
    WRITE(61,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(61,*) '    {'
    WRITE(61,*) '        type            processor;'
    WRITE(61,*) "        value           uniform 0;"
    WRITE(61,*) '    }'
        
    ENDDO
    
    
        !DOWN FACE--PROCESSOR_MI1_TO_MI2(RECTANGLE 1-3)  !CASE CYLINDER 1-1
    MJ1=1
    MJ2=PNL
    PRO2=PNI*PNJ+PNL*PNI*(MJ2-1)+PNI*(MJ1-1)+MI-1
    WRITE(21,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(21,*) '    {'
    WRITE(21,*) '        type            processor;'
    WRITE(21,*) '        inGroups        List<word> 1(processor);'
    N=CNI*CNL*2
    WRITE(21,232) "        nFaces          ",N,";"
    N=F_INTERNAL_NUM+CNJ*NK+CNJ*NK+CNI*NK+CNI*CNL*2+(PNL-2)*CNI*CNL+CNI*CNL*2+(PNL-2)*CNI*CNL*2
    WRITE(21,232) "        startFace       ",N,";"
    WRITE(21,*) '        matchTolerance  0.0001;'
    WRITE(21,232) "        myProcNo        ",PRO1,";"
    WRITE(21,232) "        neighbProcNo    ",PRO2,";"
    WRITE(21,*) '    }'
    WRITE(41,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(41,*) '    {'
    WRITE(41,*) '        type            processor;'
    WRITE(41,*) '        value           nonuniform List<scalar> '
    N=CNI*CNL*2
    WRITE(41,'(I0)') N
    WRITE(41,*) '('
    DO KP=3*NL-CNL+1,3*NL
    DO IP=1,CNI
    WRITE(41,'(E21.14)') P_DOWN_CY(IP,KP)    
    ENDDO
    ENDDO
    DO KP=3*NL+1,3*NL+CNL
    DO IP=1,CNI
    WRITE(41,'(E21.14)') P_DOWN_CY(IP,KP)    
    ENDDO
    ENDDO
    WRITE(41,*) ')'
    WRITE(41,*) ';'
    WRITE(41,*) '    }'
    WRITE(51,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(51,*) '    {'
    WRITE(51,*) '        type            processor;'
    WRITE(51,*) '        value           nonuniform List<vector> '
    N=CNI*CNL*2
    WRITE(51,'(I0)') N
    WRITE(51,*) '('
    DO KP=3*NL-CNL+1,3*NL
    DO IP=1,CNI
    UU1=U_DOWN_CY(1,IP,KP)
    UU2=U_DOWN_CY(2,IP,KP)
    UU3=U_DOWN_CY(3,IP,KP)
    WRITE(51,230) '( ',UU1,' ',UU2,' ',UU3,' )'  
    ENDDO
    ENDDO
    DO KP=3*NL+1,3*NL+CNL
    DO IP=1,CNI
    UU1=U_DOWN_CY(1,IP,KP)
    UU2=U_DOWN_CY(2,IP,KP)
    UU3=U_DOWN_CY(3,IP,KP)
    WRITE(51,230) '( ',UU1,' ',UU2,' ',UU3,' )'  
    ENDDO
    ENDDO
    WRITE(51,*) ')'
    WRITE(51,*) ';'
    WRITE(51,*) '    }'
    WRITE(61,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(61,*) '    {'
    WRITE(61,*) '        type            processor;'
    WRITE(61,*) "        value           uniform 0;"
    WRITE(61,*) '    }'
    
    
        !DOWN FACE--PROCESSOR_MI1_TO_MI2(RECTANGLE 2-3)  !CASE CYLINDER 1-1
    DO MJ1=2,PNL-1
    MJ2=PNL
    PRO2=PNI*PNJ+PNL*PNI*(MJ2-1)+PNI*(MJ1-1)+MI-1
    WRITE(21,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(21,*) '    {'
    WRITE(21,*) '        type            processor;'
    WRITE(21,*) '        inGroups        List<word> 1(processor);'
    N=CNI*CNL
    WRITE(21,232) "        nFaces          ",N,";"
    N=F_INTERNAL_NUM+CNJ*NK+CNJ*NK+CNI*NK+CNI*CNL*2+(PNL-2)*CNI*CNL+CNI*CNL*2+(PNL-2)*CNI*CNL*2+CNI*CNL*2+(MJ1-2)*CNI*CNL
    WRITE(21,232) "        startFace       ",N,";"
    WRITE(21,*) '        matchTolerance  0.0001;'
    WRITE(21,232) "        myProcNo        ",PRO1,";"
    WRITE(21,232) "        neighbProcNo    ",PRO2,";"
    WRITE(21,*) '    }'
    WRITE(41,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(41,*) '    {'
    WRITE(41,*) '        type            processor;'
    WRITE(41,*) '        value           nonuniform List<scalar> '
    N=CNI*CNL
    WRITE(41,'(I0)') N
    WRITE(41,*) '('
    DO KP=2*NL+CNL*(PNL-MJ1)+1,2*NL+CNL*(PNL-MJ1+1)
    DO IP=1,CNI
    WRITE(41,'(E21.14)') P_DOWN_CY(IP,KP)    
    ENDDO
    ENDDO
    WRITE(41,*) ')'
    WRITE(41,*) ';'
    WRITE(41,*) '    }'
    WRITE(51,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(51,*) '    {'
    WRITE(51,*) '        type            processor;'
    WRITE(51,*) '        value           nonuniform List<vector> '
    N=CNI*CNL
    WRITE(51,'(I0)') N
    WRITE(51,*) '('
    DO KP=2*NL+CNL*(PNL-MJ1)+1,2*NL+CNL*(PNL-MJ1+1)
    DO IP=1,CNI
    UU1=U_DOWN_CY(1,IP,KP)
    UU2=U_DOWN_CY(2,IP,KP)
    UU3=U_DOWN_CY(3,IP,KP)
    WRITE(51,230) '( ',UU1,' ',UU2,' ',UU3,' )'  
    ENDDO
    ENDDO
    WRITE(51,*) ')'
    WRITE(51,*) ';'
    WRITE(51,*) '    }'
    WRITE(61,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(61,*) '    {'
    WRITE(61,*) '        type            processor;'
    WRITE(61,*) "        value           uniform 0;"
    WRITE(61,*) '    }'    
    ENDDO
    
    
    
        !DOWN FACE--PROCESSOR_MI1_TO_MI2(RECTANGLE 3-3)  !CASE CYLINDER 1-1
    MJ1=PNL
    MJ2=PNL
    PRO2=PNI*PNJ+PNL*PNI*(MJ2-1)+PNI*(MJ1-1)+MI-1
    WRITE(21,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(21,*) '    {'
    WRITE(21,*) '        type            processor;'
    WRITE(21,*) '        inGroups        List<word> 1(processor);'
    N=CNI*CNL*2
    WRITE(21,232) "        nFaces          ",N,";"
    N=F_INTERNAL_NUM+CNJ*NK+CNJ*NK+CNI*NK+CNI*CNL*2+(PNL-2)*CNI*CNL+CNI*CNL*2+(PNL-2)*CNI*CNL*2+CNI*CNL*2+(PNL-2)*CNI*CNL
    WRITE(21,232) "        startFace       ",N,";"
    WRITE(21,*) '        matchTolerance  0.0001;'
    WRITE(21,232) "        myProcNo        ",PRO1,";"
    WRITE(21,232) "        neighbProcNo    ",PRO2,";"
    WRITE(21,*) '    }'
    WRITE(41,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(41,*) '    {'
    WRITE(41,*) '        type            processor;'
    WRITE(41,*) '        value           nonuniform List<scalar> '
    N=CNI*CNL*2
    WRITE(41,'(I0)') N
    WRITE(41,*) '('
    DO KP=2*NL-CNL+1,2*NL
    DO IP=1,CNI
    WRITE(41,'(E21.14)') P_DOWN_CY(IP,KP)    
    ENDDO
    ENDDO
    DO KP=2*NL+1,2*NL+CNL
    DO IP=1,CNI
    WRITE(41,'(E21.14)') P_DOWN_CY(IP,KP)    
    ENDDO
    ENDDO
    WRITE(41,*) ')'
    WRITE(41,*) ';'
    WRITE(41,*) '    }'
    WRITE(51,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(51,*) '    {'
    WRITE(51,*) '        type            processor;'
    WRITE(51,*) '        value           nonuniform List<vector> '
    N=CNI*CNL*2
    WRITE(51,'(I0)') N
    WRITE(51,*) '('
    DO KP=2*NL-CNL+1,2*NL
    DO IP=1,CNI
    UU1=U_DOWN_CY(1,IP,KP)
    UU2=U_DOWN_CY(2,IP,KP)
    UU3=U_DOWN_CY(3,IP,KP)
    WRITE(51,230) '( ',UU1,' ',UU2,' ',UU3,' )'  
    ENDDO
    ENDDO
    DO KP=2*NL+1,2*NL+CNL
    DO IP=1,CNI
    UU1=U_DOWN_CY(1,IP,KP)
    UU2=U_DOWN_CY(2,IP,KP)
    UU3=U_DOWN_CY(3,IP,KP)
    WRITE(51,230) '( ',UU1,' ',UU2,' ',UU3,' )'  
    ENDDO
    ENDDO
    WRITE(51,*) ')'
    WRITE(51,*) ';'
    WRITE(51,*) '    }'
    WRITE(61,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(61,*) '    {'
    WRITE(61,*) '        type            processor;'
    WRITE(61,*) "        value           uniform 0;"
    WRITE(61,*) '    }'
    
    
    
    
    
    
    
        
    ELSEIF (MJ.GT.1 .AND. MJ.LT.PNJ) THEN
    !CASE CYLINDER 1-2
        
        !RIGHT FACE--INLET2  !CASE CYLINDER 1-2
    IP=1
    I=IP+CNI*(MI-1)
    DO KP=1,NK
    K=KP
    DO JP=1,CNJ
    J=JP+CNJ*(MJ-1)
    F11=(CNJ+1)*(CNI+1)*(KP-1)+(CNI+1)*(JP-1)+IP-1
    F12=F11+1
    F13=F12+CNI+1
    F14=F13-1
    F21=F11+(CNJ+1)*(CNI+1)
    IF (KP.EQ.NK) F21=(CNI+1)*(JP-1)+IP-1
    F22=F21+1
    F23=F22+CNI+1
    F24=F23-1
    CALL CAL_FIJ_CY(F11,F12,F13,F14,F21,F22,F23,F24)
    WRITE(2,235) '4(',F11,' ',F21,' ',F24,' ',F14,')'
    FAD1=F_INTERNAL_NUM_ALL+NL*NL+NJ*(K-1)+J
    WRITE(3,231) FAD1
    OWN=CELL_ID_CY(IP,JP,KP)
    WRITE(4,231) OWN
    ENDDO   
    ENDDO  
    
        !DOWN FACE--PROCESSOR_MI1_TO_MI2 !CASE CYLINDER 1-2 
    JP=1
    J=JP+CNJ*(MJ-1)
    DO KP=1,NK
    K=KP
    DO IP=1,CNI
    I=IP+CNI*(MI-1)
    F11=(CNJ+1)*(CNI+1)*(KP-1)+(CNI+1)*(JP-1)+IP-1
    F12=F11+1
    F13=F12+CNI+1
    F14=F13-1
    F21=F11+(CNJ+1)*(CNI+1)
    IF (KP.EQ.NK) F21=(CNI+1)*(JP-1)+IP-1
    F22=F21+1
    F23=F22+CNI+1
    F24=F23-1
    CALL CAL_FIJ_CY(F11,F12,F13,F14,F21,F22,F23,F24)
    WRITE(2,235) '4(',F11,' ',F12,' ',F22,' ',F21,')'
    FAD1=F_DOWN_CY(IP,KP)
    WRITE(3,231) FAD1
    OWN=CELL_ID_CY(IP,JP,KP)
    WRITE(4,231) OWN
    ENDDO   
    ENDDO    
        
        !LEFT FACE--PROCESSOR_MI1_TO_MI2  !CASE CYLINDER 1-2
    IP=CNI
    I=IP+CNI*(MI-1)
    DO KP=1,NK
    K=KP
    DO JP=1,CNJ
    J=JP+CNJ*(MJ-1)
    F11=(CNJ+1)*(CNI+1)*(KP-1)+(CNI+1)*(JP-1)+IP-1
    F12=F11+1
    F13=F12+CNI+1
    F14=F13-1
    F21=F11+(CNJ+1)*(CNI+1)
    IF (KP.EQ.NK) F21=(CNI+1)*(JP-1)+IP-1
    F22=F21+1
    F23=F22+CNI+1
    F24=F23-1
    CALL CAL_FIJ_CY(F11,F12,F13,F14,F21,F22,F23,F24)
    WRITE(2,235) '4(',F12,' ',F13,' ',F23,' ',F22,')'
    FAD1=F_LEFT_CY(JP,KP)
    WRITE(3,231) FAD1
    OWN=CELL_ID_CY(IP,JP,KP)
    WRITE(4,231) OWN
    ENDDO   
    ENDDO
    
        !UP FACE--PROCESSOR_MI1_TO_MI2  !CASE CYLINDER 1-2
    JP=CNJ
    J=JP+CNJ*(MJ-1)
    DO KP=1,NK
    K=KP
    DO IP=1,CNI
    I=IP+CNI*(MI-1)
    F11=(CNJ+1)*(CNI+1)*(KP-1)+(CNI+1)*(JP-1)+IP-1
    F12=F11+1
    F13=F12+CNI+1
    F14=F13-1
    F21=F11+(CNJ+1)*(CNI+1)
    IF (KP.EQ.NK) F21=(CNI+1)*(JP-1)+IP-1
    F22=F21+1
    F23=F22+CNI+1
    F24=F23-1
    CALL CAL_FIJ_CY(F11,F12,F13,F14,F21,F22,F23,F24)
    WRITE(2,235) '4(',F14,' ',F24,' ',F23,' ',F13,')'
    FAD1=F_UP_CY(IP,KP)
    WRITE(3,231) FAD1
    OWN=CELL_ID_CY(IP,JP,KP)
    WRITE(4,231) OWN
    ENDDO   
    ENDDO    
        
    N=3
    WRITE(31,231) N+4
    WRITE(31,*) '('  
    WRITE(31,*) '0'
    WRITE(31,*) '1'
    WRITE(31,*) '2' 
    WRITE(31,*) '3' 
    DO I=1,N
    WRITE(31,*) '-1' 
    ENDDO
    WRITE(31,*) ')'
    
    WRITE(21,231) N+4
    WRITE(21,*) "("    
    WRITE(21,*) "    INLET1"
    WRITE(21,*) "    {"    
    WRITE(21,*) "        type            patch;"
    N=0
    WRITE(21,232) "        nFaces          ",N,";"
    N=F_INTERNAL_NUM
    WRITE(21,232) "        startFace       ",N,";"
    WRITE(21,*) "    }"
    WRITE(21,*) "    INLET2"
    WRITE(21,*) "    {"    
    WRITE(21,*) "        type            patch;"
    N=CNJ*NK
    WRITE(21,232) "        nFaces          ",N,";"
    N=F_INTERNAL_NUM
    WRITE(21,232) "        startFace       ",N,";"
    WRITE(21,*) "    }"
    WRITE(21,*) "    OUTLET"
    WRITE(21,*) "    {"
    WRITE(21,*) "        type            patch;"
    N=0
    WRITE(21,232) "        nFaces          ",N,";"
    N=F_INTERNAL_NUM+CNJ*NK
    WRITE(21,232) "        startFace       ",N,";"
    WRITE(21,*) "    }"    
    WRITE(21,*) "    HULL"
    WRITE(21,*) "    {"
    WRITE(21,*) "        type            wall;"
    WRITE(21,*) "        inGroups        1(wall);"
    N=0
    WRITE(21,232) "        nFaces          ",N,";"
    N=F_INTERNAL_NUM+CNJ*NK
    WRITE(21,232) "        startFace       ",N,";"
    WRITE(21,*) "    }"    
        
        !DOWN FACE--PROCESSOR_MI1_TO_MI2 !CASE CYLINDER 1-2 
    PRO1=PNI*(MJ-1)+MI-1
    PRO2=PNI*(MJ-2)+MI-1
    WRITE(21,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(21,*) '    {'
    WRITE(21,*) '        type            processor;'
    WRITE(21,*) '        inGroups        List<word> 1(processor);'
    N=CNI*NK
    WRITE(21,232) "        nFaces          ",N,";"
    N=F_INTERNAL_NUM+CNJ*NK
    WRITE(21,232) "        startFace       ",N,";"
    WRITE(21,*) '        matchTolerance  0.0001;'
    WRITE(21,232) "        myProcNo        ",PRO1,";"
    WRITE(21,232) "        neighbProcNo    ",PRO2,";"
    WRITE(21,*) '    }'
    WRITE(41,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(41,*) '    {'
    WRITE(41,*) '        type            processor;'
    WRITE(41,*) '        value           nonuniform List<scalar> '
    N=CNI*NK
    WRITE(41,'(I0)') N
    WRITE(41,*) '('
    DO KP=1,NK
    DO IP=1,CNI
    WRITE(41,'(E21.14)') P_DOWN_CY(IP,KP)    
    ENDDO
    ENDDO
    WRITE(41,*) ')'
    WRITE(41,*) ';'
    WRITE(41,*) '    }'
    WRITE(51,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(51,*) '    {'
    WRITE(51,*) '        type            processor;'
    WRITE(51,*) '        value           nonuniform List<vector> '
    N=CNI*NK
    WRITE(51,'(I0)') N
    WRITE(51,*) '('
    DO KP=1,NK
    DO IP=1,CNI
    UU1=U_DOWN_CY(1,IP,KP)
    UU2=U_DOWN_CY(2,IP,KP)
    UU3=U_DOWN_CY(3,IP,KP)
    WRITE(51,230) '( ',UU1,' ',UU2,' ',UU3,' )'  
    ENDDO
    ENDDO
    WRITE(51,*) ')'
    WRITE(51,*) ';'
    WRITE(51,*) '    }'
    WRITE(61,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(61,*) '    {'
    WRITE(61,*) '        type            processor;'
    WRITE(61,*) "        value           uniform 0;"
    WRITE(61,*) '    }'
    
        !LEFT FACE--PROCESSOR_MI1_TO_MI2  !CASE CYLINDER 1-2
    PRO2=PNI*(MJ-1)+MI
    WRITE(21,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(21,*) '    {'
    WRITE(21,*) '        type            processor;'
    WRITE(21,*) '        inGroups        List<word> 1(processor);'
    N=CNJ*NK
    WRITE(21,232) "        nFaces          ",N,";"
    N=F_INTERNAL_NUM+CNJ*NK+CNI*NK
    WRITE(21,232) "        startFace       ",N,";"
    WRITE(21,*) '        matchTolerance  0.0001;'
    WRITE(21,232) "        myProcNo        ",PRO1,";"
    WRITE(21,232) "        neighbProcNo    ",PRO2,";"
    WRITE(21,*) '    }'
    WRITE(41,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(41,*) '    {'
    WRITE(41,*) '        type            processor;'
    WRITE(41,*) '        value           nonuniform List<scalar> '
    N=CNJ*NK
    WRITE(41,'(I0)') N
    WRITE(41,*) '('
    DO KP=1,NK
    DO JP=1,CNJ
    WRITE(41,'(E21.14)') P_LEFT_CY(JP,KP)    
    ENDDO
    ENDDO
    WRITE(41,*) ')'
    WRITE(41,*) ';'
    WRITE(41,*) '    }'
    WRITE(51,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(51,*) '    {'
    WRITE(51,*) '        type            processor;'
    WRITE(51,*) '        value           nonuniform List<vector> '
    N=CNJ*NK
    WRITE(51,'(I0)') N
    WRITE(51,*) '('
    DO KP=1,NK
    DO JP=1,CNJ
    UU1=U_LEFT_CY(1,JP,KP)
    UU2=U_LEFT_CY(2,JP,KP)
    UU3=U_LEFT_CY(3,JP,KP)
    WRITE(51,230) '( ',UU1,' ',UU2,' ',UU3,' )'  
    ENDDO
    ENDDO
    WRITE(51,*) ')'
    WRITE(51,*) ';'
    WRITE(51,*) '    }'
    WRITE(61,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(61,*) '    {'
    WRITE(61,*) '        type            processor;'
    WRITE(61,*) "        value           uniform 0;"
    WRITE(61,*) '    }'
    
        !UP FACE--PROCESSOR_MI1_TO_MI2  !CASE CYLINDER 1-2
    PRO2=PNI*MJ+MI-1
    WRITE(21,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(21,*) '    {'
    WRITE(21,*) '        type            processor;'
    WRITE(21,*) '        inGroups        List<word> 1(processor);'
    N=CNI*NK
    WRITE(21,232) "        nFaces          ",N,";"
    N=F_INTERNAL_NUM+CNJ*NK+CNI*NK+CNJ*NK
    WRITE(21,232) "        startFace       ",N,";"
    WRITE(21,*) '        matchTolerance  0.0001;'
    WRITE(21,232) "        myProcNo        ",PRO1,";"
    WRITE(21,232) "        neighbProcNo    ",PRO2,";"
    WRITE(21,*) '    }'
    WRITE(41,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(41,*) '    {'
    WRITE(41,*) '        type            processor;'
    WRITE(41,*) '        value           nonuniform List<scalar> '
    N=CNI*NK
    WRITE(41,'(I0)') N
    WRITE(41,*) '('
    DO KP=1,NK
    DO IP=1,CNI
    WRITE(41,'(E21.14)') P_UP_CY(IP,KP)    
    ENDDO
    ENDDO
    WRITE(41,*) ')'
    WRITE(41,*) ';'
    WRITE(41,*) '    }'
    WRITE(51,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(51,*) '    {'
    WRITE(51,*) '        type            processor;'
    WRITE(51,*) '        value           nonuniform List<vector> '
    N=CNI*NK
    WRITE(51,'(I0)') N
    WRITE(51,*) '('
    DO KP=1,NK
    DO IP=1,CNI
    UU1=U_UP_CY(1,IP,KP)
    UU2=U_UP_CY(2,IP,KP)
    UU3=U_UP_CY(3,IP,KP)
    WRITE(51,230) '( ',UU1,' ',UU2,' ',UU3,' )'  
    ENDDO
    ENDDO
    WRITE(51,*) ')'
    WRITE(51,*) ';'
    WRITE(51,*) '    }'
    WRITE(61,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(61,*) '    {'
    WRITE(61,*) '        type            processor;'
    WRITE(61,*) "        value           uniform 0;"
    WRITE(61,*) '    }'    
        
    
    ELSEIF (MJ.EQ.PNJ) THEN
    !CASE CYLINDER 1-3
        
        !RIGHT FACE--INLET2  !CASE CYLINDER 1-3
    IP=1
    I=IP+CNI*(MI-1)
    DO KP=1,NK
    K=KP
    DO JP=1,CNJ
    J=JP+CNJ*(MJ-1)
    F11=(CNJ+1)*(CNI+1)*(KP-1)+(CNI+1)*(JP-1)+IP-1
    F12=F11+1
    F13=F12+CNI+1
    F14=F13-1
    F21=F11+(CNJ+1)*(CNI+1)
    IF (KP.EQ.NK) F21=(CNI+1)*(JP-1)+IP-1
    F22=F21+1
    F23=F22+CNI+1
    F24=F23-1
    CALL CAL_FIJ_CY(F11,F12,F13,F14,F21,F22,F23,F24)
    WRITE(2,235) '4(',F11,' ',F21,' ',F24,' ',F14,')'
    FAD1=F_INTERNAL_NUM_ALL+NL*NL+NJ*(K-1)+J
    WRITE(3,231) FAD1
    OWN=CELL_ID_CY(IP,JP,KP)
    WRITE(4,231) OWN
    ENDDO   
    ENDDO
    
        !UP FACE--OUTLET  !CASE CYLINDER 1-3
    JP=CNJ
    J=JP+CNJ*(MJ-1)
    DO IP=1,CNI
    I=IP+CNI*(MI-1)
    DO KP=1,NL
    K=KP
    F11=(CNJ+1)*(CNI+1)*(KP-1)+(CNI+1)*(JP-1)+IP-1
    F12=F11+1
    F13=F12+CNI+1
    F14=F13-1
    F21=F11+(CNJ+1)*(CNI+1)
    IF (KP.EQ.NK) F21=(CNI+1)*(JP-1)+IP-1
    F22=F21+1
    F23=F22+CNI+1
    F24=F23-1
    CALL CAL_FIJ_CY(F11,F12,F13,F14,F21,F22,F23,F24)
    WRITE(2,235) '4(',F14,' ',F24,' ',F23,' ',F13,')'
    FAD1=F_INTERNAL_NUM_ALL+NL*NL+NJ*NK+NL*(I-1)+K
    WRITE(3,231) FAD1
    OWN=CELL_ID_CY(IP,JP,KP)
    WRITE(4,231) OWN
    ENDDO   
    ENDDO
    DO IP=1,CNI
    I=IP+CNI*(MI-1)
    DO KP=1+NL,2*NL
    K=KP
    F11=(CNJ+1)*(CNI+1)*(KP-1)+(CNI+1)*(JP-1)+IP-1
    F12=F11+1
    F13=F12+CNI+1
    F14=F13-1
    F21=F11+(CNJ+1)*(CNI+1)
    IF (KP.EQ.NK) F21=(CNI+1)*(JP-1)+IP-1
    F22=F21+1
    F23=F22+CNI+1
    F24=F23-1
    CALL CAL_FIJ_CY(F11,F12,F13,F14,F21,F22,F23,F24)
    WRITE(2,235) '4(',F14,' ',F24,' ',F23,' ',F13,')'
    FAD1=F_INTERNAL_NUM_ALL+NL*NL+NJ*NK+NI*NL+NL*(I-1)+(K-NL)
    WRITE(3,231) FAD1
    OWN=CELL_ID_CY(IP,JP,KP)
    WRITE(4,231) OWN
    ENDDO   
    ENDDO
    DO IP=1,CNI
    I=IP+CNI*(MI-1)
    DO KP=1+2*NL,3*NL
    K=KP
    F11=(CNJ+1)*(CNI+1)*(KP-1)+(CNI+1)*(JP-1)+IP-1
    F12=F11+1
    F13=F12+CNI+1
    F14=F13-1
    F21=F11+(CNJ+1)*(CNI+1)
    IF (KP.EQ.NK) F21=(CNI+1)*(JP-1)+IP-1
    F22=F21+1
    F23=F22+CNI+1
    F24=F23-1
    CALL CAL_FIJ_CY(F11,F12,F13,F14,F21,F22,F23,F24)
    WRITE(2,235) '4(',F14,' ',F24,' ',F23,' ',F13,')'
    FAD1=F_INTERNAL_NUM_ALL+NL*NL+NJ*NK+NI*NL*2+NL*(I-1)+(K-2*NL)
    WRITE(3,231) FAD1
    OWN=CELL_ID_CY(IP,JP,KP)
    WRITE(4,231) OWN
    ENDDO   
    ENDDO
    DO IP=1,CNI
    I=IP+CNI*(MI-1)
    DO KP=1+3*NL,NK
    K=KP
    F11=(CNJ+1)*(CNI+1)*(KP-1)+(CNI+1)*(JP-1)+IP-1
    F12=F11+1
    F13=F12+CNI+1
    F14=F13-1
    F21=F11+(CNJ+1)*(CNI+1)
    IF (KP.EQ.NK) F21=(CNI+1)*(JP-1)+IP-1
    F22=F21+1
    F23=F22+CNI+1
    F24=F23-1
    CALL CAL_FIJ_CY(F11,F12,F13,F14,F21,F22,F23,F24)
    WRITE(2,235) '4(',F14,' ',F24,' ',F23,' ',F13,')'
    FAD1=F_INTERNAL_NUM_ALL+NL*NL+NJ*NK+NI*NL*3+NL*(I-1)+(K-3*NL)
    WRITE(3,231) FAD1
    OWN=CELL_ID_CY(IP,JP,KP)
    WRITE(4,231) OWN
    ENDDO   
    ENDDO
    
    
        !DOWN FACE--PROCESSOR_MI1_TO_MI2 !CASE CYLINDER 1-3 
    JP=1
    J=JP+CNJ*(MJ-1)
    DO KP=1,NK
    K=KP
    DO IP=1,CNI
    I=IP+CNI*(MI-1)
    F11=(CNJ+1)*(CNI+1)*(KP-1)+(CNI+1)*(JP-1)+IP-1
    F12=F11+1
    F13=F12+CNI+1
    F14=F13-1
    F21=F11+(CNJ+1)*(CNI+1)
    IF (KP.EQ.NK) F21=(CNI+1)*(JP-1)+IP-1
    F22=F21+1
    F23=F22+CNI+1
    F24=F23-1
    CALL CAL_FIJ_CY(F11,F12,F13,F14,F21,F22,F23,F24)
    WRITE(2,235) '4(',F11,' ',F12,' ',F22,' ',F21,')'
    FAD1=F_DOWN_CY(IP,KP)
    WRITE(3,231) FAD1
    OWN=CELL_ID_CY(IP,JP,KP)
    WRITE(4,231) OWN
    ENDDO   
    ENDDO 
    
        !LEFT FACE--PROCESSOR_MI1_TO_MI2  !CASE CYLINDER 1-3
    IP=CNI
    I=IP+CNI*(MI-1)
    DO KP=1,NK
    K=KP
    DO JP=1,CNJ
    J=JP+CNJ*(MJ-1)
    F11=(CNJ+1)*(CNI+1)*(KP-1)+(CNI+1)*(JP-1)+IP-1
    F12=F11+1
    F13=F12+CNI+1
    F14=F13-1
    F21=F11+(CNJ+1)*(CNI+1)
    IF (KP.EQ.NK) F21=(CNI+1)*(JP-1)+IP-1
    F22=F21+1
    F23=F22+CNI+1
    F24=F23-1
    CALL CAL_FIJ_CY(F11,F12,F13,F14,F21,F22,F23,F24)
    WRITE(2,235) '4(',F12,' ',F13,' ',F23,' ',F22,')'
    FAD1=F_LEFT_CY(JP,KP)
    WRITE(3,231) FAD1
    OWN=CELL_ID_CY(IP,JP,KP)
    WRITE(4,231) OWN
    ENDDO   
    ENDDO
    
    N=2
    WRITE(31,231) N+4
    WRITE(31,*) '('  
    WRITE(31,*) '0'
    WRITE(31,*) '1'
    WRITE(31,*) '2' 
    WRITE(31,*) '3' 
    DO I=1,N
    WRITE(31,*) '-1' 
    ENDDO
    WRITE(31,*) ')'
    
    WRITE(21,231) N+4
    WRITE(21,*) "("    
    WRITE(21,*) "    INLET1"
    WRITE(21,*) "    {"    
    WRITE(21,*) "        type            patch;"
    N=0
    WRITE(21,232) "        nFaces          ",N,";"
    N=F_INTERNAL_NUM
    WRITE(21,232) "        startFace       ",N,";"
    WRITE(21,*) "    }"
    WRITE(21,*) "    INLET2"
    WRITE(21,*) "    {"    
    WRITE(21,*) "        type            patch;"
    N=CNJ*NK
    WRITE(21,232) "        nFaces          ",N,";"
    N=F_INTERNAL_NUM
    WRITE(21,232) "        startFace       ",N,";"
    WRITE(21,*) "    }"
    WRITE(21,*) "    OUTLET"
    WRITE(21,*) "    {"
    WRITE(21,*) "        type            patch;"
    N=CNI*NK
    WRITE(21,232) "        nFaces          ",N,";"
    N=F_INTERNAL_NUM+CNJ*NK
    WRITE(21,232) "        startFace       ",N,";"
    WRITE(21,*) "    }"
    WRITE(21,*) "    HULL"
    WRITE(21,*) "    {"
    WRITE(21,*) "        type            wall;"
    WRITE(21,*) "        inGroups        1(wall);"
    N=0
    WRITE(21,232) "        nFaces          ",N,";"
    N=F_INTERNAL_NUM+CNJ*NK+CNI*NK
    WRITE(21,232) "        startFace       ",N,";"
    WRITE(21,*) "    }"
    
        !DOWN FACE--PROCESSOR_MI1_TO_MI2 !CASE CYLINDER 1-3 
    PRO1=PNI*(MJ-1)+MI-1
    PRO2=PNI*(MJ-2)+MI-1
    WRITE(21,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(21,*) '    {'
    WRITE(21,*) '        type            processor;'
    WRITE(21,*) '        inGroups        List<word> 1(processor);'
    N=CNI*NK
    WRITE(21,232) "        nFaces          ",N,";"
    N=F_INTERNAL_NUM+CNJ*NK+CNI*NK
    WRITE(21,232) "        startFace       ",N,";"
    WRITE(21,*) '        matchTolerance  0.0001;'
    WRITE(21,232) "        myProcNo        ",PRO1,";"
    WRITE(21,232) "        neighbProcNo    ",PRO2,";"
    WRITE(21,*) '    }'
    WRITE(41,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(41,*) '    {'
    WRITE(41,*) '        type            processor;'
    WRITE(41,*) '        value           nonuniform List<scalar> '
    N=CNI*NK
    WRITE(41,'(I0)') N
    WRITE(41,*) '('
    DO KP=1,NK
    DO IP=1,CNI
    WRITE(41,'(E21.14)') P_DOWN_CY(IP,KP)    
    ENDDO
    ENDDO
    WRITE(41,*) ')'
    WRITE(41,*) ';'
    WRITE(41,*) '    }'
    WRITE(51,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(51,*) '    {'
    WRITE(51,*) '        type            processor;'
    WRITE(51,*) '        value           nonuniform List<vector> '
    N=CNI*NK
    WRITE(51,'(I0)') N
    WRITE(51,*) '('
    DO KP=1,NK
    DO IP=1,CNI
    UU1=U_DOWN_CY(1,IP,KP)
    UU2=U_DOWN_CY(2,IP,KP)
    UU3=U_DOWN_CY(3,IP,KP)
    WRITE(51,230) '( ',UU1,' ',UU2,' ',UU3,' )'  
    ENDDO
    ENDDO
    WRITE(51,*) ')'
    WRITE(51,*) ';'
    WRITE(51,*) '    }'
    WRITE(61,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(61,*) '    {'
    WRITE(61,*) '        type            processor;'
    WRITE(61,*) "        value           uniform 0;"
    WRITE(61,*) '    }'
    
        !LEFT FACE--PROCESSOR_MI1_TO_MI2  !CASE CYLINDER 1-3
    PRO2=PNI*(MJ-1)+MI
    WRITE(21,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(21,*) '    {'
    WRITE(21,*) '        type            processor;'
    WRITE(21,*) '        inGroups        List<word> 1(processor);'
    N=CNJ*NK
    WRITE(21,232) "        nFaces          ",N,";"
    N=F_INTERNAL_NUM+CNJ*NK+CNI*NK+CNI*NK
    WRITE(21,232) "        startFace       ",N,";"
    WRITE(21,*) '        matchTolerance  0.0001;'
    WRITE(21,232) "        myProcNo        ",PRO1,";"
    WRITE(21,232) "        neighbProcNo    ",PRO2,";"
    WRITE(21,*) '    }'
    WRITE(41,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(41,*) '    {'
    WRITE(41,*) '        type            processor;'
    WRITE(41,*) '        value           nonuniform List<scalar> '
    N=CNJ*NK
    WRITE(41,'(I0)') N
    WRITE(41,*) '('
    DO KP=1,NK
    DO JP=1,CNJ
    WRITE(41,'(E21.14)') P_LEFT_CY(JP,KP)    
    ENDDO
    ENDDO
    WRITE(41,*) ')'
    WRITE(41,*) ';'
    WRITE(41,*) '    }'
    WRITE(51,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(51,*) '    {'
    WRITE(51,*) '        type            processor;'
    WRITE(51,*) '        value           nonuniform List<vector> '
    N=CNJ*NK
    WRITE(51,'(I0)') N
    WRITE(51,*) '('
    DO KP=1,NK
    DO JP=1,CNJ
    UU1=U_LEFT_CY(1,JP,KP)
    UU2=U_LEFT_CY(2,JP,KP)
    UU3=U_LEFT_CY(3,JP,KP)
    WRITE(51,230) '( ',UU1,' ',UU2,' ',UU3,' )'  
    ENDDO
    ENDDO
    WRITE(51,*) ')'
    WRITE(51,*) ';'
    WRITE(51,*) '    }'
    WRITE(61,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(61,*) '    {'
    WRITE(61,*) '        type            processor;'
    WRITE(61,*) "        value           uniform 0;"
    WRITE(61,*) '    }'
    
    
    ELSE
    MJ1=( (COORD(1)-PNJ)-MOD(COORD(1)-PNJ,PNL) )/PNL+1
    MJ2=MJ-PNJ-(MJ1-1)*PNL
    IF (MJ1.EQ.1) THEN
    IF (MJ2.EQ.1) THEN
    !CASE RECTANGLE 1-1-1
        
        !RIGHT FACE--INLET1  !CASE RECTANGLE 1-1-1 
    IP=1
    I=IP+CNI*(MI-1)
    DO KP=1,CNL
    K=KP+CNL*(MJ2-1)
    DO JP=1,CNL
    J=JP+CNL*(MJ1-1)
    F11=(CNL+1)*(CNI+1)*(KP-1)+(CNI+1)*(JP-1)+IP-1
    F12=F11+1
    F13=F12+CNI+1
    F14=F13-1
    F21=F11+(CNL+1)*(CNI+1)
    F22=F21+1
    F23=F22+CNI+1
    F24=F23-1
    CALL CAL_FIJ_RECT(F11,F12,F13,F14,F21,F22,F23,F24)
    WRITE(2,235) '4(',F11,' ',F21,' ',F24,' ',F14,')'
    FAD1=F_INTERNAL_NUM_ALL+NL*(K-1)+J
    WRITE(3,231) FAD1
    OWN=CELL_ID_RECT(IP,JP,KP)
    WRITE(4,231) OWN
    ENDDO   
    ENDDO
    
        !FRONT FACE--PROCESSOR_MI1_TO_MI2  !CASE RECTANGLE 1-1-1
    KP=1
    K=KP+CNL*(MJ2-1)
    DO JP=1,CNL
    J=JP+CNL*(MJ1-1)
    DO IP=1,CNI
    I=IP+CNI*(MI-1)
    F11=(CNL+1)*(CNI+1)*(KP-1)+(CNI+1)*(JP-1)+IP-1
    F12=F11+1
    F13=F12+CNI+1
    F14=F13-1
    F21=F11+(CNL+1)*(CNI+1)
    F22=F21+1
    F23=F22+CNI+1
    F24=F23-1
    CALL CAL_FIJ_RECT(F11,F12,F13,F14,F21,F22,F23,F24)
    WRITE(2,235) '4(',F11,' ',F14,' ',F13,' ',F12,')'
    FAD1=F_FRONT_RECT(IP,JP)
    WRITE(3,231) FAD1
    OWN=CELL_ID_RECT(IP,JP,KP)
    WRITE(4,231) OWN
    ENDDO   
    ENDDO
    
        !DOWN FACE--PROCESSOR_MI1_TO_MI2  !CASE RECTANGLE 1-1-1
    JP=1
    J=JP+CNL*(MJ1-1)
    DO KP=CNL,1,-1
    K=KP+CNL*(MJ2-1)
    DO IP=1,CNI
    I=IP+CNI*(MI-1)
    F11=(CNL+1)*(CNI+1)*(KP-1)+(CNI+1)*(JP-1)+IP-1
    F12=F11+1
    F13=F12+CNI+1
    F14=F13-1
    F21=F11+(CNL+1)*(CNI+1)
    F22=F21+1
    F23=F22+CNI+1
    F24=F23-1
    CALL CAL_FIJ_RECT(F11,F12,F13,F14,F21,F22,F23,F24)
    !WRITE(2,235) '4(',F11,' ',F12,' ',F22,' ',F21,')'
    WRITE(2,235) '4(',F21,' ',F11,' ',F12,' ',F22,')'
    FAD1=F_DOWN_RECT(IP,KP)
    WRITE(3,231) FAD1
    OWN=CELL_ID_RECT(IP,JP,KP)
    WRITE(4,231) OWN
    ENDDO   
    ENDDO
    
        !LEFT FACE--PROCESSOR_MI1_TO_MI2  !CASE RECTANGLE 1-1-1 
    IP=CNI
    I=IP+CNI*(MI-1)
    DO KP=1,CNL
    K=KP+CNL*(MJ2-1)
    DO JP=1,CNL
    J=JP+CNL*(MJ1-1)
    F11=(CNL+1)*(CNI+1)*(KP-1)+(CNI+1)*(JP-1)+IP-1
    F12=F11+1
    F13=F12+CNI+1
    F14=F13-1
    F21=F11+(CNL+1)*(CNI+1)
    F22=F21+1
    F23=F22+CNI+1
    F24=F23-1
    CALL CAL_FIJ_RECT(F11,F12,F13,F14,F21,F22,F23,F24)
    WRITE(2,235) '4(',F12,' ',F13,' ',F23,' ',F22,')'
    FAD1=F_LEFT_RECT(JP,KP)
    WRITE(3,231) FAD1
    OWN=CELL_ID_RECT(IP,JP,KP)
    WRITE(4,231) OWN
    ENDDO   
    ENDDO
    
        !UP FACE--PROCESSOR_MI1_TO_MI2  !CASE RECTANGLE 1-1-1
    JP=CNL
    J=JP+CNL*(MJ1-1)
    DO KP=1,CNL
    K=KP+CNL*(MJ2-1)
    DO IP=1,CNI
    I=IP+CNI*(MI-1)
    F11=(CNL+1)*(CNI+1)*(KP-1)+(CNI+1)*(JP-1)+IP-1
    F12=F11+1
    F13=F12+CNI+1
    F14=F13-1
    F21=F11+(CNL+1)*(CNI+1)
    F22=F21+1
    F23=F22+CNI+1
    F24=F23-1
    CALL CAL_FIJ_RECT(F11,F12,F13,F14,F21,F22,F23,F24)
    WRITE(2,235) '4(',F14,' ',F24,' ',F23,' ',F13,')'
    FAD1=F_UP_RECT(IP,KP)
    WRITE(3,231) FAD1
    OWN=CELL_ID_RECT(IP,JP,KP)
    WRITE(4,231) OWN
    ENDDO   
    ENDDO
    
        !BACK FACE--PROCESSOR_MI1_TO_MI2  !CASE RECTANGLE 1-1-1
    KP=CNL
    K=KP+CNL*(MJ2-1)
    DO JP=1,CNL
    J=JP+CNL*(MJ1-1)
    DO IP=1,CNI
    I=IP+CNI*(MI-1)
    F11=(CNL+1)*(CNI+1)*(KP-1)+(CNI+1)*(JP-1)+IP-1
    F12=F11+1
    F13=F12+CNI+1
    F14=F13-1
    F21=F11+(CNL+1)*(CNI+1)
    F22=F21+1
    F23=F22+CNI+1
    F24=F23-1
    CALL CAL_FIJ_RECT(F11,F12,F13,F14,F21,F22,F23,F24)
    WRITE(2,235) '4(',F21,' ',F22,' ',F23,' ',F24,')'
    FAD1=F_BACK_RECT(IP,JP)
    WRITE(3,231) FAD1
    OWN=CELL_ID_RECT(IP,JP,KP)
    WRITE(4,231) OWN
    ENDDO   
    ENDDO
    
    N=4
    WRITE(31,231) N+4
    WRITE(31,*) '('  
    WRITE(31,*) '0'
    WRITE(31,*) '1'
    WRITE(31,*) '2' 
    WRITE(31,*) '3' 
    DO I=1,N
    WRITE(31,*) '-1' 
    ENDDO
    WRITE(31,*) ')'
    
    WRITE(21,231) N+4
    WRITE(21,*) "("    
    WRITE(21,*) "    INLET1"
    WRITE(21,*) "    {"    
    WRITE(21,*) "        type            patch;"
    N=CNL*CNL
    WRITE(21,232) "        nFaces          ",N,";"
    N=F_INTERNAL_NUM
    WRITE(21,232) "        startFace       ",N,";"
    WRITE(21,*) "    }"
    WRITE(21,*) "    INLET2"
    WRITE(21,*) "    {"    
    WRITE(21,*) "        type            patch;"
    N=0
    WRITE(21,232) "        nFaces          ",N,";"
    N=F_INTERNAL_NUM+CNL*CNL
    WRITE(21,232) "        startFace       ",N,";"
    WRITE(21,*) "    }"
    WRITE(21,*) "    OUTLET"
    WRITE(21,*) "    {"
    WRITE(21,*) "        type            patch;"
    N=0
    WRITE(21,232) "        nFaces          ",N,";"
    N=F_INTERNAL_NUM+CNL*CNL
    WRITE(21,232) "        startFace       ",N,";"
    WRITE(21,*) "    }"
    WRITE(21,*) "    HULL"
    WRITE(21,*) "    {"
    WRITE(21,*) "        type            wall;"
    WRITE(21,*) "        inGroups        1(wall);"
    N=0
    WRITE(21,232) "        nFaces          ",N,";"
    N=F_INTERNAL_NUM+CNL*CNL
    WRITE(21,232) "        startFace       ",N,";"
    WRITE(21,*) "    }"
    
        !FRONT/DOWN FACE--PROCESSOR_MI1_TO_MI2  !CASE RECTANGLE 1-1-1 
    PRO1=PNI*PNJ+PNL*PNI*(MJ2-1)+PNI*(MJ1-1)+MI-1
    MJ=1
    PRO2=PNI*(MJ-1)+MI-1
    WRITE(21,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(21,*) '    {'
    WRITE(21,*) '        type            processor;'
    WRITE(21,*) '        inGroups        List<word> 1(processor);'
    N=CNI*CNL*2
    WRITE(21,232) "        nFaces          ",N,";"
    N=F_INTERNAL_NUM+CNL*CNL
    WRITE(21,232) "        startFace       ",N,";"
    WRITE(21,*) '        matchTolerance  0.0001;'
    WRITE(21,232) "        myProcNo        ",PRO1,";"
    WRITE(21,232) "        neighbProcNo    ",PRO2,";"
    WRITE(21,*) '    }'
    WRITE(41,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(41,*) '    {'
    WRITE(41,*) '        type            processor;'
    WRITE(41,*) '        value           nonuniform List<scalar> '
    N=CNI*CNL*2
    WRITE(41,'(I0)') N
    WRITE(41,*) '('
    DO JP=1,CNL
    DO IP=1,CNI
    WRITE(41,'(E21.14)') P_FRONT_RECT(IP,JP)    
    ENDDO
    ENDDO
    DO KP=1,CNL
    DO IP=1,CNI
    WRITE(41,'(E21.14)') P_DOWN_RECT(IP,KP)    
    ENDDO
    ENDDO
    WRITE(41,*) ')'
    WRITE(41,*) ';'
    WRITE(41,*) '    }'
    WRITE(51,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(51,*) '    {'
    WRITE(51,*) '        type            processor;'
    WRITE(51,*) '        value           nonuniform List<vector> '
    N=CNI*CNL*2
    WRITE(51,'(I0)') N
    WRITE(51,*) '('
    DO JP=1,CNL
    DO IP=1,CNI
    UU1=U_FRONT_RECT(1,IP,JP)
    UU2=U_FRONT_RECT(2,IP,JP)
    UU3=U_FRONT_RECT(3,IP,JP)
    WRITE(51,230) '( ',UU1,' ',UU2,' ',UU3,' )'  
    ENDDO
    ENDDO
    DO KP=1,CNL
    DO IP=1,CNI
    UU1=U_DOWN_RECT(1,IP,KP)
    UU2=U_DOWN_RECT(2,IP,KP)
    UU3=U_DOWN_RECT(3,IP,KP)
    WRITE(51,230) '( ',UU1,' ',UU2,' ',UU3,' )'  
    ENDDO
    ENDDO
    WRITE(51,*) ')'
    WRITE(51,*) ';'
    WRITE(51,*) '    }'
    WRITE(61,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(61,*) '    {'
    WRITE(61,*) '        type            processor;'
    WRITE(61,*) "        value           uniform 0;"
    WRITE(61,*) '    }'
   
    
        !LEFT FACE--PROCESSOR_MI1_TO_MI2  !CASE RECTANGLE 1-1-1 
    PRO2=PNI*PNJ+PNL*PNI*(MJ2-1)+PNI*(MJ1-1)+MI
    WRITE(21,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(21,*) '    {'
    WRITE(21,*) '        type            processor;'
    WRITE(21,*) '        inGroups        List<word> 1(processor);'
    N=CNL*CNL
    WRITE(21,232) "        nFaces          ",N,";"
    N=F_INTERNAL_NUM+CNL*CNL+CNI*CNL*2
    WRITE(21,232) "        startFace       ",N,";"
    WRITE(21,*) '        matchTolerance  0.0001;'
    WRITE(21,232) "        myProcNo        ",PRO1,";"
    WRITE(21,232) "        neighbProcNo    ",PRO2,";"
    WRITE(21,*) '    }'
    WRITE(41,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(41,*) '    {'
    WRITE(41,*) '        type            processor;'
    WRITE(41,*) '        value           nonuniform List<scalar> '
    N=CNL*CNL
    WRITE(41,'(I0)') N
    WRITE(41,*) '('
    DO KP=1,CNL
    DO JP=1,CNL
    WRITE(41,'(E21.14)') P_LEFT_RECT(JP,KP)    
    ENDDO
    ENDDO
    WRITE(41,*) ')'
    WRITE(41,*) ';'
    WRITE(41,*) '    }'
    WRITE(51,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(51,*) '    {'
    WRITE(51,*) '        type            processor;'
    WRITE(51,*) '        value           nonuniform List<vector> '
    N=CNL*CNL
    WRITE(51,'(I0)') N
    WRITE(51,*) '('
    DO KP=1,CNL
    DO JP=1,CNL
    UU1=U_LEFT_RECT(1,JP,KP)
    UU2=U_LEFT_RECT(2,JP,KP)
    UU3=U_LEFT_RECT(3,JP,KP)
    WRITE(51,230) '( ',UU1,' ',UU2,' ',UU3,' )'  
    ENDDO
    ENDDO
    WRITE(51,*) ')'
    WRITE(51,*) ';'
    WRITE(51,*) '    }'
    WRITE(61,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(61,*) '    {'
    WRITE(61,*) '        type            processor;'
    WRITE(61,*) "        value           uniform 0;"
    WRITE(61,*) '    }' 
    
        !UP FACE--PROCESSOR_MI1_TO_MI2  !CASE RECTANGLE 1-1-1
    PRO2=PNI*PNJ+PNL*PNI*(MJ2-1)+PNI*MJ1+MI-1
    WRITE(21,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(21,*) '    {'
    WRITE(21,*) '        type            processor;'
    WRITE(21,*) '        inGroups        List<word> 1(processor);'
    N=CNI*CNL
    WRITE(21,232) "        nFaces          ",N,";"
    N=F_INTERNAL_NUM+CNL*CNL+CNI*CNL*2+CNL*CNL
    WRITE(21,232) "        startFace       ",N,";"
    WRITE(21,*) '        matchTolerance  0.0001;'
    WRITE(21,232) "        myProcNo        ",PRO1,";"
    WRITE(21,232) "        neighbProcNo    ",PRO2,";"
    WRITE(21,*) '    }'
    WRITE(41,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(41,*) '    {'
    WRITE(41,*) '        type            processor;'
    WRITE(41,*) '        value           nonuniform List<scalar> '
    N=CNI*CNL
    WRITE(41,'(I0)') N
    WRITE(41,*) '('
    DO KP=1,CNL
    DO IP=1,CNI
    WRITE(41,'(E21.14)') P_UP_RECT(IP,KP)    
    ENDDO
    ENDDO
    WRITE(41,*) ')'
    WRITE(41,*) ';'
    WRITE(41,*) '    }'
    WRITE(51,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(51,*) '    {'
    WRITE(51,*) '        type            processor;'
    WRITE(51,*) '        value           nonuniform List<vector> '
    N=CNI*CNL
    WRITE(51,'(I0)') N
    WRITE(51,*) '('
    DO KP=1,CNL
    DO IP=1,CNI
    UU1=U_UP_RECT(1,IP,KP)
    UU2=U_UP_RECT(2,IP,KP)
    UU3=U_UP_RECT(3,IP,KP)
    WRITE(51,230) '( ',UU1,' ',UU2,' ',UU3,' )'  
    ENDDO
    ENDDO
    WRITE(51,*) ')'
    WRITE(51,*) ';'
    WRITE(51,*) '    }'
    WRITE(61,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(61,*) '    {'
    WRITE(61,*) '        type            processor;'
    WRITE(61,*) "        value           uniform 0;"
    WRITE(61,*) '    }'
    
        !BACK FACE--PROCESSOR_MI1_TO_MI2  !CASE RECTANGLE 1-1-1
    PRO2=PNI*PNJ+PNL*PNI*MJ2+PNI*(MJ1-1)+MI-1
    WRITE(21,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(21,*) '    {'
    WRITE(21,*) '        type            processor;'
    WRITE(21,*) '        inGroups        List<word> 1(processor);'
    N=CNI*CNL
    WRITE(21,232) "        nFaces          ",N,";"
    N=F_INTERNAL_NUM+CNL*CNL+CNI*CNL*2+CNL*CNL+CNI*CNL
    WRITE(21,232) "        startFace       ",N,";"
    WRITE(21,*) '        matchTolerance  0.0001;'
    WRITE(21,232) "        myProcNo        ",PRO1,";"
    WRITE(21,232) "        neighbProcNo    ",PRO2,";"
    WRITE(21,*) '    }'
    WRITE(41,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(41,*) '    {'
    WRITE(41,*) '        type            processor;'
    WRITE(41,*) '        value           nonuniform List<scalar> '
    N=CNI*CNL
    WRITE(41,'(I0)') N
    WRITE(41,*) '('
    DO JP=1,CNL
    DO IP=1,CNI
    WRITE(41,'(E21.14)') P_BACK_RECT(IP,JP)    
    ENDDO
    ENDDO
    WRITE(41,*) ')'
    WRITE(41,*) ';'
    WRITE(41,*) '    }'
    WRITE(51,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(51,*) '    {'
    WRITE(51,*) '        type            processor;'
    WRITE(51,*) '        value           nonuniform List<vector> '
    N=CNI*CNL
    WRITE(51,'(I0)') N
    WRITE(51,*) '('
    DO JP=1,CNL
    DO IP=1,CNI
    UU1=U_BACK_RECT(1,IP,JP)
    UU2=U_BACK_RECT(2,IP,JP)
    UU3=U_BACK_RECT(3,IP,JP)
    WRITE(51,230) '( ',UU1,' ',UU2,' ',UU3,' )'  
    ENDDO
    ENDDO
    WRITE(51,*) ')'
    WRITE(51,*) ';'
    WRITE(51,*) '    }'
    WRITE(61,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(61,*) '    {'
    WRITE(61,*) '        type            processor;'
    WRITE(61,*) "        value           uniform 0;"
    WRITE(61,*) '    }'
    
    
    
        
    ELSEIF (MJ2.EQ.PNL) THEN
    !CASE RECTANGLE 1-1-3
    
        !RIGHT FACE--INLET1  !CASE RECTANGLE 1-1-3 
    IP=1
    I=IP+CNI*(MI-1)
    DO KP=1,CNL
    K=KP+CNL*(MJ2-1)
    DO JP=1,CNL
    J=JP+CNL*(MJ1-1)
    F11=(CNL+1)*(CNI+1)*(KP-1)+(CNI+1)*(JP-1)+IP-1
    F12=F11+1
    F13=F12+CNI+1
    F14=F13-1
    F21=F11+(CNL+1)*(CNI+1)
    F22=F21+1
    F23=F22+CNI+1
    F24=F23-1
    CALL CAL_FIJ_RECT(F11,F12,F13,F14,F21,F22,F23,F24)
    WRITE(2,235) '4(',F11,' ',F21,' ',F24,' ',F14,')'
    FAD1=F_INTERNAL_NUM_ALL+NL*(K-1)+J
    WRITE(3,231) FAD1
    OWN=CELL_ID_RECT(IP,JP,KP)
    WRITE(4,231) OWN
    ENDDO   
    ENDDO
    
        !BACK FACE--PROCESSOR_MI1_TO_MI2  !CASE RECTANGLE 1-1-3
    KP=CNL
    K=KP+CNL*(MJ2-1)
    DO JP=CNL,1,-1
    J=JP+CNL*(MJ1-1)
    DO IP=1,CNI
    I=IP+CNI*(MI-1)
    F11=(CNL+1)*(CNI+1)*(KP-1)+(CNI+1)*(JP-1)+IP-1
    F12=F11+1
    F13=F12+CNI+1
    F14=F13-1
    F21=F11+(CNL+1)*(CNI+1)
    F22=F21+1
    F23=F22+CNI+1
    F24=F23-1
    CALL CAL_FIJ_RECT(F11,F12,F13,F14,F21,F22,F23,F24)
    !WRITE(2,235) '4(',F21,' ',F22,' ',F23,' ',F24,')'
    WRITE(2,235) '4(',F24,' ',F21,' ',F22,' ',F23,')'
    FAD1=F_BACK_RECT(IP,JP)
    WRITE(3,231) FAD1
    OWN=CELL_ID_RECT(IP,JP,KP)
    WRITE(4,231) OWN
    ENDDO   
    ENDDO
    
        !DOWN FACE--PROCESSOR_MI1_TO_MI2  !CASE RECTANGLE 1-1-3
    JP=1
    J=JP+CNL*(MJ1-1)
    DO KP=CNL,1,-1
    K=KP+CNL*(MJ2-1)
    DO IP=1,CNI
    I=IP+CNI*(MI-1)
    F11=(CNL+1)*(CNI+1)*(KP-1)+(CNI+1)*(JP-1)+IP-1
    F12=F11+1
    F13=F12+CNI+1
    F14=F13-1
    F21=F11+(CNL+1)*(CNI+1)
    F22=F21+1
    F23=F22+CNI+1
    F24=F23-1
    CALL CAL_FIJ_RECT(F11,F12,F13,F14,F21,F22,F23,F24)
    !WRITE(2,235) '4(',F11,' ',F12,' ',F22,' ',F21,')'
    WRITE(2,235) '4(',F21,' ',F11,' ',F12,' ',F22,')'
    FAD1=F_DOWN_RECT(IP,KP)
    WRITE(3,231) FAD1
    OWN=CELL_ID_RECT(IP,JP,KP)
    WRITE(4,231) OWN
    ENDDO   
    ENDDO
    
        !FRONT FACE--PROCESSOR_MI1_TO_MI2  !CASE RECTANGLE 1-1-3
    KP=1
    K=KP+CNL*(MJ2-1)
    DO JP=1,CNL
    J=JP+CNL*(MJ1-1)
    DO IP=1,CNI
    I=IP+CNI*(MI-1)
    F11=(CNL+1)*(CNI+1)*(KP-1)+(CNI+1)*(JP-1)+IP-1
    F12=F11+1
    F13=F12+CNI+1
    F14=F13-1
    F21=F11+(CNL+1)*(CNI+1)
    F22=F21+1
    F23=F22+CNI+1
    F24=F23-1
    CALL CAL_FIJ_RECT(F11,F12,F13,F14,F21,F22,F23,F24)
    WRITE(2,235) '4(',F11,' ',F14,' ',F13,' ',F12,')'
    FAD1=F_FRONT_RECT(IP,JP)
    WRITE(3,231) FAD1
    OWN=CELL_ID_RECT(IP,JP,KP)
    WRITE(4,231) OWN
    ENDDO   
    ENDDO
    
        !LEFT FACE--PROCESSOR_MI1_TO_MI2  !CASE RECTANGLE 1-1-3 
    IP=CNI
    I=IP+CNI*(MI-1)
    DO KP=1,CNL
    K=KP+CNL*(MJ2-1)
    DO JP=1,CNL
    J=JP+CNL*(MJ1-1)
    F11=(CNL+1)*(CNI+1)*(KP-1)+(CNI+1)*(JP-1)+IP-1
    F12=F11+1
    F13=F12+CNI+1
    F14=F13-1
    F21=F11+(CNL+1)*(CNI+1)
    F22=F21+1
    F23=F22+CNI+1
    F24=F23-1
    CALL CAL_FIJ_RECT(F11,F12,F13,F14,F21,F22,F23,F24)
    WRITE(2,235) '4(',F12,' ',F13,' ',F23,' ',F22,')'
    FAD1=F_LEFT_RECT(JP,KP)
    WRITE(3,231) FAD1
    OWN=CELL_ID_RECT(IP,JP,KP)
    WRITE(4,231) OWN
    ENDDO   
    ENDDO
    
        !UP FACE--PROCESSOR_MI1_TO_MI2  !CASE RECTANGLE 1-1-3
    JP=CNL
    J=JP+CNL*(MJ1-1)
    DO KP=1,CNL
    K=KP+CNL*(MJ2-1)
    DO IP=1,CNI
    I=IP+CNI*(MI-1)
    F11=(CNL+1)*(CNI+1)*(KP-1)+(CNI+1)*(JP-1)+IP-1
    F12=F11+1
    F13=F12+CNI+1
    F14=F13-1
    F21=F11+(CNL+1)*(CNI+1)
    F22=F21+1
    F23=F22+CNI+1
    F24=F23-1
    CALL CAL_FIJ_RECT(F11,F12,F13,F14,F21,F22,F23,F24)
    WRITE(2,235) '4(',F14,' ',F24,' ',F23,' ',F13,')'
    FAD1=F_UP_RECT(IP,KP)
    WRITE(3,231) FAD1
    OWN=CELL_ID_RECT(IP,JP,KP)
    WRITE(4,231) OWN
    ENDDO   
    ENDDO
    
    N=4
    WRITE(31,231) N+4
    WRITE(31,*) '('  
    WRITE(31,*) '0'
    WRITE(31,*) '1'
    WRITE(31,*) '2' 
    WRITE(31,*) '3' 
    DO I=1,N
    WRITE(31,*) '-1' 
    ENDDO
    WRITE(31,*) ')'
    
    WRITE(21,231) N+4
    WRITE(21,*) "("    
    WRITE(21,*) "    INLET1"
    WRITE(21,*) "    {"    
    WRITE(21,*) "        type            patch;"
    N=CNL*CNL
    WRITE(21,232) "        nFaces          ",N,";"
    N=F_INTERNAL_NUM
    WRITE(21,232) "        startFace       ",N,";"
    WRITE(21,*) "    }"
    WRITE(21,*) "    INLET2"
    WRITE(21,*) "    {"    
    WRITE(21,*) "        type            patch;"
    N=0
    WRITE(21,232) "        nFaces          ",N,";"
    N=F_INTERNAL_NUM+CNL*CNL
    WRITE(21,232) "        startFace       ",N,";"
    WRITE(21,*) "    }"
    WRITE(21,*) "    OUTLET"
    WRITE(21,*) "    {"
    WRITE(21,*) "        type            patch;"
    N=0
    WRITE(21,232) "        nFaces          ",N,";"
    N=F_INTERNAL_NUM+CNL*CNL
    WRITE(21,232) "        startFace       ",N,";"
    WRITE(21,*) "    }"
    WRITE(21,*) "    HULL"
    WRITE(21,*) "    {"
    WRITE(21,*) "        type            wall;"
    WRITE(21,*) "        inGroups        1(wall);"
    N=0
    WRITE(21,232) "        nFaces          ",N,";"
    N=F_INTERNAL_NUM+CNL*CNL
    WRITE(21,232) "        startFace       ",N,";"
    WRITE(21,*) "    }"
    
        !BACK/DOWN FACE--PROCESSOR_MI1_TO_MI2  !CASE RECTANGLE 1-1-3
    PRO1=PNI*PNJ+PNL*PNI*(MJ2-1)+PNI*(MJ1-1)+MI-1
    MJ=1
    PRO2=PNI*(MJ-1)+MI-1
    WRITE(21,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(21,*) '    {'
    WRITE(21,*) '        type            processor;'
    WRITE(21,*) '        inGroups        List<word> 1(processor);'
    N=CNI*CNL*2
    WRITE(21,232) "        nFaces          ",N,";"
    N=F_INTERNAL_NUM+CNL*CNL
    WRITE(21,232) "        startFace       ",N,";"
    WRITE(21,*) '        matchTolerance  0.0001;'
    WRITE(21,232) "        myProcNo        ",PRO1,";"
    WRITE(21,232) "        neighbProcNo    ",PRO2,";"
    WRITE(21,*) '    }'
    WRITE(41,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(41,*) '    {'
    WRITE(41,*) '        type            processor;'
    WRITE(41,*) '        value           nonuniform List<scalar> '
    N=CNI*CNL*2
    WRITE(41,'(I0)') N
    WRITE(41,*) '('
    DO JP=1,CNL
    DO IP=1,CNI
    WRITE(41,'(E21.14)') P_BACK_RECT(IP,JP)    
    ENDDO
    ENDDO
    DO KP=1,CNL
    DO IP=1,CNI
    WRITE(41,'(E21.14)') P_DOWN_RECT(IP,KP)    
    ENDDO
    ENDDO
    WRITE(41,*) ')'
    WRITE(41,*) ';'
    WRITE(41,*) '    }'
    WRITE(51,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(51,*) '    {'
    WRITE(51,*) '        type            processor;'
    WRITE(51,*) '        value           nonuniform List<vector> '
    N=CNI*CNL*2
    WRITE(51,'(I0)') N
    WRITE(51,*) '('
    DO JP=1,CNL
    DO IP=1,CNI
    UU1=U_BACK_RECT(1,IP,JP)
    UU2=U_BACK_RECT(2,IP,JP)
    UU3=U_BACK_RECT(3,IP,JP)
    WRITE(51,230) '( ',UU1,' ',UU2,' ',UU3,' )'  
    ENDDO
    ENDDO
    DO KP=1,CNL
    DO IP=1,CNI
    UU1=U_DOWN_RECT(1,IP,KP)
    UU2=U_DOWN_RECT(2,IP,KP)
    UU3=U_DOWN_RECT(3,IP,KP)
    WRITE(51,230) '( ',UU1,' ',UU2,' ',UU3,' )'  
    ENDDO
    ENDDO
    WRITE(51,*) ')'
    WRITE(51,*) ';'
    WRITE(51,*) '    }'
    WRITE(61,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(61,*) '    {'
    WRITE(61,*) '        type            processor;'
    WRITE(61,*) "        value           uniform 0;"
    WRITE(61,*) '    }'
    
        !FRONT FACE--PROCESSOR_MI1_TO_MI2  !CASE RECTANGLE 1-1-3
    PRO2=PNI*PNJ+PNL*PNI*(MJ2-2)+PNI*(MJ1-1)+MI-1
    WRITE(21,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(21,*) '    {'
    WRITE(21,*) '        type            processor;'
    WRITE(21,*) '        inGroups        List<word> 1(processor);'
    N=CNI*CNL
    WRITE(21,232) "        nFaces          ",N,";"
    N=F_INTERNAL_NUM+CNL*CNL+CNI*CNL*2
    WRITE(21,232) "        startFace       ",N,";"
    WRITE(21,*) '        matchTolerance  0.0001;'
    WRITE(21,232) "        myProcNo        ",PRO1,";"
    WRITE(21,232) "        neighbProcNo    ",PRO2,";"
    WRITE(21,*) '    }'
    WRITE(41,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(41,*) '    {'
    WRITE(41,*) '        type            processor;'
    WRITE(41,*) '        value           nonuniform List<scalar> '
    N=CNI*CNL
    WRITE(41,'(I0)') N
    WRITE(41,*) '('
    DO JP=1,CNL
    DO IP=1,CNI
    WRITE(41,'(E21.14)') P_FRONT_RECT(IP,JP)    
    ENDDO
    ENDDO
    WRITE(41,*) ')'
    WRITE(41,*) ';'
    WRITE(41,*) '    }'
    WRITE(51,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(51,*) '    {'
    WRITE(51,*) '        type            processor;'
    WRITE(51,*) '        value           nonuniform List<vector> '
    N=CNI*CNL
    WRITE(51,'(I0)') N
    WRITE(51,*) '('
    DO JP=1,CNL
    DO IP=1,CNI
    UU1=U_FRONT_RECT(1,IP,JP)
    UU2=U_FRONT_RECT(2,IP,JP)
    UU3=U_FRONT_RECT(3,IP,JP)
    WRITE(51,230) '( ',UU1,' ',UU2,' ',UU3,' )'  
    ENDDO
    ENDDO
    WRITE(51,*) ')'
    WRITE(51,*) ';'
    WRITE(51,*) '    }'
    WRITE(61,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(61,*) '    {'
    WRITE(61,*) '        type            processor;'
    WRITE(61,*) "        value           uniform 0;"
    WRITE(61,*) '    }'
    
        !LEFT FACE--PROCESSOR_MI1_TO_MI2  !CASE RECTANGLE 1-1-3 
    PRO2=PNI*PNJ+PNL*PNI*(MJ2-1)+PNI*(MJ1-1)+MI
    WRITE(21,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(21,*) '    {'
    WRITE(21,*) '        type            processor;'
    WRITE(21,*) '        inGroups        List<word> 1(processor);'
    N=CNL*CNL
    WRITE(21,232) "        nFaces          ",N,";"
    N=F_INTERNAL_NUM+CNL*CNL+CNI*CNL*2+CNI*CNL
    WRITE(21,232) "        startFace       ",N,";"
    WRITE(21,*) '        matchTolerance  0.0001;'
    WRITE(21,232) "        myProcNo        ",PRO1,";"
    WRITE(21,232) "        neighbProcNo    ",PRO2,";"
    WRITE(21,*) '    }'
    WRITE(41,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(41,*) '    {'
    WRITE(41,*) '        type            processor;'
    WRITE(41,*) '        value           nonuniform List<scalar> '
    N=CNL*CNL
    WRITE(41,'(I0)') N
    WRITE(41,*) '('
    DO KP=1,CNL
    DO JP=1,CNL
    WRITE(41,'(E21.14)') P_LEFT_RECT(JP,KP)    
    ENDDO
    ENDDO
    WRITE(41,*) ')'
    WRITE(41,*) ';'
    WRITE(41,*) '    }'
    WRITE(51,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(51,*) '    {'
    WRITE(51,*) '        type            processor;'
    WRITE(51,*) '        value           nonuniform List<vector> '
    N=CNL*CNL
    WRITE(51,'(I0)') N
    WRITE(51,*) '('
    DO KP=1,CNL
    DO JP=1,CNL
    UU1=U_LEFT_RECT(1,JP,KP)
    UU2=U_LEFT_RECT(2,JP,KP)
    UU3=U_LEFT_RECT(3,JP,KP)
    WRITE(51,230) '( ',UU1,' ',UU2,' ',UU3,' )'  
    ENDDO
    ENDDO
    WRITE(51,*) ')'
    WRITE(51,*) ';'
    WRITE(51,*) '    }'
    WRITE(61,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(61,*) '    {'
    WRITE(61,*) '        type            processor;'
    WRITE(61,*) "        value           uniform 0;"
    WRITE(61,*) '    }'
    
        !UP FACE--PROCESSOR_MI1_TO_MI2  !CASE RECTANGLE 1-1-3
    PRO2=PNI*PNJ+PNL*PNI*(MJ2-1)+PNI*MJ1+MI-1
    WRITE(21,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(21,*) '    {'
    WRITE(21,*) '        type            processor;'
    WRITE(21,*) '        inGroups        List<word> 1(processor);'
    N=CNI*CNL
    WRITE(21,232) "        nFaces          ",N,";"
    N=F_INTERNAL_NUM+CNL*CNL+CNI*CNL*2+CNI*CNL+CNL*CNL
    WRITE(21,232) "        startFace       ",N,";"
    WRITE(21,*) '        matchTolerance  0.0001;'
    WRITE(21,232) "        myProcNo        ",PRO1,";"
    WRITE(21,232) "        neighbProcNo    ",PRO2,";"
    WRITE(21,*) '    }'
    WRITE(41,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(41,*) '    {'
    WRITE(41,*) '        type            processor;'
    WRITE(41,*) '        value           nonuniform List<scalar> '
    N=CNI*CNL
    WRITE(41,'(I0)') N
    WRITE(41,*) '('
    DO KP=1,CNL
    DO IP=1,CNI
    WRITE(41,'(E21.14)') P_UP_RECT(IP,KP)    
    ENDDO
    ENDDO
    WRITE(41,*) ')'
    WRITE(41,*) ';'
    WRITE(41,*) '    }'
    WRITE(51,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(51,*) '    {'
    WRITE(51,*) '        type            processor;'
    WRITE(51,*) '        value           nonuniform List<vector> '
    N=CNI*CNL
    WRITE(51,'(I0)') N
    WRITE(51,*) '('
    DO KP=1,CNL
    DO IP=1,CNI
    UU1=U_UP_RECT(1,IP,KP)
    UU2=U_UP_RECT(2,IP,KP)
    UU3=U_UP_RECT(3,IP,KP)
    WRITE(51,230) '( ',UU1,' ',UU2,' ',UU3,' )'  
    ENDDO
    ENDDO
    WRITE(51,*) ')'
    WRITE(51,*) ';'
    WRITE(51,*) '    }'
    WRITE(61,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(61,*) '    {'
    WRITE(61,*) '        type            processor;'
    WRITE(61,*) "        value           uniform 0;"
    WRITE(61,*) '    }'
    
    
    
    
        
    ELSE
    !CASE RECTANGLE 1-1-2
    
        !RIGHT FACE--INLET1  !CASE RECTANGLE 1-1-2 
    IP=1
    I=IP+CNI*(MI-1)
    DO KP=1,CNL
    K=KP+CNL*(MJ2-1)
    DO JP=1,CNL
    J=JP+CNL*(MJ1-1)
    F11=(CNL+1)*(CNI+1)*(KP-1)+(CNI+1)*(JP-1)+IP-1
    F12=F11+1
    F13=F12+CNI+1
    F14=F13-1
    F21=F11+(CNL+1)*(CNI+1)
    F22=F21+1
    F23=F22+CNI+1
    F24=F23-1
    CALL CAL_FIJ_RECT(F11,F12,F13,F14,F21,F22,F23,F24)
    WRITE(2,235) '4(',F11,' ',F21,' ',F24,' ',F14,')'
    FAD1=F_INTERNAL_NUM_ALL+NL*(K-1)+J
    WRITE(3,231) FAD1
    OWN=CELL_ID_RECT(IP,JP,KP)
    WRITE(4,231) OWN
    ENDDO   
    ENDDO    
        
        !DOWN FACE--PROCESSOR_MI1_TO_MI2  !CASE RECTANGLE 1-1-2
    JP=1
    J=JP+CNL*(MJ1-1)
    DO KP=CNL,1,-1
    K=KP+CNL*(MJ2-1)
    DO IP=1,CNI
    I=IP+CNI*(MI-1)
    F11=(CNL+1)*(CNI+1)*(KP-1)+(CNI+1)*(JP-1)+IP-1
    F12=F11+1
    F13=F12+CNI+1
    F14=F13-1
    F21=F11+(CNL+1)*(CNI+1)
    F22=F21+1
    F23=F22+CNI+1
    F24=F23-1
    CALL CAL_FIJ_RECT(F11,F12,F13,F14,F21,F22,F23,F24)
    !WRITE(2,235) '4(',F11,' ',F12,' ',F22,' ',F21,')'
    WRITE(2,235) '4(',F21,' ',F11,' ',F12,' ',F22,')'
    FAD1=F_DOWN_RECT(IP,KP)
    WRITE(3,231) FAD1
    OWN=CELL_ID_RECT(IP,JP,KP)
    WRITE(4,231) OWN
    ENDDO   
    ENDDO
    
        !FRONT FACE--PROCESSOR_MI1_TO_MI2  !CASE RECTANGLE 1-1-2
    KP=1
    K=KP+CNL*(MJ2-1)
    DO JP=1,CNL
    J=JP+CNL*(MJ1-1)
    DO IP=1,CNI
    I=IP+CNI*(MI-1)
    F11=(CNL+1)*(CNI+1)*(KP-1)+(CNI+1)*(JP-1)+IP-1
    F12=F11+1
    F13=F12+CNI+1
    F14=F13-1
    F21=F11+(CNL+1)*(CNI+1)
    F22=F21+1
    F23=F22+CNI+1
    F24=F23-1
    CALL CAL_FIJ_RECT(F11,F12,F13,F14,F21,F22,F23,F24)
    WRITE(2,235) '4(',F11,' ',F14,' ',F13,' ',F12,')'
    FAD1=F_FRONT_RECT(IP,JP)
    WRITE(3,231) FAD1
    OWN=CELL_ID_RECT(IP,JP,KP)
    WRITE(4,231) OWN
    ENDDO   
    ENDDO
    
        !LEFT FACE--PROCESSOR_MI1_TO_MI2  !CASE RECTANGLE 1-1-2 
    IP=CNI
    I=IP+CNI*(MI-1)
    DO KP=1,CNL
    K=KP+CNL*(MJ2-1)
    DO JP=1,CNL
    J=JP+CNL*(MJ1-1)
    F11=(CNL+1)*(CNI+1)*(KP-1)+(CNI+1)*(JP-1)+IP-1
    F12=F11+1
    F13=F12+CNI+1
    F14=F13-1
    F21=F11+(CNL+1)*(CNI+1)
    F22=F21+1
    F23=F22+CNI+1
    F24=F23-1
    CALL CAL_FIJ_RECT(F11,F12,F13,F14,F21,F22,F23,F24)
    WRITE(2,235) '4(',F12,' ',F13,' ',F23,' ',F22,')'
    FAD1=F_LEFT_RECT(JP,KP)
    WRITE(3,231) FAD1
    OWN=CELL_ID_RECT(IP,JP,KP)
    WRITE(4,231) OWN
    ENDDO   
    ENDDO
    
        !UP FACE--PROCESSOR_MI1_TO_MI2  !CASE RECTANGLE 1-1-2
    JP=CNL
    J=JP+CNL*(MJ1-1)
    DO KP=1,CNL
    K=KP+CNL*(MJ2-1)
    DO IP=1,CNI
    I=IP+CNI*(MI-1)
    F11=(CNL+1)*(CNI+1)*(KP-1)+(CNI+1)*(JP-1)+IP-1
    F12=F11+1
    F13=F12+CNI+1
    F14=F13-1
    F21=F11+(CNL+1)*(CNI+1)
    F22=F21+1
    F23=F22+CNI+1
    F24=F23-1
    CALL CAL_FIJ_RECT(F11,F12,F13,F14,F21,F22,F23,F24)
    WRITE(2,235) '4(',F14,' ',F24,' ',F23,' ',F13,')'
    FAD1=F_UP_RECT(IP,KP)
    WRITE(3,231) FAD1
    OWN=CELL_ID_RECT(IP,JP,KP)
    WRITE(4,231) OWN
    ENDDO   
    ENDDO
    
        !BACK FACE--PROCESSOR_MI1_TO_MI2  !CASE RECTANGLE 1-1-2
    KP=CNL
    K=KP+CNL*(MJ2-1)
    DO JP=1,CNL
    J=JP+CNL*(MJ1-1)
    DO IP=1,CNI
    I=IP+CNI*(MI-1)
    F11=(CNL+1)*(CNI+1)*(KP-1)+(CNI+1)*(JP-1)+IP-1
    F12=F11+1
    F13=F12+CNI+1
    F14=F13-1
    F21=F11+(CNL+1)*(CNI+1)
    F22=F21+1
    F23=F22+CNI+1
    F24=F23-1
    CALL CAL_FIJ_RECT(F11,F12,F13,F14,F21,F22,F23,F24)
    WRITE(2,235) '4(',F21,' ',F22,' ',F23,' ',F24,')'
    FAD1=F_BACK_RECT(IP,JP)
    WRITE(3,231) FAD1
    OWN=CELL_ID_RECT(IP,JP,KP)
    WRITE(4,231) OWN
    ENDDO   
    ENDDO
    
    N=5
    WRITE(31,231) N+4
    WRITE(31,*) '('  
    WRITE(31,*) '0'
    WRITE(31,*) '1'
    WRITE(31,*) '2' 
    WRITE(31,*) '3' 
    DO I=1,N
    WRITE(31,*) '-1' 
    ENDDO
    WRITE(31,*) ')'
    
    WRITE(21,231) N+4
    WRITE(21,*) "("    
    WRITE(21,*) "    INLET1"
    WRITE(21,*) "    {"    
    WRITE(21,*) "        type            patch;"
    N=CNL*CNL
    WRITE(21,232) "        nFaces          ",N,";"
    N=F_INTERNAL_NUM
    WRITE(21,232) "        startFace       ",N,";"
    WRITE(21,*) "    }"
    WRITE(21,*) "    INLET2"
    WRITE(21,*) "    {"    
    WRITE(21,*) "        type            patch;"
    N=0
    WRITE(21,232) "        nFaces          ",N,";"
    N=F_INTERNAL_NUM+CNL*CNL
    WRITE(21,232) "        startFace       ",N,";"
    WRITE(21,*) "    }"
    WRITE(21,*) "    OUTLET"
    WRITE(21,*) "    {"
    WRITE(21,*) "        type            patch;"
    N=0
    WRITE(21,232) "        nFaces          ",N,";"
    N=F_INTERNAL_NUM+CNL*CNL
    WRITE(21,232) "        startFace       ",N,";"
    WRITE(21,*) "    }"
    WRITE(21,*) "    HULL"
    WRITE(21,*) "    {"
    WRITE(21,*) "        type            wall;"
    WRITE(21,*) "        inGroups        1(wall);"
    N=0
    WRITE(21,232) "        nFaces          ",N,";"
    N=F_INTERNAL_NUM+CNL*CNL
    WRITE(21,232) "        startFace       ",N,";"
    WRITE(21,*) "    }"
    
        !DOWN FACE--PROCESSOR_MI1_TO_MI2  !CASE RECTANGLE 1-1-2
    PRO1=PNI*PNJ+PNL*PNI*(MJ2-1)+PNI*(MJ1-1)+MI-1
    MJ=1
    PRO2=PNI*(MJ-1)+MI-1
    WRITE(21,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(21,*) '    {'
    WRITE(21,*) '        type            processor;'
    WRITE(21,*) '        inGroups        List<word> 1(processor);'
    N=CNI*CNL
    WRITE(21,232) "        nFaces          ",N,";"
    N=F_INTERNAL_NUM+CNL*CNL
    WRITE(21,232) "        startFace       ",N,";"
    WRITE(21,*) '        matchTolerance  0.0001;'
    WRITE(21,232) "        myProcNo        ",PRO1,";"
    WRITE(21,232) "        neighbProcNo    ",PRO2,";"
    WRITE(21,*) '    }'
    WRITE(41,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(41,*) '    {'
    WRITE(41,*) '        type            processor;'
    WRITE(41,*) '        value           nonuniform List<scalar> '
    N=CNI*CNL
    WRITE(41,'(I0)') N
    WRITE(41,*) '('
    DO KP=1,CNL
    DO IP=1,CNI
    WRITE(41,'(E21.14)') P_DOWN_RECT(IP,KP)    
    ENDDO
    ENDDO
    WRITE(41,*) ')'
    WRITE(41,*) ';'
    WRITE(41,*) '    }'
    WRITE(51,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(51,*) '    {'
    WRITE(51,*) '        type            processor;'
    WRITE(51,*) '        value           nonuniform List<vector> '
    N=CNI*CNL
    WRITE(51,'(I0)') N
    WRITE(51,*) '('
    DO KP=1,CNL
    DO IP=1,CNI
    UU1=U_DOWN_RECT(1,IP,KP)
    UU2=U_DOWN_RECT(2,IP,KP)
    UU3=U_DOWN_RECT(3,IP,KP)
    WRITE(51,230) '( ',UU1,' ',UU2,' ',UU3,' )'  
    ENDDO
    ENDDO
    WRITE(51,*) ')'
    WRITE(51,*) ';'
    WRITE(51,*) '    }'
    WRITE(61,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(61,*) '    {'
    WRITE(61,*) '        type            processor;'
    WRITE(61,*) "        value           uniform 0;"
    WRITE(61,*) '    }'
    
        !FRONT FACE--PROCESSOR_MI1_TO_MI2  !CASE RECTANGLE 1-1-2
    PRO2=PNI*PNJ+PNL*PNI*(MJ2-2)+PNI*(MJ1-1)+MI-1
    WRITE(21,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(21,*) '    {'
    WRITE(21,*) '        type            processor;'
    WRITE(21,*) '        inGroups        List<word> 1(processor);'
    N=CNI*CNL
    WRITE(21,232) "        nFaces          ",N,";"
    N=F_INTERNAL_NUM+CNL*CNL+CNI*CNL
    WRITE(21,232) "        startFace       ",N,";"
    WRITE(21,*) '        matchTolerance  0.0001;'
    WRITE(21,232) "        myProcNo        ",PRO1,";"
    WRITE(21,232) "        neighbProcNo    ",PRO2,";"
    WRITE(21,*) '    }'
    WRITE(41,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(41,*) '    {'
    WRITE(41,*) '        type            processor;'
    WRITE(41,*) '        value           nonuniform List<scalar> '
    N=CNI*CNL
    WRITE(41,'(I0)') N
    WRITE(41,*) '('
    DO JP=1,CNL
    DO IP=1,CNI
    WRITE(41,'(E21.14)') P_FRONT_RECT(IP,JP)    
    ENDDO
    ENDDO
    WRITE(41,*) ')'
    WRITE(41,*) ';'
    WRITE(41,*) '    }'
    WRITE(51,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(51,*) '    {'
    WRITE(51,*) '        type            processor;'
    WRITE(51,*) '        value           nonuniform List<vector> '
    N=CNI*CNL
    WRITE(51,'(I0)') N
    WRITE(51,*) '('
    DO JP=1,CNL
    DO IP=1,CNI
    UU1=U_FRONT_RECT(1,IP,JP)
    UU2=U_FRONT_RECT(2,IP,JP)
    UU3=U_FRONT_RECT(3,IP,JP)
    WRITE(51,230) '( ',UU1,' ',UU2,' ',UU3,' )'  
    ENDDO
    ENDDO
    WRITE(51,*) ')'
    WRITE(51,*) ';'
    WRITE(51,*) '    }'
    WRITE(61,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(61,*) '    {'
    WRITE(61,*) '        type            processor;'
    WRITE(61,*) "        value           uniform 0;"
    WRITE(61,*) '    }'
    
        !LEFT FACE--PROCESSOR_MI1_TO_MI2  !CASE RECTANGLE 1-1-2
    PRO2=PNI*PNJ+PNL*PNI*(MJ2-1)+PNI*(MJ1-1)+MI
    WRITE(21,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(21,*) '    {'
    WRITE(21,*) '        type            processor;'
    WRITE(21,*) '        inGroups        List<word> 1(processor);'
    N=CNL*CNL
    WRITE(21,232) "        nFaces          ",N,";"
    N=F_INTERNAL_NUM+CNL*CNL+CNI*CNL+CNI*CNL
    WRITE(21,232) "        startFace       ",N,";"
    WRITE(21,*) '        matchTolerance  0.0001;'
    WRITE(21,232) "        myProcNo        ",PRO1,";"
    WRITE(21,232) "        neighbProcNo    ",PRO2,";"
    WRITE(21,*) '    }'
    WRITE(41,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(41,*) '    {'
    WRITE(41,*) '        type            processor;'
    WRITE(41,*) '        value           nonuniform List<scalar> '
    N=CNL*CNL
    WRITE(41,'(I0)') N
    WRITE(41,*) '('
    DO KP=1,CNL
    DO JP=1,CNL
    WRITE(41,'(E21.14)') P_LEFT_RECT(JP,KP)    
    ENDDO
    ENDDO
    WRITE(41,*) ')'
    WRITE(41,*) ';'
    WRITE(41,*) '    }'
    WRITE(51,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(51,*) '    {'
    WRITE(51,*) '        type            processor;'
    WRITE(51,*) '        value           nonuniform List<vector> '
    N=CNL*CNL
    WRITE(51,'(I0)') N
    WRITE(51,*) '('
    DO KP=1,CNL
    DO JP=1,CNL
    UU1=U_LEFT_RECT(1,JP,KP)
    UU2=U_LEFT_RECT(2,JP,KP)
    UU3=U_LEFT_RECT(3,JP,KP)
    WRITE(51,230) '( ',UU1,' ',UU2,' ',UU3,' )'  
    ENDDO
    ENDDO
    WRITE(51,*) ')'
    WRITE(51,*) ';'
    WRITE(51,*) '    }'
    WRITE(61,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(61,*) '    {'
    WRITE(61,*) '        type            processor;'
    WRITE(61,*) "        value           uniform 0;"
    WRITE(61,*) '    }'
    
        !UP FACE--PROCESSOR_MI1_TO_MI2  !CASE RECTANGLE 1-1-2
    PRO2=PNI*PNJ+PNL*PNI*(MJ2-1)+PNI*MJ1+MI-1
    WRITE(21,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(21,*) '    {'
    WRITE(21,*) '        type            processor;'
    WRITE(21,*) '        inGroups        List<word> 1(processor);'
    N=CNI*CNL
    WRITE(21,232) "        nFaces          ",N,";"
    N=F_INTERNAL_NUM+CNL*CNL+CNI*CNL+CNI*CNL+CNL*CNL
    WRITE(21,232) "        startFace       ",N,";"
    WRITE(21,*) '        matchTolerance  0.0001;'
    WRITE(21,232) "        myProcNo        ",PRO1,";"
    WRITE(21,232) "        neighbProcNo    ",PRO2,";"
    WRITE(21,*) '    }'
    WRITE(41,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(41,*) '    {'
    WRITE(41,*) '        type            processor;'
    WRITE(41,*) '        value           nonuniform List<scalar> '
    N=CNI*CNL
    WRITE(41,'(I0)') N
    WRITE(41,*) '('
    DO KP=1,CNL
    DO IP=1,CNI
    WRITE(41,'(E21.14)') P_UP_RECT(IP,KP)    
    ENDDO
    ENDDO
    WRITE(41,*) ')'
    WRITE(41,*) ';'
    WRITE(41,*) '    }'
    WRITE(51,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(51,*) '    {'
    WRITE(51,*) '        type            processor;'
    WRITE(51,*) '        value           nonuniform List<vector> '
    N=CNI*CNL
    WRITE(51,'(I0)') N
    WRITE(51,*) '('
    DO KP=1,CNL
    DO IP=1,CNI
    UU1=U_UP_RECT(1,IP,KP)
    UU2=U_UP_RECT(2,IP,KP)
    UU3=U_UP_RECT(3,IP,KP)
    WRITE(51,230) '( ',UU1,' ',UU2,' ',UU3,' )'  
    ENDDO
    ENDDO
    WRITE(51,*) ')'
    WRITE(51,*) ';'
    WRITE(51,*) '    }'
    WRITE(61,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(61,*) '    {'
    WRITE(61,*) '        type            processor;'
    WRITE(61,*) "        value           uniform 0;"
    WRITE(61,*) '    }'
    
        !BACK FACE--PROCESSOR_MI1_TO_MI2  !CASE RECTANGLE 1-1-2
    PRO2=PNI*PNJ+PNL*PNI*MJ2+PNI*(MJ1-1)+MI-1
    WRITE(21,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(21,*) '    {'
    WRITE(21,*) '        type            processor;'
    WRITE(21,*) '        inGroups        List<word> 1(processor);'
    N=CNI*CNL
    WRITE(21,232) "        nFaces          ",N,";"
    N=F_INTERNAL_NUM+CNL*CNL+CNI*CNL+CNI*CNL+CNL*CNL+CNI*CNL
    WRITE(21,232) "        startFace       ",N,";"
    WRITE(21,*) '        matchTolerance  0.0001;'
    WRITE(21,232) "        myProcNo        ",PRO1,";"
    WRITE(21,232) "        neighbProcNo    ",PRO2,";"
    WRITE(21,*) '    }'
    WRITE(41,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(41,*) '    {'
    WRITE(41,*) '        type            processor;'
    WRITE(41,*) '        value           nonuniform List<scalar> '
    N=CNI*CNL
    WRITE(41,'(I0)') N
    WRITE(41,*) '('
    DO JP=1,CNL
    DO IP=1,CNI
    WRITE(41,'(E21.14)') P_BACK_RECT(IP,JP)    
    ENDDO
    ENDDO
    WRITE(41,*) ')'
    WRITE(41,*) ';'
    WRITE(41,*) '    }'
    WRITE(51,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(51,*) '    {'
    WRITE(51,*) '        type            processor;'
    WRITE(51,*) '        value           nonuniform List<vector> '
    N=CNI*CNL
    WRITE(51,'(I0)') N
    WRITE(51,*) '('
    DO JP=1,CNL
    DO IP=1,CNI
    UU1=U_BACK_RECT(1,IP,JP)
    UU2=U_BACK_RECT(2,IP,JP)
    UU3=U_BACK_RECT(3,IP,JP)
    WRITE(51,230) '( ',UU1,' ',UU2,' ',UU3,' )'  
    ENDDO
    ENDDO
    WRITE(51,*) ')'
    WRITE(51,*) ';'
    WRITE(51,*) '    }'
    WRITE(61,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(61,*) '    {'
    WRITE(61,*) '        type            processor;'
    WRITE(61,*) "        value           uniform 0;"
    WRITE(61,*) '    }'    
        
        
        
        
        
        
        
        
        
        
        
    ENDIF    
    ELSEIF (MJ1.EQ.PNL) THEN
    IF (MJ2.EQ.1) THEN
    !CASE RECTANGLE 1-3-1
    
        !RIGHT FACE--INLET1  !CASE RECTANGLE 1-3-1 
    IP=1
    I=IP+CNI*(MI-1)
    DO KP=1,CNL
    K=KP+CNL*(MJ2-1)
    DO JP=1,CNL
    J=JP+CNL*(MJ1-1)
    F11=(CNL+1)*(CNI+1)*(KP-1)+(CNI+1)*(JP-1)+IP-1
    F12=F11+1
    F13=F12+CNI+1
    F14=F13-1
    F21=F11+(CNL+1)*(CNI+1)
    F22=F21+1
    F23=F22+CNI+1
    F24=F23-1
    CALL CAL_FIJ_RECT(F11,F12,F13,F14,F21,F22,F23,F24)
    WRITE(2,235) '4(',F11,' ',F21,' ',F24,' ',F14,')'
    FAD1=F_INTERNAL_NUM_ALL+NL*(K-1)+J
    WRITE(3,231) FAD1
    OWN=CELL_ID_RECT(IP,JP,KP)
    WRITE(4,231) OWN
    ENDDO   
    ENDDO    
        
        !FRONT FACE--PROCESSOR_MI1_TO_MI2  !CASE RECTANGLE 1-3-1
    KP=1
    K=KP+CNL*(MJ2-1)
    DO JP=1,CNL
    J=JP+CNL*(MJ1-1)
    DO IP=1,CNI
    I=IP+CNI*(MI-1)
    F11=(CNL+1)*(CNI+1)*(KP-1)+(CNI+1)*(JP-1)+IP-1
    F12=F11+1
    F13=F12+CNI+1
    F14=F13-1
    F21=F11+(CNL+1)*(CNI+1)
    F22=F21+1
    F23=F22+CNI+1
    F24=F23-1
    CALL CAL_FIJ_RECT(F11,F12,F13,F14,F21,F22,F23,F24)
    WRITE(2,235) '4(',F11,' ',F14,' ',F13,' ',F12,')'
    FAD1=F_FRONT_RECT(IP,JP)
    WRITE(3,231) FAD1
    OWN=CELL_ID_RECT(IP,JP,KP)
    WRITE(4,231) OWN
    ENDDO   
    ENDDO    
        
        !UP FACE--PROCESSOR_MI1_TO_MI2  !CASE RECTANGLE 1-3-1
    JP=CNL
    J=JP+CNL*(MJ1-1)
    DO KP=1,CNL
    K=KP+CNL*(MJ2-1)
    DO IP=1,CNI
    I=IP+CNI*(MI-1)
    F11=(CNL+1)*(CNI+1)*(KP-1)+(CNI+1)*(JP-1)+IP-1
    F12=F11+1
    F13=F12+CNI+1
    F14=F13-1
    F21=F11+(CNL+1)*(CNI+1)
    F22=F21+1
    F23=F22+CNI+1
    F24=F23-1
    CALL CAL_FIJ_RECT(F11,F12,F13,F14,F21,F22,F23,F24)
    WRITE(2,235) '4(',F14,' ',F24,' ',F23,' ',F13,')'
    FAD1=F_UP_RECT(IP,KP)
    WRITE(3,231) FAD1
    OWN=CELL_ID_RECT(IP,JP,KP)
    WRITE(4,231) OWN
    ENDDO   
    ENDDO    
        
        !DOWN FACE--PROCESSOR_MI1_TO_MI2  !CASE RECTANGLE 1-3-1
    JP=1
    J=JP+CNL*(MJ1-1)
    DO KP=1,CNL
    K=KP+CNL*(MJ2-1)
    DO IP=1,CNI
    I=IP+CNI*(MI-1)
    F11=(CNL+1)*(CNI+1)*(KP-1)+(CNI+1)*(JP-1)+IP-1
    F12=F11+1
    F13=F12+CNI+1
    F14=F13-1
    F21=F11+(CNL+1)*(CNI+1)
    F22=F21+1
    F23=F22+CNI+1
    F24=F23-1
    CALL CAL_FIJ_RECT(F11,F12,F13,F14,F21,F22,F23,F24)
    WRITE(2,235) '4(',F11,' ',F12,' ',F22,' ',F21,')'
    FAD1=F_DOWN_RECT(IP,KP)
    WRITE(3,231) FAD1
    OWN=CELL_ID_RECT(IP,JP,KP)
    WRITE(4,231) OWN
    ENDDO   
    ENDDO    
        
        !LEFT FACE--PROCESSOR_MI1_TO_MI2  !CASE RECTANGLE 1-3-1 
    IP=CNI
    I=IP+CNI*(MI-1)
    DO KP=1,CNL
    K=KP+CNL*(MJ2-1)
    DO JP=1,CNL
    J=JP+CNL*(MJ1-1)
    F11=(CNL+1)*(CNI+1)*(KP-1)+(CNI+1)*(JP-1)+IP-1
    F12=F11+1
    F13=F12+CNI+1
    F14=F13-1
    F21=F11+(CNL+1)*(CNI+1)
    F22=F21+1
    F23=F22+CNI+1
    F24=F23-1
    CALL CAL_FIJ_RECT(F11,F12,F13,F14,F21,F22,F23,F24)
    WRITE(2,235) '4(',F12,' ',F13,' ',F23,' ',F22,')'
    FAD1=F_LEFT_RECT(JP,KP)
    WRITE(3,231) FAD1
    OWN=CELL_ID_RECT(IP,JP,KP)
    WRITE(4,231) OWN
    ENDDO   
    ENDDO    
        
        !BACK FACE--PROCESSOR_MI1_TO_MI2  !CASE RECTANGLE 1-3-1
    KP=CNL
    K=KP+CNL*(MJ2-1)
    DO JP=1,CNL
    J=JP+CNL*(MJ1-1)
    DO IP=1,CNI
    I=IP+CNI*(MI-1)
    F11=(CNL+1)*(CNI+1)*(KP-1)+(CNI+1)*(JP-1)+IP-1
    F12=F11+1
    F13=F12+CNI+1
    F14=F13-1
    F21=F11+(CNL+1)*(CNI+1)
    F22=F21+1
    F23=F22+CNI+1
    F24=F23-1
    CALL CAL_FIJ_RECT(F11,F12,F13,F14,F21,F22,F23,F24)
    WRITE(2,235) '4(',F21,' ',F22,' ',F23,' ',F24,')'
    FAD1=F_BACK_RECT(IP,JP)
    WRITE(3,231) FAD1
    OWN=CELL_ID_RECT(IP,JP,KP)
    WRITE(4,231) OWN
    ENDDO   
    ENDDO    
        
    N=4
    WRITE(31,231) N+4
    WRITE(31,*) '('  
    WRITE(31,*) '0'
    WRITE(31,*) '1'
    WRITE(31,*) '2' 
    WRITE(31,*) '3' 
    DO I=1,N
    WRITE(31,*) '-1' 
    ENDDO
    WRITE(31,*) ')'
    
    WRITE(21,231) N+4
    WRITE(21,*) "("    
    WRITE(21,*) "    INLET1"
    WRITE(21,*) "    {"    
    WRITE(21,*) "        type            patch;"
    N=CNL*CNL
    WRITE(21,232) "        nFaces          ",N,";"
    N=F_INTERNAL_NUM
    WRITE(21,232) "        startFace       ",N,";"
    WRITE(21,*) "    }"
    WRITE(21,*) "    INLET2"
    WRITE(21,*) "    {"    
    WRITE(21,*) "        type            patch;"
    N=0
    WRITE(21,232) "        nFaces          ",N,";"
    N=F_INTERNAL_NUM+CNL*CNL
    WRITE(21,232) "        startFace       ",N,";"
    WRITE(21,*) "    }"
    WRITE(21,*) "    OUTLET"
    WRITE(21,*) "    {"
    WRITE(21,*) "        type            patch;"
    N=0
    WRITE(21,232) "        nFaces          ",N,";"
    N=F_INTERNAL_NUM+CNL*CNL
    WRITE(21,232) "        startFace       ",N,";"
    WRITE(21,*) "    }" 
    WRITE(21,*) "    HULL"
    WRITE(21,*) "    {"
    WRITE(21,*) "        type            wall;"
    WRITE(21,*) "        inGroups        1(wall);"
    N=0
    WRITE(21,232) "        nFaces          ",N,";"
    N=F_INTERNAL_NUM+CNL*CNL
    WRITE(21,232) "        startFace       ",N,";"
    WRITE(21,*) "    }"
       
        !FRONT/UP FACE--PROCESSOR_MI1_TO_MI2  !CASE RECTANGLE 1-3-1
    PRO1=PNI*PNJ+PNL*PNI*(MJ2-1)+PNI*(MJ1-1)+MI-1
    MJ=1
    PRO2=PNI*(MJ-1)+MI-1
    WRITE(21,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(21,*) '    {'
    WRITE(21,*) '        type            processor;'
    WRITE(21,*) '        inGroups        List<word> 1(processor);'
    N=CNI*CNL*2
    WRITE(21,232) "        nFaces          ",N,";"
    N=F_INTERNAL_NUM+CNL*CNL
    WRITE(21,232) "        startFace       ",N,";"
    WRITE(21,*) '        matchTolerance  0.0001;'
    WRITE(21,232) "        myProcNo        ",PRO1,";"
    WRITE(21,232) "        neighbProcNo    ",PRO2,";"
    WRITE(21,*) '    }'
    WRITE(41,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(41,*) '    {'
    WRITE(41,*) '        type            processor;'
    WRITE(41,*) '        value           nonuniform List<scalar> '
    N=CNI*CNL*2
    WRITE(41,'(I0)') N
    WRITE(41,*) '('
    DO JP=1,CNL
    DO IP=1,CNI
    WRITE(41,'(E21.14)') P_FRONT_RECT(IP,JP)    
    ENDDO
    ENDDO
    DO KP=1,CNL
    DO IP=1,CNI
    WRITE(41,'(E21.14)') P_UP_RECT(IP,KP)    
    ENDDO
    ENDDO
    WRITE(41,*) ')'
    WRITE(41,*) ';'
    WRITE(41,*) '    }'
    WRITE(51,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(51,*) '    {'
    WRITE(51,*) '        type            processor;'
    WRITE(51,*) '        value           nonuniform List<vector> '
    N=CNI*CNL*2
    WRITE(51,'(I0)') N
    WRITE(51,*) '('
    DO JP=1,CNL
    DO IP=1,CNI
    UU1=U_FRONT_RECT(1,IP,JP)
    UU2=U_FRONT_RECT(2,IP,JP)
    UU3=U_FRONT_RECT(3,IP,JP)
    WRITE(51,230) '( ',UU1,' ',UU2,' ',UU3,' )'  
    ENDDO
    ENDDO
    DO KP=1,CNL
    DO IP=1,CNI
    UU1=U_UP_RECT(1,IP,KP)
    UU2=U_UP_RECT(2,IP,KP)
    UU3=U_UP_RECT(3,IP,KP)
    WRITE(51,230) '( ',UU1,' ',UU2,' ',UU3,' )'  
    ENDDO
    ENDDO
    WRITE(51,*) ')'
    WRITE(51,*) ';'
    WRITE(51,*) '    }'
    WRITE(61,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(61,*) '    {'
    WRITE(61,*) '        type            processor;'
    WRITE(61,*) "        value           uniform 0;"
    WRITE(61,*) '    }'
    
        !DOWN FACE--PROCESSOR_MI1_TO_MI2  !CASE RECTANGLE 1-3-1
    PRO2=PNI*PNJ+PNL*PNI*(MJ2-1)+PNI*(MJ1-2)+MI-1
    WRITE(21,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(21,*) '    {'
    WRITE(21,*) '        type            processor;'
    WRITE(21,*) '        inGroups        List<word> 1(processor);'
    N=CNI*CNL
    WRITE(21,232) "        nFaces          ",N,";"
    N=F_INTERNAL_NUM+CNL*CNL+CNI*CNL*2
    WRITE(21,232) "        startFace       ",N,";"
    WRITE(21,*) '        matchTolerance  0.0001;'
    WRITE(21,232) "        myProcNo        ",PRO1,";"
    WRITE(21,232) "        neighbProcNo    ",PRO2,";"
    WRITE(21,*) '    }'
    WRITE(41,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(41,*) '    {'
    WRITE(41,*) '        type            processor;'
    WRITE(41,*) '        value           nonuniform List<scalar> '
    N=CNI*CNL
    WRITE(41,'(I0)') N
    WRITE(41,*) '('
    DO KP=1,CNL
    DO IP=1,CNI
    WRITE(41,'(E21.14)') P_DOWN_RECT(IP,KP)    
    ENDDO
    ENDDO
    WRITE(41,*) ')'
    WRITE(41,*) ';'
    WRITE(41,*) '    }'
    WRITE(51,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(51,*) '    {'
    WRITE(51,*) '        type            processor;'
    WRITE(51,*) '        value           nonuniform List<vector> '
    N=CNI*CNL
    WRITE(51,'(I0)') N
    WRITE(51,*) '('
    DO KP=1,CNL
    DO IP=1,CNI
    UU1=U_DOWN_RECT(1,IP,KP)
    UU2=U_DOWN_RECT(2,IP,KP)
    UU3=U_DOWN_RECT(3,IP,KP)
    WRITE(51,230) '( ',UU1,' ',UU2,' ',UU3,' )'  
    ENDDO
    ENDDO
    WRITE(51,*) ')'
    WRITE(51,*) ';'
    WRITE(51,*) '    }'
    WRITE(61,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(61,*) '    {'
    WRITE(61,*) '        type            processor;'
    WRITE(61,*) "        value           uniform 0;"
    WRITE(61,*) '    }'
    
        !LEFT FACE--PROCESSOR_MI1_TO_MI2  !CASE RECTANGLE 1-3-1
    PRO2=PNI*PNJ+PNL*PNI*(MJ2-1)+PNI*(MJ1-1)+MI
    WRITE(21,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(21,*) '    {'
    WRITE(21,*) '        type            processor;'
    WRITE(21,*) '        inGroups        List<word> 1(processor);'
    N=CNL*CNL
    WRITE(21,232) "        nFaces          ",N,";"
    N=F_INTERNAL_NUM+CNL*CNL+CNI*CNL*2+CNI*CNL
    WRITE(21,232) "        startFace       ",N,";"
    WRITE(21,*) '        matchTolerance  0.0001;'
    WRITE(21,232) "        myProcNo        ",PRO1,";"
    WRITE(21,232) "        neighbProcNo    ",PRO2,";"
    WRITE(21,*) '    }'
    WRITE(41,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(41,*) '    {'
    WRITE(41,*) '        type            processor;'
    WRITE(41,*) '        value           nonuniform List<scalar> '
    N=CNL*CNL
    WRITE(41,'(I0)') N
    WRITE(41,*) '('
    DO KP=1,CNL
    DO JP=1,CNL
    WRITE(41,'(E21.14)') P_LEFT_RECT(JP,KP)    
    ENDDO
    ENDDO
    WRITE(41,*) ')'
    WRITE(41,*) ';'
    WRITE(41,*) '    }'
    WRITE(51,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(51,*) '    {'
    WRITE(51,*) '        type            processor;'
    WRITE(51,*) '        value           nonuniform List<vector> '
    N=CNL*CNL
    WRITE(51,'(I0)') N
    WRITE(51,*) '('
    DO KP=1,CNL
    DO JP=1,CNL
    UU1=U_LEFT_RECT(1,JP,KP)
    UU2=U_LEFT_RECT(2,JP,KP)
    UU3=U_LEFT_RECT(3,JP,KP)
    WRITE(51,230) '( ',UU1,' ',UU2,' ',UU3,' )'  
    ENDDO
    ENDDO
    WRITE(51,*) ')'
    WRITE(51,*) ';'
    WRITE(51,*) '    }'
    WRITE(61,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(61,*) '    {'
    WRITE(61,*) '        type            processor;'
    WRITE(61,*) "        value           uniform 0;"
    WRITE(61,*) '    }'
    
        !BACK FACE--PROCESSOR_MI1_TO_MI2  !CASE RECTANGLE 1-3-1
    PRO2=PNI*PNJ+PNL*PNI*MJ2+PNI*(MJ1-1)+MI-1
    WRITE(21,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(21,*) '    {'
    WRITE(21,*) '        type            processor;'
    WRITE(21,*) '        inGroups        List<word> 1(processor);'
    N=CNI*CNL
    WRITE(21,232) "        nFaces          ",N,";"
    N=F_INTERNAL_NUM+CNL*CNL+CNI*CNL*2+CNI*CNL+CNL*CNL
    WRITE(21,232) "        startFace       ",N,";"
    WRITE(21,*) '        matchTolerance  0.0001;'
    WRITE(21,232) "        myProcNo        ",PRO1,";"
    WRITE(21,232) "        neighbProcNo    ",PRO2,";"
    WRITE(21,*) '    }'
    WRITE(41,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(41,*) '    {'
    WRITE(41,*) '        type            processor;'
    WRITE(41,*) '        value           nonuniform List<scalar> '
    N=CNI*CNL
    WRITE(41,'(I0)') N
    WRITE(41,*) '('
    DO JP=1,CNL
    DO IP=1,CNI
    WRITE(41,'(E21.14)') P_BACK_RECT(IP,JP)    
    ENDDO
    ENDDO
    WRITE(41,*) ')'
    WRITE(41,*) ';'
    WRITE(41,*) '    }'
    WRITE(51,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(51,*) '    {'
    WRITE(51,*) '        type            processor;'
    WRITE(51,*) '        value           nonuniform List<vector> '
    N=CNI*CNL
    WRITE(51,'(I0)') N
    WRITE(51,*) '('
    DO JP=1,CNL
    DO IP=1,CNI
    UU1=U_BACK_RECT(1,IP,JP)
    UU2=U_BACK_RECT(2,IP,JP)
    UU3=U_BACK_RECT(3,IP,JP)
    WRITE(51,230) '( ',UU1,' ',UU2,' ',UU3,' )'  
    ENDDO
    ENDDO
    WRITE(51,*) ')'
    WRITE(51,*) ';'
    WRITE(51,*) '    }'
    WRITE(61,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(61,*) '    {'
    WRITE(61,*) '        type            processor;'
    WRITE(61,*) "        value           uniform 0;"
    WRITE(61,*) '    }'    
        
        
        
        
        
        
        
        
        
        
        
        
        
    ELSEIF (MJ2.EQ.PNL) THEN
    !CASE RECTANGLE 1-3-3
    
        !RIGHT FACE--INLET1  !CASE RECTANGLE 1-3-3 
    IP=1
    I=IP+CNI*(MI-1)
    DO KP=1,CNL
    K=KP+CNL*(MJ2-1)
    DO JP=1,CNL
    J=JP+CNL*(MJ1-1)
    F11=(CNL+1)*(CNI+1)*(KP-1)+(CNI+1)*(JP-1)+IP-1
    F12=F11+1
    F13=F12+CNI+1
    F14=F13-1
    F21=F11+(CNL+1)*(CNI+1)
    F22=F21+1
    F23=F22+CNI+1
    F24=F23-1
    CALL CAL_FIJ_RECT(F11,F12,F13,F14,F21,F22,F23,F24)
    WRITE(2,235) '4(',F11,' ',F21,' ',F24,' ',F14,')'
    FAD1=F_INTERNAL_NUM_ALL+NL*(K-1)+J
    WRITE(3,231) FAD1
    OWN=CELL_ID_RECT(IP,JP,KP)
    WRITE(4,231) OWN
    ENDDO   
    ENDDO    
        
        !UP FACE--PROCESSOR_MI1_TO_MI2  !CASE RECTANGLE 1-3-3
    JP=CNL
    J=JP+CNL*(MJ1-1)
    DO KP=1,CNL
    K=KP+CNL*(MJ2-1)
    DO IP=1,CNI
    I=IP+CNI*(MI-1)
    F11=(CNL+1)*(CNI+1)*(KP-1)+(CNI+1)*(JP-1)+IP-1
    F12=F11+1
    F13=F12+CNI+1
    F14=F13-1
    F21=F11+(CNL+1)*(CNI+1)
    F22=F21+1
    F23=F22+CNI+1
    F24=F23-1
    CALL CAL_FIJ_RECT(F11,F12,F13,F14,F21,F22,F23,F24)
    WRITE(2,235) '4(',F14,' ',F24,' ',F23,' ',F13,')'
    FAD1=F_UP_RECT(IP,KP)
    WRITE(3,231) FAD1
    OWN=CELL_ID_RECT(IP,JP,KP)
    WRITE(4,231) OWN
    ENDDO   
    ENDDO
    
        !BACK FACE--PROCESSOR_MI1_TO_MI2  !CASE RECTANGLE 1-3-3
    KP=CNL
    K=KP+CNL*(MJ2-1)
    DO JP=CNL,1,-1
    J=JP+CNL*(MJ1-1)
    DO IP=1,CNI
    I=IP+CNI*(MI-1)
    F11=(CNL+1)*(CNI+1)*(KP-1)+(CNI+1)*(JP-1)+IP-1
    F12=F11+1
    F13=F12+CNI+1
    F14=F13-1
    F21=F11+(CNL+1)*(CNI+1)
    F22=F21+1
    F23=F22+CNI+1
    F24=F23-1
    CALL CAL_FIJ_RECT(F11,F12,F13,F14,F21,F22,F23,F24)
    !WRITE(2,235) '4(',F21,' ',F22,' ',F23,' ',F24,')'
    WRITE(2,235) '4(',F24,' ',F21,' ',F22,' ',F23,')'
    FAD1=F_BACK_RECT(IP,JP)
    WRITE(3,231) FAD1
    OWN=CELL_ID_RECT(IP,JP,KP)
    WRITE(4,231) OWN
    ENDDO   
    ENDDO   
        
        !FRONT FACE--PROCESSOR_MI1_TO_MI2  !CASE RECTANGLE 1-3-3
    KP=1
    K=KP+CNL*(MJ2-1)
    DO JP=1,CNL
    J=JP+CNL*(MJ1-1)
    DO IP=1,CNI
    I=IP+CNI*(MI-1)
    F11=(CNL+1)*(CNI+1)*(KP-1)+(CNI+1)*(JP-1)+IP-1
    F12=F11+1
    F13=F12+CNI+1
    F14=F13-1
    F21=F11+(CNL+1)*(CNI+1)
    F22=F21+1
    F23=F22+CNI+1
    F24=F23-1
    CALL CAL_FIJ_RECT(F11,F12,F13,F14,F21,F22,F23,F24)
    WRITE(2,235) '4(',F11,' ',F14,' ',F13,' ',F12,')'
    FAD1=F_FRONT_RECT(IP,JP)
    WRITE(3,231) FAD1
    OWN=CELL_ID_RECT(IP,JP,KP)
    WRITE(4,231) OWN
    ENDDO   
    ENDDO
    
        !DOWN FACE--PROCESSOR_MI1_TO_MI2  !CASE RECTANGLE 1-3-3
    JP=1
    J=JP+CNL*(MJ1-1)
    DO KP=1,CNL
    K=KP+CNL*(MJ2-1)
    DO IP=1,CNI
    I=IP+CNI*(MI-1)
    F11=(CNL+1)*(CNI+1)*(KP-1)+(CNI+1)*(JP-1)+IP-1
    F12=F11+1
    F13=F12+CNI+1
    F14=F13-1
    F21=F11+(CNL+1)*(CNI+1)
    F22=F21+1
    F23=F22+CNI+1
    F24=F23-1
    CALL CAL_FIJ_RECT(F11,F12,F13,F14,F21,F22,F23,F24)
    WRITE(2,235) '4(',F11,' ',F12,' ',F22,' ',F21,')'
    FAD1=F_DOWN_RECT(IP,KP)
    WRITE(3,231) FAD1
    OWN=CELL_ID_RECT(IP,JP,KP)
    WRITE(4,231) OWN
    ENDDO   
    ENDDO
    
        !LEFT FACE--PROCESSOR_MI1_TO_MI2  !CASE RECTANGLE 1-3-3 
    IP=CNI
    I=IP+CNI*(MI-1)
    DO KP=1,CNL
    K=KP+CNL*(MJ2-1)
    DO JP=1,CNL
    J=JP+CNL*(MJ1-1)
    F11=(CNL+1)*(CNI+1)*(KP-1)+(CNI+1)*(JP-1)+IP-1
    F12=F11+1
    F13=F12+CNI+1
    F14=F13-1
    F21=F11+(CNL+1)*(CNI+1)
    F22=F21+1
    F23=F22+CNI+1
    F24=F23-1
    CALL CAL_FIJ_RECT(F11,F12,F13,F14,F21,F22,F23,F24)
    WRITE(2,235) '4(',F12,' ',F13,' ',F23,' ',F22,')'
    FAD1=F_LEFT_RECT(JP,KP)
    WRITE(3,231) FAD1
    OWN=CELL_ID_RECT(IP,JP,KP)
    WRITE(4,231) OWN
    ENDDO   
    ENDDO    
        
    N=4
    WRITE(31,231) N+4
    WRITE(31,*) '('  
    WRITE(31,*) '0'
    WRITE(31,*) '1'
    WRITE(31,*) '2' 
    WRITE(31,*) '3' 
    DO I=1,N
    WRITE(31,*) '-1' 
    ENDDO
    WRITE(31,*) ')'
    
    WRITE(21,231) N+4
    WRITE(21,*) "("    
    WRITE(21,*) "    INLET1"
    WRITE(21,*) "    {"    
    WRITE(21,*) "        type            patch;"
    N=CNL*CNL
    WRITE(21,232) "        nFaces          ",N,";"
    N=F_INTERNAL_NUM
    WRITE(21,232) "        startFace       ",N,";"
    WRITE(21,*) "    }"
    WRITE(21,*) "    INLET2"
    WRITE(21,*) "    {"    
    WRITE(21,*) "        type            patch;"
    N=0
    WRITE(21,232) "        nFaces          ",N,";"
    N=F_INTERNAL_NUM+CNL*CNL
    WRITE(21,232) "        startFace       ",N,";"
    WRITE(21,*) "    }"
    WRITE(21,*) "    OUTLET"
    WRITE(21,*) "    {"
    WRITE(21,*) "        type            patch;"
    N=0
    WRITE(21,232) "        nFaces          ",N,";"
    N=F_INTERNAL_NUM+CNL*CNL
    WRITE(21,232) "        startFace       ",N,";"
    WRITE(21,*) "    }" 
    WRITE(21,*) "    HULL"
    WRITE(21,*) "    {"
    WRITE(21,*) "        type            wall;"
    WRITE(21,*) "        inGroups        1(wall);"
    N=0
    WRITE(21,232) "        nFaces          ",N,";"
    N=F_INTERNAL_NUM+CNL*CNL
    WRITE(21,232) "        startFace       ",N,";"
    WRITE(21,*) "    }"
       
        !UP/BACK FACE--PROCESSOR_MI1_TO_MI2  !CASE RECTANGLE 1-3-3
    PRO1=PNI*PNJ+PNL*PNI*(MJ2-1)+PNI*(MJ1-1)+MI-1
    MJ=1
    PRO2=PNI*(MJ-1)+MI-1
    WRITE(21,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(21,*) '    {'
    WRITE(21,*) '        type            processor;'
    WRITE(21,*) '        inGroups        List<word> 1(processor);'
    N=CNI*CNL*2
    WRITE(21,232) "        nFaces          ",N,";"
    N=F_INTERNAL_NUM+CNL*CNL
    WRITE(21,232) "        startFace       ",N,";"
    WRITE(21,*) '        matchTolerance  0.0001;'
    WRITE(21,232) "        myProcNo        ",PRO1,";"
    WRITE(21,232) "        neighbProcNo    ",PRO2,";"
    WRITE(21,*) '    }'
    WRITE(41,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(41,*) '    {'
    WRITE(41,*) '        type            processor;'
    WRITE(41,*) '        value           nonuniform List<scalar> '
    N=CNI*CNL*2
    WRITE(41,'(I0)') N
    WRITE(41,*) '('
    DO KP=1,CNL
    DO IP=1,CNI
    WRITE(41,'(E21.14)') P_UP_RECT(IP,KP)    
    ENDDO
    ENDDO
    DO JP=1,CNL
    DO IP=1,CNI
    WRITE(41,'(E21.14)') P_BACK_RECT(IP,JP)    
    ENDDO
    ENDDO
    WRITE(41,*) ')'
    WRITE(41,*) ';'
    WRITE(41,*) '    }'
    WRITE(51,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(51,*) '    {'
    WRITE(51,*) '        type            processor;'
    WRITE(51,*) '        value           nonuniform List<vector> '
    N=CNI*CNL*2
    WRITE(51,'(I0)') N
    WRITE(51,*) '('
    DO KP=1,CNL
    DO IP=1,CNI
    UU1=U_UP_RECT(1,IP,KP)
    UU2=U_UP_RECT(2,IP,KP)
    UU3=U_UP_RECT(3,IP,KP)
    WRITE(51,230) '( ',UU1,' ',UU2,' ',UU3,' )'  
    ENDDO
    ENDDO
    DO JP=1,CNL
    DO IP=1,CNI
    UU1=U_BACK_RECT(1,IP,JP)
    UU2=U_BACK_RECT(2,IP,JP)
    UU3=U_BACK_RECT(3,IP,JP)
    WRITE(51,230) '( ',UU1,' ',UU2,' ',UU3,' )'  
    ENDDO
    ENDDO
    WRITE(51,*) ')'
    WRITE(51,*) ';'
    WRITE(51,*) '    }'
    WRITE(61,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(61,*) '    {'
    WRITE(61,*) '        type            processor;'
    WRITE(61,*) "        value           uniform 0;"
    WRITE(61,*) '    }'
    !����ɾ��������Ӧ���Ƕ����
    !DO JP=1,CNL
    !DO IP=1,CNI
    !WRITE(41,'(E21.14)') P_BACK_RECT(IP,JP)    
    !ENDDO
    !ENDDO
    
        !FRONT FACE--PROCESSOR_MI1_TO_MI2  !CASE RECTANGLE 1-3-3
    PRO2=PNI*PNJ+PNL*PNI*(MJ2-2)+PNI*(MJ1-1)+MI-1
    WRITE(21,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(21,*) '    {'
    WRITE(21,*) '        type            processor;'
    WRITE(21,*) '        inGroups        List<word> 1(processor);'
    N=CNI*CNL
    WRITE(21,232) "        nFaces          ",N,";"
    N=F_INTERNAL_NUM+CNL*CNL+CNI*CNL*2
    WRITE(21,232) "        startFace       ",N,";"
    WRITE(21,*) '        matchTolerance  0.0001;'
    WRITE(21,232) "        myProcNo        ",PRO1,";"
    WRITE(21,232) "        neighbProcNo    ",PRO2,";"
    WRITE(21,*) '    }'
    WRITE(41,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(41,*) '    {'
    WRITE(41,*) '        type            processor;'
    WRITE(41,*) '        value           nonuniform List<scalar> '
    N=CNI*CNL
    WRITE(41,'(I0)') N
    WRITE(41,*) '('
    DO JP=1,CNL
    DO IP=1,CNI
    WRITE(41,'(E21.14)') P_FRONT_RECT(IP,JP)    
    ENDDO
    ENDDO
    WRITE(41,*) ')'
    WRITE(41,*) ';'
    WRITE(41,*) '    }'
    WRITE(51,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(51,*) '    {'
    WRITE(51,*) '        type            processor;'
    WRITE(51,*) '        value           nonuniform List<vector> '
    N=CNI*CNL
    WRITE(51,'(I0)') N
    WRITE(51,*) '('
    DO JP=1,CNL
    DO IP=1,CNI
    UU1=U_FRONT_RECT(1,IP,JP)
    UU2=U_FRONT_RECT(2,IP,JP)
    UU3=U_FRONT_RECT(3,IP,JP)
    WRITE(51,230) '( ',UU1,' ',UU2,' ',UU3,' )'  
    ENDDO
    ENDDO
    WRITE(51,*) ')'
    WRITE(51,*) ';'
    WRITE(51,*) '    }'
    WRITE(61,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(61,*) '    {'
    WRITE(61,*) '        type            processor;'
    WRITE(61,*) "        value           uniform 0;"
    WRITE(61,*) '    }'
    
        !DOWN FACE--PROCESSOR_MI1_TO_MI2  !CASE RECTANGLE 1-3-3
    PRO2=PNI*PNJ+PNL*PNI*(MJ2-1)+PNI*(MJ1-2)+MI-1
    WRITE(21,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(21,*) '    {'
    WRITE(21,*) '        type            processor;'
    WRITE(21,*) '        inGroups        List<word> 1(processor);'
    N=CNI*CNL
    WRITE(21,232) "        nFaces          ",N,";"
    N=F_INTERNAL_NUM+CNL*CNL+CNI*CNL*2+CNI*CNL
    WRITE(21,232) "        startFace       ",N,";"
    WRITE(21,*) '        matchTolerance  0.0001;'
    WRITE(21,232) "        myProcNo        ",PRO1,";"
    WRITE(21,232) "        neighbProcNo    ",PRO2,";"
    WRITE(21,*) '    }'
    WRITE(41,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(41,*) '    {'
    WRITE(41,*) '        type            processor;'
    WRITE(41,*) '        value           nonuniform List<scalar> '
    N=CNI*CNL
    WRITE(41,'(I0)') N
    WRITE(41,*) '('
    DO KP=1,CNL
    DO IP=1,CNI
    WRITE(41,'(E21.14)') P_DOWN_RECT(IP,KP)    
    ENDDO
    ENDDO
    WRITE(41,*) ')'
    WRITE(41,*) ';'
    WRITE(41,*) '    }'
    WRITE(51,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(51,*) '    {'
    WRITE(51,*) '        type            processor;'
    WRITE(51,*) '        value           nonuniform List<vector> '
    N=CNI*CNL
    WRITE(51,'(I0)') N
    WRITE(51,*) '('
    DO KP=1,CNL
    DO IP=1,CNI
    UU1=U_DOWN_RECT(1,IP,KP)
    UU2=U_DOWN_RECT(2,IP,KP)
    UU3=U_DOWN_RECT(3,IP,KP)
    WRITE(51,230) '( ',UU1,' ',UU2,' ',UU3,' )'  
    ENDDO
    ENDDO
    WRITE(51,*) ')'
    WRITE(51,*) ';'
    WRITE(51,*) '    }'
    WRITE(61,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(61,*) '    {'
    WRITE(61,*) '        type            processor;'
    WRITE(61,*) "        value           uniform 0;"
    WRITE(61,*) '    }'
    
        !LEFT FACE--PROCESSOR_MI1_TO_MI2  !CASE RECTANGLE 1-3-3
    PRO2=PNI*PNJ+PNL*PNI*(MJ2-1)+PNI*(MJ1-1)+MI
    WRITE(21,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(21,*) '    {'
    WRITE(21,*) '        type            processor;'
    WRITE(21,*) '        inGroups        List<word> 1(processor);'
    N=CNL*CNL
    WRITE(21,232) "        nFaces          ",N,";"
    N=F_INTERNAL_NUM+CNL*CNL+CNI*CNL*2+CNI*CNL+CNI*CNL
    WRITE(21,232) "        startFace       ",N,";"
    WRITE(21,*) '        matchTolerance  0.0001;'
    WRITE(21,232) "        myProcNo        ",PRO1,";"
    WRITE(21,232) "        neighbProcNo    ",PRO2,";"
    WRITE(21,*) '    }'
    WRITE(41,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(41,*) '    {'
    WRITE(41,*) '        type            processor;'
    WRITE(41,*) '        value           nonuniform List<scalar> '
    N=CNL*CNL
    WRITE(41,'(I0)') N
    WRITE(41,*) '('
    DO KP=1,CNL
    DO JP=1,CNL
    WRITE(41,'(E21.14)') P_LEFT_RECT(JP,KP)    
    ENDDO
    ENDDO
    WRITE(41,*) ')'
    WRITE(41,*) ';'
    WRITE(41,*) '    }'
    WRITE(51,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(51,*) '    {'
    WRITE(51,*) '        type            processor;'
    WRITE(51,*) '        value           nonuniform List<vector> '
    N=CNL*CNL
    WRITE(51,'(I0)') N
    WRITE(51,*) '('
    DO KP=1,CNL
    DO JP=1,CNL
    UU1=U_LEFT_RECT(1,JP,KP)
    UU2=U_LEFT_RECT(2,JP,KP)
    UU3=U_LEFT_RECT(3,JP,KP)
    WRITE(51,230) '( ',UU1,' ',UU2,' ',UU3,' )'  
    ENDDO
    ENDDO
    WRITE(51,*) ')'
    WRITE(51,*) ';'
    WRITE(51,*) '    }'
    WRITE(61,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(61,*) '    {'
    WRITE(61,*) '        type            processor;'
    WRITE(61,*) "        value           uniform 0;"
    WRITE(61,*) '    }'    
        
        
        
        
        
        
        
        
        
        
        
        
        
    ELSE
    !CASE RECTANGLE 1-3-2
    
        !RIGHT FACE--INLET1  !CASE RECTANGLE 1-3-2 
    IP=1
    I=IP+CNI*(MI-1)
    DO KP=1,CNL
    K=KP+CNL*(MJ2-1)
    DO JP=1,CNL
    J=JP+CNL*(MJ1-1)
    F11=(CNL+1)*(CNI+1)*(KP-1)+(CNI+1)*(JP-1)+IP-1
    F12=F11+1
    F13=F12+CNI+1
    F14=F13-1
    F21=F11+(CNL+1)*(CNI+1)
    F22=F21+1
    F23=F22+CNI+1
    F24=F23-1
    CALL CAL_FIJ_RECT(F11,F12,F13,F14,F21,F22,F23,F24)
    WRITE(2,235) '4(',F11,' ',F21,' ',F24,' ',F14,')'
    FAD1=F_INTERNAL_NUM_ALL+NL*(K-1)+J
    WRITE(3,231) FAD1
    OWN=CELL_ID_RECT(IP,JP,KP)
    WRITE(4,231) OWN
    ENDDO   
    ENDDO    
        
        !UP FACE--PROCESSOR_MI1_TO_MI2  !CASE RECTANGLE 1-3-2
    JP=CNL
    J=JP+CNL*(MJ1-1)
    DO KP=1,CNL
    K=KP+CNL*(MJ2-1)
    DO IP=1,CNI
    I=IP+CNI*(MI-1)
    F11=(CNL+1)*(CNI+1)*(KP-1)+(CNI+1)*(JP-1)+IP-1
    F12=F11+1
    F13=F12+CNI+1
    F14=F13-1
    F21=F11+(CNL+1)*(CNI+1)
    F22=F21+1
    F23=F22+CNI+1
    F24=F23-1
    CALL CAL_FIJ_RECT(F11,F12,F13,F14,F21,F22,F23,F24)
    WRITE(2,235) '4(',F14,' ',F24,' ',F23,' ',F13,')'
    FAD1=F_UP_RECT(IP,KP)
    WRITE(3,231) FAD1
    OWN=CELL_ID_RECT(IP,JP,KP)
    WRITE(4,231) OWN
    ENDDO   
    ENDDO    
        
        !FRONT FACE--PROCESSOR_MI1_TO_MI2  !CASE RECTANGLE 1-3-2
    KP=1
    K=KP+CNL*(MJ2-1)
    DO JP=1,CNL
    J=JP+CNL*(MJ1-1)
    DO IP=1,CNI
    I=IP+CNI*(MI-1)
    F11=(CNL+1)*(CNI+1)*(KP-1)+(CNI+1)*(JP-1)+IP-1
    F12=F11+1
    F13=F12+CNI+1
    F14=F13-1
    F21=F11+(CNL+1)*(CNI+1)
    F22=F21+1
    F23=F22+CNI+1
    F24=F23-1
    CALL CAL_FIJ_RECT(F11,F12,F13,F14,F21,F22,F23,F24)
    WRITE(2,235) '4(',F11,' ',F14,' ',F13,' ',F12,')'
    FAD1=F_FRONT_RECT(IP,JP)
    WRITE(3,231) FAD1
    OWN=CELL_ID_RECT(IP,JP,KP)
    WRITE(4,231) OWN
    ENDDO   
    ENDDO
    
        !DOWN FACE--PROCESSOR_MI1_TO_MI2  !CASE RECTANGLE 1-3-2
    JP=1
    J=JP+CNL*(MJ1-1)
    DO KP=1,CNL
    K=KP+CNL*(MJ2-1)
    DO IP=1,CNI
    I=IP+CNI*(MI-1)
    F11=(CNL+1)*(CNI+1)*(KP-1)+(CNI+1)*(JP-1)+IP-1
    F12=F11+1
    F13=F12+CNI+1
    F14=F13-1
    F21=F11+(CNL+1)*(CNI+1)
    F22=F21+1
    F23=F22+CNI+1
    F24=F23-1
    CALL CAL_FIJ_RECT(F11,F12,F13,F14,F21,F22,F23,F24)
    WRITE(2,235) '4(',F11,' ',F12,' ',F22,' ',F21,')'
    FAD1=F_DOWN_RECT(IP,KP)
    WRITE(3,231) FAD1
    OWN=CELL_ID_RECT(IP,JP,KP)
    WRITE(4,231) OWN
    ENDDO   
    ENDDO
    
        !LEFT FACE--PROCESSOR_MI1_TO_MI2  !CASE RECTANGLE 1-3-2 
    IP=CNI
    I=IP+CNI*(MI-1)
    DO KP=1,CNL
    K=KP+CNL*(MJ2-1)
    DO JP=1,CNL
    J=JP+CNL*(MJ1-1)
    F11=(CNL+1)*(CNI+1)*(KP-1)+(CNI+1)*(JP-1)+IP-1
    F12=F11+1
    F13=F12+CNI+1
    F14=F13-1
    F21=F11+(CNL+1)*(CNI+1)
    F22=F21+1
    F23=F22+CNI+1
    F24=F23-1
    CALL CAL_FIJ_RECT(F11,F12,F13,F14,F21,F22,F23,F24)
    WRITE(2,235) '4(',F12,' ',F13,' ',F23,' ',F22,')'
    FAD1=F_LEFT_RECT(JP,KP)
    WRITE(3,231) FAD1
    OWN=CELL_ID_RECT(IP,JP,KP)
    WRITE(4,231) OWN
    ENDDO   
    ENDDO
    
        !BACK FACE--PROCESSOR_MI1_TO_MI2  !CASE RECTANGLE 1-3-2
    KP=CNL
    K=KP+CNL*(MJ2-1)
    DO JP=1,CNL
    J=JP+CNL*(MJ1-1)
    DO IP=1,CNI
    I=IP+CNI*(MI-1)
    F11=(CNL+1)*(CNI+1)*(KP-1)+(CNI+1)*(JP-1)+IP-1
    F12=F11+1
    F13=F12+CNI+1
    F14=F13-1
    F21=F11+(CNL+1)*(CNI+1)
    F22=F21+1
    F23=F22+CNI+1
    F24=F23-1
    CALL CAL_FIJ_RECT(F11,F12,F13,F14,F21,F22,F23,F24)
    WRITE(2,235) '4(',F21,' ',F22,' ',F23,' ',F24,')'
    FAD1=F_BACK_RECT(IP,JP)
    WRITE(3,231) FAD1
    OWN=CELL_ID_RECT(IP,JP,KP)
    WRITE(4,231) OWN
    ENDDO   
    ENDDO    
        
    N=5
    WRITE(31,231) N+4
    WRITE(31,*) '('  
    WRITE(31,*) '0'
    WRITE(31,*) '1'
    WRITE(31,*) '2' 
    WRITE(31,*) '3' 
    DO I=1,N
    WRITE(31,*) '-1' 
    ENDDO
    WRITE(31,*) ')'
    
    WRITE(21,231) N+4
    WRITE(21,*) "("    
    WRITE(21,*) "    INLET1"
    WRITE(21,*) "    {"    
    WRITE(21,*) "        type            patch;"
    N=CNL*CNL
    WRITE(21,232) "        nFaces          ",N,";"
    N=F_INTERNAL_NUM
    WRITE(21,232) "        startFace       ",N,";"
    WRITE(21,*) "    }"
    WRITE(21,*) "    INLET2"
    WRITE(21,*) "    {"    
    WRITE(21,*) "        type            patch;"
    N=0
    WRITE(21,232) "        nFaces          ",N,";"
    N=F_INTERNAL_NUM+CNL*CNL
    WRITE(21,232) "        startFace       ",N,";"
    WRITE(21,*) "    }"
    WRITE(21,*) "    OUTLET"
    WRITE(21,*) "    {"
    WRITE(21,*) "        type            patch;"
    N=0
    WRITE(21,232) "        nFaces          ",N,";"
    N=F_INTERNAL_NUM+CNL*CNL
    WRITE(21,232) "        startFace       ",N,";"
    WRITE(21,*) "    }" 
    WRITE(21,*) "    HULL"
    WRITE(21,*) "    {"
    WRITE(21,*) "        type            wall;"
    WRITE(21,*) "        inGroups        1(wall);"
    N=0
    WRITE(21,232) "        nFaces          ",N,";"
    N=F_INTERNAL_NUM+CNL*CNL
    WRITE(21,232) "        startFace       ",N,";"
    WRITE(21,*) "    }"
       
        !UP FACE--PROCESSOR_MI1_TO_MI2  !CASE RECTANGLE 1-3-2
    PRO1=PNI*PNJ+PNL*PNI*(MJ2-1)+PNI*(MJ1-1)+MI-1
    MJ=1
    PRO2=PNI*(MJ-1)+MI-1
    WRITE(21,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(21,*) '    {'
    WRITE(21,*) '        type            processor;'
    WRITE(21,*) '        inGroups        List<word> 1(processor);'
    N=CNI*CNL
    WRITE(21,232) "        nFaces          ",N,";"
    N=F_INTERNAL_NUM+CNL*CNL
    WRITE(21,232) "        startFace       ",N,";"
    WRITE(21,*) '        matchTolerance  0.0001;'
    WRITE(21,232) "        myProcNo        ",PRO1,";"
    WRITE(21,232) "        neighbProcNo    ",PRO2,";"
    WRITE(21,*) '    }'
    WRITE(41,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(41,*) '    {'
    WRITE(41,*) '        type            processor;'
    WRITE(41,*) '        value           nonuniform List<scalar> '
    N=CNI*CNL
    WRITE(41,'(I0)') N
    WRITE(41,*) '('
    DO KP=1,CNL
    DO IP=1,CNI
    WRITE(41,'(E21.14)') P_UP_RECT(IP,KP)    
    ENDDO
    ENDDO
    WRITE(41,*) ')'
    WRITE(41,*) ';'
    WRITE(41,*) '    }'
    WRITE(51,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(51,*) '    {'
    WRITE(51,*) '        type            processor;'
    WRITE(51,*) '        value           nonuniform List<vector> '
    N=CNI*CNL
    WRITE(51,'(I0)') N
    WRITE(51,*) '('
    DO KP=1,CNL
    DO IP=1,CNI
    UU1=U_UP_RECT(1,IP,KP)
    UU2=U_UP_RECT(2,IP,KP)
    UU3=U_UP_RECT(3,IP,KP)
    WRITE(51,230) '( ',UU1,' ',UU2,' ',UU3,' )'  
    ENDDO
    ENDDO
    WRITE(51,*) ')'
    WRITE(51,*) ';'
    WRITE(51,*) '    }'
    WRITE(61,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(61,*) '    {'
    WRITE(61,*) '        type            processor;'
    WRITE(61,*) "        value           uniform 0;"
    WRITE(61,*) '    }'
    
        !FRONT FACE--PROCESSOR_MI1_TO_MI2  !CASE RECTANGLE 1-3-2
    PRO2=PNI*PNJ+PNL*PNI*(MJ2-2)+PNI*(MJ1-1)+MI-1
    WRITE(21,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(21,*) '    {'
    WRITE(21,*) '        type            processor;'
    WRITE(21,*) '        inGroups        List<word> 1(processor);'
    N=CNI*CNL
    WRITE(21,232) "        nFaces          ",N,";"
    N=F_INTERNAL_NUM+CNL*CNL+CNI*CNL
    WRITE(21,232) "        startFace       ",N,";"
    WRITE(21,*) '        matchTolerance  0.0001;'
    WRITE(21,232) "        myProcNo        ",PRO1,";"
    WRITE(21,232) "        neighbProcNo    ",PRO2,";"
    WRITE(21,*) '    }'
    WRITE(41,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(41,*) '    {'
    WRITE(41,*) '        type            processor;'
    WRITE(41,*) '        value           nonuniform List<scalar> '
    N=CNI*CNL
    WRITE(41,'(I0)') N
    WRITE(41,*) '('
    DO JP=1,CNL
    DO IP=1,CNI
    WRITE(41,'(E21.14)') P_FRONT_RECT(IP,JP)    
    ENDDO
    ENDDO
    WRITE(41,*) ')'
    WRITE(41,*) ';'
    WRITE(41,*) '    }'
    WRITE(51,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(51,*) '    {'
    WRITE(51,*) '        type            processor;'
    WRITE(51,*) '        value           nonuniform List<vector> '
    N=CNI*CNL
    WRITE(51,'(I0)') N
    WRITE(51,*) '('
    DO JP=1,CNL
    DO IP=1,CNI
    UU1=U_FRONT_RECT(1,IP,JP)
    UU2=U_FRONT_RECT(2,IP,JP)
    UU3=U_FRONT_RECT(3,IP,JP)
    WRITE(51,230) '( ',UU1,' ',UU2,' ',UU3,' )'  
    ENDDO
    ENDDO
    WRITE(51,*) ')'
    WRITE(51,*) ';'
    WRITE(51,*) '    }'
    WRITE(61,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(61,*) '    {'
    WRITE(61,*) '        type            processor;'
    WRITE(61,*) "        value           uniform 0;"
    WRITE(61,*) '    }'
    
        !DOWN FACE--PROCESSOR_MI1_TO_MI2  !CASE RECTANGLE 1-3-2
    PRO2=PNI*PNJ+PNL*PNI*(MJ2-1)+PNI*(MJ1-2)+MI-1
    WRITE(21,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(21,*) '    {'
    WRITE(21,*) '        type            processor;'
    WRITE(21,*) '        inGroups        List<word> 1(processor);'
    N=CNI*CNL
    WRITE(21,232) "        nFaces          ",N,";"
    N=F_INTERNAL_NUM+CNL*CNL+CNI*CNL+CNI*CNL
    WRITE(21,232) "        startFace       ",N,";"
    WRITE(21,*) '        matchTolerance  0.0001;'
    WRITE(21,232) "        myProcNo        ",PRO1,";"
    WRITE(21,232) "        neighbProcNo    ",PRO2,";"
    WRITE(21,*) '    }'
    WRITE(41,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(41,*) '    {'
    WRITE(41,*) '        type            processor;'
    WRITE(41,*) '        value           nonuniform List<scalar> '
    N=CNI*CNL
    WRITE(41,'(I0)') N
    WRITE(41,*) '('
    DO KP=1,CNL
    DO IP=1,CNI
    WRITE(41,'(E21.14)') P_DOWN_RECT(IP,KP)    
    ENDDO
    ENDDO
    WRITE(41,*) ')'
    WRITE(41,*) ';'
    WRITE(41,*) '    }'
    WRITE(51,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(51,*) '    {'
    WRITE(51,*) '        type            processor;'
    WRITE(51,*) '        value           nonuniform List<vector> '
    N=CNI*CNL
    WRITE(51,'(I0)') N
    WRITE(51,*) '('
    DO KP=1,CNL
    DO IP=1,CNI
    UU1=U_DOWN_RECT(1,IP,KP)
    UU2=U_DOWN_RECT(2,IP,KP)
    UU3=U_DOWN_RECT(3,IP,KP)
    WRITE(51,230) '( ',UU1,' ',UU2,' ',UU3,' )'  
    ENDDO
    ENDDO
    WRITE(51,*) ')'
    WRITE(51,*) ';'
    WRITE(51,*) '    }'
    WRITE(61,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(61,*) '    {'
    WRITE(61,*) '        type            processor;'
    WRITE(61,*) "        value           uniform 0;"
    WRITE(61,*) '    }'
    
        !LEFT FACE--PROCESSOR_MI1_TO_MI2  !CASE RECTANGLE 1-3-2
    PRO2=PNI*PNJ+PNL*PNI*(MJ2-1)+PNI*(MJ1-1)+MI
    WRITE(21,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(21,*) '    {'
    WRITE(21,*) '        type            processor;'
    WRITE(21,*) '        inGroups        List<word> 1(processor);'
    N=CNL*CNL
    WRITE(21,232) "        nFaces          ",N,";"
    N=F_INTERNAL_NUM+CNL*CNL+CNI*CNL+CNI*CNL+CNI*CNL
    WRITE(21,232) "        startFace       ",N,";"
    WRITE(21,*) '        matchTolerance  0.0001;'
    WRITE(21,232) "        myProcNo        ",PRO1,";"
    WRITE(21,232) "        neighbProcNo    ",PRO2,";"
    WRITE(21,*) '    }'
    WRITE(41,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(41,*) '    {'
    WRITE(41,*) '        type            processor;'
    WRITE(41,*) '        value           nonuniform List<scalar> '
    N=CNL*CNL
    WRITE(41,'(I0)') N
    WRITE(41,*) '('
    DO KP=1,CNL
    DO JP=1,CNL
    WRITE(41,'(E21.14)') P_LEFT_RECT(JP,KP)    
    ENDDO
    ENDDO
    WRITE(41,*) ')'
    WRITE(41,*) ';'
    WRITE(41,*) '    }'
    WRITE(51,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(51,*) '    {'
    WRITE(51,*) '        type            processor;'
    WRITE(51,*) '        value           nonuniform List<vector> '
    N=CNL*CNL
    WRITE(51,'(I0)') N
    WRITE(51,*) '('
    DO KP=1,CNL
    DO JP=1,CNL
    UU1=U_LEFT_RECT(1,JP,KP)
    UU2=U_LEFT_RECT(2,JP,KP)
    UU3=U_LEFT_RECT(3,JP,KP)
    WRITE(51,230) '( ',UU1,' ',UU2,' ',UU3,' )'  
    ENDDO
    ENDDO
    WRITE(51,*) ')'
    WRITE(51,*) ';'
    WRITE(51,*) '    }'
    WRITE(61,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(61,*) '    {'
    WRITE(61,*) '        type            processor;'
    WRITE(61,*) "        value           uniform 0;"
    WRITE(61,*) '    }'
    
        !BACK FACE--PROCESSOR_MI1_TO_MI2  !CASE RECTANGLE 1-3-2
    PRO2=PNI*PNJ+PNL*PNI*MJ2+PNI*(MJ1-1)+MI-1
    WRITE(21,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(21,*) '    {'
    WRITE(21,*) '        type            processor;'
    WRITE(21,*) '        inGroups        List<word> 1(processor);'
    N=CNI*CNL
    WRITE(21,232) "        nFaces          ",N,";"
    N=F_INTERNAL_NUM+CNL*CNL+CNI*CNL+CNI*CNL+CNI*CNL+CNL*CNL
    WRITE(21,232) "        startFace       ",N,";"
    WRITE(21,*) '        matchTolerance  0.0001;'
    WRITE(21,232) "        myProcNo        ",PRO1,";"
    WRITE(21,232) "        neighbProcNo    ",PRO2,";"
    WRITE(21,*) '    }'
    WRITE(41,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(41,*) '    {'
    WRITE(41,*) '        type            processor;'
    WRITE(41,*) '        value           nonuniform List<scalar> '
    N=CNI*CNL
    WRITE(41,'(I0)') N
    WRITE(41,*) '('
    DO JP=1,CNL
    DO IP=1,CNI
    WRITE(41,'(E21.14)') P_BACK_RECT(IP,JP)    
    ENDDO
    ENDDO
    WRITE(41,*) ')'
    WRITE(41,*) ';'
    WRITE(41,*) '    }'
    WRITE(51,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(51,*) '    {'
    WRITE(51,*) '        type            processor;'
    WRITE(51,*) '        value           nonuniform List<vector> '
    N=CNI*CNL
    WRITE(51,'(I0)') N
    WRITE(51,*) '('
    DO JP=1,CNL
    DO IP=1,CNI
    UU1=U_BACK_RECT(1,IP,JP)
    UU2=U_BACK_RECT(2,IP,JP)
    UU3=U_BACK_RECT(3,IP,JP)
    WRITE(51,230) '( ',UU1,' ',UU2,' ',UU3,' )'  
    ENDDO
    ENDDO
    WRITE(51,*) ')'
    WRITE(51,*) ';'
    WRITE(51,*) '    }'
    WRITE(61,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(61,*) '    {'
    WRITE(61,*) '        type            processor;'
    WRITE(61,*) "        value           uniform 0;"
    WRITE(61,*) '    }'    
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
    ENDIF 
    ELSE
    IF (MJ2.EQ.1) THEN
    !CASE RECTANGLE 1-2-1
    
        !RIGHT FACE--INLET1  !CASE RECTANGLE 1-2-1 
    IP=1
    I=IP+CNI*(MI-1)
    DO KP=1,CNL
    K=KP+CNL*(MJ2-1)
    DO JP=1,CNL
    J=JP+CNL*(MJ1-1)
    F11=(CNL+1)*(CNI+1)*(KP-1)+(CNI+1)*(JP-1)+IP-1
    F12=F11+1
    F13=F12+CNI+1
    F14=F13-1
    F21=F11+(CNL+1)*(CNI+1)
    F22=F21+1
    F23=F22+CNI+1
    F24=F23-1
    CALL CAL_FIJ_RECT(F11,F12,F13,F14,F21,F22,F23,F24)
    WRITE(2,235) '4(',F11,' ',F21,' ',F24,' ',F14,')'
    FAD1=F_INTERNAL_NUM_ALL+NL*(K-1)+J
    WRITE(3,231) FAD1
    OWN=CELL_ID_RECT(IP,JP,KP)
    WRITE(4,231) OWN
    ENDDO   
    ENDDO    
        
        !FRONT FACE--PROCESSOR_MI1_TO_MI2  !CASE RECTANGLE 1-2-1
    KP=1
    K=KP+CNL*(MJ2-1)
    DO JP=1,CNL
    J=JP+CNL*(MJ1-1)
    DO IP=1,CNI
    I=IP+CNI*(MI-1)
    F11=(CNL+1)*(CNI+1)*(KP-1)+(CNI+1)*(JP-1)+IP-1
    F12=F11+1
    F13=F12+CNI+1
    F14=F13-1
    F21=F11+(CNL+1)*(CNI+1)
    F22=F21+1
    F23=F22+CNI+1
    F24=F23-1
    CALL CAL_FIJ_RECT(F11,F12,F13,F14,F21,F22,F23,F24)
    WRITE(2,235) '4(',F11,' ',F14,' ',F13,' ',F12,')'
    FAD1=F_FRONT_RECT(IP,JP)
    WRITE(3,231) FAD1
    OWN=CELL_ID_RECT(IP,JP,KP)
    WRITE(4,231) OWN
    ENDDO   
    ENDDO    
        
        !DOWN FACE--PROCESSOR_MI1_TO_MI2  !CASE RECTANGLE 1-2-1
    JP=1
    J=JP+CNL*(MJ1-1)
    DO KP=1,CNL
    K=KP+CNL*(MJ2-1)
    DO IP=1,CNI
    I=IP+CNI*(MI-1)
    F11=(CNL+1)*(CNI+1)*(KP-1)+(CNI+1)*(JP-1)+IP-1
    F12=F11+1
    F13=F12+CNI+1
    F14=F13-1
    F21=F11+(CNL+1)*(CNI+1)
    F22=F21+1
    F23=F22+CNI+1
    F24=F23-1
    CALL CAL_FIJ_RECT(F11,F12,F13,F14,F21,F22,F23,F24)
    WRITE(2,235) '4(',F11,' ',F12,' ',F22,' ',F21,')'
    FAD1=F_DOWN_RECT(IP,KP)
    WRITE(3,231) FAD1
    OWN=CELL_ID_RECT(IP,JP,KP)
    WRITE(4,231) OWN
    ENDDO   
    ENDDO
    
        !LEFT FACE--PROCESSOR_MI1_TO_MI2  !CASE RECTANGLE 1-2-1 
    IP=CNI
    I=IP+CNI*(MI-1)
    DO KP=1,CNL
    K=KP+CNL*(MJ2-1)
    DO JP=1,CNL
    J=JP+CNL*(MJ1-1)
    F11=(CNL+1)*(CNI+1)*(KP-1)+(CNI+1)*(JP-1)+IP-1
    F12=F11+1
    F13=F12+CNI+1
    F14=F13-1
    F21=F11+(CNL+1)*(CNI+1)
    F22=F21+1
    F23=F22+CNI+1
    F24=F23-1
    CALL CAL_FIJ_RECT(F11,F12,F13,F14,F21,F22,F23,F24)
    WRITE(2,235) '4(',F12,' ',F13,' ',F23,' ',F22,')'
    FAD1=F_LEFT_RECT(JP,KP)
    WRITE(3,231) FAD1
    OWN=CELL_ID_RECT(IP,JP,KP)
    WRITE(4,231) OWN
    ENDDO   
    ENDDO
    
        !UP FACE--PROCESSOR_MI1_TO_MI2  !CASE RECTANGLE 1-2-1
    JP=CNL
    J=JP+CNL*(MJ1-1)
    DO KP=1,CNL
    K=KP+CNL*(MJ2-1)
    DO IP=1,CNI
    I=IP+CNI*(MI-1)
    F11=(CNL+1)*(CNI+1)*(KP-1)+(CNI+1)*(JP-1)+IP-1
    F12=F11+1
    F13=F12+CNI+1
    F14=F13-1
    F21=F11+(CNL+1)*(CNI+1)
    F22=F21+1
    F23=F22+CNI+1
    F24=F23-1
    CALL CAL_FIJ_RECT(F11,F12,F13,F14,F21,F22,F23,F24)
    WRITE(2,235) '4(',F14,' ',F24,' ',F23,' ',F13,')'
    FAD1=F_UP_RECT(IP,KP)
    WRITE(3,231) FAD1
    OWN=CELL_ID_RECT(IP,JP,KP)
    WRITE(4,231) OWN
    ENDDO   
    ENDDO
    
        !BACK FACE--PROCESSOR_MI1_TO_MI2  !CASE RECTANGLE 1-2-1
    KP=CNL
    K=KP+CNL*(MJ2-1)
    DO JP=1,CNL
    J=JP+CNL*(MJ1-1)
    DO IP=1,CNI
    I=IP+CNI*(MI-1)
    F11=(CNL+1)*(CNI+1)*(KP-1)+(CNI+1)*(JP-1)+IP-1
    F12=F11+1
    F13=F12+CNI+1
    F14=F13-1
    F21=F11+(CNL+1)*(CNI+1)
    F22=F21+1
    F23=F22+CNI+1
    F24=F23-1
    CALL CAL_FIJ_RECT(F11,F12,F13,F14,F21,F22,F23,F24)
    WRITE(2,235) '4(',F21,' ',F22,' ',F23,' ',F24,')'
    FAD1=F_BACK_RECT(IP,JP)
    WRITE(3,231) FAD1
    OWN=CELL_ID_RECT(IP,JP,KP)
    WRITE(4,231) OWN
    ENDDO   
    ENDDO    
        
    N=5
    WRITE(31,231) N+4
    WRITE(31,*) '('  
    WRITE(31,*) '0'
    WRITE(31,*) '1'
    WRITE(31,*) '2' 
    WRITE(31,*) '3' 
    DO I=1,N
    WRITE(31,*) '-1' 
    ENDDO
    WRITE(31,*) ')'
    
    WRITE(21,231) N+4
    WRITE(21,*) "("    
    WRITE(21,*) "    INLET1"
    WRITE(21,*) "    {"    
    WRITE(21,*) "        type            patch;"
    N=CNL*CNL
    WRITE(21,232) "        nFaces          ",N,";"
    N=F_INTERNAL_NUM
    WRITE(21,232) "        startFace       ",N,";"
    WRITE(21,*) "    }"
    WRITE(21,*) "    INLET2"
    WRITE(21,*) "    {"    
    WRITE(21,*) "        type            patch;"
    N=0
    WRITE(21,232) "        nFaces          ",N,";"
    N=F_INTERNAL_NUM+CNL*CNL
    WRITE(21,232) "        startFace       ",N,";"
    WRITE(21,*) "    }"
    WRITE(21,*) "    OUTLET"
    WRITE(21,*) "    {"
    WRITE(21,*) "        type            patch;"
    N=0
    WRITE(21,232) "        nFaces          ",N,";"
    N=F_INTERNAL_NUM+CNL*CNL
    WRITE(21,232) "        startFace       ",N,";"
    WRITE(21,*) "    }" 
    WRITE(21,*) "    HULL"
    WRITE(21,*) "    {"
    WRITE(21,*) "        type            wall;"
    WRITE(21,*) "        inGroups        1(wall);"
    N=0
    WRITE(21,232) "        nFaces          ",N,";"
    N=F_INTERNAL_NUM+CNL*CNL
    WRITE(21,232) "        startFace       ",N,";"
    WRITE(21,*) "    }"
       
        !FRONT FACE--PROCESSOR_MI1_TO_MI2  !CASE RECTANGLE 1-2-1
    PRO1=PNI*PNJ+PNL*PNI*(MJ2-1)+PNI*(MJ1-1)+MI-1
    MJ=1
    PRO2=PNI*(MJ-1)+MI-1
    WRITE(21,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(21,*) '    {'
    WRITE(21,*) '        type            processor;'
    WRITE(21,*) '        inGroups        List<word> 1(processor);'
    N=CNI*CNL
    WRITE(21,232) "        nFaces          ",N,";"
    N=F_INTERNAL_NUM+CNL*CNL
    WRITE(21,232) "        startFace       ",N,";"
    WRITE(21,*) '        matchTolerance  0.0001;'
    WRITE(21,232) "        myProcNo        ",PRO1,";"
    WRITE(21,232) "        neighbProcNo    ",PRO2,";"
    WRITE(21,*) '    }'
    WRITE(41,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(41,*) '    {'
    WRITE(41,*) '        type            processor;'
    WRITE(41,*) '        value           nonuniform List<scalar> '
    N=CNI*CNL
    WRITE(41,'(I0)') N
    WRITE(41,*) '('
    DO JP=1,CNL
    DO IP=1,CNI
    WRITE(41,'(E21.14)') P_FRONT_RECT(IP,JP)    
    ENDDO
    ENDDO
    WRITE(41,*) ')'
    WRITE(41,*) ';'
    WRITE(41,*) '    }'
    WRITE(51,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(51,*) '    {'
    WRITE(51,*) '        type            processor;'
    WRITE(51,*) '        value           nonuniform List<vector> '
    N=CNI*CNL
    WRITE(51,'(I0)') N
    WRITE(51,*) '('
    DO JP=1,CNL
    DO IP=1,CNI
    UU1=U_FRONT_RECT(1,IP,JP)
    UU2=U_FRONT_RECT(2,IP,JP)
    UU3=U_FRONT_RECT(3,IP,JP)
    WRITE(51,230) '( ',UU1,' ',UU2,' ',UU3,' )'  
    ENDDO
    ENDDO
    WRITE(51,*) ')'
    WRITE(51,*) ';'
    WRITE(51,*) '    }'
    WRITE(61,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(61,*) '    {'
    WRITE(61,*) '        type            processor;'
    WRITE(61,*) "        value           uniform 0;"
    WRITE(61,*) '    }'
    
        !DOWN FACE--PROCESSOR_MI1_TO_MI2  !CASE RECTANGLE 1-2-1
    PRO2=PNI*PNJ+PNL*PNI*(MJ2-1)+PNI*(MJ1-2)+MI-1
    WRITE(21,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(21,*) '    {'
    WRITE(21,*) '        type            processor;'
    WRITE(21,*) '        inGroups        List<word> 1(processor);'
    N=CNI*CNL
    WRITE(21,232) "        nFaces          ",N,";"
    N=F_INTERNAL_NUM+CNL*CNL+CNI*CNL
    WRITE(21,232) "        startFace       ",N,";"
    WRITE(21,*) '        matchTolerance  0.0001;'
    WRITE(21,232) "        myProcNo        ",PRO1,";"
    WRITE(21,232) "        neighbProcNo    ",PRO2,";"
    WRITE(21,*) '    }'
    WRITE(41,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(41,*) '    {'
    WRITE(41,*) '        type            processor;'
    WRITE(41,*) '        value           nonuniform List<scalar> '
    N=CNI*CNL
    WRITE(41,'(I0)') N
    WRITE(41,*) '('
    DO KP=1,CNL
    DO IP=1,CNI
    WRITE(41,'(E21.14)') P_DOWN_RECT(IP,KP)    
    ENDDO
    ENDDO
    WRITE(41,*) ')'
    WRITE(41,*) ';'
    WRITE(41,*) '    }'
    WRITE(51,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(51,*) '    {'
    WRITE(51,*) '        type            processor;'
    WRITE(51,*) '        value           nonuniform List<vector> '
    N=CNI*CNL
    WRITE(51,'(I0)') N
    WRITE(51,*) '('
    DO KP=1,CNL
    DO IP=1,CNI
    UU1=U_DOWN_RECT(1,IP,KP)
    UU2=U_DOWN_RECT(2,IP,KP)
    UU3=U_DOWN_RECT(3,IP,KP)
    WRITE(51,230) '( ',UU1,' ',UU2,' ',UU3,' )'  
    ENDDO
    ENDDO
    WRITE(51,*) ')'
    WRITE(51,*) ';'
    WRITE(51,*) '    }'
    WRITE(61,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(61,*) '    {'
    WRITE(61,*) '        type            processor;'
    WRITE(61,*) "        value           uniform 0;"
    WRITE(61,*) '    }'
    
        !LEFT FACE--PROCESSOR_MI1_TO_MI2  !CASE RECTANGLE 1-2-1 
    PRO2=PNI*PNJ+PNL*PNI*(MJ2-1)+PNI*(MJ1-1)+MI
    WRITE(21,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(21,*) '    {'
    WRITE(21,*) '        type            processor;'
    WRITE(21,*) '        inGroups        List<word> 1(processor);'
    N=CNL*CNL
    WRITE(21,232) "        nFaces          ",N,";"
    N=F_INTERNAL_NUM+CNL*CNL+CNI*CNL+CNI*CNL
    WRITE(21,232) "        startFace       ",N,";"
    WRITE(21,*) '        matchTolerance  0.0001;'
    WRITE(21,232) "        myProcNo        ",PRO1,";"
    WRITE(21,232) "        neighbProcNo    ",PRO2,";"
    WRITE(21,*) '    }'
    WRITE(41,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(41,*) '    {'
    WRITE(41,*) '        type            processor;'
    WRITE(41,*) '        value           nonuniform List<scalar> '
    N=CNL*CNL
    WRITE(41,'(I0)') N
    WRITE(41,*) '('
    DO KP=1,CNL
    DO JP=1,CNL
    WRITE(41,'(E21.14)') P_LEFT_RECT(JP,KP)    
    ENDDO
    ENDDO
    WRITE(41,*) ')'
    WRITE(41,*) ';'
    WRITE(41,*) '    }'
    WRITE(51,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(51,*) '    {'
    WRITE(51,*) '        type            processor;'
    WRITE(51,*) '        value           nonuniform List<vector> '
    N=CNL*CNL
    WRITE(51,'(I0)') N
    WRITE(51,*) '('
    DO KP=1,CNL
    DO JP=1,CNL
    UU1=U_LEFT_RECT(1,JP,KP)
    UU2=U_LEFT_RECT(2,JP,KP)
    UU3=U_LEFT_RECT(3,JP,KP)
    WRITE(51,230) '( ',UU1,' ',UU2,' ',UU3,' )'  
    ENDDO
    ENDDO
    WRITE(51,*) ')'
    WRITE(51,*) ';'
    WRITE(51,*) '    }'
    WRITE(61,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(61,*) '    {'
    WRITE(61,*) '        type            processor;'
    WRITE(61,*) "        value           uniform 0;"
    WRITE(61,*) '    }' 
    
        !UP FACE--PROCESSOR_MI1_TO_MI2  !CASE RECTANGLE 1-2-1
    PRO2=PNI*PNJ+PNL*PNI*(MJ2-1)+PNI*MJ1+MI-1
    WRITE(21,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(21,*) '    {'
    WRITE(21,*) '        type            processor;'
    WRITE(21,*) '        inGroups        List<word> 1(processor);'
    N=CNI*CNL
    WRITE(21,232) "        nFaces          ",N,";"
    N=F_INTERNAL_NUM+CNL*CNL+CNI*CNL+CNI*CNL+CNL*CNL
    WRITE(21,232) "        startFace       ",N,";"
    WRITE(21,*) '        matchTolerance  0.0001;'
    WRITE(21,232) "        myProcNo        ",PRO1,";"
    WRITE(21,232) "        neighbProcNo    ",PRO2,";"
    WRITE(21,*) '    }'
    WRITE(41,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(41,*) '    {'
    WRITE(41,*) '        type            processor;'
    WRITE(41,*) '        value           nonuniform List<scalar> '
    N=CNI*CNL
    WRITE(41,'(I0)') N
    WRITE(41,*) '('
    DO KP=1,CNL
    DO IP=1,CNI
    WRITE(41,'(E21.14)') P_UP_RECT(IP,KP)    
    ENDDO
    ENDDO
    WRITE(41,*) ')'
    WRITE(41,*) ';'
    WRITE(41,*) '    }'
    WRITE(51,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(51,*) '    {'
    WRITE(51,*) '        type            processor;'
    WRITE(51,*) '        value           nonuniform List<vector> '
    N=CNI*CNL
    WRITE(51,'(I0)') N
    WRITE(51,*) '('
    DO KP=1,CNL
    DO IP=1,CNI
    UU1=U_UP_RECT(1,IP,KP)
    UU2=U_UP_RECT(2,IP,KP)
    UU3=U_UP_RECT(3,IP,KP)
    WRITE(51,230) '( ',UU1,' ',UU2,' ',UU3,' )'  
    ENDDO
    ENDDO
    WRITE(51,*) ')'
    WRITE(51,*) ';'
    WRITE(51,*) '    }'
    WRITE(61,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(61,*) '    {'
    WRITE(61,*) '        type            processor;'
    WRITE(61,*) "        value           uniform 0;"
    WRITE(61,*) '    }'
    
        !BACK FACE--PROCESSOR_MI1_TO_MI2  !CASE RECTANGLE 1-2-1
    PRO2=PNI*PNJ+PNL*PNI*MJ2+PNI*(MJ1-1)+MI-1
    WRITE(21,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(21,*) '    {'
    WRITE(21,*) '        type            processor;'
    WRITE(21,*) '        inGroups        List<word> 1(processor);'
    N=CNI*CNL
    WRITE(21,232) "        nFaces          ",N,";"
    N=F_INTERNAL_NUM+CNL*CNL+CNI*CNL+CNI*CNL+CNL*CNL+CNI*CNL
    WRITE(21,232) "        startFace       ",N,";"
    WRITE(21,*) '        matchTolerance  0.0001;'
    WRITE(21,232) "        myProcNo        ",PRO1,";"
    WRITE(21,232) "        neighbProcNo    ",PRO2,";"
    WRITE(21,*) '    }'
    WRITE(41,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(41,*) '    {'
    WRITE(41,*) '        type            processor;'
    WRITE(41,*) '        value           nonuniform List<scalar> '
    N=CNI*CNL
    WRITE(41,'(I0)') N
    WRITE(41,*) '('
    DO JP=1,CNL
    DO IP=1,CNI
    WRITE(41,'(E21.14)') P_BACK_RECT(IP,JP)    
    ENDDO
    ENDDO
    WRITE(41,*) ')'
    WRITE(41,*) ';'
    WRITE(41,*) '    }'
    WRITE(51,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(51,*) '    {'
    WRITE(51,*) '        type            processor;'
    WRITE(51,*) '        value           nonuniform List<vector> '
    N=CNI*CNL
    WRITE(51,'(I0)') N
    WRITE(51,*) '('
    DO JP=1,CNL
    DO IP=1,CNI
    UU1=U_BACK_RECT(1,IP,JP)
    UU2=U_BACK_RECT(2,IP,JP)
    UU3=U_BACK_RECT(3,IP,JP)
    WRITE(51,230) '( ',UU1,' ',UU2,' ',UU3,' )'  
    ENDDO
    ENDDO
    WRITE(51,*) ')'
    WRITE(51,*) ';'
    WRITE(51,*) '    }'
    WRITE(61,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(61,*) '    {'
    WRITE(61,*) '        type            processor;'
    WRITE(61,*) "        value           uniform 0;"
    WRITE(61,*) '    }'    
        
        
        
        
        
        
        
        
        
        
        
        
        
    ELSEIF (MJ2.EQ.PNL) THEN
    !CASE RECTANGLE 1-2-3
    
        !RIGHT FACE--INLET1  !CASE RECTANGLE 1-2-3 
    IP=1
    I=IP+CNI*(MI-1)
    DO KP=1,CNL
    K=KP+CNL*(MJ2-1)
    DO JP=1,CNL
    J=JP+CNL*(MJ1-1)
    F11=(CNL+1)*(CNI+1)*(KP-1)+(CNI+1)*(JP-1)+IP-1
    F12=F11+1
    F13=F12+CNI+1
    F14=F13-1
    F21=F11+(CNL+1)*(CNI+1)
    F22=F21+1
    F23=F22+CNI+1
    F24=F23-1
    CALL CAL_FIJ_RECT(F11,F12,F13,F14,F21,F22,F23,F24)
    WRITE(2,235) '4(',F11,' ',F21,' ',F24,' ',F14,')'
    FAD1=F_INTERNAL_NUM_ALL+NL*(K-1)+J
    WRITE(3,231) FAD1
    OWN=CELL_ID_RECT(IP,JP,KP)
    WRITE(4,231) OWN
    ENDDO   
    ENDDO    
        
        !BACK FACE--PROCESSOR_MI1_TO_MI2  !CASE RECTANGLE 1-2-3
    KP=CNL
    K=KP+CNL*(MJ2-1)
    DO JP=CNL,1,-1
    J=JP+CNL*(MJ1-1)
    DO IP=1,CNI
    I=IP+CNI*(MI-1)
    F11=(CNL+1)*(CNI+1)*(KP-1)+(CNI+1)*(JP-1)+IP-1
    F12=F11+1
    F13=F12+CNI+1
    F14=F13-1
    F21=F11+(CNL+1)*(CNI+1)
    F22=F21+1
    F23=F22+CNI+1
    F24=F23-1
    CALL CAL_FIJ_RECT(F11,F12,F13,F14,F21,F22,F23,F24)
    !WRITE(2,235) '4(',F21,' ',F22,' ',F23,' ',F24,')'
    WRITE(2,235) '4(',F24,' ',F21,' ',F22,' ',F23,')'
    FAD1=F_BACK_RECT(IP,JP)
    WRITE(3,231) FAD1
    OWN=CELL_ID_RECT(IP,JP,KP)
    WRITE(4,231) OWN
    ENDDO   
    ENDDO    
        
        !FRONT FACE--PROCESSOR_MI1_TO_MI2  !CASE RECTANGLE 1-2-3
    KP=1
    K=KP+CNL*(MJ2-1)
    DO JP=1,CNL
    J=JP+CNL*(MJ1-1)
    DO IP=1,CNI
    I=IP+CNI*(MI-1)
    F11=(CNL+1)*(CNI+1)*(KP-1)+(CNI+1)*(JP-1)+IP-1
    F12=F11+1
    F13=F12+CNI+1
    F14=F13-1
    F21=F11+(CNL+1)*(CNI+1)
    F22=F21+1
    F23=F22+CNI+1
    F24=F23-1
    CALL CAL_FIJ_RECT(F11,F12,F13,F14,F21,F22,F23,F24)
    WRITE(2,235) '4(',F11,' ',F14,' ',F13,' ',F12,')'
    FAD1=F_FRONT_RECT(IP,JP)
    WRITE(3,231) FAD1
    OWN=CELL_ID_RECT(IP,JP,KP)
    WRITE(4,231) OWN
    ENDDO   
    ENDDO
    
        !DOWN FACE--PROCESSOR_MI1_TO_MI2  !CASE RECTANGLE 1-2-3
    JP=1
    J=JP+CNL*(MJ1-1)
    DO KP=1,CNL
    K=KP+CNL*(MJ2-1)
    DO IP=1,CNI
    I=IP+CNI*(MI-1)
    F11=(CNL+1)*(CNI+1)*(KP-1)+(CNI+1)*(JP-1)+IP-1
    F12=F11+1
    F13=F12+CNI+1
    F14=F13-1
    F21=F11+(CNL+1)*(CNI+1)
    F22=F21+1
    F23=F22+CNI+1
    F24=F23-1
    CALL CAL_FIJ_RECT(F11,F12,F13,F14,F21,F22,F23,F24)
    WRITE(2,235) '4(',F11,' ',F12,' ',F22,' ',F21,')'
    FAD1=F_DOWN_RECT(IP,KP)
    WRITE(3,231) FAD1
    OWN=CELL_ID_RECT(IP,JP,KP)
    WRITE(4,231) OWN
    ENDDO   
    ENDDO
    
        !LEFT FACE--PROCESSOR_MI1_TO_MI2  !CASE RECTANGLE 1-2-3 
    IP=CNI
    I=IP+CNI*(MI-1)
    DO KP=1,CNL
    K=KP+CNL*(MJ2-1)
    DO JP=1,CNL
    J=JP+CNL*(MJ1-1)
    F11=(CNL+1)*(CNI+1)*(KP-1)+(CNI+1)*(JP-1)+IP-1
    F12=F11+1
    F13=F12+CNI+1
    F14=F13-1
    F21=F11+(CNL+1)*(CNI+1)
    F22=F21+1
    F23=F22+CNI+1
    F24=F23-1
    CALL CAL_FIJ_RECT(F11,F12,F13,F14,F21,F22,F23,F24)
    WRITE(2,235) '4(',F12,' ',F13,' ',F23,' ',F22,')'
    FAD1=F_LEFT_RECT(JP,KP)
    WRITE(3,231) FAD1
    OWN=CELL_ID_RECT(IP,JP,KP)
    WRITE(4,231) OWN
    ENDDO   
    ENDDO
    
        !UP FACE--PROCESSOR_MI1_TO_MI2  !CASE RECTANGLE 1-2-3
    JP=CNL
    J=JP+CNL*(MJ1-1)
    DO KP=1,CNL
    K=KP+CNL*(MJ2-1)
    DO IP=1,CNI
    I=IP+CNI*(MI-1)
    F11=(CNL+1)*(CNI+1)*(KP-1)+(CNI+1)*(JP-1)+IP-1
    F12=F11+1
    F13=F12+CNI+1
    F14=F13-1
    F21=F11+(CNL+1)*(CNI+1)
    F22=F21+1
    F23=F22+CNI+1
    F24=F23-1
    CALL CAL_FIJ_RECT(F11,F12,F13,F14,F21,F22,F23,F24)
    WRITE(2,235) '4(',F14,' ',F24,' ',F23,' ',F13,')'
    FAD1=F_UP_RECT(IP,KP)
    WRITE(3,231) FAD1
    OWN=CELL_ID_RECT(IP,JP,KP)
    WRITE(4,231) OWN
    ENDDO   
    ENDDO    
        
    N=5
    WRITE(31,231) N+4
    WRITE(31,*) '('  
    WRITE(31,*) '0'
    WRITE(31,*) '1'
    WRITE(31,*) '2' 
    WRITE(31,*) '3' 
    DO I=1,N
    WRITE(31,*) '-1' 
    ENDDO
    WRITE(31,*) ')'
    
    WRITE(21,231) N+4
    WRITE(21,*) "("    
    WRITE(21,*) "    INLET1"
    WRITE(21,*) "    {"    
    WRITE(21,*) "        type            patch;"
    N=CNL*CNL
    WRITE(21,232) "        nFaces          ",N,";"
    N=F_INTERNAL_NUM
    WRITE(21,232) "        startFace       ",N,";"
    WRITE(21,*) "    }"
    WRITE(21,*) "    INLET2"
    WRITE(21,*) "    {"    
    WRITE(21,*) "        type            patch;"
    N=0
    WRITE(21,232) "        nFaces          ",N,";"
    N=F_INTERNAL_NUM+CNL*CNL
    WRITE(21,232) "        startFace       ",N,";"
    WRITE(21,*) "    }"
    WRITE(21,*) "    OUTLET"
    WRITE(21,*) "    {"
    WRITE(21,*) "        type            patch;"
    N=0
    WRITE(21,232) "        nFaces          ",N,";"
    N=F_INTERNAL_NUM+CNL*CNL
    WRITE(21,232) "        startFace       ",N,";"
    WRITE(21,*) "    }" 
    WRITE(21,*) "    HULL"
    WRITE(21,*) "    {"
    WRITE(21,*) "        type            wall;"
    WRITE(21,*) "        inGroups        1(wall);"
    N=0
    WRITE(21,232) "        nFaces          ",N,";"
    N=F_INTERNAL_NUM+CNL*CNL
    WRITE(21,232) "        startFace       ",N,";"
    WRITE(21,*) "    }"
       
        !BACK FACE--PROCESSOR_MI1_TO_MI2  !CASE RECTANGLE 1-2-3
    PRO1=PNI*PNJ+PNL*PNI*(MJ2-1)+PNI*(MJ1-1)+MI-1
    MJ=1
    PRO2=PNI*(MJ-1)+MI-1
    WRITE(21,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(21,*) '    {'
    WRITE(21,*) '        type            processor;'
    WRITE(21,*) '        inGroups        List<word> 1(processor);'
    N=CNI*CNL
    WRITE(21,232) "        nFaces          ",N,";"
    N=F_INTERNAL_NUM+CNL*CNL
    WRITE(21,232) "        startFace       ",N,";"
    WRITE(21,*) '        matchTolerance  0.0001;'
    WRITE(21,232) "        myProcNo        ",PRO1,";"
    WRITE(21,232) "        neighbProcNo    ",PRO2,";"
    WRITE(21,*) '    }'
    WRITE(41,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(41,*) '    {'
    WRITE(41,*) '        type            processor;'
    WRITE(41,*) '        value           nonuniform List<scalar> '
    N=CNI*CNL
    WRITE(41,'(I0)') N
    WRITE(41,*) '('
    DO JP=1,CNL
    DO IP=1,CNI
    WRITE(41,'(E21.14)') P_BACK_RECT(IP,JP)    
    ENDDO
    ENDDO
    WRITE(41,*) ')'
    WRITE(41,*) ';'
    WRITE(41,*) '    }'
    WRITE(51,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(51,*) '    {'
    WRITE(51,*) '        type            processor;'
    WRITE(51,*) '        value           nonuniform List<vector> '
    N=CNI*CNL
    WRITE(51,'(I0)') N
    WRITE(51,*) '('
    DO JP=1,CNL
    DO IP=1,CNI
    UU1=U_BACK_RECT(1,IP,JP)
    UU2=U_BACK_RECT(2,IP,JP)
    UU3=U_BACK_RECT(3,IP,JP)
    WRITE(51,230) '( ',UU1,' ',UU2,' ',UU3,' )'  
    ENDDO
    ENDDO
    WRITE(51,*) ')'
    WRITE(51,*) ';'
    WRITE(51,*) '    }'
    WRITE(61,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(61,*) '    {'
    WRITE(61,*) '        type            processor;'
    WRITE(61,*) "        value           uniform 0;"
    WRITE(61,*) '    }'
    
        !FRONT FACE--PROCESSOR_MI1_TO_MI2  !CASE RECTANGLE 1-2-3
    PRO2=PNI*PNJ+PNL*PNI*(MJ2-2)+PNI*(MJ1-1)+MI-1
    WRITE(21,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(21,*) '    {'
    WRITE(21,*) '        type            processor;'
    WRITE(21,*) '        inGroups        List<word> 1(processor);'
    N=CNI*CNL
    WRITE(21,232) "        nFaces          ",N,";"
    N=F_INTERNAL_NUM+CNL*CNL+CNI*CNL
    WRITE(21,232) "        startFace       ",N,";"
    WRITE(21,*) '        matchTolerance  0.0001;'
    WRITE(21,232) "        myProcNo        ",PRO1,";"
    WRITE(21,232) "        neighbProcNo    ",PRO2,";"
    WRITE(21,*) '    }'
    WRITE(41,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(41,*) '    {'
    WRITE(41,*) '        type            processor;'
    WRITE(41,*) '        value           nonuniform List<scalar> '
    N=CNI*CNL
    WRITE(41,'(I0)') N
    WRITE(41,*) '('
    DO JP=1,CNL
    DO IP=1,CNI
    WRITE(41,'(E21.14)') P_FRONT_RECT(IP,JP)    
    ENDDO
    ENDDO
    WRITE(41,*) ')'
    WRITE(41,*) ';'
    WRITE(41,*) '    }'
    WRITE(51,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(51,*) '    {'
    WRITE(51,*) '        type            processor;'
    WRITE(51,*) '        value           nonuniform List<vector> '
    N=CNI*CNL
    WRITE(51,'(I0)') N
    WRITE(51,*) '('
    DO JP=1,CNL
    DO IP=1,CNI
    UU1=U_FRONT_RECT(1,IP,JP)
    UU2=U_FRONT_RECT(2,IP,JP)
    UU3=U_FRONT_RECT(3,IP,JP)
    WRITE(51,230) '( ',UU1,' ',UU2,' ',UU3,' )'  
    ENDDO
    ENDDO
    WRITE(51,*) ')'
    WRITE(51,*) ';'
    WRITE(51,*) '    }'
    WRITE(61,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(61,*) '    {'
    WRITE(61,*) '        type            processor;'
    WRITE(61,*) "        value           uniform 0;"
    WRITE(61,*) '    }'
    
        !DOWN FACE--PROCESSOR_MI1_TO_MI2  !CASE RECTANGLE 1-2-3
    PRO2=PNI*PNJ+PNL*PNI*(MJ2-1)+PNI*(MJ1-2)+MI-1
    WRITE(21,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(21,*) '    {'
    WRITE(21,*) '        type            processor;'
    WRITE(21,*) '        inGroups        List<word> 1(processor);'
    N=CNI*CNL
    WRITE(21,232) "        nFaces          ",N,";"
    N=F_INTERNAL_NUM+CNL*CNL+CNI*CNL+CNI*CNL
    WRITE(21,232) "        startFace       ",N,";"
    WRITE(21,*) '        matchTolerance  0.0001;'
    WRITE(21,232) "        myProcNo        ",PRO1,";"
    WRITE(21,232) "        neighbProcNo    ",PRO2,";"
    WRITE(21,*) '    }'
    WRITE(41,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(41,*) '    {'
    WRITE(41,*) '        type            processor;'
    WRITE(41,*) '        value           nonuniform List<scalar> '
    N=CNI*CNL
    WRITE(41,'(I0)') N
    WRITE(41,*) '('
    DO KP=1,CNL
    DO IP=1,CNI
    WRITE(41,'(E21.14)') P_DOWN_RECT(IP,KP)    
    ENDDO
    ENDDO
    WRITE(41,*) ')'
    WRITE(41,*) ';'
    WRITE(41,*) '    }'
    WRITE(51,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(51,*) '    {'
    WRITE(51,*) '        type            processor;'
    WRITE(51,*) '        value           nonuniform List<vector> '
    N=CNI*CNL
    WRITE(51,'(I0)') N
    WRITE(51,*) '('
    DO KP=1,CNL
    DO IP=1,CNI
    UU1=U_DOWN_RECT(1,IP,KP)
    UU2=U_DOWN_RECT(2,IP,KP)
    UU3=U_DOWN_RECT(3,IP,KP)
    WRITE(51,230) '( ',UU1,' ',UU2,' ',UU3,' )'  
    ENDDO
    ENDDO
    WRITE(51,*) ')'
    WRITE(51,*) ';'
    WRITE(51,*) '    }'
    WRITE(61,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(61,*) '    {'
    WRITE(61,*) '        type            processor;'
    WRITE(61,*) "        value           uniform 0;"
    WRITE(61,*) '    }'
    
        !LEFT FACE--PROCESSOR_MI1_TO_MI2  !CASE RECTANGLE 1-2-3
    PRO2=PNI*PNJ+PNL*PNI*(MJ2-1)+PNI*(MJ1-1)+MI
    WRITE(21,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(21,*) '    {'
    WRITE(21,*) '        type            processor;'
    WRITE(21,*) '        inGroups        List<word> 1(processor);'
    N=CNL*CNL
    WRITE(21,232) "        nFaces          ",N,";"
    N=F_INTERNAL_NUM+CNL*CNL+CNI*CNL+CNI*CNL+CNI*CNL
    WRITE(21,232) "        startFace       ",N,";"
    WRITE(21,*) '        matchTolerance  0.0001;'
    WRITE(21,232) "        myProcNo        ",PRO1,";"
    WRITE(21,232) "        neighbProcNo    ",PRO2,";"
    WRITE(21,*) '    }'
    WRITE(41,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(41,*) '    {'
    WRITE(41,*) '        type            processor;'
    WRITE(41,*) '        value           nonuniform List<scalar> '
    N=CNL*CNL
    WRITE(41,'(I0)') N
    WRITE(41,*) '('
    DO KP=1,CNL
    DO JP=1,CNL
    WRITE(41,'(E21.14)') P_LEFT_RECT(JP,KP)    
    ENDDO
    ENDDO
    WRITE(41,*) ')'
    WRITE(41,*) ';'
    WRITE(41,*) '    }'
    WRITE(51,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(51,*) '    {'
    WRITE(51,*) '        type            processor;'
    WRITE(51,*) '        value           nonuniform List<vector> '
    N=CNL*CNL
    WRITE(51,'(I0)') N
    WRITE(51,*) '('
    DO KP=1,CNL
    DO JP=1,CNL
    UU1=U_LEFT_RECT(1,JP,KP)
    UU2=U_LEFT_RECT(2,JP,KP)
    UU3=U_LEFT_RECT(3,JP,KP)
    WRITE(51,230) '( ',UU1,' ',UU2,' ',UU3,' )'  
    ENDDO
    ENDDO
    WRITE(51,*) ')'
    WRITE(51,*) ';'
    WRITE(51,*) '    }'
    WRITE(61,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(61,*) '    {'
    WRITE(61,*) '        type            processor;'
    WRITE(61,*) "        value           uniform 0;"
    WRITE(61,*) '    }'
    
        !UP FACE--PROCESSOR_MI1_TO_MI2  !CASE RECTANGLE 1-2-3
    PRO2=PNI*PNJ+PNL*PNI*(MJ2-1)+PNI*MJ1+MI-1
    WRITE(21,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(21,*) '    {'
    WRITE(21,*) '        type            processor;'
    WRITE(21,*) '        inGroups        List<word> 1(processor);'
    N=CNI*CNL
    WRITE(21,232) "        nFaces          ",N,";"
    N=F_INTERNAL_NUM+CNL*CNL+CNI*CNL+CNI*CNL+CNI*CNL+CNL*CNL
    WRITE(21,232) "        startFace       ",N,";"
    WRITE(21,*) '        matchTolerance  0.0001;'
    WRITE(21,232) "        myProcNo        ",PRO1,";"
    WRITE(21,232) "        neighbProcNo    ",PRO2,";"
    WRITE(21,*) '    }'
    WRITE(41,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(41,*) '    {'
    WRITE(41,*) '        type            processor;'
    WRITE(41,*) '        value           nonuniform List<scalar> '
    N=CNI*CNL
    WRITE(41,'(I0)') N
    WRITE(41,*) '('
    DO KP=1,CNL
    DO IP=1,CNI
    WRITE(41,'(E21.14)') P_UP_RECT(IP,KP)    
    ENDDO
    ENDDO
    WRITE(41,*) ')'
    WRITE(41,*) ';'
    WRITE(41,*) '    }'
    WRITE(51,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(51,*) '    {'
    WRITE(51,*) '        type            processor;'
    WRITE(51,*) '        value           nonuniform List<vector> '
    N=CNI*CNL
    WRITE(51,'(I0)') N
    WRITE(51,*) '('
    DO KP=1,CNL
    DO IP=1,CNI
    UU1=U_UP_RECT(1,IP,KP)
    UU2=U_UP_RECT(2,IP,KP)
    UU3=U_UP_RECT(3,IP,KP)
    WRITE(51,230) '( ',UU1,' ',UU2,' ',UU3,' )'  
    ENDDO
    ENDDO
    WRITE(51,*) ')'
    WRITE(51,*) ';'
    WRITE(51,*) '    }'
    WRITE(61,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(61,*) '    {'
    WRITE(61,*) '        type            processor;'
    WRITE(61,*) "        value           uniform 0;"
    WRITE(61,*) '    }'    
  
        
    ELSE
    !CASE RECTANGLE 1-2-2
    
        !RIGHT FACE--INLET1  !CASE RECTANGLE 1-2-2 
    IP=1
    I=IP+CNI*(MI-1)
    DO KP=1,CNL
    K=KP+CNL*(MJ2-1)
    DO JP=1,CNL
    J=JP+CNL*(MJ1-1)
    F11=(CNL+1)*(CNI+1)*(KP-1)+(CNI+1)*(JP-1)+IP-1
    F12=F11+1
    F13=F12+CNI+1
    F14=F13-1
    F21=F11+(CNL+1)*(CNI+1)
    F22=F21+1
    F23=F22+CNI+1
    F24=F23-1
    CALL CAL_FIJ_RECT(F11,F12,F13,F14,F21,F22,F23,F24)
    WRITE(2,235) '4(',F11,' ',F21,' ',F24,' ',F14,')'
    FAD1=F_INTERNAL_NUM_ALL+NL*(K-1)+J
    WRITE(3,231) FAD1
    OWN=CELL_ID_RECT(IP,JP,KP)
    WRITE(4,231) OWN
    ENDDO   
    ENDDO    
        
        !FRONT FACE--PROCESSOR_MI1_TO_MI2  !CASE RECTANGLE 1-2-2
    KP=1
    K=KP+CNL*(MJ2-1)
    DO JP=1,CNL
    J=JP+CNL*(MJ1-1)
    DO IP=1,CNI
    I=IP+CNI*(MI-1)
    F11=(CNL+1)*(CNI+1)*(KP-1)+(CNI+1)*(JP-1)+IP-1
    F12=F11+1
    F13=F12+CNI+1
    F14=F13-1
    F21=F11+(CNL+1)*(CNI+1)
    F22=F21+1
    F23=F22+CNI+1
    F24=F23-1
    CALL CAL_FIJ_RECT(F11,F12,F13,F14,F21,F22,F23,F24)
    WRITE(2,235) '4(',F11,' ',F14,' ',F13,' ',F12,')'
    FAD1=F_FRONT_RECT(IP,JP)
    WRITE(3,231) FAD1
    OWN=CELL_ID_RECT(IP,JP,KP)
    WRITE(4,231) OWN
    ENDDO   
    ENDDO
    
        !DOWN FACE--PROCESSOR_MI1_TO_MI2  !CASE RECTANGLE 1-2-2
    JP=1
    J=JP+CNL*(MJ1-1)
    DO KP=1,CNL
    K=KP+CNL*(MJ2-1)
    DO IP=1,CNI
    I=IP+CNI*(MI-1)
    F11=(CNL+1)*(CNI+1)*(KP-1)+(CNI+1)*(JP-1)+IP-1
    F12=F11+1
    F13=F12+CNI+1
    F14=F13-1
    F21=F11+(CNL+1)*(CNI+1)
    F22=F21+1
    F23=F22+CNI+1
    F24=F23-1
    CALL CAL_FIJ_RECT(F11,F12,F13,F14,F21,F22,F23,F24)
    WRITE(2,235) '4(',F11,' ',F12,' ',F22,' ',F21,')'
    FAD1=F_DOWN_RECT(IP,KP)
    WRITE(3,231) FAD1
    OWN=CELL_ID_RECT(IP,JP,KP)
    WRITE(4,231) OWN
    ENDDO   
    ENDDO
    
        !LEFT FACE--PROCESSOR_MI1_TO_MI2  !CASE RECTANGLE 1-2-2 
    IP=CNI
    I=IP+CNI*(MI-1)
    DO KP=1,CNL
    K=KP+CNL*(MJ2-1)
    DO JP=1,CNL
    J=JP+CNL*(MJ1-1)
    F11=(CNL+1)*(CNI+1)*(KP-1)+(CNI+1)*(JP-1)+IP-1
    F12=F11+1
    F13=F12+CNI+1
    F14=F13-1
    F21=F11+(CNL+1)*(CNI+1)
    F22=F21+1
    F23=F22+CNI+1
    F24=F23-1
    CALL CAL_FIJ_RECT(F11,F12,F13,F14,F21,F22,F23,F24)
    WRITE(2,235) '4(',F12,' ',F13,' ',F23,' ',F22,')'
    FAD1=F_LEFT_RECT(JP,KP)
    WRITE(3,231) FAD1
    OWN=CELL_ID_RECT(IP,JP,KP)
    WRITE(4,231) OWN
    ENDDO   
    ENDDO
    
        !UP FACE--PROCESSOR_MI1_TO_MI2  !CASE RECTANGLE 1-2-2
    JP=CNL
    J=JP+CNL*(MJ1-1)
    DO KP=1,CNL
    K=KP+CNL*(MJ2-1)
    DO IP=1,CNI
    I=IP+CNI*(MI-1)
    F11=(CNL+1)*(CNI+1)*(KP-1)+(CNI+1)*(JP-1)+IP-1
    F12=F11+1
    F13=F12+CNI+1
    F14=F13-1
    F21=F11+(CNL+1)*(CNI+1)
    F22=F21+1
    F23=F22+CNI+1
    F24=F23-1
    CALL CAL_FIJ_RECT(F11,F12,F13,F14,F21,F22,F23,F24)
    WRITE(2,235) '4(',F14,' ',F24,' ',F23,' ',F13,')'
    FAD1=F_UP_RECT(IP,KP)
    WRITE(3,231) FAD1
    OWN=CELL_ID_RECT(IP,JP,KP)
    WRITE(4,231) OWN
    ENDDO   
    ENDDO
    
        !BACK FACE--PROCESSOR_MI1_TO_MI2  !CASE RECTANGLE 1-2-2
    KP=CNL
    K=KP+CNL*(MJ2-1)
    DO JP=1,CNL
    J=JP+CNL*(MJ1-1)
    DO IP=1,CNI
    I=IP+CNI*(MI-1)
    F11=(CNL+1)*(CNI+1)*(KP-1)+(CNI+1)*(JP-1)+IP-1
    F12=F11+1
    F13=F12+CNI+1
    F14=F13-1
    F21=F11+(CNL+1)*(CNI+1)
    F22=F21+1
    F23=F22+CNI+1
    F24=F23-1
    CALL CAL_FIJ_RECT(F11,F12,F13,F14,F21,F22,F23,F24)
    WRITE(2,235) '4(',F21,' ',F22,' ',F23,' ',F24,')'
    FAD1=F_BACK_RECT(IP,JP)
    WRITE(3,231) FAD1
    OWN=CELL_ID_RECT(IP,JP,KP)
    WRITE(4,231) OWN
    ENDDO   
    ENDDO    
        
    N=5
    WRITE(31,231) N+4
    WRITE(31,*) '('  
    WRITE(31,*) '0'
    WRITE(31,*) '1'
    WRITE(31,*) '2' 
    WRITE(31,*) '3' 
    DO I=1,N
    WRITE(31,*) '-1' 
    ENDDO
    WRITE(31,*) ')'
    
    WRITE(21,231) N+4
    WRITE(21,*) "("    
    WRITE(21,*) "    INLET1"
    WRITE(21,*) "    {"    
    WRITE(21,*) "        type            patch;"
    N=CNL*CNL
    WRITE(21,232) "        nFaces          ",N,";"
    N=F_INTERNAL_NUM
    WRITE(21,232) "        startFace       ",N,";"
    WRITE(21,*) "    }"
    WRITE(21,*) "    INLET2"
    WRITE(21,*) "    {"    
    WRITE(21,*) "        type            patch;"
    N=0
    WRITE(21,232) "        nFaces          ",N,";"
    N=F_INTERNAL_NUM+CNL*CNL
    WRITE(21,232) "        startFace       ",N,";"
    WRITE(21,*) "    }"
    WRITE(21,*) "    OUTLET"
    WRITE(21,*) "    {"
    WRITE(21,*) "        type            patch;"
    N=0
    WRITE(21,232) "        nFaces          ",N,";"
    N=F_INTERNAL_NUM+CNL*CNL
    WRITE(21,232) "        startFace       ",N,";"
    WRITE(21,*) "    }" 
    WRITE(21,*) "    HULL"
    WRITE(21,*) "    {"
    WRITE(21,*) "        type            wall;"
    WRITE(21,*) "        inGroups        1(wall);"
    N=0
    WRITE(21,232) "        nFaces          ",N,";"
    N=F_INTERNAL_NUM+CNL*CNL
    WRITE(21,232) "        startFace       ",N,";"
    WRITE(21,*) "    }"
       
        !FRONT FACE--PROCESSOR_MI1_TO_MI2  !CASE RECTANGLE 1-2-2
    PRO1=PNI*PNJ+PNL*PNI*(MJ2-1)+PNI*(MJ1-1)+MI-1
    PRO2=PNI*PNJ+PNL*PNI*(MJ2-2)+PNI*(MJ1-1)+MI-1
    WRITE(21,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(21,*) '    {'
    WRITE(21,*) '        type            processor;'
    WRITE(21,*) '        inGroups        List<word> 1(processor);'
    N=CNI*CNL
    WRITE(21,232) "        nFaces          ",N,";"
    N=F_INTERNAL_NUM+CNL*CNL
    WRITE(21,232) "        startFace       ",N,";"
    WRITE(21,*) '        matchTolerance  0.0001;'
    WRITE(21,232) "        myProcNo        ",PRO1,";"
    WRITE(21,232) "        neighbProcNo    ",PRO2,";"
    WRITE(21,*) '    }'
    WRITE(41,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(41,*) '    {'
    WRITE(41,*) '        type            processor;'
    WRITE(41,*) '        value           nonuniform List<scalar> '
    N=CNI*CNL
    WRITE(41,'(I0)') N
    WRITE(41,*) '('
    DO JP=1,CNL
    DO IP=1,CNI
    WRITE(41,'(E21.14)') P_FRONT_RECT(IP,JP)    
    ENDDO
    ENDDO
    WRITE(41,*) ')'
    WRITE(41,*) ';'
    WRITE(41,*) '    }'
    WRITE(51,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(51,*) '    {'
    WRITE(51,*) '        type            processor;'
    WRITE(51,*) '        value           nonuniform List<vector> '
    N=CNI*CNL
    WRITE(51,'(I0)') N
    WRITE(51,*) '('
    DO JP=1,CNL
    DO IP=1,CNI
    UU1=U_FRONT_RECT(1,IP,JP)
    UU2=U_FRONT_RECT(2,IP,JP)
    UU3=U_FRONT_RECT(3,IP,JP)
    WRITE(51,230) '( ',UU1,' ',UU2,' ',UU3,' )'  
    ENDDO
    ENDDO
    WRITE(51,*) ')'
    WRITE(51,*) ';'
    WRITE(51,*) '    }'
    WRITE(61,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(61,*) '    {'
    WRITE(61,*) '        type            processor;'
    WRITE(61,*) "        value           uniform 0;"
    WRITE(61,*) '    }'
    
        !DOWN FACE--PROCESSOR_MI1_TO_MI2  !CASE RECTANGLE 1-2-2
    PRO2=PNI*PNJ+PNL*PNI*(MJ2-1)+PNI*(MJ1-2)+MI-1
    WRITE(21,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(21,*) '    {'
    WRITE(21,*) '        type            processor;'
    WRITE(21,*) '        inGroups        List<word> 1(processor);'
    N=CNI*CNL
    WRITE(21,232) "        nFaces          ",N,";"
    N=F_INTERNAL_NUM+CNL*CNL+CNI*CNL
    WRITE(21,232) "        startFace       ",N,";"
    WRITE(21,*) '        matchTolerance  0.0001;'
    WRITE(21,232) "        myProcNo        ",PRO1,";"
    WRITE(21,232) "        neighbProcNo    ",PRO2,";"
    WRITE(21,*) '    }'
    WRITE(41,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(41,*) '    {'
    WRITE(41,*) '        type            processor;'
    WRITE(41,*) '        value           nonuniform List<scalar> '
    N=CNI*CNL
    WRITE(41,'(I0)') N
    WRITE(41,*) '('
    DO KP=1,CNL
    DO IP=1,CNI
    WRITE(41,'(E21.14)') P_DOWN_RECT(IP,KP)    
    ENDDO
    ENDDO
    WRITE(41,*) ')'
    WRITE(41,*) ';'
    WRITE(41,*) '    }'
    WRITE(51,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(51,*) '    {'
    WRITE(51,*) '        type            processor;'
    WRITE(51,*) '        value           nonuniform List<vector> '
    N=CNI*CNL
    WRITE(51,'(I0)') N
    WRITE(51,*) '('
    DO KP=1,CNL
    DO IP=1,CNI
    UU1=U_DOWN_RECT(1,IP,KP)
    UU2=U_DOWN_RECT(2,IP,KP)
    UU3=U_DOWN_RECT(3,IP,KP)
    WRITE(51,230) '( ',UU1,' ',UU2,' ',UU3,' )'  
    ENDDO
    ENDDO
    WRITE(51,*) ')'
    WRITE(51,*) ';'
    WRITE(51,*) '    }'
    WRITE(61,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(61,*) '    {'
    WRITE(61,*) '        type            processor;'
    WRITE(61,*) "        value           uniform 0;"
    WRITE(61,*) '    }' 
    
        !LEFT FACE--PROCESSOR_MI1_TO_MI2  !CASE RECTANGLE 1-2-2
    PRO2=PNI*PNJ+PNL*PNI*(MJ2-1)+PNI*(MJ1-1)+MI
    WRITE(21,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(21,*) '    {'
    WRITE(21,*) '        type            processor;'
    WRITE(21,*) '        inGroups        List<word> 1(processor);'
    N=CNL*CNL
    WRITE(21,232) "        nFaces          ",N,";"
    N=F_INTERNAL_NUM+CNL*CNL+CNI*CNL+CNI*CNL
    WRITE(21,232) "        startFace       ",N,";"
    WRITE(21,*) '        matchTolerance  0.0001;'
    WRITE(21,232) "        myProcNo        ",PRO1,";"
    WRITE(21,232) "        neighbProcNo    ",PRO2,";"
    WRITE(21,*) '    }'
    WRITE(41,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(41,*) '    {'
    WRITE(41,*) '        type            processor;'
    WRITE(41,*) '        value           nonuniform List<scalar> '
    N=CNL*CNL
    WRITE(41,'(I0)') N
    WRITE(41,*) '('
    DO KP=1,CNL
    DO JP=1,CNL
    WRITE(41,'(E21.14)') P_LEFT_RECT(JP,KP)    
    ENDDO
    ENDDO
    WRITE(41,*) ')'
    WRITE(41,*) ';'
    WRITE(41,*) '    }'
    WRITE(51,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(51,*) '    {'
    WRITE(51,*) '        type            processor;'
    WRITE(51,*) '        value           nonuniform List<vector> '
    N=CNL*CNL
    WRITE(51,'(I0)') N
    WRITE(51,*) '('
    DO KP=1,CNL
    DO JP=1,CNL
    UU1=U_LEFT_RECT(1,JP,KP)
    UU2=U_LEFT_RECT(2,JP,KP)
    UU3=U_LEFT_RECT(3,JP,KP)
    WRITE(51,230) '( ',UU1,' ',UU2,' ',UU3,' )'  
    ENDDO
    ENDDO
    WRITE(51,*) ')'
    WRITE(51,*) ';'
    WRITE(51,*) '    }'
    WRITE(61,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(61,*) '    {'
    WRITE(61,*) '        type            processor;'
    WRITE(61,*) "        value           uniform 0;"
    WRITE(61,*) '    }' 
    
        !UP FACE--PROCESSOR_MI1_TO_MI2  !CASE RECTANGLE 1-2-2
    PRO2=PNI*PNJ+PNL*PNI*(MJ2-1)+PNI*MJ1+MI-1
    WRITE(21,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(21,*) '    {'
    WRITE(21,*) '        type            processor;'
    WRITE(21,*) '        inGroups        List<word> 1(processor);'
    N=CNI*CNL
    WRITE(21,232) "        nFaces          ",N,";"
    N=F_INTERNAL_NUM+CNL*CNL+CNI*CNL+CNI*CNL+CNL*CNL
    WRITE(21,232) "        startFace       ",N,";"
    WRITE(21,*) '        matchTolerance  0.0001;'
    WRITE(21,232) "        myProcNo        ",PRO1,";"
    WRITE(21,232) "        neighbProcNo    ",PRO2,";"
    WRITE(21,*) '    }'
    WRITE(41,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(41,*) '    {'
    WRITE(41,*) '        type            processor;'
    WRITE(41,*) '        value           nonuniform List<scalar> '
    N=CNI*CNL
    WRITE(41,'(I0)') N
    WRITE(41,*) '('
    DO KP=1,CNL
    DO IP=1,CNI
    WRITE(41,'(E21.14)') P_UP_RECT(IP,KP)    
    ENDDO
    ENDDO
    WRITE(41,*) ')'
    WRITE(41,*) ';'
    WRITE(41,*) '    }'
    WRITE(51,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(51,*) '    {'
    WRITE(51,*) '        type            processor;'
    WRITE(51,*) '        value           nonuniform List<vector> '
    N=CNI*CNL
    WRITE(51,'(I0)') N
    WRITE(51,*) '('
    DO KP=1,CNL
    DO IP=1,CNI
    UU1=U_UP_RECT(1,IP,KP)
    UU2=U_UP_RECT(2,IP,KP)
    UU3=U_UP_RECT(3,IP,KP)
    WRITE(51,230) '( ',UU1,' ',UU2,' ',UU3,' )'  
    ENDDO
    ENDDO
    WRITE(51,*) ')'
    WRITE(51,*) ';'
    WRITE(51,*) '    }'
    WRITE(61,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(61,*) '    {'
    WRITE(61,*) '        type            processor;'
    WRITE(61,*) "        value           uniform 0;"
    WRITE(61,*) '    }'
    
        !BACK FACE--PROCESSOR_MI1_TO_MI2  !CASE RECTANGLE 1-2-2
    PRO2=PNI*PNJ+PNL*PNI*MJ2+PNI*(MJ1-1)+MI-1
    WRITE(21,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(21,*) '    {'
    WRITE(21,*) '        type            processor;'
    WRITE(21,*) '        inGroups        List<word> 1(processor);'
    N=CNI*CNL
    WRITE(21,232) "        nFaces          ",N,";"
    N=F_INTERNAL_NUM+CNL*CNL+CNI*CNL+CNI*CNL+CNL*CNL+CNI*CNL
    WRITE(21,232) "        startFace       ",N,";"
    WRITE(21,*) '        matchTolerance  0.0001;'
    WRITE(21,232) "        myProcNo        ",PRO1,";"
    WRITE(21,232) "        neighbProcNo    ",PRO2,";"
    WRITE(21,*) '    }'
    WRITE(41,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(41,*) '    {'
    WRITE(41,*) '        type            processor;'
    WRITE(41,*) '        value           nonuniform List<scalar> '
    N=CNI*CNL
    WRITE(41,'(I0)') N
    WRITE(41,*) '('
    DO JP=1,CNL
    DO IP=1,CNI
    WRITE(41,'(E21.14)') P_BACK_RECT(IP,JP)    
    ENDDO
    ENDDO
    WRITE(41,*) ')'
    WRITE(41,*) ';'
    WRITE(41,*) '    }'
    WRITE(51,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(51,*) '    {'
    WRITE(51,*) '        type            processor;'
    WRITE(51,*) '        value           nonuniform List<vector> '
    N=CNI*CNL
    WRITE(51,'(I0)') N
    WRITE(51,*) '('
    DO JP=1,CNL
    DO IP=1,CNI
    UU1=U_BACK_RECT(1,IP,JP)
    UU2=U_BACK_RECT(2,IP,JP)
    UU3=U_BACK_RECT(3,IP,JP)
    WRITE(51,230) '( ',UU1,' ',UU2,' ',UU3,' )'  
    ENDDO
    ENDDO
    WRITE(51,*) ')'
    WRITE(51,*) ';'
    WRITE(51,*) '    }'
    WRITE(61,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(61,*) '    {'
    WRITE(61,*) '        type            processor;'
    WRITE(61,*) "        value           uniform 0;"
    WRITE(61,*) '    }'     
        
        
        
        
        
        
    ENDIF 
    ENDIF             
    ENDIF 
        
    ELSEIF (MI.EQ.NPR) THEN
!****************************************NI-END*******************************************
    IF (MJ.EQ.1) THEN
    !CASE CYLINDER 3-1
    
        !LEFT FACE--HULL  !CASE CYLINDER 3-1
    IP=CNI
    I=IP+CNI*(MI-1)
    DO KP=1,NK
    K=KP
    DO JP=1,CNJ
    J=JP+CNJ*(MJ-1)
    F11=(CNJ+1)*(CNI+1)*(KP-1)+(CNI+1)*(JP-1)+IP-1
    F12=F11+1
    F13=F12+CNI+1
    F14=F13-1
    F21=F11+(CNJ+1)*(CNI+1)
    IF (KP.EQ.NK) F21=(CNI+1)*(JP-1)+IP-1
    F22=F21+1
    F23=F22+CNI+1
    F24=F23-1
    CALL CAL_FIJ_CY(F11,F12,F13,F14,F21,F22,F23,F24)
    WRITE(2,235) '4(',F12,' ',F13,' ',F23,' ',F22,')'
    FAD1=F_INTERNAL_NUM_ALL+NL*NL+NJ*NK+NI*NK+NJ*(K-1)+J
    WRITE(3,231) FAD1
    OWN=CELL_ID_CY(IP,JP,KP)
    WRITE(4,231) OWN
    ENDDO   
    ENDDO
    
        !RIGHT FACE--PROCESSOR_MI1_TO_MI2  !CASE CYLINDER 3-1
    IP=1
    I=IP+CNI*(MI-1)
    DO KP=1,NK
    K=KP
    DO JP=1,CNJ
    J=JP+CNJ*(MJ-1)
    F11=(CNJ+1)*(CNI+1)*(KP-1)+(CNI+1)*(JP-1)+IP-1
    F12=F11+1
    F13=F12+CNI+1
    F14=F13-1
    F21=F11+(CNJ+1)*(CNI+1)
    IF (KP.EQ.NK) F21=(CNI+1)*(JP-1)+IP-1
    F22=F21+1
    F23=F22+CNI+1
    F24=F23-1
    CALL CAL_FIJ_CY(F11,F12,F13,F14,F21,F22,F23,F24)
    WRITE(2,235) '4(',F11,' ',F21,' ',F24,' ',F14,')'
    FAD1=F_RIGHT_CY(JP,KP)
    WRITE(3,231) FAD1
    OWN=CELL_ID_CY(IP,JP,KP)
    WRITE(4,231) OWN
    ENDDO   
    ENDDO
    
        !UP FACE--PROCESSOR_MI1_TO_MI2  !CASE CYLINDER 3-1
    JP=CNJ
    J=JP+CNJ*(MJ-1)
    DO KP=1,NK
    K=KP
    DO IP=1,CNI
    I=IP+CNI*(MI-1)
    F11=(CNJ+1)*(CNI+1)*(KP-1)+(CNI+1)*(JP-1)+IP-1
    F12=F11+1
    F13=F12+CNI+1
    F14=F13-1
    F21=F11+(CNJ+1)*(CNI+1)
    IF (KP.EQ.NK) F21=(CNI+1)*(JP-1)+IP-1
    F22=F21+1
    F23=F22+CNI+1
    F24=F23-1
    CALL CAL_FIJ_CY(F11,F12,F13,F14,F21,F22,F23,F24)
    WRITE(2,235) '4(',F14,' ',F24,' ',F23,' ',F13,')'
    FAD1=F_UP_CY(IP,KP)
    WRITE(3,231) FAD1
    OWN=CELL_ID_CY(IP,JP,KP)
    WRITE(4,231) OWN
    ENDDO   
    ENDDO
    
        !DOWN FACE--PROCESSOR_MI1_TO_MI2(RECTANGLE 1-1)  !CASE CYLINDER 3-1
    JP=1
    J=JP+CNJ*(MJ-1)
    DO KP=1,CNL
    K=KP
    DO IP=1,CNI
    I=IP+CNI*(MI-1)
    F11=(CNJ+1)*(CNI+1)*(KP-1)+(CNI+1)*(JP-1)+IP-1
    F12=F11+1
    F13=F12+CNI+1
    F14=F13-1
    F21=F11+(CNJ+1)*(CNI+1)
    IF (KP.EQ.NK) F21=(CNI+1)*(JP-1)+IP-1
    F22=F21+1
    F23=F22+CNI+1
    F24=F23-1
    CALL CAL_FIJ_CY(F11,F12,F13,F14,F21,F22,F23,F24)
    WRITE(2,235) '4(',F11,' ',F12,' ',F22,' ',F21,')'
    FAD1=F_DOWN_CY(IP,KP)
    WRITE(3,231) FAD1
    OWN=CELL_ID_CY(IP,JP,KP)
    WRITE(4,231) OWN
    ENDDO   
    ENDDO
    DO KP=NK-CNL+1,NK
    K=KP
    DO IP=1,CNI
    I=IP+CNI*(MI-1)
    F11=(CNJ+1)*(CNI+1)*(KP-1)+(CNI+1)*(JP-1)+IP-1
    F12=F11+1
    F13=F12+CNI+1
    F14=F13-1
    F21=F11+(CNJ+1)*(CNI+1)
    IF (KP.EQ.NK) F21=(CNI+1)*(JP-1)+IP-1
    F22=F21+1
    F23=F22+CNI+1
    F24=F23-1
    CALL CAL_FIJ_CY(F11,F12,F13,F14,F21,F22,F23,F24)
    WRITE(2,235) '4(',F11,' ',F12,' ',F22,' ',F21,')'
    FAD1=F_DOWN_CY(IP,KP)
    WRITE(3,231) FAD1
    OWN=CELL_ID_CY(IP,JP,KP)
    WRITE(4,231) OWN
    ENDDO   
    ENDDO
    
        !DOWN FACE--PROCESSOR_MI1_TO_MI2(RECTANGLE 2-1)  !CASE CYLINDER 3-1
    JP=1
    J=JP+CNJ*(MJ-1)
    DO NN=2,PNL-1
    DO KP=CNL*(NN-1)+1,CNL*NN
    K=KP
    DO IP=1,CNI
    I=IP+CNI*(MI-1)
    F11=(CNJ+1)*(CNI+1)*(KP-1)+(CNI+1)*(JP-1)+IP-1
    F12=F11+1
    F13=F12+CNI+1
    F14=F13-1
    F21=F11+(CNJ+1)*(CNI+1)
    IF (KP.EQ.NK) F21=(CNI+1)*(JP-1)+IP-1
    F22=F21+1
    F23=F22+CNI+1
    F24=F23-1
    CALL CAL_FIJ_CY(F11,F12,F13,F14,F21,F22,F23,F24)
    WRITE(2,235) '4(',F11,' ',F12,' ',F22,' ',F21,')'
    FAD1=F_DOWN_CY(IP,KP)
    WRITE(3,231) FAD1
    OWN=CELL_ID_CY(IP,JP,KP)
    WRITE(4,231) OWN
    ENDDO   
    ENDDO
    ENDDO
    
        !DOWN FACE--PROCESSOR_MI1_TO_MI2(RECTANGLE 3-1)  !CASE CYLINDER 3-1 
    JP=1
    J=JP+CNJ*(MJ-1)
    DO KP=NL-CNL+1,NL
    K=KP
    DO IP=1,CNI
    I=IP+CNI*(MI-1)
    F11=(CNJ+1)*(CNI+1)*(KP-1)+(CNI+1)*(JP-1)+IP-1
    F12=F11+1
    F13=F12+CNI+1
    F14=F13-1
    F21=F11+(CNJ+1)*(CNI+1)
    IF (KP.EQ.NK) F21=(CNI+1)*(JP-1)+IP-1
    F22=F21+1
    F23=F22+CNI+1
    F24=F23-1
    CALL CAL_FIJ_CY(F11,F12,F13,F14,F21,F22,F23,F24)
    WRITE(2,235) '4(',F11,' ',F12,' ',F22,' ',F21,')'
    FAD1=F_DOWN_CY(IP,KP)
    WRITE(3,231) FAD1
    OWN=CELL_ID_CY(IP,JP,KP)
    WRITE(4,231) OWN
    ENDDO   
    ENDDO
    DO KP=NL+1,NL+CNL
    K=KP
    DO IP=1,CNI
    I=IP+CNI*(MI-1)
    F11=(CNJ+1)*(CNI+1)*(KP-1)+(CNI+1)*(JP-1)+IP-1
    F12=F11+1
    F13=F12+CNI+1
    F14=F13-1
    F21=F11+(CNJ+1)*(CNI+1)
    IF (KP.EQ.NK) F21=(CNI+1)*(JP-1)+IP-1
    F22=F21+1
    F23=F22+CNI+1
    F24=F23-1
    CALL CAL_FIJ_CY(F11,F12,F13,F14,F21,F22,F23,F24)
    WRITE(2,235) '4(',F11,' ',F12,' ',F22,' ',F21,')'
    FAD1=F_DOWN_CY(IP,KP)
    WRITE(3,231) FAD1
    OWN=CELL_ID_CY(IP,JP,KP)
    WRITE(4,231) OWN
    ENDDO   
    ENDDO
        
        !DOWN FACE--PROCESSOR_MI1_TO_MI2(RECTANGLE 1-2 / 3-2)  !CASE CYLINDER 3-1  
    JP=1
    J=JP+CNJ*(MJ-1)
    DO NN=2,PNL-1
        
    DO KP=3*NL+CNL*(PNL-NN)+1,3*NL+CNL*(PNL-NN+1)
    K=KP
    DO IP=1,CNI
    I=IP+CNI*(MI-1)
    F11=(CNJ+1)*(CNI+1)*(KP-1)+(CNI+1)*(JP-1)+IP-1
    F12=F11+1
    F13=F12+CNI+1
    F14=F13-1
    F21=F11+(CNJ+1)*(CNI+1)
    IF (KP.EQ.NK) F21=(CNI+1)*(JP-1)+IP-1
    F22=F21+1
    F23=F22+CNI+1
    F24=F23-1
    CALL CAL_FIJ_CY(F11,F12,F13,F14,F21,F22,F23,F24)
    WRITE(2,235) '4(',F11,' ',F12,' ',F22,' ',F21,')'
    FAD1=F_DOWN_CY(IP,KP)
    WRITE(3,231) FAD1
    OWN=CELL_ID_CY(IP,JP,KP)
    WRITE(4,231) OWN
    ENDDO   
    ENDDO
    
    DO KP=NL+CNL*(NN-1)+1,NL+CNL*NN
    K=KP
    DO IP=1,CNI
    I=IP+CNI*(MI-1)
    F11=(CNJ+1)*(CNI+1)*(KP-1)+(CNI+1)*(JP-1)+IP-1
    F12=F11+1
    F13=F12+CNI+1
    F14=F13-1
    F21=F11+(CNJ+1)*(CNI+1)
    IF (KP.EQ.NK) F21=(CNI+1)*(JP-1)+IP-1
    F22=F21+1
    F23=F22+CNI+1
    F24=F23-1
    CALL CAL_FIJ_CY(F11,F12,F13,F14,F21,F22,F23,F24)
    WRITE(2,235) '4(',F11,' ',F12,' ',F22,' ',F21,')'
    FAD1=F_DOWN_CY(IP,KP)
    WRITE(3,231) FAD1
    OWN=CELL_ID_CY(IP,JP,KP)
    WRITE(4,231) OWN
    ENDDO   
    ENDDO
    
    ENDDO    
        
        !DOWN FACE--PROCESSOR_MI1_TO_MI2(RECTANGLE 1-3)  !CASE CYLINDER 3-1 
    JP=1
    J=JP+CNJ*(MJ-1)
    DO KP=3*NL-CNL+1,3*NL
    K=KP
    DO IP=1,CNI
    I=IP+CNI*(MI-1)
    F11=(CNJ+1)*(CNI+1)*(KP-1)+(CNI+1)*(JP-1)+IP-1
    F12=F11+1
    F13=F12+CNI+1
    F14=F13-1
    F21=F11+(CNJ+1)*(CNI+1)
    IF (KP.EQ.NK) F21=(CNI+1)*(JP-1)+IP-1
    F22=F21+1
    F23=F22+CNI+1
    F24=F23-1
    CALL CAL_FIJ_CY(F11,F12,F13,F14,F21,F22,F23,F24)
    WRITE(2,235) '4(',F11,' ',F12,' ',F22,' ',F21,')'
    FAD1=F_DOWN_CY(IP,KP)
    WRITE(3,231) FAD1
    OWN=CELL_ID_CY(IP,JP,KP)
    WRITE(4,231) OWN
    ENDDO   
    ENDDO
    DO KP=3*NL+1,3*NL+CNL
    K=KP
    DO IP=1,CNI
    I=IP+CNI*(MI-1)
    F11=(CNJ+1)*(CNI+1)*(KP-1)+(CNI+1)*(JP-1)+IP-1
    F12=F11+1
    F13=F12+CNI+1
    F14=F13-1
    F21=F11+(CNJ+1)*(CNI+1)
    IF (KP.EQ.NK) F21=(CNI+1)*(JP-1)+IP-1
    F22=F21+1
    F23=F22+CNI+1
    F24=F23-1
    CALL CAL_FIJ_CY(F11,F12,F13,F14,F21,F22,F23,F24)
    WRITE(2,235) '4(',F11,' ',F12,' ',F22,' ',F21,')'
    FAD1=F_DOWN_CY(IP,KP)
    WRITE(3,231) FAD1
    OWN=CELL_ID_CY(IP,JP,KP)
    WRITE(4,231) OWN
    ENDDO   
    ENDDO
    
        !DOWN FACE--PROCESSOR_MI1_TO_MI2(RECTANGLE 2-3)  !CASE CYLINDER 3-1  
    JP=1
    J=JP+CNJ*(MJ-1)
    DO NN=2,PNL-1
    DO KP=2*NL+CNL*(PNL-NN)+1,2*NL+CNL*(PNL-NN+1)
    K=KP
    DO IP=1,CNI
    I=IP+CNI*(MI-1)
    F11=(CNJ+1)*(CNI+1)*(KP-1)+(CNI+1)*(JP-1)+IP-1
    F12=F11+1
    F13=F12+CNI+1
    F14=F13-1
    F21=F11+(CNJ+1)*(CNI+1)
    IF (KP.EQ.NK) F21=(CNI+1)*(JP-1)+IP-1
    F22=F21+1
    F23=F22+CNI+1
    F24=F23-1
    CALL CAL_FIJ_CY(F11,F12,F13,F14,F21,F22,F23,F24)
    WRITE(2,235) '4(',F11,' ',F12,' ',F22,' ',F21,')'
    FAD1=F_DOWN_CY(IP,KP)
    WRITE(3,231) FAD1
    OWN=CELL_ID_CY(IP,JP,KP)
    WRITE(4,231) OWN
    ENDDO   
    ENDDO
    ENDDO
    
        !DOWN FACE--PROCESSOR_MI1_TO_MI2(RECTANGLE 3-3)  !CASE CYLINDER 3-1 
    JP=1
    J=JP+CNJ*(MJ-1)
    DO KP=2*NL-CNL+1,2*NL
    K=KP
    DO IP=1,CNI
    I=IP+CNI*(MI-1)
    F11=(CNJ+1)*(CNI+1)*(KP-1)+(CNI+1)*(JP-1)+IP-1
    F12=F11+1
    F13=F12+CNI+1
    F14=F13-1
    F21=F11+(CNJ+1)*(CNI+1)
    IF (KP.EQ.NK) F21=(CNI+1)*(JP-1)+IP-1
    F22=F21+1
    F23=F22+CNI+1
    F24=F23-1
    CALL CAL_FIJ_CY(F11,F12,F13,F14,F21,F22,F23,F24)
    WRITE(2,235) '4(',F11,' ',F12,' ',F22,' ',F21,')'
    FAD1=F_DOWN_CY(IP,KP)
    WRITE(3,231) FAD1
    OWN=CELL_ID_CY(IP,JP,KP)
    WRITE(4,231) OWN
    ENDDO   
    ENDDO
    DO KP=2*NL+1,2*NL+CNL
    K=KP
    DO IP=1,CNI
    I=IP+CNI*(MI-1)
    F11=(CNJ+1)*(CNI+1)*(KP-1)+(CNI+1)*(JP-1)+IP-1
    F12=F11+1
    F13=F12+CNI+1
    F14=F13-1
    F21=F11+(CNJ+1)*(CNI+1)
    IF (KP.EQ.NK) F21=(CNI+1)*(JP-1)+IP-1
    F22=F21+1
    F23=F22+CNI+1
    F24=F23-1
    CALL CAL_FIJ_CY(F11,F12,F13,F14,F21,F22,F23,F24)
    WRITE(2,235) '4(',F11,' ',F12,' ',F22,' ',F21,')'
    FAD1=F_DOWN_CY(IP,KP)
    WRITE(3,231) FAD1
    OWN=CELL_ID_CY(IP,JP,KP)
    WRITE(4,231) OWN
    ENDDO   
    ENDDO    
        
    N=2+4*(PNL-1)
    WRITE(31,231) N+4
    WRITE(31,*) '('  
    WRITE(31,*) '0'
    WRITE(31,*) '1'
    WRITE(31,*) '2' 
    WRITE(31,*) '3' 
    DO I=1,N
    WRITE(31,*) '-1' 
    ENDDO
    WRITE(31,*) ')'
    
    WRITE(21,231) N+4
    WRITE(21,*) "("    
    WRITE(21,*) "    INLET1"
    WRITE(21,*) "    {"    
    WRITE(21,*) "        type            patch;"
    N=0
    WRITE(21,232) "        nFaces          ",N,";"
    N=F_INTERNAL_NUM
    WRITE(21,232) "        startFace       ",N,";"
    WRITE(21,*) "    }"
    WRITE(21,*) "    INLET2"
    WRITE(21,*) "    {"    
    WRITE(21,*) "        type            patch;"
    N=0
    WRITE(21,232) "        nFaces          ",N,";"
    N=F_INTERNAL_NUM
    WRITE(21,232) "        startFace       ",N,";"
    WRITE(21,*) "    }"
    WRITE(21,*) "    OUTLET"
    WRITE(21,*) "    {"
    WRITE(21,*) "        type            patch;"
    N=0
    WRITE(21,232) "        nFaces          ",N,";"
    N=F_INTERNAL_NUM
    WRITE(21,232) "        startFace       ",N,";"
    WRITE(21,*) "    }" 
    WRITE(21,*) "    HULL"
    WRITE(21,*) "    {"
    WRITE(21,*) "        type            wall;"
    WRITE(21,*) "        inGroups        1(wall);"
    N=CNJ*NK
    WRITE(21,232) "        nFaces          ",N,";"
    N=F_INTERNAL_NUM
    WRITE(21,232) "        startFace       ",N,";"
    WRITE(21,*) "    }"
       
        !LEFT FACE--PROCESSOR_MI1_TO_MI2  !CASE CYLINDER 3-1
    PRO1=PNI*(MJ-1)+MI-1
    PRO2=PNI*(MJ-1)+MI-2
    WRITE(21,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(21,*) '    {'
    WRITE(21,*) '        type            processor;'
    WRITE(21,*) '        inGroups        List<word> 1(processor);'
    N=CNJ*NK
    WRITE(21,232) "        nFaces          ",N,";"
    N=F_INTERNAL_NUM+CNJ*NK
    WRITE(21,232) "        startFace       ",N,";"
    WRITE(21,*) '        matchTolerance  0.0001;'
    WRITE(21,232) "        myProcNo        ",PRO1,";"
    WRITE(21,232) "        neighbProcNo    ",PRO2,";"
    WRITE(21,*) '    }'
    WRITE(41,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(41,*) '    {'
    WRITE(41,*) '        type            processor;'
    WRITE(41,*) '        value           nonuniform List<scalar> '
    N=CNJ*NK
    WRITE(41,'(I0)') N
    WRITE(41,*) '('
    DO KP=1,NK
    DO JP=1,CNJ
    WRITE(41,'(E21.14)') P_RIGHT_CY(JP,KP)    
    ENDDO
    ENDDO
    WRITE(41,*) ')'
    WRITE(41,*) ';'
    WRITE(41,*) '    }'
    WRITE(51,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(51,*) '    {'
    WRITE(51,*) '        type            processor;'
    WRITE(51,*) '        value           nonuniform List<vector> '
    N=CNJ*NK
    WRITE(51,'(I0)') N
    WRITE(51,*) '('
    DO KP=1,NK
    DO JP=1,CNJ
    UU1=U_RIGHT_CY(1,JP,KP)
    UU2=U_RIGHT_CY(2,JP,KP)
    UU3=U_RIGHT_CY(3,JP,KP)
    WRITE(51,230) '( ',UU1,' ',UU2,' ',UU3,' )'  
    ENDDO
    ENDDO
    WRITE(51,*) ')'
    WRITE(51,*) ';'
    WRITE(51,*) '    }'
    WRITE(61,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(61,*) '    {'
    WRITE(61,*) '        type            processor;'
    WRITE(61,*) "        value           uniform 0;"
    WRITE(61,*) '    }'
        !UP FACE--PROCESSOR_MI1_TO_MI2  !CASE CYLINDER 3-1
    PRO2=PNI*MJ+MI-1
    WRITE(21,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(21,*) '    {'
    WRITE(21,*) '        type            processor;'
    WRITE(21,*) '        inGroups        List<word> 1(processor);'
    N=CNI*NK
    WRITE(21,232) "        nFaces          ",N,";"
    N=F_INTERNAL_NUM+CNJ*NK+CNJ*NK
    WRITE(21,232) "        startFace       ",N,";"
    WRITE(21,*) '        matchTolerance  0.0001;'
    WRITE(21,232) "        myProcNo        ",PRO1,";"
    WRITE(21,232) "        neighbProcNo    ",PRO2,";"
    WRITE(21,*) '    }'
    WRITE(41,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(41,*) '    {'
    WRITE(41,*) '        type            processor;'
    WRITE(41,*) '        value           nonuniform List<scalar> '
    N=CNI*NK
    WRITE(41,'(I0)') N
    WRITE(41,*) '('
    DO KP=1,NK
    DO IP=1,CNI
    WRITE(41,'(E21.14)') P_LEFT_CY(IP,KP)    
    ENDDO
    ENDDO
    WRITE(41,*) ')'
    WRITE(41,*) ';'
    WRITE(41,*) '    }'
    WRITE(51,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(51,*) '    {'
    WRITE(51,*) '        type            processor;'
    WRITE(51,*) '        value           nonuniform List<vector> '
    N=CNI*NK
    WRITE(51,'(I0)') N
    WRITE(51,*) '('
    DO KP=1,NK
    DO IP=1,CNI
    UU1=U_UP_CY(1,IP,KP)
    UU2=U_UP_CY(2,IP,KP)
    UU3=U_UP_CY(3,IP,KP)
    WRITE(51,230) '( ',UU1,' ',UU2,' ',UU3,' )'  
    ENDDO
    ENDDO
    WRITE(51,*) ')'
    WRITE(51,*) ';'
    WRITE(51,*) '    }'
    WRITE(61,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(61,*) '    {'
    WRITE(61,*) '        type            processor;'
    WRITE(61,*) "        value           uniform 0;"
    WRITE(61,*) '    }'
        !DOWN FACE--PROCESSOR_MI1_TO_MI2(RECTANGLE 1-1)  !CASE CYLINDER 3-1
    MJ1=1
    MJ2=1
    PRO2=PNI*PNJ+PNL*PNI*(MJ2-1)+PNI*(MJ1-1)+MI-1
    WRITE(21,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(21,*) '    {'
    WRITE(21,*) '        type            processor;'
    WRITE(21,*) '        inGroups        List<word> 1(processor);'
    N=CNI*CNL*2
    WRITE(21,232) "        nFaces          ",N,";"
    N=F_INTERNAL_NUM+CNJ*NK+CNJ*NK+CNI*NK
    WRITE(21,232) "        startFace       ",N,";"
    WRITE(21,*) '        matchTolerance  0.0001;'
    WRITE(21,232) "        myProcNo        ",PRO1,";"
    WRITE(21,232) "        neighbProcNo    ",PRO2,";"
    WRITE(21,*) '    }'
    WRITE(41,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(41,*) '    {'
    WRITE(41,*) '        type            processor;'
    WRITE(41,*) '        value           nonuniform List<scalar> '
    N=CNI*CNL*2
    WRITE(41,'(I0)') N
    WRITE(41,*) '('
    DO KP=1,CNL
    DO IP=1,CNI
    WRITE(41,'(E21.14)') P_DOWN_CY(IP,KP)    
    ENDDO
    ENDDO
    DO KP=NK-CNL+1,NK
    DO IP=1,CNI
    WRITE(41,'(E21.14)') P_DOWN_CY(IP,KP)    
    ENDDO
    ENDDO
    WRITE(41,*) ')'
    WRITE(41,*) ';'
    WRITE(41,*) '    }'
    WRITE(51,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(51,*) '    {'
    WRITE(51,*) '        type            processor;'
    WRITE(51,*) '        value           nonuniform List<vector> '
    N=CNI*CNL*2
    WRITE(51,'(I0)') N
    WRITE(51,*) '('
    DO KP=1,CNL
    DO IP=1,CNI
    UU1=U_DOWN_CY(1,IP,KP)
    UU2=U_DOWN_CY(2,IP,KP)
    UU3=U_DOWN_CY(3,IP,KP)
    WRITE(51,230) '( ',UU1,' ',UU2,' ',UU3,' )'  
    ENDDO
    ENDDO
    DO KP=NK-CNL+1,NK
    DO IP=1,CNI
    UU1=U_DOWN_CY(1,IP,KP)
    UU2=U_DOWN_CY(2,IP,KP)
    UU3=U_DOWN_CY(3,IP,KP)
    WRITE(51,230) '( ',UU1,' ',UU2,' ',UU3,' )'  
    ENDDO
    ENDDO
    WRITE(51,*) ')'
    WRITE(51,*) ';'
    WRITE(51,*) '    }'
    WRITE(61,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(61,*) '    {'
    WRITE(61,*) '        type            processor;'
    WRITE(61,*) "        value           uniform 0;"
    WRITE(61,*) '    }'
    
        !DOWN FACE--PROCESSOR_MI1_TO_MI2(RECTANGLE 2-1)  !CASE CYLINDER 3-1
    DO MJ1=2,PNL-1
    MJ2=1
    PRO2=PNI*PNJ+PNL*PNI*(MJ2-1)+PNI*(MJ1-1)+MI-1
    WRITE(21,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(21,*) '    {'
    WRITE(21,*) '        type            processor;'
    WRITE(21,*) '        inGroups        List<word> 1(processor);'
    N=CNI*CNL
    WRITE(21,232) "        nFaces          ",N,";"
    N=F_INTERNAL_NUM+CNJ*NK+CNJ*NK+CNI*NK+CNI*CNL*2+(MJ1-2)*CNI*CNL
    WRITE(21,232) "        startFace       ",N,";"
    WRITE(21,*) '        matchTolerance  0.0001;'
    WRITE(21,232) "        myProcNo        ",PRO1,";"
    WRITE(21,232) "        neighbProcNo    ",PRO2,";"
    WRITE(21,*) '    }'
    WRITE(41,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(41,*) '    {'
    WRITE(41,*) '        type            processor;'
    WRITE(41,*) '        value           nonuniform List<scalar> '
    N=CNI*CNL
    WRITE(41,'(I0)') N
    WRITE(41,*) '('
    DO KP=CNL*(MJ1-1)+1,CNL*MJ1
    DO IP=1,CNI
    WRITE(41,'(E21.14)') P_DOWN_CY(IP,KP)    
    ENDDO
    ENDDO
    WRITE(41,*) ')'
    WRITE(41,*) ';'
    WRITE(41,*) '    }'
    WRITE(51,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(51,*) '    {'
    WRITE(51,*) '        type            processor;'
    WRITE(51,*) '        value           nonuniform List<vector> '
    N=CNI*CNL
    WRITE(51,'(I0)') N
    WRITE(51,*) '('
    DO KP=CNL*(MJ1-1)+1,CNL*MJ1
    DO IP=1,CNI
    UU1=U_DOWN_CY(1,IP,KP)
    UU2=U_DOWN_CY(2,IP,KP)
    UU3=U_DOWN_CY(3,IP,KP)
    WRITE(51,230) '( ',UU1,' ',UU2,' ',UU3,' )'  
    ENDDO
    ENDDO
    WRITE(51,*) ')'
    WRITE(51,*) ';'
    WRITE(51,*) '    }'
    WRITE(61,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(61,*) '    {'
    WRITE(61,*) '        type            processor;'
    WRITE(61,*) "        value           uniform 0;"
    WRITE(61,*) '    }'    
    ENDDO
    
        !DOWN FACE--PROCESSOR_MI1_TO_MI2(RECTANGLE 3-1)  !CASE CYLINDER 3-1
    MJ1=PNL
    MJ2=1
    PRO2=PNI*PNJ+PNL*PNI*(MJ2-1)+PNI*(MJ1-1)+MI-1
    WRITE(21,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(21,*) '    {'
    WRITE(21,*) '        type            processor;'
    WRITE(21,*) '        inGroups        List<word> 1(processor);'
    N=CNI*CNL*2
    WRITE(21,232) "        nFaces          ",N,";"
    N=F_INTERNAL_NUM+CNJ*NK+CNJ*NK+CNI*NK+CNI*CNL*2+(PNL-2)*CNI*CNL
    WRITE(21,232) "        startFace       ",N,";"
    WRITE(21,*) '        matchTolerance  0.0001;'
    WRITE(21,232) "        myProcNo        ",PRO1,";"
    WRITE(21,232) "        neighbProcNo    ",PRO2,";"
    WRITE(21,*) '    }'
    WRITE(41,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(41,*) '    {'
    WRITE(41,*) '        type            processor;'
    WRITE(41,*) '        value           nonuniform List<scalar> '
    N=CNI*CNL*2
    WRITE(41,'(I0)') N
    WRITE(41,*) '('
    DO KP=NL-CNL+1,NL
    DO IP=1,CNI
    WRITE(41,'(E21.14)') P_DOWN_CY(IP,KP)    
    ENDDO
    ENDDO
    DO KP=NL+1,NL+CNL
    DO IP=1,CNI
    WRITE(41,'(E21.14)') P_DOWN_CY(IP,KP)    
    ENDDO
    ENDDO
    WRITE(41,*) ')'
    WRITE(41,*) ';'
    WRITE(41,*) '    }'
    WRITE(51,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(51,*) '    {'
    WRITE(51,*) '        type            processor;'
    WRITE(51,*) '        value           nonuniform List<vector> '
    N=CNI*CNL*2
    WRITE(51,'(I0)') N
    WRITE(51,*) '('
    DO KP=NL-CNL+1,NL
    DO IP=1,CNI
    UU1=U_DOWN_CY(1,IP,KP)
    UU2=U_DOWN_CY(2,IP,KP)
    UU3=U_DOWN_CY(3,IP,KP)
    WRITE(51,230) '( ',UU1,' ',UU2,' ',UU3,' )'  
    ENDDO
    ENDDO
    DO KP=NL+1,NL+CNL
    DO IP=1,CNI
    UU1=U_DOWN_CY(1,IP,KP)
    UU2=U_DOWN_CY(2,IP,KP)
    UU3=U_DOWN_CY(3,IP,KP)
    WRITE(51,230) '( ',UU1,' ',UU2,' ',UU3,' )'  
    ENDDO
    ENDDO
    WRITE(51,*) ')'
    WRITE(51,*) ';'
    WRITE(51,*) '    }'
    WRITE(61,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(61,*) '    {'
    WRITE(61,*) '        type            processor;'
    WRITE(61,*) "        value           uniform 0;"
    WRITE(61,*) '    }'
    
        !DOWN FACE--PROCESSOR_MI1_TO_MI2(RECTANGLE 1-2 / 3-2)  !CASE CYLINDER 3-1
    DO MJ2=2,PNL-1
    MJ1=1
    PRO2=PNI*PNJ+PNL*PNI*(MJ2-1)+PNI*(MJ1-1)+MI-1
    WRITE(21,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(21,*) '    {'
    WRITE(21,*) '        type            processor;'
    WRITE(21,*) '        inGroups        List<word> 1(processor);'
    N=CNI*CNL
    WRITE(21,232) "        nFaces          ",N,";"
    N=F_INTERNAL_NUM+CNJ*NK+CNJ*NK+CNI*NK+CNI*CNL*2+(PNL-2)*CNI*CNL+CNI*CNL*2+(MJ2-2)*CNI*CNL*2
    WRITE(21,232) "        startFace       ",N,";"
    WRITE(21,*) '        matchTolerance  0.0001;'
    WRITE(21,232) "        myProcNo        ",PRO1,";"
    WRITE(21,232) "        neighbProcNo    ",PRO2,";"
    WRITE(21,*) '    }'
    WRITE(41,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(41,*) '    {'
    WRITE(41,*) '        type            processor;'
    WRITE(41,*) '        value           nonuniform List<scalar> '
    N=CNI*CNL
    WRITE(41,'(I0)') N
    WRITE(41,*) '('
    DO KP=3*NL+CNL*(PNL-MJ2)+1,3*NL+CNL*(PNL-MJ2+1)
    DO IP=1,CNI
    WRITE(41,'(E21.14)') P_DOWN_CY(IP,KP)    
    ENDDO
    ENDDO
    WRITE(41,*) ')'
    WRITE(41,*) ';'
    WRITE(41,*) '    }'
    WRITE(51,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(51,*) '    {'
    WRITE(51,*) '        type            processor;'
    WRITE(51,*) '        value           nonuniform List<vector> '
    N=CNI*CNL
    WRITE(51,'(I0)') N
    WRITE(51,*) '('
    DO KP=3*NL+CNL*(PNL-MJ2)+1,3*NL+CNL*(PNL-MJ2+1)
    DO IP=1,CNI
    UU1=U_DOWN_CY(1,IP,KP)
    UU2=U_DOWN_CY(2,IP,KP)
    UU3=U_DOWN_CY(3,IP,KP)
    WRITE(51,230) '( ',UU1,' ',UU2,' ',UU3,' )'  
    ENDDO
    ENDDO
    WRITE(51,*) ')'
    WRITE(51,*) ';'
    WRITE(51,*) '    }'
    WRITE(61,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(61,*) '    {'
    WRITE(61,*) '        type            processor;'
    WRITE(61,*) "        value           uniform 0;"
    WRITE(61,*) '    }' 
    
    MJ1=PNL
    PRO2=PNI*PNJ+PNL*PNI*(MJ2-1)+PNI*(MJ1-1)+MI-1
    WRITE(21,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(21,*) '    {'
    WRITE(21,*) '        type            processor;'
    WRITE(21,*) '        inGroups        List<word> 1(processor);'
    N=CNI*CNL
    WRITE(21,232) "        nFaces          ",N,";"
    N=F_INTERNAL_NUM+CNJ*NK+CNJ*NK+CNI*NK+CNI*CNL*2+(PNL-2)*CNI*CNL+CNI*CNL*2+(MJ2-2)*CNI*CNL*2+CNI*CNL
    WRITE(21,232) "        startFace       ",N,";"
    WRITE(21,*) '        matchTolerance  0.0001;'
    WRITE(21,232) "        myProcNo        ",PRO1,";"
    WRITE(21,232) "        neighbProcNo    ",PRO2,";"
    WRITE(21,*) '    }'
    WRITE(41,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(41,*) '    {'
    WRITE(41,*) '        type            processor;'
    WRITE(41,*) '        value           nonuniform List<scalar> '
    N=CNI*CNL
    WRITE(41,'(I0)') N
    WRITE(41,*) '('
    DO KP=NL+CNL*(MJ2-1)+1,NL+CNL*MJ2
    DO IP=1,CNI
    WRITE(41,'(E21.14)') P_DOWN_CY(IP,KP)    
    ENDDO
    ENDDO
    WRITE(41,*) ')'
    WRITE(41,*) ';'
    WRITE(41,*) '    }'
    WRITE(51,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(51,*) '    {'
    WRITE(51,*) '        type            processor;'
    WRITE(51,*) '        value           nonuniform List<vector> '
    N=CNI*CNL
    WRITE(51,'(I0)') N
    WRITE(51,*) '('
    DO KP=NL+CNL*(MJ2-1)+1,NL+CNL*MJ2
    DO IP=1,CNI
    UU1=U_DOWN_CY(1,IP,KP)
    UU2=U_DOWN_CY(2,IP,KP)
    UU3=U_DOWN_CY(3,IP,KP)
    WRITE(51,230) '( ',UU1,' ',UU2,' ',UU3,' )'  
    ENDDO
    ENDDO
    WRITE(51,*) ')'
    WRITE(51,*) ';'
    WRITE(51,*) '    }'
    WRITE(61,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(61,*) '    {'
    WRITE(61,*) '        type            processor;'
    WRITE(61,*) "        value           uniform 0;"
    WRITE(61,*) '    }'
        
    ENDDO
    
    
        !DOWN FACE--PROCESSOR_MI1_TO_MI2(RECTANGLE 1-3)  !CASE CYLINDER 3-1
    MJ1=1
    MJ2=PNL
    PRO2=PNI*PNJ+PNL*PNI*(MJ2-1)+PNI*(MJ1-1)+MI-1
    WRITE(21,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(21,*) '    {'
    WRITE(21,*) '        type            processor;'
    WRITE(21,*) '        inGroups        List<word> 1(processor);'
    N=CNI*CNL*2
    WRITE(21,232) "        nFaces          ",N,";"
    N=F_INTERNAL_NUM+CNJ*NK+CNJ*NK+CNI*NK+CNI*CNL*2+(PNL-2)*CNI*CNL+CNI*CNL*2+(PNL-2)*CNI*CNL*2
    WRITE(21,232) "        startFace       ",N,";"
    WRITE(21,*) '        matchTolerance  0.0001;'
    WRITE(21,232) "        myProcNo        ",PRO1,";"
    WRITE(21,232) "        neighbProcNo    ",PRO2,";"
    WRITE(21,*) '    }'
    WRITE(41,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(41,*) '    {'
    WRITE(41,*) '        type            processor;'
    WRITE(41,*) '        value           nonuniform List<scalar> '
    N=CNI*CNL*2
    WRITE(41,'(I0)') N
    WRITE(41,*) '('
    DO KP=3*NL-CNL+1,3*NL
    DO IP=1,CNI
    WRITE(41,'(E21.14)') P_DOWN_CY(IP,KP)    
    ENDDO
    ENDDO
    DO KP=3*NL+1,3*NL+CNL
    DO IP=1,CNI
    WRITE(41,'(E21.14)') P_DOWN_CY(IP,KP)    
    ENDDO
    ENDDO
    WRITE(41,*) ')'
    WRITE(41,*) ';'
    WRITE(41,*) '    }'
    WRITE(51,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(51,*) '    {'
    WRITE(51,*) '        type            processor;'
    WRITE(51,*) '        value           nonuniform List<vector> '
    N=CNI*CNL*2
    WRITE(51,'(I0)') N
    WRITE(51,*) '('
    DO KP=3*NL-CNL+1,3*NL
    DO IP=1,CNI
    UU1=U_DOWN_CY(1,IP,KP)
    UU2=U_DOWN_CY(2,IP,KP)
    UU3=U_DOWN_CY(3,IP,KP)
    WRITE(51,230) '( ',UU1,' ',UU2,' ',UU3,' )'  
    ENDDO
    ENDDO
    DO KP=3*NL+1,3*NL+CNL
    DO IP=1,CNI
    UU1=U_DOWN_CY(1,IP,KP)
    UU2=U_DOWN_CY(2,IP,KP)
    UU3=U_DOWN_CY(3,IP,KP)
    WRITE(51,230) '( ',UU1,' ',UU2,' ',UU3,' )'  
    ENDDO
    ENDDO
    WRITE(51,*) ')'
    WRITE(51,*) ';'
    WRITE(51,*) '    }'
    WRITE(61,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(61,*) '    {'
    WRITE(61,*) '        type            processor;'
    WRITE(61,*) "        value           uniform 0;"
    WRITE(61,*) '    }'
    
    
        !DOWN FACE--PROCESSOR_MI1_TO_MI2(RECTANGLE 2-3)  !CASE CYLINDER 3-1
    DO MJ1=2,PNL-1
    MJ2=PNL
    PRO2=PNI*PNJ+PNL*PNI*(MJ2-1)+PNI*(MJ1-1)+MI-1
    WRITE(21,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(21,*) '    {'
    WRITE(21,*) '        type            processor;'
    WRITE(21,*) '        inGroups        List<word> 1(processor);'
    N=CNI*CNL
    WRITE(21,232) "        nFaces          ",N,";"
    N=F_INTERNAL_NUM+CNJ*NK+CNJ*NK+CNI*NK+CNI*CNL*2+(PNL-2)*CNI*CNL+CNI*CNL*2+(PNL-2)*CNI*CNL*2+CNI*CNL*2+(MJ1-2)*CNI*CNL
    WRITE(21,232) "        startFace       ",N,";"
    WRITE(21,*) '        matchTolerance  0.0001;'
    WRITE(21,232) "        myProcNo        ",PRO1,";"
    WRITE(21,232) "        neighbProcNo    ",PRO2,";"
    WRITE(21,*) '    }'
    WRITE(41,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(41,*) '    {'
    WRITE(41,*) '        type            processor;'
    WRITE(41,*) '        value           nonuniform List<scalar> '
    N=CNI*CNL
    WRITE(41,'(I0)') N
    WRITE(41,*) '('
    DO KP=2*NL+CNL*(PNL-MJ1)+1,2*NL+CNL*(PNL-MJ1+1)
    DO IP=1,CNI
    WRITE(41,'(E21.14)') P_DOWN_CY(IP,KP)    
    ENDDO
    ENDDO
    WRITE(41,*) ')'
    WRITE(41,*) ';'
    WRITE(41,*) '    }'
    WRITE(51,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(51,*) '    {'
    WRITE(51,*) '        type            processor;'
    WRITE(51,*) '        value           nonuniform List<vector> '
    N=CNI*CNL
    WRITE(51,'(I0)') N
    WRITE(51,*) '('
    DO KP=2*NL+CNL*(PNL-MJ1)+1,2*NL+CNL*(PNL-MJ1+1)
    DO IP=1,CNI
    UU1=U_DOWN_CY(1,IP,KP)
    UU2=U_DOWN_CY(2,IP,KP)
    UU3=U_DOWN_CY(3,IP,KP)
    WRITE(51,230) '( ',UU1,' ',UU2,' ',UU3,' )'  
    ENDDO
    ENDDO
    WRITE(51,*) ')'
    WRITE(51,*) ';'
    WRITE(51,*) '    }'
    WRITE(61,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(61,*) '    {'
    WRITE(61,*) '        type            processor;'
    WRITE(61,*) "        value           uniform 0;"
    WRITE(61,*) '    }'    
    ENDDO
    
    
    
        !DOWN FACE--PROCESSOR_MI1_TO_MI2(RECTANGLE 3-3)  !CASE CYLINDER 3-1
    MJ1=PNL
    MJ2=PNL
    PRO2=PNI*PNJ+PNL*PNI*(MJ2-1)+PNI*(MJ1-1)+MI-1
    WRITE(21,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(21,*) '    {'
    WRITE(21,*) '        type            processor;'
    WRITE(21,*) '        inGroups        List<word> 1(processor);'
    N=CNI*CNL*2
    WRITE(21,232) "        nFaces          ",N,";"
    N=F_INTERNAL_NUM+CNJ*NK+CNJ*NK+CNI*NK+CNI*CNL*2+(PNL-2)*CNI*CNL+CNI*CNL*2+(PNL-2)*CNI*CNL*2+CNI*CNL*2+(PNL-2)*CNI*CNL
    WRITE(21,232) "        startFace       ",N,";"
    WRITE(21,*) '        matchTolerance  0.0001;'
    WRITE(21,232) "        myProcNo        ",PRO1,";"
    WRITE(21,232) "        neighbProcNo    ",PRO2,";"
    WRITE(21,*) '    }'
    WRITE(41,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(41,*) '    {'
    WRITE(41,*) '        type            processor;'
    WRITE(41,*) '        value           nonuniform List<scalar> '
    N=CNI*CNL*2
    WRITE(41,'(I0)') N
    WRITE(41,*) '('
    DO KP=2*NL-CNL+1,2*NL
    DO IP=1,CNI
    WRITE(41,'(E21.14)') P_DOWN_CY(IP,KP)    
    ENDDO
    ENDDO
    DO KP=2*NL+1,2*NL+CNL
    DO IP=1,CNI
    WRITE(41,'(E21.14)') P_DOWN_CY(IP,KP)    
    ENDDO
    ENDDO
    WRITE(41,*) ')'
    WRITE(41,*) ';'
    WRITE(41,*) '    }'
    WRITE(51,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(51,*) '    {'
    WRITE(51,*) '        type            processor;'
    WRITE(51,*) '        value           nonuniform List<vector> '
    N=CNI*CNL*2
    WRITE(51,'(I0)') N
    WRITE(51,*) '('
    DO KP=2*NL-CNL+1,2*NL
    DO IP=1,CNI
    UU1=U_DOWN_CY(1,IP,KP)
    UU2=U_DOWN_CY(2,IP,KP)
    UU3=U_DOWN_CY(3,IP,KP)
    WRITE(51,230) '( ',UU1,' ',UU2,' ',UU3,' )'  
    ENDDO
    ENDDO
    DO KP=2*NL+1,2*NL+CNL
    DO IP=1,CNI
    UU1=U_DOWN_CY(1,IP,KP)
    UU2=U_DOWN_CY(2,IP,KP)
    UU3=U_DOWN_CY(3,IP,KP)
    WRITE(51,230) '( ',UU1,' ',UU2,' ',UU3,' )'  
    ENDDO
    ENDDO
    WRITE(51,*) ')'
    WRITE(51,*) ';'
    WRITE(51,*) '    }'
    WRITE(61,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(61,*) '    {'
    WRITE(61,*) '        type            processor;'
    WRITE(61,*) "        value           uniform 0;"
    WRITE(61,*) '    }'    
        
        
        
        
        
        
        
        
        
        
        
        
    
    ELSEIF (MJ.GT.1 .AND. MJ.LT.PNJ) THEN
    !CASE CYLINDER 3-2
        
        !LEFT FACE--HULL  !CASE CYLINDER 3-2
    IP=CNI
    I=IP+CNI*(MI-1)
    DO KP=1,NK
    K=KP
    DO JP=1,CNJ
    J=JP+CNJ*(MJ-1)
    F11=(CNJ+1)*(CNI+1)*(KP-1)+(CNI+1)*(JP-1)+IP-1
    F12=F11+1
    F13=F12+CNI+1
    F14=F13-1
    F21=F11+(CNJ+1)*(CNI+1)
    IF (KP.EQ.NK) F21=(CNI+1)*(JP-1)+IP-1
    F22=F21+1
    F23=F22+CNI+1
    F24=F23-1
    CALL CAL_FIJ_CY(F11,F12,F13,F14,F21,F22,F23,F24)
    WRITE(2,235) '4(',F12,' ',F13,' ',F23,' ',F22,')'
    FAD1=F_INTERNAL_NUM_ALL+NL*NL+NJ*NK+NI*NK+NJ*(K-1)+J
    WRITE(3,231) FAD1
    OWN=CELL_ID_CY(IP,JP,KP)
    WRITE(4,231) OWN
    ENDDO   
    ENDDO  
    
        !DOWN FACE--PROCESSOR_MI1_TO_MI2 !CASE CYLINDER 3-2 
    JP=1
    J=JP+CNJ*(MJ-1)
    DO KP=1,NK
    K=KP
    DO IP=1,CNI
    I=IP+CNI*(MI-1)
    F11=(CNJ+1)*(CNI+1)*(KP-1)+(CNI+1)*(JP-1)+IP-1
    F12=F11+1
    F13=F12+CNI+1
    F14=F13-1
    F21=F11+(CNJ+1)*(CNI+1)
    IF (KP.EQ.NK) F21=(CNI+1)*(JP-1)+IP-1
    F22=F21+1
    F23=F22+CNI+1
    F24=F23-1
    CALL CAL_FIJ_CY(F11,F12,F13,F14,F21,F22,F23,F24)
    WRITE(2,235) '4(',F11,' ',F12,' ',F22,' ',F21,')'
    FAD1=F_DOWN_CY(IP,KP)
    WRITE(3,231) FAD1
    OWN=CELL_ID_CY(IP,JP,KP)
    WRITE(4,231) OWN
    ENDDO   
    ENDDO   
    
        !RIGHT FACE--PROCESSOR_MI1_TO_MI2  !CASE CYLINDER 3-2
    IP=1
    I=IP+CNI*(MI-1)
    DO KP=1,NK
    K=KP
    DO JP=1,CNJ
    J=JP+CNJ*(MJ-1)
    F11=(CNJ+1)*(CNI+1)*(KP-1)+(CNI+1)*(JP-1)+IP-1
    F12=F11+1
    F13=F12+CNI+1
    F14=F13-1
    F21=F11+(CNJ+1)*(CNI+1)
    IF (KP.EQ.NK) F21=(CNI+1)*(JP-1)+IP-1
    F22=F21+1
    F23=F22+CNI+1
    F24=F23-1
    CALL CAL_FIJ_CY(F11,F12,F13,F14,F21,F22,F23,F24)
    WRITE(2,235) '4(',F11,' ',F21,' ',F24,' ',F14,')'
    FAD1=F_RIGHT_CY(JP,KP)
    WRITE(3,231) FAD1
    OWN=CELL_ID_CY(IP,JP,KP)
    WRITE(4,231) OWN
    ENDDO   
    ENDDO
    
        !UP FACE--PROCESSOR_MI1_TO_MI2  !CASE CYLINDER 3-2
    JP=CNJ
    J=JP+CNJ*(MJ-1)
    DO KP=1,NK
    K=KP
    DO IP=1,CNI
    I=IP+CNI*(MI-1)
    F11=(CNJ+1)*(CNI+1)*(KP-1)+(CNI+1)*(JP-1)+IP-1
    F12=F11+1
    F13=F12+CNI+1
    F14=F13-1
    F21=F11+(CNJ+1)*(CNI+1)
    IF (KP.EQ.NK) F21=(CNI+1)*(JP-1)+IP-1
    F22=F21+1
    F23=F22+CNI+1
    F24=F23-1
    CALL CAL_FIJ_CY(F11,F12,F13,F14,F21,F22,F23,F24)
    WRITE(2,235) '4(',F14,' ',F24,' ',F23,' ',F13,')'
    FAD1=F_UP_CY(IP,KP)
    WRITE(3,231) FAD1
    OWN=CELL_ID_CY(IP,JP,KP)
    WRITE(4,231) OWN
    ENDDO   
    ENDDO
     
    N=3
    WRITE(31,231) N+4
    WRITE(31,*) '('  
    WRITE(31,*) '0'
    WRITE(31,*) '1'
    WRITE(31,*) '2' 
    WRITE(31,*) '3' 
    DO I=1,N
    WRITE(31,*) '-1' 
    ENDDO
    WRITE(31,*) ')'
    
    WRITE(21,231) N+4
    WRITE(21,*) "("    
    WRITE(21,*) "    INLET1"
    WRITE(21,*) "    {"    
    WRITE(21,*) "        type            patch;"
    N=0
    WRITE(21,232) "        nFaces          ",N,";"
    N=F_INTERNAL_NUM
    WRITE(21,232) "        startFace       ",N,";"
    WRITE(21,*) "    }"
    WRITE(21,*) "    INLET2"
    WRITE(21,*) "    {"    
    WRITE(21,*) "        type            patch;"
    N=0
    WRITE(21,232) "        nFaces          ",N,";"
    N=F_INTERNAL_NUM
    WRITE(21,232) "        startFace       ",N,";"
    WRITE(21,*) "    }"
    WRITE(21,*) "    OUTLET"
    WRITE(21,*) "    {"
    WRITE(21,*) "        type            patch;"
    N=0
    WRITE(21,232) "        nFaces          ",N,";"
    N=F_INTERNAL_NUM
    WRITE(21,232) "        startFace       ",N,";"
    WRITE(21,*) "    }"
    WRITE(21,*) "    HULL"
    WRITE(21,*) "    {"
    WRITE(21,*) "        type            wall;"
    WRITE(21,*) "        inGroups        1(wall);"
    N=CNJ*NK
    WRITE(21,232) "        nFaces          ",N,";"
    N=F_INTERNAL_NUM
    WRITE(21,232) "        startFace       ",N,";"
    WRITE(21,*) "    }"
    
        !DOWN FACE--PROCESSOR_MI1_TO_MI2 !CASE CYLINDER 3-2 
    PRO1=PNI*(MJ-1)+MI-1
    PRO2=PNI*(MJ-2)+MI-1
    WRITE(21,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(21,*) '    {'
    WRITE(21,*) '        type            processor;'
    WRITE(21,*) '        inGroups        List<word> 1(processor);'
    N=CNI*NK
    WRITE(21,232) "        nFaces          ",N,";"
    N=F_INTERNAL_NUM+CNJ*NK
    WRITE(21,232) "        startFace       ",N,";"
    WRITE(21,*) '        matchTolerance  0.0001;'
    WRITE(21,232) "        myProcNo        ",PRO1,";"
    WRITE(21,232) "        neighbProcNo    ",PRO2,";"
    WRITE(21,*) '    }'
    WRITE(41,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(41,*) '    {'
    WRITE(41,*) '        type            processor;'
    WRITE(41,*) '        value           nonuniform List<scalar> '
    N=CNI*NK
    WRITE(41,'(I0)') N
    WRITE(41,*) '('
    DO KP=1,NK
    DO IP=1,CNI
    WRITE(41,'(E21.14)') P_DOWN_CY(IP,KP)    
    ENDDO
    ENDDO
    WRITE(41,*) ')'
    WRITE(41,*) ';'
    WRITE(41,*) '    }'
    WRITE(51,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(51,*) '    {'
    WRITE(51,*) '        type            processor;'
    WRITE(51,*) '        value           nonuniform List<vector> '
    N=CNI*NK
    WRITE(51,'(I0)') N
    WRITE(51,*) '('
    DO KP=1,NK
    DO IP=1,CNI
    UU1=U_DOWN_CY(1,IP,KP)
    UU2=U_DOWN_CY(2,IP,KP)
    UU3=U_DOWN_CY(3,IP,KP)
    WRITE(51,230) '( ',UU1,' ',UU2,' ',UU3,' )'  
    ENDDO
    ENDDO
    WRITE(51,*) ')'
    WRITE(51,*) ';'
    WRITE(51,*) '    }'
    WRITE(61,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(61,*) '    {'
    WRITE(61,*) '        type            processor;'
    WRITE(61,*) "        value           uniform 0;"
    WRITE(61,*) '    }'
    
        !LEFT FACE--PROCESSOR_MI1_TO_MI2  !CASE CYLINDER 3-2
    PRO2=PNI*(MJ-1)+MI-2
    WRITE(21,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(21,*) '    {'
    WRITE(21,*) '        type            processor;'
    WRITE(21,*) '        inGroups        List<word> 1(processor);'
    N=CNJ*NK
    WRITE(21,232) "        nFaces          ",N,";"
    N=F_INTERNAL_NUM+CNJ*NK+CNI*NK
    WRITE(21,232) "        startFace       ",N,";"
    WRITE(21,*) '        matchTolerance  0.0001;'
    WRITE(21,232) "        myProcNo        ",PRO1,";"
    WRITE(21,232) "        neighbProcNo    ",PRO2,";"
    WRITE(21,*) '    }'
    WRITE(41,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(41,*) '    {'
    WRITE(41,*) '        type            processor;'
    WRITE(41,*) '        value           nonuniform List<scalar> '
    N=CNJ*NK
    WRITE(41,'(I0)') N
    WRITE(41,*) '('
    DO KP=1,NK
    DO JP=1,CNJ
    WRITE(41,'(E21.14)') P_RIGHT_CY(JP,KP)    
    ENDDO
    ENDDO
    WRITE(41,*) ')'
    WRITE(41,*) ';'
    WRITE(41,*) '    }'
    WRITE(51,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(51,*) '    {'
    WRITE(51,*) '        type            processor;'
    WRITE(51,*) '        value           nonuniform List<vector> '
    N=CNJ*NK
    WRITE(51,'(I0)') N
    WRITE(51,*) '('
    DO KP=1,NK
    DO JP=1,CNJ
    UU1=U_RIGHT_CY(1,JP,KP)
    UU2=U_RIGHT_CY(2,JP,KP)
    UU3=U_RIGHT_CY(3,JP,KP)
    WRITE(51,230) '( ',UU1,' ',UU2,' ',UU3,' )'  
    ENDDO
    ENDDO
    WRITE(51,*) ')'
    WRITE(51,*) ';'
    WRITE(51,*) '    }'
    WRITE(61,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(61,*) '    {'
    WRITE(61,*) '        type            processor;'
    WRITE(61,*) "        value           uniform 0;"
    WRITE(61,*) '    }'
    
        !UP FACE--PROCESSOR_MI1_TO_MI2  !CASE CYLINDER 3-2
    PRO2=PNI*MJ+MI-1
    WRITE(21,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(21,*) '    {'
    WRITE(21,*) '        type            processor;'
    WRITE(21,*) '        inGroups        List<word> 1(processor);'
    N=CNI*NK
    WRITE(21,232) "        nFaces          ",N,";"
    N=F_INTERNAL_NUM+CNJ*NK+CNI*NK+CNJ*NK
    WRITE(21,232) "        startFace       ",N,";"
    WRITE(21,*) '        matchTolerance  0.0001;'
    WRITE(21,232) "        myProcNo        ",PRO1,";"
    WRITE(21,232) "        neighbProcNo    ",PRO2,";"
    WRITE(21,*) '    }'
    WRITE(41,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(41,*) '    {'
    WRITE(41,*) '        type            processor;'
    WRITE(41,*) '        value           nonuniform List<scalar> '
    N=CNI*NK
    WRITE(41,'(I0)') N
    WRITE(41,*) '('
    DO KP=1,NK
    DO IP=1,CNI
    WRITE(41,'(E21.14)') P_UP_CY(IP,KP)    
    ENDDO
    ENDDO
    WRITE(41,*) ')'
    WRITE(41,*) ';'
    WRITE(41,*) '    }'
    WRITE(51,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(51,*) '    {'
    WRITE(51,*) '        type            processor;'
    WRITE(51,*) '        value           nonuniform List<vector> '
    N=CNI*NK
    WRITE(51,'(I0)') N
    WRITE(51,*) '('
    DO KP=1,NK
    DO IP=1,CNI
    UU1=U_UP_CY(1,IP,KP)
    UU2=U_UP_CY(2,IP,KP)
    UU3=U_UP_CY(3,IP,KP)
    WRITE(51,230) '( ',UU1,' ',UU2,' ',UU3,' )'  
    ENDDO
    ENDDO
    WRITE(51,*) ')'
    WRITE(51,*) ';'
    WRITE(51,*) '    }'
    WRITE(61,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(61,*) '    {'
    WRITE(61,*) '        type            processor;'
    WRITE(61,*) "        value           uniform 0;"
    WRITE(61,*) '    }'
    
    
    
    
    
    
    
    
    
    
    
    ELSEIF (MJ.EQ.PNJ) THEN
    !CASE CYLINDER 3-3
        
        !UP FACE--OUTLET  !CASE CYLINDER 3-3
    JP=CNJ
    J=JP+CNJ*(MJ-1)
    DO IP=1,CNI
    I=IP+CNI*(MI-1)
    DO KP=1,NL
    K=KP
    F11=(CNJ+1)*(CNI+1)*(KP-1)+(CNI+1)*(JP-1)+IP-1
    F12=F11+1
    F13=F12+CNI+1
    F14=F13-1
    F21=F11+(CNJ+1)*(CNI+1)
    IF (KP.EQ.NK) F21=(CNI+1)*(JP-1)+IP-1
    F22=F21+1
    F23=F22+CNI+1
    F24=F23-1
    CALL CAL_FIJ_CY(F11,F12,F13,F14,F21,F22,F23,F24)
    WRITE(2,235) '4(',F14,' ',F24,' ',F23,' ',F13,')'
    FAD1=F_INTERNAL_NUM_ALL+NL*NL+NJ*NK+NL*(I-1)+K
    WRITE(3,231) FAD1
    OWN=CELL_ID_CY(IP,JP,KP)
    WRITE(4,231) OWN
    ENDDO   
    ENDDO
    DO IP=1,CNI
    I=IP+CNI*(MI-1)
    DO KP=1+NL,2*NL
    K=KP
    F11=(CNJ+1)*(CNI+1)*(KP-1)+(CNI+1)*(JP-1)+IP-1
    F12=F11+1
    F13=F12+CNI+1
    F14=F13-1
    F21=F11+(CNJ+1)*(CNI+1)
    IF (KP.EQ.NK) F21=(CNI+1)*(JP-1)+IP-1
    F22=F21+1
    F23=F22+CNI+1
    F24=F23-1
    CALL CAL_FIJ_CY(F11,F12,F13,F14,F21,F22,F23,F24)
    WRITE(2,235) '4(',F14,' ',F24,' ',F23,' ',F13,')'
    FAD1=F_INTERNAL_NUM_ALL+NL*NL+NJ*NK+NI*NL+NL*(I-1)+(K-NL)
    WRITE(3,231) FAD1
    OWN=CELL_ID_CY(IP,JP,KP)
    WRITE(4,231) OWN
    ENDDO   
    ENDDO
    DO IP=1,CNI
    I=IP+CNI*(MI-1)
    DO KP=1+2*NL,3*NL
    K=KP
    F11=(CNJ+1)*(CNI+1)*(KP-1)+(CNI+1)*(JP-1)+IP-1
    F12=F11+1
    F13=F12+CNI+1
    F14=F13-1
    F21=F11+(CNJ+1)*(CNI+1)
    IF (KP.EQ.NK) F21=(CNI+1)*(JP-1)+IP-1
    F22=F21+1
    F23=F22+CNI+1
    F24=F23-1
    CALL CAL_FIJ_CY(F11,F12,F13,F14,F21,F22,F23,F24)
    WRITE(2,235) '4(',F14,' ',F24,' ',F23,' ',F13,')'
    FAD1=F_INTERNAL_NUM_ALL+NL*NL+NJ*NK+NI*NL*2+NL*(I-1)+(K-2*NL)
    WRITE(3,231) FAD1
    OWN=CELL_ID_CY(IP,JP,KP)
    WRITE(4,231) OWN
    ENDDO   
    ENDDO
    DO IP=1,CNI
    I=IP+CNI*(MI-1)
    DO KP=1+3*NL,NK
    K=KP
    F11=(CNJ+1)*(CNI+1)*(KP-1)+(CNI+1)*(JP-1)+IP-1
    F12=F11+1
    F13=F12+CNI+1
    F14=F13-1
    F21=F11+(CNJ+1)*(CNI+1)
    IF (KP.EQ.NK) F21=(CNI+1)*(JP-1)+IP-1
    F22=F21+1
    F23=F22+CNI+1
    F24=F23-1
    CALL CAL_FIJ_CY(F11,F12,F13,F14,F21,F22,F23,F24)
    WRITE(2,235) '4(',F14,' ',F24,' ',F23,' ',F13,')'
    FAD1=F_INTERNAL_NUM_ALL+NL*NL+NJ*NK+NI*NL*3+NL*(I-1)+(K-3*NL)
    WRITE(3,231) FAD1
    OWN=CELL_ID_CY(IP,JP,KP)
    WRITE(4,231) OWN
    ENDDO   
    ENDDO
    
        !LEFT FACE--HULL  !CASE CYLINDER 3-3
    IP=CNI
    I=IP+CNI*(MI-1)
    DO KP=1,NK
    K=KP
    DO JP=1,CNJ
    J=JP+CNJ*(MJ-1)
    F11=(CNJ+1)*(CNI+1)*(KP-1)+(CNI+1)*(JP-1)+IP-1
    F12=F11+1
    F13=F12+CNI+1
    F14=F13-1
    F21=F11+(CNJ+1)*(CNI+1)
    IF (KP.EQ.NK) F21=(CNI+1)*(JP-1)+IP-1
    F22=F21+1
    F23=F22+CNI+1
    F24=F23-1
    CALL CAL_FIJ_CY(F11,F12,F13,F14,F21,F22,F23,F24)
    WRITE(2,235) '4(',F12,' ',F13,' ',F23,' ',F22,')'
    FAD1=F_INTERNAL_NUM_ALL+NL*NL+NJ*NK+NI*NK+NJ*(K-1)+J
    WRITE(3,231) FAD1
    OWN=CELL_ID_CY(IP,JP,KP)
    WRITE(4,231) OWN
    ENDDO   
    ENDDO
    
        !DOWN FACE--PROCESSOR_MI1_TO_MI2 !CASE CYLINDER 3-3 
    JP=1
    J=JP+CNJ*(MJ-1)
    DO KP=1,NK
    K=KP
    DO IP=1,CNI
    I=IP+CNI*(MI-1)
    F11=(CNJ+1)*(CNI+1)*(KP-1)+(CNI+1)*(JP-1)+IP-1
    F12=F11+1
    F13=F12+CNI+1
    F14=F13-1
    F21=F11+(CNJ+1)*(CNI+1)
    IF (KP.EQ.NK) F21=(CNI+1)*(JP-1)+IP-1
    F22=F21+1
    F23=F22+CNI+1
    F24=F23-1
    CALL CAL_FIJ_CY(F11,F12,F13,F14,F21,F22,F23,F24)
    WRITE(2,235) '4(',F11,' ',F12,' ',F22,' ',F21,')'
    FAD1=F_DOWN_CY(IP,KP)
    WRITE(3,231) FAD1
    OWN=CELL_ID_CY(IP,JP,KP)
    WRITE(4,231) OWN
    ENDDO   
    ENDDO   
    
        !RIGHT FACE--PROCESSOR_MI1_TO_MI2  !CASE CYLINDER 3-3
    IP=1
    I=IP+CNI*(MI-1)
    DO KP=1,NK
    K=KP
    DO JP=1,CNJ
    J=JP+CNJ*(MJ-1)
    F11=(CNJ+1)*(CNI+1)*(KP-1)+(CNI+1)*(JP-1)+IP-1
    F12=F11+1
    F13=F12+CNI+1
    F14=F13-1
    F21=F11+(CNJ+1)*(CNI+1)
    IF (KP.EQ.NK) F21=(CNI+1)*(JP-1)+IP-1
    F22=F21+1
    F23=F22+CNI+1
    F24=F23-1
    CALL CAL_FIJ_CY(F11,F12,F13,F14,F21,F22,F23,F24)
    WRITE(2,235) '4(',F11,' ',F21,' ',F24,' ',F14,')'
    FAD1=F_RIGHT_CY(JP,KP)
    WRITE(3,231) FAD1
    OWN=CELL_ID_CY(IP,JP,KP)
    WRITE(4,231) OWN
    ENDDO   
    ENDDO
    
    N=2
    WRITE(31,231) N+4
    WRITE(31,*) '('  
    WRITE(31,*) '0'
    WRITE(31,*) '1'
    WRITE(31,*) '2' 
    WRITE(31,*) '3' 
    DO I=1,N
    WRITE(31,*) '-1' 
    ENDDO
    WRITE(31,*) ')'
    
    WRITE(21,231) N+4
    WRITE(21,*) "("    
    WRITE(21,*) "    INLET1"
    WRITE(21,*) "    {"    
    WRITE(21,*) "        type            patch;"
    N=0
    WRITE(21,232) "        nFaces          ",N,";"
    N=F_INTERNAL_NUM
    WRITE(21,232) "        startFace       ",N,";"
    WRITE(21,*) "    }"
    WRITE(21,*) "    INLET2"
    WRITE(21,*) "    {"    
    WRITE(21,*) "        type            patch;"
    N=0
    WRITE(21,232) "        nFaces          ",N,";"
    N=F_INTERNAL_NUM
    WRITE(21,232) "        startFace       ",N,";"
    WRITE(21,*) "    }"
    WRITE(21,*) "    OUTLET"
    WRITE(21,*) "    {"
    WRITE(21,*) "        type            patch;"
    N=CNI*NK
    WRITE(21,232) "        nFaces          ",N,";"
    N=F_INTERNAL_NUM
    WRITE(21,232) "        startFace       ",N,";"
    WRITE(21,*) "    }"
    WRITE(21,*) "    HULL"
    WRITE(21,*) "    {"
    WRITE(21,*) "        type            wall;"
    WRITE(21,*) "        inGroups        1(wall);"
    N=CNJ*NK
    WRITE(21,232) "        nFaces          ",N,";"
    N=F_INTERNAL_NUM+CNI*NK
    WRITE(21,232) "        startFace       ",N,";"
    WRITE(21,*) "    }"
    
        !DOWN FACE--PROCESSOR_MI1_TO_MI2 !CASE CYLINDER 3-3 
    PRO1=PNI*(MJ-1)+MI-1
    PRO2=PNI*(MJ-2)+MI-1
    WRITE(21,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(21,*) '    {'
    WRITE(21,*) '        type            processor;'
    WRITE(21,*) '        inGroups        List<word> 1(processor);'
    N=CNI*NK
    WRITE(21,232) "        nFaces          ",N,";"
    N=F_INTERNAL_NUM+CNJ*NK+CNI*NK
    WRITE(21,232) "        startFace       ",N,";"
    WRITE(21,*) '        matchTolerance  0.0001;'
    WRITE(21,232) "        myProcNo        ",PRO1,";"
    WRITE(21,232) "        neighbProcNo    ",PRO2,";"
    WRITE(21,*) '    }'
    WRITE(41,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(41,*) '    {'
    WRITE(41,*) '        type            processor;'
    WRITE(41,*) '        value           nonuniform List<scalar> '
    N=CNI*NK
    WRITE(41,'(I0)') N
    WRITE(41,*) '('
    DO KP=1,NK
    DO IP=1,CNI
    WRITE(41,'(E21.14)') P_DOWN_CY(IP,KP)    
    ENDDO
    ENDDO
    WRITE(41,*) ')'
    WRITE(41,*) ';'
    WRITE(41,*) '    }'
    WRITE(51,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(51,*) '    {'
    WRITE(51,*) '        type            processor;'
    WRITE(51,*) '        value           nonuniform List<vector> '
    N=CNI*NK
    WRITE(51,'(I0)') N
    WRITE(51,*) '('
    DO KP=1,NK
    DO IP=1,CNI
    UU1=U_DOWN_CY(1,IP,KP)
    UU2=U_DOWN_CY(2,IP,KP)
    UU3=U_DOWN_CY(3,IP,KP)
    WRITE(51,230) '( ',UU1,' ',UU2,' ',UU3,' )'  
    ENDDO
    ENDDO
    WRITE(51,*) ')'
    WRITE(51,*) ';'
    WRITE(51,*) '    }'
    WRITE(61,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(61,*) '    {'
    WRITE(61,*) '        type            processor;'
    WRITE(61,*) "        value           uniform 0;"
    WRITE(61,*) '    }'
    
        !RIGHT FACE--PROCESSOR_MI1_TO_MI2  !CASE CYLINDER 3-3 
    PRO2=PNI*(MJ-1)+MI-2
    WRITE(21,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(21,*) '    {'
    WRITE(21,*) '        type            processor;'
    WRITE(21,*) '        inGroups        List<word> 1(processor);'
    N=CNJ*NK
    WRITE(21,232) "        nFaces          ",N,";"
    N=F_INTERNAL_NUM+CNJ*NK+CNI*NK+CNI*NK
    WRITE(21,232) "        startFace       ",N,";"
    WRITE(21,*) '        matchTolerance  0.0001;'
    WRITE(21,232) "        myProcNo        ",PRO1,";"
    WRITE(21,232) "        neighbProcNo    ",PRO2,";"
    WRITE(21,*) '    }'
    WRITE(41,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(41,*) '    {'
    WRITE(41,*) '        type            processor;'
    WRITE(41,*) '        value           nonuniform List<scalar> '
    N=CNJ*NK
    WRITE(41,'(I0)') N
    WRITE(41,*) '('
    DO KP=1,NK
    DO JP=1,CNJ
    WRITE(41,'(E21.14)') P_RIGHT_CY(JP,KP)    
    ENDDO
    ENDDO
    WRITE(41,*) ')'
    WRITE(41,*) ';'
    WRITE(41,*) '    }'
    WRITE(51,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(51,*) '    {'
    WRITE(51,*) '        type            processor;'
    WRITE(51,*) '        value           nonuniform List<vector> '
    N=CNJ*NK
    WRITE(51,'(I0)') N
    WRITE(51,*) '('
    DO KP=1,NK
    DO JP=1,CNJ
    UU1=U_RIGHT_CY(1,JP,KP)
    UU2=U_RIGHT_CY(2,JP,KP)
    UU3=U_RIGHT_CY(3,JP,KP)
    WRITE(51,230) '( ',UU1,' ',UU2,' ',UU3,' )'  
    ENDDO
    ENDDO
    WRITE(51,*) ')'
    WRITE(51,*) ';'
    WRITE(51,*) '    }'
    WRITE(61,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(61,*) '    {'
    WRITE(61,*) '        type            processor;'
    WRITE(61,*) "        value           uniform 0;"
    WRITE(61,*) '    }'
    
    
    
    
    
    
    
    
    
    
    
    
    ELSE
    MJ1=( (COORD(1)-PNJ)-MOD(COORD(1)-PNJ,PNL) )/PNL+1
    MJ2=MJ-PNJ-(MJ1-1)*PNL
    IF (MJ1.EQ.1) THEN
    IF (MJ2.EQ.1) THEN
    !CASE RECTANGLE 3-1-1
        
        !LEFT FACE--HULL  !CASE RECTANGLE 3-1-1 
    IP=CNI
    I=IP+CNI*(MI-1)
    DO KP=1,CNL
    K=KP+CNL*(MJ2-1)
    DO JP=1,CNL
    J=JP+CNL*(MJ1-1)
    F11=(CNL+1)*(CNI+1)*(KP-1)+(CNI+1)*(JP-1)+IP-1
    F12=F11+1
    F13=F12+CNI+1
    F14=F13-1
    F21=F11+(CNL+1)*(CNI+1)
    F22=F21+1
    F23=F22+CNI+1
    F24=F23-1
    CALL CAL_FIJ_RECT(F11,F12,F13,F14,F21,F22,F23,F24)
    WRITE(2,235) '4(',F12,' ',F13,' ',F23,' ',F22,')'
    FAD1=F_INTERNAL_NUM_ALL+NL*NL+NJ*NK+NI*NK+NJ*NK+NL*(K-1)+J
    WRITE(3,231) FAD1
    OWN=CELL_ID_RECT(IP,JP,KP)
    WRITE(4,231) OWN
    ENDDO   
    ENDDO
    
        !FRONT FACE--PROCESSOR_MI1_TO_MI2  !CASE RECTANGLE 3-1-1
    KP=1
    K=KP+CNL*(MJ2-1)
    DO JP=1,CNL
    J=JP+CNL*(MJ1-1)
    DO IP=1,CNI
    I=IP+CNI*(MI-1)
    F11=(CNL+1)*(CNI+1)*(KP-1)+(CNI+1)*(JP-1)+IP-1
    F12=F11+1
    F13=F12+CNI+1
    F14=F13-1
    F21=F11+(CNL+1)*(CNI+1)
    F22=F21+1
    F23=F22+CNI+1
    F24=F23-1
    CALL CAL_FIJ_RECT(F11,F12,F13,F14,F21,F22,F23,F24)
    WRITE(2,235) '4(',F11,' ',F14,' ',F13,' ',F12,')'
    FAD1=F_FRONT_RECT(IP,JP)
    WRITE(3,231) FAD1
    OWN=CELL_ID_RECT(IP,JP,KP)
    WRITE(4,231) OWN
    ENDDO   
    ENDDO
    
        !DOWN FACE--PROCESSOR_MI1_TO_MI2  !CASE RECTANGLE 3-1-1
    JP=1
    J=JP+CNL*(MJ1-1)
    DO KP=CNL,1,-1
    K=KP+CNL*(MJ2-1)
    DO IP=1,CNI
    I=IP+CNI*(MI-1)
    F11=(CNL+1)*(CNI+1)*(KP-1)+(CNI+1)*(JP-1)+IP-1
    F12=F11+1
    F13=F12+CNI+1
    F14=F13-1
    F21=F11+(CNL+1)*(CNI+1)
    F22=F21+1
    F23=F22+CNI+1
    F24=F23-1
    CALL CAL_FIJ_RECT(F11,F12,F13,F14,F21,F22,F23,F24)
    !WRITE(2,235) '4(',F11,' ',F12,' ',F22,' ',F21,')'
    WRITE(2,235) '4(',F21,' ',F11,' ',F12,' ',F22,')'
    FAD1=F_DOWN_RECT(IP,KP)
    WRITE(3,231) FAD1
    OWN=CELL_ID_RECT(IP,JP,KP)
    WRITE(4,231) OWN
    ENDDO   
    ENDDO
    
        !RIGHT FACE--PROCESSOR_MI1_TO_MI2  !CASE RECTANGLE 3-1-1 
    IP=1
    I=IP+CNI*(MI-1)
    DO KP=1,CNL
    K=KP+CNL*(MJ2-1)
    DO JP=1,CNL
    J=JP+CNL*(MJ1-1)
    F11=(CNL+1)*(CNI+1)*(KP-1)+(CNI+1)*(JP-1)+IP-1
    F12=F11+1
    F13=F12+CNI+1
    F14=F13-1
    F21=F11+(CNL+1)*(CNI+1)
    F22=F21+1
    F23=F22+CNI+1
    F24=F23-1
    CALL CAL_FIJ_RECT(F11,F12,F13,F14,F21,F22,F23,F24)
    WRITE(2,235) '4(',F11,' ',F21,' ',F24,' ',F14,')'
    FAD1=F_RIGHT_RECT(JP,KP)
    WRITE(3,231) FAD1
    OWN=CELL_ID_RECT(IP,JP,KP)
    WRITE(4,231) OWN
    ENDDO   
    ENDDO
    
        !UP FACE--PROCESSOR_MI1_TO_MI2  !CASE RECTANGLE 3-1-1
    JP=CNL
    J=JP+CNL*(MJ1-1)
    DO KP=1,CNL
    K=KP+CNL*(MJ2-1)
    DO IP=1,CNI
    I=IP+CNI*(MI-1)
    F11=(CNL+1)*(CNI+1)*(KP-1)+(CNI+1)*(JP-1)+IP-1
    F12=F11+1
    F13=F12+CNI+1
    F14=F13-1
    F21=F11+(CNL+1)*(CNI+1)
    F22=F21+1
    F23=F22+CNI+1
    F24=F23-1
    CALL CAL_FIJ_RECT(F11,F12,F13,F14,F21,F22,F23,F24)
    WRITE(2,235) '4(',F14,' ',F24,' ',F23,' ',F13,')'
    FAD1=F_UP_RECT(IP,KP)
    WRITE(3,231) FAD1
    OWN=CELL_ID_RECT(IP,JP,KP)
    WRITE(4,231) OWN
    ENDDO   
    ENDDO
    
        !BACK FACE--PROCESSOR_MI1_TO_MI2  !CASE RECTANGLE 3-1-1
    KP=CNL
    K=KP+CNL*(MJ2-1)
    DO JP=1,CNL
    J=JP+CNL*(MJ1-1)
    DO IP=1,CNI
    I=IP+CNI*(MI-1)
    F11=(CNL+1)*(CNI+1)*(KP-1)+(CNI+1)*(JP-1)+IP-1
    F12=F11+1
    F13=F12+CNI+1
    F14=F13-1
    F21=F11+(CNL+1)*(CNI+1)
    F22=F21+1
    F23=F22+CNI+1
    F24=F23-1
    CALL CAL_FIJ_RECT(F11,F12,F13,F14,F21,F22,F23,F24)
    WRITE(2,235) '4(',F21,' ',F22,' ',F23,' ',F24,')'
    FAD1=F_BACK_RECT(IP,JP)
    WRITE(3,231) FAD1
    OWN=CELL_ID_RECT(IP,JP,KP)
    WRITE(4,231) OWN
    ENDDO   
    ENDDO
    
    N=4
    WRITE(31,231) N+4
    WRITE(31,*) '('  
    WRITE(31,*) '0'
    WRITE(31,*) '1'
    WRITE(31,*) '2' 
    WRITE(31,*) '3' 
    DO I=1,N
    WRITE(31,*) '-1' 
    ENDDO
    WRITE(31,*) ')'
    
    WRITE(21,231) N+4
    WRITE(21,*) "("    
    WRITE(21,*) "    INLET1"
    WRITE(21,*) "    {"    
    WRITE(21,*) "        type            patch;"
    N=0
    WRITE(21,232) "        nFaces          ",N,";"
    N=F_INTERNAL_NUM
    WRITE(21,232) "        startFace       ",N,";"
    WRITE(21,*) "    }"
    WRITE(21,*) "    INLET2"
    WRITE(21,*) "    {"    
    WRITE(21,*) "        type            patch;"
    N=0
    WRITE(21,232) "        nFaces          ",N,";"
    N=F_INTERNAL_NUM
    WRITE(21,232) "        startFace       ",N,";"
    WRITE(21,*) "    }"
    WRITE(21,*) "    OUTLET"
    WRITE(21,*) "    {"
    WRITE(21,*) "        type            patch;"
    N=0
    WRITE(21,232) "        nFaces          ",N,";"
    N=F_INTERNAL_NUM
    WRITE(21,232) "        startFace       ",N,";"
    WRITE(21,*) "    }"
    WRITE(21,*) "    HULL"
    WRITE(21,*) "    {"
    WRITE(21,*) "        type            wall;"
    WRITE(21,*) "        inGroups        1(wall);"
    N=CNL*CNL
    WRITE(21,232) "        nFaces          ",N,";"
    N=F_INTERNAL_NUM
    WRITE(21,232) "        startFace       ",N,";"
    WRITE(21,*) "    }"
    
        !FRONT/DOWN FACE--PROCESSOR_MI1_TO_MI2  !CASE RECTANGLE 3-1-1
    PRO1=PNI*PNJ+PNL*PNI*(MJ2-1)+PNI*(MJ1-1)+MI-1
    MJ=1
    PRO2=PNI*(MJ-1)+MI-1
    WRITE(21,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(21,*) '    {'
    WRITE(21,*) '        type            processor;'
    WRITE(21,*) '        inGroups        List<word> 1(processor);'
    N=CNI*CNL*2
    WRITE(21,232) "        nFaces          ",N,";"
    N=F_INTERNAL_NUM+CNL*CNL
    WRITE(21,232) "        startFace       ",N,";"
    WRITE(21,*) '        matchTolerance  0.0001;'
    WRITE(21,232) "        myProcNo        ",PRO1,";"
    WRITE(21,232) "        neighbProcNo    ",PRO2,";"
    WRITE(21,*) '    }'
    WRITE(41,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(41,*) '    {'
    WRITE(41,*) '        type            processor;'
    WRITE(41,*) '        value           nonuniform List<scalar> '
    N=CNI*CNL*2
    WRITE(41,'(I0)') N
    WRITE(41,*) '('
    DO JP=1,CNL
    DO IP=1,CNI
    WRITE(41,'(E21.14)') P_FRONT_RECT(IP,JP)    
    ENDDO
    ENDDO
    DO KP=1,CNL
    DO IP=1,CNI
    WRITE(41,'(E21.14)') P_DOWN_RECT(IP,KP)    
    ENDDO
    ENDDO
    WRITE(41,*) ')'
    WRITE(41,*) ';'
    WRITE(41,*) '    }'
    WRITE(51,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(51,*) '    {'
    WRITE(51,*) '        type            processor;'
    WRITE(51,*) '        value           nonuniform List<vector> '
    N=CNI*CNL*2
    WRITE(51,'(I0)') N
    WRITE(51,*) '('
    DO JP=1,CNL
    DO IP=1,CNI
    UU1=U_FRONT_RECT(1,IP,JP)
    UU2=U_FRONT_RECT(2,IP,JP)
    UU3=U_FRONT_RECT(3,IP,JP)
    WRITE(51,230) '( ',UU1,' ',UU2,' ',UU3,' )'  
    ENDDO
    ENDDO
    DO KP=1,CNL
    DO IP=1,CNI
    UU1=U_DOWN_RECT(1,IP,KP)
    UU2=U_DOWN_RECT(2,IP,KP)
    UU3=U_DOWN_RECT(3,IP,KP)
    WRITE(51,230) '( ',UU1,' ',UU2,' ',UU3,' )'  
    ENDDO
    ENDDO
    WRITE(51,*) ')'
    WRITE(51,*) ';'
    WRITE(51,*) '    }'
    WRITE(61,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(61,*) '    {'
    WRITE(61,*) '        type            processor;'
    WRITE(61,*) "        value           uniform 0;"
    WRITE(61,*) '    }'
    
        !RIGHT FACE--PROCESSOR_MI1_TO_MI2  !CASE RECTANGLE 3-1-1 
    PRO2=PNI*PNJ+PNL*PNI*(MJ2-1)+PNI*(MJ1-1)+MI-2
    WRITE(21,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(21,*) '    {'
    WRITE(21,*) '        type            processor;'
    WRITE(21,*) '        inGroups        List<word> 1(processor);'
    N=CNL*CNL
    WRITE(21,232) "        nFaces          ",N,";"
    N=F_INTERNAL_NUM+CNL*CNL+CNI*CNL*2
    WRITE(21,232) "        startFace       ",N,";"
    WRITE(21,*) '        matchTolerance  0.0001;'
    WRITE(21,232) "        myProcNo        ",PRO1,";"
    WRITE(21,232) "        neighbProcNo    ",PRO2,";"
    WRITE(21,*) '    }'
    WRITE(41,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(41,*) '    {'
    WRITE(41,*) '        type            processor;'
    WRITE(41,*) '        value           nonuniform List<scalar> '
    N=CNL*CNL
    WRITE(41,'(I0)') N
    WRITE(41,*) '('
    DO KP=1,CNL
    DO JP=1,CNL
    WRITE(41,'(E21.14)') P_RIGHT_RECT(JP,KP)    
    ENDDO
    ENDDO
    WRITE(41,*) ')'
    WRITE(41,*) ';'
    WRITE(41,*) '    }'
    WRITE(51,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(51,*) '    {'
    WRITE(51,*) '        type            processor;'
    WRITE(51,*) '        value           nonuniform List<vector> '
    N=CNL*CNL
    WRITE(51,'(I0)') N
    WRITE(51,*) '('
    DO KP=1,CNL
    DO JP=1,CNL
    UU1=U_RIGHT_RECT(1,JP,KP)
    UU2=U_RIGHT_RECT(2,JP,KP)
    UU3=U_RIGHT_RECT(3,JP,KP)
    WRITE(51,230) '( ',UU1,' ',UU2,' ',UU3,' )'  
    ENDDO
    ENDDO
    WRITE(51,*) ')'
    WRITE(51,*) ';'
    WRITE(51,*) '    }'
    WRITE(61,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(61,*) '    {'
    WRITE(61,*) '        type            processor;'
    WRITE(61,*) "        value           uniform 0;"
    WRITE(61,*) '    }'
    
        !UP FACE--PROCESSOR_MI1_TO_MI2  !CASE RECTANGLE 3-1-1
    PRO2=PNI*PNJ+PNL*PNI*(MJ2-1)+PNI*MJ1+MI-1
    WRITE(21,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(21,*) '    {'
    WRITE(21,*) '        type            processor;'
    WRITE(21,*) '        inGroups        List<word> 1(processor);'
    N=CNI*CNL
    WRITE(21,232) "        nFaces          ",N,";"
    N=F_INTERNAL_NUM+CNL*CNL+CNI*CNL*2+CNL*CNL
    WRITE(21,232) "        startFace       ",N,";"
    WRITE(21,*) '        matchTolerance  0.0001;'
    WRITE(21,232) "        myProcNo        ",PRO1,";"
    WRITE(21,232) "        neighbProcNo    ",PRO2,";"
    WRITE(21,*) '    }'
    WRITE(41,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(41,*) '    {'
    WRITE(41,*) '        type            processor;'
    WRITE(41,*) '        value           nonuniform List<scalar> '
    N=CNI*CNL
    WRITE(41,'(I0)') N
    WRITE(41,*) '('
    DO KP=1,CNL
    DO IP=1,CNI
    WRITE(41,'(E21.14)') P_UP_RECT(IP,KP)    
    ENDDO
    ENDDO
    WRITE(41,*) ')'
    WRITE(41,*) ';'
    WRITE(41,*) '    }'
    WRITE(51,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(51,*) '    {'
    WRITE(51,*) '        type            processor;'
    WRITE(51,*) '        value           nonuniform List<vector> '
    N=CNI*CNL
    WRITE(51,'(I0)') N
    WRITE(51,*) '('
    DO KP=1,CNL
    DO IP=1,CNI
    UU1=U_UP_RECT(1,IP,KP)
    UU2=U_UP_RECT(2,IP,KP)
    UU3=U_UP_RECT(3,IP,KP)
    WRITE(51,230) '( ',UU1,' ',UU2,' ',UU3,' )'  
    ENDDO
    ENDDO
    WRITE(51,*) ')'
    WRITE(51,*) ';'
    WRITE(51,*) '    }'
    WRITE(61,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(61,*) '    {'
    WRITE(61,*) '        type            processor;'
    WRITE(61,*) "        value           uniform 0;"
    WRITE(61,*) '    }'
    
        !BACK FACE--PROCESSOR_MI1_TO_MI2  !CASE RECTANGLE 3-1-1
    PRO2=PNI*PNJ+PNL*PNI*MJ2+PNI*(MJ1-1)+MI-1
    WRITE(21,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(21,*) '    {'
    WRITE(21,*) '        type            processor;'
    WRITE(21,*) '        inGroups        List<word> 1(processor);'
    N=CNI*CNL
    WRITE(21,232) "        nFaces          ",N,";"
    N=F_INTERNAL_NUM+CNL*CNL+CNI*CNL*2+CNL*CNL+CNI*CNL
    WRITE(21,232) "        startFace       ",N,";"
    WRITE(21,*) '        matchTolerance  0.0001;'
    WRITE(21,232) "        myProcNo        ",PRO1,";"
    WRITE(21,232) "        neighbProcNo    ",PRO2,";"
    WRITE(21,*) '    }'
    WRITE(41,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(41,*) '    {'
    WRITE(41,*) '        type            processor;'
    WRITE(41,*) '        value           nonuniform List<scalar> '
    N=CNI*CNL
    WRITE(41,'(I0)') N
    WRITE(41,*) '('
    DO JP=1,CNL
    DO IP=1,CNI
    WRITE(41,'(E21.14)') P_BACK_RECT(IP,JP)    
    ENDDO
    ENDDO
    WRITE(41,*) ')'
    WRITE(41,*) ';'
    WRITE(41,*) '    }'
    WRITE(51,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(51,*) '    {'
    WRITE(51,*) '        type            processor;'
    WRITE(51,*) '        value           nonuniform List<vector> '
    N=CNI*CNL
    WRITE(51,'(I0)') N
    WRITE(51,*) '('
    DO JP=1,CNL
    DO IP=1,CNI
    UU1=U_BACK_RECT(1,IP,JP)
    UU2=U_BACK_RECT(2,IP,JP)
    UU3=U_BACK_RECT(3,IP,JP)
    WRITE(51,230) '( ',UU1,' ',UU2,' ',UU3,' )'  
    ENDDO
    ENDDO
    WRITE(51,*) ')'
    WRITE(51,*) ';'
    WRITE(51,*) '    }'
    WRITE(61,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(61,*) '    {'
    WRITE(61,*) '        type            processor;'
    WRITE(61,*) "        value           uniform 0;"
    WRITE(61,*) '    }'
    
    
    
    
    
    
    
    
    
    ELSEIF (MJ2.EQ.PNL) THEN
    !CASE RECTANGLE 3-1-3
        
        !LEFT FACE--HULL  !CASE RECTANGLE 3-1-3 
    IP=CNI
    I=IP+CNI*(MI-1)
    DO KP=1,CNL
    K=KP+CNL*(MJ2-1)
    DO JP=1,CNL
    J=JP+CNL*(MJ1-1)
    F11=(CNL+1)*(CNI+1)*(KP-1)+(CNI+1)*(JP-1)+IP-1
    F12=F11+1
    F13=F12+CNI+1
    F14=F13-1
    F21=F11+(CNL+1)*(CNI+1)
    F22=F21+1
    F23=F22+CNI+1
    F24=F23-1
    CALL CAL_FIJ_RECT(F11,F12,F13,F14,F21,F22,F23,F24)
    WRITE(2,235) '4(',F12,' ',F13,' ',F23,' ',F22,')'
    FAD1=F_INTERNAL_NUM_ALL+NL*NL+NJ*NK+NI*NK+NJ*NK+NL*(K-1)+J
    WRITE(3,231) FAD1
    OWN=CELL_ID_RECT(IP,JP,KP)
    WRITE(4,231) OWN
    ENDDO   
    ENDDO
    
        !BACK FACE--PROCESSOR_MI1_TO_MI2  !CASE RECTANGLE 3-1-3
    KP=CNL
    K=KP+CNL*(MJ2-1)
    DO JP=CNL,1,-1
    J=JP+CNL*(MJ1-1)
    DO IP=1,CNI
    I=IP+CNI*(MI-1)
    F11=(CNL+1)*(CNI+1)*(KP-1)+(CNI+1)*(JP-1)+IP-1
    F12=F11+1
    F13=F12+CNI+1
    F14=F13-1
    F21=F11+(CNL+1)*(CNI+1)
    F22=F21+1
    F23=F22+CNI+1
    F24=F23-1
    CALL CAL_FIJ_RECT(F11,F12,F13,F14,F21,F22,F23,F24)
    !WRITE(2,235) '4(',F21,' ',F22,' ',F23,' ',F24,')'
    WRITE(2,235) '4(',F24,' ',F21,' ',F22,' ',F23,')'
    FAD1=F_BACK_RECT(IP,JP)
    WRITE(3,231) FAD1
    OWN=CELL_ID_RECT(IP,JP,KP)
    WRITE(4,231) OWN
    ENDDO   
    ENDDO
    
        !DOWN FACE--PROCESSOR_MI1_TO_MI2  !CASE RECTANGLE 3-1-3
    JP=1
    J=JP+CNL*(MJ1-1)
    DO KP=CNL,1,-1
    K=KP+CNL*(MJ2-1)
    DO IP=1,CNI
    I=IP+CNI*(MI-1)
    F11=(CNL+1)*(CNI+1)*(KP-1)+(CNI+1)*(JP-1)+IP-1
    F12=F11+1
    F13=F12+CNI+1
    F14=F13-1
    F21=F11+(CNL+1)*(CNI+1)
    F22=F21+1
    F23=F22+CNI+1
    F24=F23-1
    CALL CAL_FIJ_RECT(F11,F12,F13,F14,F21,F22,F23,F24)
    !WRITE(2,235) '4(',F11,' ',F12,' ',F22,' ',F21,')'
    WRITE(2,235) '4(',F21,' ',F11,' ',F12,' ',F22,')'
    FAD1=F_DOWN_RECT(IP,KP)
    WRITE(3,231) FAD1
    OWN=CELL_ID_RECT(IP,JP,KP)
    WRITE(4,231) OWN
    ENDDO   
    ENDDO
    
        !FRONT FACE--PROCESSOR_MI1_TO_MI2  !CASE RECTANGLE 3-1-3
    KP=1
    K=KP+CNL*(MJ2-1)
    DO JP=1,CNL
    J=JP+CNL*(MJ1-1)
    DO IP=1,CNI
    I=IP+CNI*(MI-1)
    F11=(CNL+1)*(CNI+1)*(KP-1)+(CNI+1)*(JP-1)+IP-1
    F12=F11+1
    F13=F12+CNI+1
    F14=F13-1
    F21=F11+(CNL+1)*(CNI+1)
    F22=F21+1
    F23=F22+CNI+1
    F24=F23-1
    CALL CAL_FIJ_RECT(F11,F12,F13,F14,F21,F22,F23,F24)
    WRITE(2,235) '4(',F11,' ',F14,' ',F13,' ',F12,')'
    FAD1=F_FRONT_RECT(IP,JP)
    WRITE(3,231) FAD1
    OWN=CELL_ID_RECT(IP,JP,KP)
    WRITE(4,231) OWN
    ENDDO   
    ENDDO
    
        !RIGHT FACE--PROCESSOR_MI1_TO_MI2  !CASE RECTANGLE 3-1-3 
    IP=1
    I=IP+CNI*(MI-1)
    DO KP=1,CNL
    K=KP+CNL*(MJ2-1)
    DO JP=1,CNL
    J=JP+CNL*(MJ1-1)
    F11=(CNL+1)*(CNI+1)*(KP-1)+(CNI+1)*(JP-1)+IP-1
    F12=F11+1
    F13=F12+CNI+1
    F14=F13-1
    F21=F11+(CNL+1)*(CNI+1)
    F22=F21+1
    F23=F22+CNI+1
    F24=F23-1
    CALL CAL_FIJ_RECT(F11,F12,F13,F14,F21,F22,F23,F24)
    WRITE(2,235) '4(',F11,' ',F21,' ',F24,' ',F14,')'
    FAD1=F_RIGHT_RECT(JP,KP)
    WRITE(3,231) FAD1
    OWN=CELL_ID_RECT(IP,JP,KP)
    WRITE(4,231) OWN
    ENDDO   
    ENDDO
    
        !UP FACE--PROCESSOR_MI1_TO_MI2  !CASE RECTANGLE 3-1-3
    JP=CNL
    J=JP+CNL*(MJ1-1)
    DO KP=1,CNL
    K=KP+CNL*(MJ2-1)
    DO IP=1,CNI
    I=IP+CNI*(MI-1)
    F11=(CNL+1)*(CNI+1)*(KP-1)+(CNI+1)*(JP-1)+IP-1
    F12=F11+1
    F13=F12+CNI+1
    F14=F13-1
    F21=F11+(CNL+1)*(CNI+1)
    F22=F21+1
    F23=F22+CNI+1
    F24=F23-1
    CALL CAL_FIJ_RECT(F11,F12,F13,F14,F21,F22,F23,F24)
    WRITE(2,235) '4(',F14,' ',F24,' ',F23,' ',F13,')'
    FAD1=F_UP_RECT(IP,KP)
    WRITE(3,231) FAD1
    OWN=CELL_ID_RECT(IP,JP,KP)
    WRITE(4,231) OWN
    ENDDO   
    ENDDO
    
    N=4
    WRITE(31,231) N+4
    WRITE(31,*) '('  
    WRITE(31,*) '0'
    WRITE(31,*) '1'
    WRITE(31,*) '2' 
    WRITE(31,*) '3' 
    DO I=1,N
    WRITE(31,*) '-1' 
    ENDDO
    WRITE(31,*) ')'
    
    WRITE(21,231) N+4
    WRITE(21,*) "("    
    WRITE(21,*) "    INLET1"
    WRITE(21,*) "    {"    
    WRITE(21,*) "        type            patch;"
    N=0
    WRITE(21,232) "        nFaces          ",N,";"
    N=F_INTERNAL_NUM
    WRITE(21,232) "        startFace       ",N,";"
    WRITE(21,*) "    }"
    WRITE(21,*) "    INLET2"
    WRITE(21,*) "    {"    
    WRITE(21,*) "        type            patch;"
    N=0
    WRITE(21,232) "        nFaces          ",N,";"
    N=F_INTERNAL_NUM
    WRITE(21,232) "        startFace       ",N,";"
    WRITE(21,*) "    }"
    WRITE(21,*) "    OUTLET"
    WRITE(21,*) "    {"
    WRITE(21,*) "        type            patch;"
    N=0
    WRITE(21,232) "        nFaces          ",N,";"
    N=F_INTERNAL_NUM
    WRITE(21,232) "        startFace       ",N,";"
    WRITE(21,*) "    }"
    WRITE(21,*) "    HULL"
    WRITE(21,*) "    {"
    WRITE(21,*) "        type            wall;"
    WRITE(21,*) "        inGroups        1(wall);"
    N=CNL*CNL
    WRITE(21,232) "        nFaces          ",N,";"
    N=F_INTERNAL_NUM
    WRITE(21,232) "        startFace       ",N,";"
    WRITE(21,*) "    }"
    
        !BACK/DOWN FACE--PROCESSOR_MI1_TO_MI2  !CASE RECTANGLE 3-1-3
    PRO1=PNI*PNJ+PNL*PNI*(MJ2-1)+PNI*(MJ1-1)+MI-1
    MJ=1
    PRO2=PNI*(MJ-1)+MI-1
    WRITE(21,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(21,*) '    {'
    WRITE(21,*) '        type            processor;'
    WRITE(21,*) '        inGroups        List<word> 1(processor);'
    N=CNI*CNL*2
    WRITE(21,232) "        nFaces          ",N,";"
    N=F_INTERNAL_NUM+CNL*CNL
    WRITE(21,232) "        startFace       ",N,";"
    WRITE(21,*) '        matchTolerance  0.0001;'
    WRITE(21,232) "        myProcNo        ",PRO1,";"
    WRITE(21,232) "        neighbProcNo    ",PRO2,";"
    WRITE(21,*) '    }'
    WRITE(41,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(41,*) '    {'
    WRITE(41,*) '        type            processor;'
    WRITE(41,*) '        value           nonuniform List<scalar> '
    N=CNI*CNL*2
    WRITE(41,'(I0)') N
    WRITE(41,*) '('
    DO JP=1,CNL
    DO IP=1,CNI
    WRITE(41,'(E21.14)') P_BACK_RECT(IP,JP)    
    ENDDO
    ENDDO
    DO KP=1,CNL
    DO IP=1,CNI
    WRITE(41,'(E21.14)') P_DOWN_RECT(IP,KP)    
    ENDDO
    ENDDO
    WRITE(41,*) ')'
    WRITE(41,*) ';'
    WRITE(41,*) '    }'
    WRITE(51,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(51,*) '    {'
    WRITE(51,*) '        type            processor;'
    WRITE(51,*) '        value           nonuniform List<vector> '
    N=CNI*CNL*2
    WRITE(51,'(I0)') N
    WRITE(51,*) '('
    DO JP=1,CNL
    DO IP=1,CNI
    UU1=U_BACK_RECT(1,IP,JP)
    UU2=U_BACK_RECT(2,IP,JP)
    UU3=U_BACK_RECT(3,IP,JP)
    WRITE(51,230) '( ',UU1,' ',UU2,' ',UU3,' )'  
    ENDDO
    ENDDO
    DO KP=1,CNL
    DO IP=1,CNI
    UU1=U_DOWN_RECT(1,IP,KP)
    UU2=U_DOWN_RECT(2,IP,KP)
    UU3=U_DOWN_RECT(3,IP,KP)
    WRITE(51,230) '( ',UU1,' ',UU2,' ',UU3,' )'  
    ENDDO
    ENDDO
    WRITE(51,*) ')'
    WRITE(51,*) ';'
    WRITE(51,*) '    }'
    WRITE(61,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(61,*) '    {'
    WRITE(61,*) '        type            processor;'
    WRITE(61,*) "        value           uniform 0;"
    WRITE(61,*) '    }'
    
        !FRONT FACE--PROCESSOR_MI1_TO_MI2  !CASE RECTANGLE 3-1-3
    PRO2=PNI*PNJ+PNL*PNI*(MJ2-2)+PNI*(MJ1-1)+MI-1
    WRITE(21,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(21,*) '    {'
    WRITE(21,*) '        type            processor;'
    WRITE(21,*) '        inGroups        List<word> 1(processor);'
    N=CNI*CNL
    WRITE(21,232) "        nFaces          ",N,";"
    N=F_INTERNAL_NUM+CNL*CNL+CNI*CNL*2
    WRITE(21,232) "        startFace       ",N,";"
    WRITE(21,*) '        matchTolerance  0.0001;'
    WRITE(21,232) "        myProcNo        ",PRO1,";"
    WRITE(21,232) "        neighbProcNo    ",PRO2,";"
    WRITE(21,*) '    }'
    WRITE(41,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(41,*) '    {'
    WRITE(41,*) '        type            processor;'
    WRITE(41,*) '        value           nonuniform List<scalar> '
    N=CNI*CNL
    WRITE(41,'(I0)') N
    WRITE(41,*) '('
    DO JP=1,CNL
    DO IP=1,CNI
    WRITE(41,'(E21.14)') P_FRONT_RECT(IP,JP)    
    ENDDO
    ENDDO
    WRITE(41,*) ')'
    WRITE(41,*) ';'
    WRITE(41,*) '    }'
    WRITE(51,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(51,*) '    {'
    WRITE(51,*) '        type            processor;'
    WRITE(51,*) '        value           nonuniform List<vector> '
    N=CNI*CNL
    WRITE(51,'(I0)') N
    WRITE(51,*) '('
    DO JP=1,CNL
    DO IP=1,CNI
    UU1=U_FRONT_RECT(1,IP,JP)
    UU2=U_FRONT_RECT(2,IP,JP)
    UU3=U_FRONT_RECT(3,IP,JP)
    WRITE(51,230) '( ',UU1,' ',UU2,' ',UU3,' )'  
    ENDDO
    ENDDO
    WRITE(51,*) ')'
    WRITE(51,*) ';'
    WRITE(51,*) '    }'
    WRITE(61,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(61,*) '    {'
    WRITE(61,*) '        type            processor;'
    WRITE(61,*) "        value           uniform 0;"
    WRITE(61,*) '    }'
    
        !RIGHT FACE--PROCESSOR_MI1_TO_MI2  !CASE RECTANGLE 3-1-3
    PRO2=PNI*PNJ+PNL*PNI*(MJ2-1)+PNI*(MJ1-1)+MI-2
    WRITE(21,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(21,*) '    {'
    WRITE(21,*) '        type            processor;'
    WRITE(21,*) '        inGroups        List<word> 1(processor);'
    N=CNL*CNL
    WRITE(21,232) "        nFaces          ",N,";"
    N=F_INTERNAL_NUM+CNL*CNL+CNI*CNL*2+CNI*CNL
    WRITE(21,232) "        startFace       ",N,";"
    WRITE(21,*) '        matchTolerance  0.0001;'
    WRITE(21,232) "        myProcNo        ",PRO1,";"
    WRITE(21,232) "        neighbProcNo    ",PRO2,";"
    WRITE(21,*) '    }'
    WRITE(41,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(41,*) '    {'
    WRITE(41,*) '        type            processor;'
    WRITE(41,*) '        value           nonuniform List<scalar> '
    N=CNL*CNL
    WRITE(41,'(I0)') N
    WRITE(41,*) '('
    DO KP=1,CNL
    DO JP=1,CNL
    WRITE(41,'(E21.14)') P_RIGHT_RECT(JP,KP)    
    ENDDO
    ENDDO
    WRITE(41,*) ')'
    WRITE(41,*) ';'
    WRITE(41,*) '    }'
    WRITE(51,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(51,*) '    {'
    WRITE(51,*) '        type            processor;'
    WRITE(51,*) '        value           nonuniform List<vector> '
    N=CNL*CNL
    WRITE(51,'(I0)') N
    WRITE(51,*) '('
    DO KP=1,CNL
    DO JP=1,CNL
    UU1=U_RIGHT_RECT(1,JP,KP)
    UU2=U_RIGHT_RECT(2,JP,KP)
    UU3=U_RIGHT_RECT(3,JP,KP)
    WRITE(51,230) '( ',UU1,' ',UU2,' ',UU3,' )'  
    ENDDO
    ENDDO
    WRITE(51,*) ')'
    WRITE(51,*) ';'
    WRITE(51,*) '    }'
    WRITE(61,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(61,*) '    {'
    WRITE(61,*) '        type            processor;'
    WRITE(61,*) "        value           uniform 0;"
    WRITE(61,*) '    }'
    
        !UP FACE--PROCESSOR_MI1_TO_MI2  !CASE RECTANGLE 3-1-3
    PRO2=PNI*PNJ+PNL*PNI*(MJ2-1)+PNI*MJ1+MI-1
    WRITE(21,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(21,*) '    {'
    WRITE(21,*) '        type            processor;'
    WRITE(21,*) '        inGroups        List<word> 1(processor);'
    N=CNI*CNL
    WRITE(21,232) "        nFaces          ",N,";"
    N=F_INTERNAL_NUM+CNL*CNL+CNI*CNL*2+CNI*CNL+CNL*CNL
    WRITE(21,232) "        startFace       ",N,";"
    WRITE(21,*) '        matchTolerance  0.0001;'
    WRITE(21,232) "        myProcNo        ",PRO1,";"
    WRITE(21,232) "        neighbProcNo    ",PRO2,";"
    WRITE(21,*) '    }'
    WRITE(41,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(41,*) '    {'
    WRITE(41,*) '        type            processor;'
    WRITE(41,*) '        value           nonuniform List<scalar> '
    N=CNI*CNL
    WRITE(41,'(I0)') N
    WRITE(41,*) '('
    DO KP=1,CNL
    DO IP=1,CNI
    WRITE(41,'(E21.14)') P_UP_RECT(IP,KP)    
    ENDDO
    ENDDO
    WRITE(41,*) ')'
    WRITE(41,*) ';'
    WRITE(41,*) '    }'
    WRITE(51,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(51,*) '    {'
    WRITE(51,*) '        type            processor;'
    WRITE(51,*) '        value           nonuniform List<vector> '
    N=CNI*CNL
    WRITE(51,'(I0)') N
    WRITE(51,*) '('
    DO KP=1,CNL
    DO IP=1,CNI
    UU1=U_UP_RECT(1,IP,KP)
    UU2=U_UP_RECT(2,IP,KP)
    UU3=U_UP_RECT(3,IP,KP)
    WRITE(51,230) '( ',UU1,' ',UU2,' ',UU3,' )'  
    ENDDO
    ENDDO
    WRITE(51,*) ')'
    WRITE(51,*) ';'
    WRITE(51,*) '    }'
    WRITE(61,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(61,*) '    {'
    WRITE(61,*) '        type            processor;'
    WRITE(61,*) "        value           uniform 0;"
    WRITE(61,*) '    }'
    
    
    ELSE
    !CASE RECTANGLE 3-1-2
        
        !LEFT FACE--HULL  !CASE RECTANGLE 3-1-2 
    IP=CNI
    I=IP+CNI*(MI-1)
    DO KP=1,CNL
    K=KP+CNL*(MJ2-1)
    DO JP=1,CNL
    J=JP+CNL*(MJ1-1)
    F11=(CNL+1)*(CNI+1)*(KP-1)+(CNI+1)*(JP-1)+IP-1
    F12=F11+1
    F13=F12+CNI+1
    F14=F13-1
    F21=F11+(CNL+1)*(CNI+1)
    F22=F21+1
    F23=F22+CNI+1
    F24=F23-1
    CALL CAL_FIJ_RECT(F11,F12,F13,F14,F21,F22,F23,F24)
    WRITE(2,235) '4(',F12,' ',F13,' ',F23,' ',F22,')'
    FAD1=F_INTERNAL_NUM_ALL+NL*NL+NJ*NK+NI*NK+NJ*NK+NL*(K-1)+J
    WRITE(3,231) FAD1
    OWN=CELL_ID_RECT(IP,JP,KP)
    WRITE(4,231) OWN
    ENDDO   
    ENDDO
    
        !DOWN FACE--PROCESSOR_MI1_TO_MI2  !CASE RECTANGLE 3-1-2
    JP=1
    J=JP+CNL*(MJ1-1)
    DO KP=CNL,1,-1
    K=KP+CNL*(MJ2-1)
    DO IP=1,CNI
    I=IP+CNI*(MI-1)
    F11=(CNL+1)*(CNI+1)*(KP-1)+(CNI+1)*(JP-1)+IP-1
    F12=F11+1
    F13=F12+CNI+1
    F14=F13-1
    F21=F11+(CNL+1)*(CNI+1)
    F22=F21+1
    F23=F22+CNI+1
    F24=F23-1
    CALL CAL_FIJ_RECT(F11,F12,F13,F14,F21,F22,F23,F24)
    !WRITE(2,235) '4(',F11,' ',F12,' ',F22,' ',F21,')'
    WRITE(2,235) '4(',F21,' ',F11,' ',F12,' ',F22,')'
    FAD1=F_DOWN_RECT(IP,KP)
    WRITE(3,231) FAD1
    OWN=CELL_ID_RECT(IP,JP,KP)
    WRITE(4,231) OWN
    ENDDO   
    ENDDO
    
        !FRONT FACE--PROCESSOR_MI1_TO_MI2  !CASE RECTANGLE 3-1-2
    KP=1
    K=KP+CNL*(MJ2-1)
    DO JP=1,CNL
    J=JP+CNL*(MJ1-1)
    DO IP=1,CNI
    I=IP+CNI*(MI-1)
    F11=(CNL+1)*(CNI+1)*(KP-1)+(CNI+1)*(JP-1)+IP-1
    F12=F11+1
    F13=F12+CNI+1
    F14=F13-1
    F21=F11+(CNL+1)*(CNI+1)
    F22=F21+1
    F23=F22+CNI+1
    F24=F23-1
    CALL CAL_FIJ_RECT(F11,F12,F13,F14,F21,F22,F23,F24)
    WRITE(2,235) '4(',F11,' ',F14,' ',F13,' ',F12,')'
    FAD1=F_FRONT_RECT(IP,JP)
    WRITE(3,231) FAD1
    OWN=CELL_ID_RECT(IP,JP,KP)
    WRITE(4,231) OWN
    ENDDO   
    ENDDO
    
        !RIGHT FACE--PROCESSOR_MI1_TO_MI2  !CASE RECTANGLE 3-1-2 
    IP=1
    I=IP+CNI*(MI-1)
    DO KP=1,CNL
    K=KP+CNL*(MJ2-1)
    DO JP=1,CNL
    J=JP+CNL*(MJ1-1)
    F11=(CNL+1)*(CNI+1)*(KP-1)+(CNI+1)*(JP-1)+IP-1
    F12=F11+1
    F13=F12+CNI+1
    F14=F13-1
    F21=F11+(CNL+1)*(CNI+1)
    F22=F21+1
    F23=F22+CNI+1
    F24=F23-1
    CALL CAL_FIJ_RECT(F11,F12,F13,F14,F21,F22,F23,F24)
    WRITE(2,235) '4(',F11,' ',F21,' ',F24,' ',F14,')'
    FAD1=F_RIGHT_RECT(JP,KP)
    WRITE(3,231) FAD1
    OWN=CELL_ID_RECT(IP,JP,KP)
    WRITE(4,231) OWN
    ENDDO   
    ENDDO
    
        !UP FACE--PROCESSOR_MI1_TO_MI2  !CASE RECTANGLE 3-1-2
    JP=CNL
    J=JP+CNL*(MJ1-1)
    DO KP=1,CNL
    K=KP+CNL*(MJ2-1)
    DO IP=1,CNI
    I=IP+CNI*(MI-1)
    F11=(CNL+1)*(CNI+1)*(KP-1)+(CNI+1)*(JP-1)+IP-1
    F12=F11+1
    F13=F12+CNI+1
    F14=F13-1
    F21=F11+(CNL+1)*(CNI+1)
    F22=F21+1
    F23=F22+CNI+1
    F24=F23-1
    CALL CAL_FIJ_RECT(F11,F12,F13,F14,F21,F22,F23,F24)
    WRITE(2,235) '4(',F14,' ',F24,' ',F23,' ',F13,')'
    FAD1=F_UP_RECT(IP,KP)
    WRITE(3,231) FAD1
    OWN=CELL_ID_RECT(IP,JP,KP)
    WRITE(4,231) OWN
    ENDDO   
    ENDDO
    
        !BACK FACE--PROCESSOR_MI1_TO_MI2  !CASE RECTANGLE 3-1-2
    KP=CNL
    K=KP+CNL*(MJ2-1)
    DO JP=1,CNL
    J=JP+CNL*(MJ1-1)
    DO IP=1,CNI
    I=IP+CNI*(MI-1)
    F11=(CNL+1)*(CNI+1)*(KP-1)+(CNI+1)*(JP-1)+IP-1
    F12=F11+1
    F13=F12+CNI+1
    F14=F13-1
    F21=F11+(CNL+1)*(CNI+1)
    F22=F21+1
    F23=F22+CNI+1
    F24=F23-1
    CALL CAL_FIJ_RECT(F11,F12,F13,F14,F21,F22,F23,F24)
    WRITE(2,235) '4(',F21,' ',F22,' ',F23,' ',F24,')'
    FAD1=F_BACK_RECT(IP,JP)
    WRITE(3,231) FAD1
    OWN=CELL_ID_RECT(IP,JP,KP)
    WRITE(4,231) OWN
    ENDDO   
    ENDDO
    
    N=5
    WRITE(31,231) N+4
    WRITE(31,*) '('  
    WRITE(31,*) '0'
    WRITE(31,*) '1'
    WRITE(31,*) '2' 
    WRITE(31,*) '3' 
    DO I=1,N
    WRITE(31,*) '-1' 
    ENDDO
    WRITE(31,*) ')'
    
    WRITE(21,231) N+4
    WRITE(21,*) "("    
    WRITE(21,*) "    INLET1"
    WRITE(21,*) "    {"    
    WRITE(21,*) "        type            patch;"
    N=0
    WRITE(21,232) "        nFaces          ",N,";"
    N=F_INTERNAL_NUM
    WRITE(21,232) "        startFace       ",N,";"
    WRITE(21,*) "    }"
    WRITE(21,*) "    INLET2"
    WRITE(21,*) "    {"    
    WRITE(21,*) "        type            patch;"
    N=0
    WRITE(21,232) "        nFaces          ",N,";"
    N=F_INTERNAL_NUM
    WRITE(21,232) "        startFace       ",N,";"
    WRITE(21,*) "    }"
    WRITE(21,*) "    OUTLET"
    WRITE(21,*) "    {"
    WRITE(21,*) "        type            patch;"
    N=0
    WRITE(21,232) "        nFaces          ",N,";"
    N=F_INTERNAL_NUM
    WRITE(21,232) "        startFace       ",N,";"
    WRITE(21,*) "    }"
    WRITE(21,*) "    HULL"
    WRITE(21,*) "    {"
    WRITE(21,*) "        type            wall;"
    WRITE(21,*) "        inGroups        1(wall);"
    N=CNL*CNL
    WRITE(21,232) "        nFaces          ",N,";"
    N=F_INTERNAL_NUM
    WRITE(21,232) "        startFace       ",N,";"
    WRITE(21,*) "    }"
    
        !DOWN FACE--PROCESSOR_MI1_TO_MI2  !CASE RECTANGLE 3-1-2
    PRO1=PNI*PNJ+PNL*PNI*(MJ2-1)+PNI*(MJ1-1)+MI-1
    MJ=1
    PRO2=PNI*(MJ-1)+MI-1
    WRITE(21,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(21,*) '    {'
    WRITE(21,*) '        type            processor;'
    WRITE(21,*) '        inGroups        List<word> 1(processor);'
    N=CNI*CNL
    WRITE(21,232) "        nFaces          ",N,";"
    N=F_INTERNAL_NUM+CNL*CNL
    WRITE(21,232) "        startFace       ",N,";"
    WRITE(21,*) '        matchTolerance  0.0001;'
    WRITE(21,232) "        myProcNo        ",PRO1,";"
    WRITE(21,232) "        neighbProcNo    ",PRO2,";"
    WRITE(21,*) '    }'
    WRITE(41,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(41,*) '    {'
    WRITE(41,*) '        type            processor;'
    WRITE(41,*) '        value           nonuniform List<scalar> '
    N=CNI*CNL
    WRITE(41,'(I0)') N
    WRITE(41,*) '('
    DO KP=1,CNL
    DO IP=1,CNI
    WRITE(41,'(E21.14)') P_DOWN_RECT(IP,KP)    
    ENDDO
    ENDDO
    WRITE(41,*) ')'
    WRITE(41,*) ';'
    WRITE(41,*) '    }'
    WRITE(51,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(51,*) '    {'
    WRITE(51,*) '        type            processor;'
    WRITE(51,*) '        value           nonuniform List<vector> '
    N=CNI*CNL
    WRITE(51,'(I0)') N
    WRITE(51,*) '('
    DO KP=1,CNL
    DO IP=1,CNI
    UU1=U_DOWN_RECT(1,IP,KP)
    UU2=U_DOWN_RECT(2,IP,KP)
    UU3=U_DOWN_RECT(3,IP,KP)
    WRITE(51,230) '( ',UU1,' ',UU2,' ',UU3,' )'  
    ENDDO
    ENDDO
    WRITE(51,*) ')'
    WRITE(51,*) ';'
    WRITE(51,*) '    }'
    WRITE(61,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(61,*) '    {'
    WRITE(61,*) '        type            processor;'
    WRITE(61,*) "        value           uniform 0;"
    WRITE(61,*) '    }'
    
        !FRONT FACE--PROCESSOR_MI1_TO_MI2  !CASE RECTANGLE 3-1-2
    PRO2=PNI*PNJ+PNL*PNI*(MJ2-2)+PNI*(MJ1-1)+MI-1
    WRITE(21,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(21,*) '    {'
    WRITE(21,*) '        type            processor;'
    WRITE(21,*) '        inGroups        List<word> 1(processor);'
    N=CNI*CNL
    WRITE(21,232) "        nFaces          ",N,";"
    N=F_INTERNAL_NUM+CNL*CNL+CNI*CNL
    WRITE(21,232) "        startFace       ",N,";"
    WRITE(21,*) '        matchTolerance  0.0001;'
    WRITE(21,232) "        myProcNo        ",PRO1,";"
    WRITE(21,232) "        neighbProcNo    ",PRO2,";"
    WRITE(21,*) '    }'
    WRITE(41,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(41,*) '    {'
    WRITE(41,*) '        type            processor;'
    WRITE(41,*) '        value           nonuniform List<scalar> '
    N=CNI*CNL
    WRITE(41,'(I0)') N
    WRITE(41,*) '('
    DO JP=1,CNL
    DO IP=1,CNI
    WRITE(41,'(E21.14)') P_FRONT_RECT(IP,JP)    
    ENDDO
    ENDDO
    WRITE(41,*) ')'
    WRITE(41,*) ';'
    WRITE(41,*) '    }'
    WRITE(51,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(51,*) '    {'
    WRITE(51,*) '        type            processor;'
    WRITE(51,*) '        value           nonuniform List<vector> '
    N=CNI*CNL
    WRITE(51,'(I0)') N
    WRITE(51,*) '('
    DO JP=1,CNL
    DO IP=1,CNI
    UU1=U_FRONT_RECT(1,IP,JP)
    UU2=U_FRONT_RECT(2,IP,JP)
    UU3=U_FRONT_RECT(3,IP,JP)
    WRITE(51,230) '( ',UU1,' ',UU2,' ',UU3,' )'  
    ENDDO
    ENDDO
    WRITE(51,*) ')'
    WRITE(51,*) ';'
    WRITE(51,*) '    }'
    WRITE(61,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(61,*) '    {'
    WRITE(61,*) '        type            processor;'
    WRITE(61,*) "        value           uniform 0;"
    WRITE(61,*) '    }'
    
        !RIGHT FACE--PROCESSOR_MI1_TO_MI2  !CASE RECTANGLE 3-1-2
    PRO2=PNI*PNJ+PNL*PNI*(MJ2-1)+PNI*(MJ1-1)+MI-2
    WRITE(21,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(21,*) '    {'
    WRITE(21,*) '        type            processor;'
    WRITE(21,*) '        inGroups        List<word> 1(processor);'
    N=CNL*CNL
    WRITE(21,232) "        nFaces          ",N,";"
    N=F_INTERNAL_NUM+CNL*CNL+CNI*CNL+CNI*CNL
    WRITE(21,232) "        startFace       ",N,";"
    WRITE(21,*) '        matchTolerance  0.0001;'
    WRITE(21,232) "        myProcNo        ",PRO1,";"
    WRITE(21,232) "        neighbProcNo    ",PRO2,";"
    WRITE(21,*) '    }'
    WRITE(41,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(41,*) '    {'
    WRITE(41,*) '        type            processor;'
    WRITE(41,*) '        value           nonuniform List<scalar> '
    N=CNL*CNL
    WRITE(41,'(I0)') N
    WRITE(41,*) '('
    DO KP=1,CNL
    DO JP=1,CNL
    WRITE(41,'(E21.14)') P_RIGHT_RECT(JP,KP)    
    ENDDO
    ENDDO
    WRITE(41,*) ')'
    WRITE(41,*) ';'
    WRITE(41,*) '    }'
    WRITE(51,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(51,*) '    {'
    WRITE(51,*) '        type            processor;'
    WRITE(51,*) '        value           nonuniform List<vector> '
    N=CNL*CNL
    WRITE(51,'(I0)') N
    WRITE(51,*) '('
    DO KP=1,CNL
    DO JP=1,CNL
    UU1=U_RIGHT_RECT(1,JP,KP)
    UU2=U_RIGHT_RECT(2,JP,KP)
    UU3=U_RIGHT_RECT(3,JP,KP)
    WRITE(51,230) '( ',UU1,' ',UU2,' ',UU3,' )'  
    ENDDO
    ENDDO
    WRITE(51,*) ')'
    WRITE(51,*) ';'
    WRITE(51,*) '    }'
    WRITE(61,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(61,*) '    {'
    WRITE(61,*) '        type            processor;'
    WRITE(61,*) "        value           uniform 0;"
    WRITE(61,*) '    }'
    
        !UP FACE--PROCESSOR_MI1_TO_MI2  !CASE RECTANGLE 3-1-2
    PRO2=PNI*PNJ+PNL*PNI*(MJ2-1)+PNI*MJ1+MI-1
    WRITE(21,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(21,*) '    {'
    WRITE(21,*) '        type            processor;'
    WRITE(21,*) '        inGroups        List<word> 1(processor);'
    N=CNI*CNL
    WRITE(21,232) "        nFaces          ",N,";"
    N=F_INTERNAL_NUM+CNL*CNL+CNI*CNL+CNI*CNL+CNL*CNL
    WRITE(21,232) "        startFace       ",N,";"
    WRITE(21,*) '        matchTolerance  0.0001;'
    WRITE(21,232) "        myProcNo        ",PRO1,";"
    WRITE(21,232) "        neighbProcNo    ",PRO2,";"
    WRITE(21,*) '    }'
    WRITE(41,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(41,*) '    {'
    WRITE(41,*) '        type            processor;'
    WRITE(41,*) '        value           nonuniform List<scalar> '
    N=CNI*CNL
    WRITE(41,'(I0)') N
    WRITE(41,*) '('
    DO KP=1,CNL
    DO IP=1,CNI
    WRITE(41,'(E21.14)') P_UP_RECT(IP,KP)    
    ENDDO
    ENDDO
    WRITE(41,*) ')'
    WRITE(41,*) ';'
    WRITE(41,*) '    }'
    WRITE(51,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(51,*) '    {'
    WRITE(51,*) '        type            processor;'
    WRITE(51,*) '        value           nonuniform List<vector> '
    N=CNI*CNL
    WRITE(51,'(I0)') N
    WRITE(51,*) '('
    DO KP=1,CNL
    DO IP=1,CNI
    UU1=U_UP_RECT(1,IP,KP)
    UU2=U_UP_RECT(2,IP,KP)
    UU3=U_UP_RECT(3,IP,KP)
    WRITE(51,230) '( ',UU1,' ',UU2,' ',UU3,' )'  
    ENDDO
    ENDDO
    WRITE(51,*) ')'
    WRITE(51,*) ';'
    WRITE(51,*) '    }'
    WRITE(61,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(61,*) '    {'
    WRITE(61,*) '        type            processor;'
    WRITE(61,*) "        value           uniform 0;"
    WRITE(61,*) '    }'
    
        !BACK FACE--PROCESSOR_MI1_TO_MI2  !CASE RECTANGLE 3-1-2
    PRO2=PNI*PNJ+PNL*PNI*MJ2+PNI*(MJ1-1)+MI-1
    WRITE(21,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(21,*) '    {'
    WRITE(21,*) '        type            processor;'
    WRITE(21,*) '        inGroups        List<word> 1(processor);'
    N=CNI*CNL
    WRITE(21,232) "        nFaces          ",N,";"
    N=F_INTERNAL_NUM+CNL*CNL+CNI*CNL+CNI*CNL+CNL*CNL+CNI*CNL
    WRITE(21,232) "        startFace       ",N,";"
    WRITE(21,*) '        matchTolerance  0.0001;'
    WRITE(21,232) "        myProcNo        ",PRO1,";"
    WRITE(21,232) "        neighbProcNo    ",PRO2,";"
    WRITE(21,*) '    }'
    WRITE(41,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(41,*) '    {'
    WRITE(41,*) '        type            processor;'
    WRITE(41,*) '        value           nonuniform List<scalar> '
    N=CNI*CNL
    WRITE(41,'(I0)') N
    WRITE(41,*) '('
    DO JP=1,CNL
    DO IP=1,CNI
    WRITE(41,'(E21.14)') P_BACK_RECT(IP,JP)    
    ENDDO
    ENDDO
    WRITE(41,*) ')'
    WRITE(41,*) ';'
    WRITE(41,*) '    }'
    WRITE(51,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(51,*) '    {'
    WRITE(51,*) '        type            processor;'
    WRITE(51,*) '        value           nonuniform List<vector> '
    N=CNI*CNL
    WRITE(51,'(I0)') N
    WRITE(51,*) '('
    DO JP=1,CNL
    DO IP=1,CNI
    UU1=U_BACK_RECT(1,IP,JP)
    UU2=U_BACK_RECT(2,IP,JP)
    UU3=U_BACK_RECT(3,IP,JP)
    WRITE(51,230) '( ',UU1,' ',UU2,' ',UU3,' )'  
    ENDDO
    ENDDO
    WRITE(51,*) ')'
    WRITE(51,*) ';'
    WRITE(51,*) '    }'
    WRITE(61,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(61,*) '    {'
    WRITE(61,*) '        type            processor;'
    WRITE(61,*) "        value           uniform 0;"
    WRITE(61,*) '    }'
    
    
    
    
    
    
    
    ENDIF    
    ELSEIF (MJ1.EQ.PNL) THEN
    IF (MJ2.EQ.1) THEN
    !CASE RECTANGLE 3-3-1
        
        !LEFT FACE--HULL  !CASE RECTANGLE 3-3-1 
    IP=CNI
    I=IP+CNI*(MI-1)
    DO KP=1,CNL
    K=KP+CNL*(MJ2-1)
    DO JP=1,CNL
    J=JP+CNL*(MJ1-1)
    F11=(CNL+1)*(CNI+1)*(KP-1)+(CNI+1)*(JP-1)+IP-1
    F12=F11+1
    F13=F12+CNI+1
    F14=F13-1
    F21=F11+(CNL+1)*(CNI+1)
    F22=F21+1
    F23=F22+CNI+1
    F24=F23-1
    CALL CAL_FIJ_RECT(F11,F12,F13,F14,F21,F22,F23,F24)
    WRITE(2,235) '4(',F12,' ',F13,' ',F23,' ',F22,')'
    FAD1=F_INTERNAL_NUM_ALL+NL*NL+NJ*NK+NI*NK+NJ*NK+NL*(K-1)+J
    WRITE(3,231) FAD1
    OWN=CELL_ID_RECT(IP,JP,KP)
    WRITE(4,231) OWN
    ENDDO   
    ENDDO
    
        !FRONT FACE--PROCESSOR_MI1_TO_MI2  !CASE RECTANGLE 3-3-1
    KP=1
    K=KP+CNL*(MJ2-1)
    DO JP=1,CNL
    J=JP+CNL*(MJ1-1)
    DO IP=1,CNI
    I=IP+CNI*(MI-1)
    F11=(CNL+1)*(CNI+1)*(KP-1)+(CNI+1)*(JP-1)+IP-1
    F12=F11+1
    F13=F12+CNI+1
    F14=F13-1
    F21=F11+(CNL+1)*(CNI+1)
    F22=F21+1
    F23=F22+CNI+1
    F24=F23-1
    CALL CAL_FIJ_RECT(F11,F12,F13,F14,F21,F22,F23,F24)
    WRITE(2,235) '4(',F11,' ',F14,' ',F13,' ',F12,')'
    FAD1=F_FRONT_RECT(IP,JP)
    WRITE(3,231) FAD1
    OWN=CELL_ID_RECT(IP,JP,KP)
    WRITE(4,231) OWN
    ENDDO   
    ENDDO
    
        !UP FACE--PROCESSOR_MI1_TO_MI2  !CASE RECTANGLE 3-3-1
    JP=CNL
    J=JP+CNL*(MJ1-1)
    DO KP=1,CNL
    K=KP+CNL*(MJ2-1)
    DO IP=1,CNI
    I=IP+CNI*(MI-1)
    F11=(CNL+1)*(CNI+1)*(KP-1)+(CNI+1)*(JP-1)+IP-1
    F12=F11+1
    F13=F12+CNI+1
    F14=F13-1
    F21=F11+(CNL+1)*(CNI+1)
    F22=F21+1
    F23=F22+CNI+1
    F24=F23-1
    CALL CAL_FIJ_RECT(F11,F12,F13,F14,F21,F22,F23,F24)
    WRITE(2,235) '4(',F14,' ',F24,' ',F23,' ',F13,')'
    FAD1=F_UP_RECT(IP,KP)
    WRITE(3,231) FAD1
    OWN=CELL_ID_RECT(IP,JP,KP)
    WRITE(4,231) OWN
    ENDDO   
    ENDDO
    
        !DOWN FACE--PROCESSOR_MI1_TO_MI2  !CASE RECTANGLE 3-3-1
    JP=1
    J=JP+CNL*(MJ1-1)
    DO KP=1,CNL
    K=KP+CNL*(MJ2-1)
    DO IP=1,CNI
    I=IP+CNI*(MI-1)
    F11=(CNL+1)*(CNI+1)*(KP-1)+(CNI+1)*(JP-1)+IP-1
    F12=F11+1
    F13=F12+CNI+1
    F14=F13-1
    F21=F11+(CNL+1)*(CNI+1)
    F22=F21+1
    F23=F22+CNI+1
    F24=F23-1
    CALL CAL_FIJ_RECT(F11,F12,F13,F14,F21,F22,F23,F24)
    WRITE(2,235) '4(',F11,' ',F12,' ',F22,' ',F21,')'
    FAD1=F_DOWN_RECT(IP,KP)
    WRITE(3,231) FAD1
    OWN=CELL_ID_RECT(IP,JP,KP)
    WRITE(4,231) OWN
    ENDDO   
    ENDDO
    
        !RIGHT FACE--PROCESSOR_MI1_TO_MI2  !CASE RECTANGLE 3-3-1 
    IP=1
    I=IP+CNI*(MI-1)
    DO KP=1,CNL
    K=KP+CNL*(MJ2-1)
    DO JP=1,CNL
    J=JP+CNL*(MJ1-1)
    F11=(CNL+1)*(CNI+1)*(KP-1)+(CNI+1)*(JP-1)+IP-1
    F12=F11+1
    F13=F12+CNI+1
    F14=F13-1
    F21=F11+(CNL+1)*(CNI+1)
    F22=F21+1
    F23=F22+CNI+1
    F24=F23-1
    CALL CAL_FIJ_RECT(F11,F12,F13,F14,F21,F22,F23,F24)
    WRITE(2,235) '4(',F11,' ',F21,' ',F24,' ',F14,')'
    FAD1=F_RIGHT_RECT(JP,KP)
    WRITE(3,231) FAD1
    OWN=CELL_ID_RECT(IP,JP,KP)
    WRITE(4,231) OWN
    ENDDO   
    ENDDO
    
        !BACK FACE--PROCESSOR_MI1_TO_MI2  !CASE RECTANGLE 3-3-1
    KP=CNL
    K=KP+CNL*(MJ2-1)
    DO JP=1,CNL
    J=JP+CNL*(MJ1-1)
    DO IP=1,CNI
    I=IP+CNI*(MI-1)
    F11=(CNL+1)*(CNI+1)*(KP-1)+(CNI+1)*(JP-1)+IP-1
    F12=F11+1
    F13=F12+CNI+1
    F14=F13-1
    F21=F11+(CNL+1)*(CNI+1)
    F22=F21+1
    F23=F22+CNI+1
    F24=F23-1
    CALL CAL_FIJ_RECT(F11,F12,F13,F14,F21,F22,F23,F24)
    WRITE(2,235) '4(',F21,' ',F22,' ',F23,' ',F24,')'
    FAD1=F_BACK_RECT(IP,JP)
    WRITE(3,231) FAD1
    OWN=CELL_ID_RECT(IP,JP,KP)
    WRITE(4,231) OWN
    ENDDO   
    ENDDO
    
    N=4
    WRITE(31,231) N+4
    WRITE(31,*) '('  
    WRITE(31,*) '0'
    WRITE(31,*) '1'
    WRITE(31,*) '2' 
    WRITE(31,*) '3' 
    DO I=1,N
    WRITE(31,*) '-1' 
    ENDDO
    WRITE(31,*) ')'
    
    WRITE(21,231) N+4
    WRITE(21,*) "("    
    WRITE(21,*) "    INLET1"
    WRITE(21,*) "    {"    
    WRITE(21,*) "        type            patch;"
    N=0
    WRITE(21,232) "        nFaces          ",N,";"
    N=F_INTERNAL_NUM
    WRITE(21,232) "        startFace       ",N,";"
    WRITE(21,*) "    }"
    WRITE(21,*) "    INLET2"
    WRITE(21,*) "    {"    
    WRITE(21,*) "        type            patch;"
    N=0
    WRITE(21,232) "        nFaces          ",N,";"
    N=F_INTERNAL_NUM
    WRITE(21,232) "        startFace       ",N,";"
    WRITE(21,*) "    }"
    WRITE(21,*) "    OUTLET"
    WRITE(21,*) "    {"
    WRITE(21,*) "        type            patch;"
    N=0
    WRITE(21,232) "        nFaces          ",N,";"
    N=F_INTERNAL_NUM
    WRITE(21,232) "        startFace       ",N,";"
    WRITE(21,*) "    }"
    WRITE(21,*) "    HULL"
    WRITE(21,*) "    {"
    WRITE(21,*) "        type            wall;"
    WRITE(21,*) "        inGroups        1(wall);"
    N=CNL*CNL
    WRITE(21,232) "        nFaces          ",N,";"
    N=F_INTERNAL_NUM
    WRITE(21,232) "        startFace       ",N,";"
    WRITE(21,*) "    }"
    
        !FRONT/UP FACE--PROCESSOR_MI1_TO_MI2  !CASE RECTANGLE 3-3-1
    PRO1=PNI*PNJ+PNL*PNI*(MJ2-1)+PNI*(MJ1-1)+MI-1
    MJ=1
    PRO2=PNI*(MJ-1)+MI-1
    WRITE(21,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(21,*) '    {'
    WRITE(21,*) '        type            processor;'
    WRITE(21,*) '        inGroups        List<word> 1(processor);'
    N=CNI*CNL*2
    WRITE(21,232) "        nFaces          ",N,";"
    N=F_INTERNAL_NUM+CNL*CNL
    WRITE(21,232) "        startFace       ",N,";"
    WRITE(21,*) '        matchTolerance  0.0001;'
    WRITE(21,232) "        myProcNo        ",PRO1,";"
    WRITE(21,232) "        neighbProcNo    ",PRO2,";"
    WRITE(21,*) '    }'
    WRITE(41,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(41,*) '    {'
    WRITE(41,*) '        type            processor;'
    WRITE(41,*) '        value           nonuniform List<scalar> '
    N=CNI*CNL*2
    WRITE(41,'(I0)') N
    WRITE(41,*) '('
    DO JP=1,CNL
    DO IP=1,CNI
    WRITE(41,'(E21.14)') P_FRONT_RECT(IP,JP)    
    ENDDO
    ENDDO
    DO KP=1,CNL
    DO IP=1,CNI
    WRITE(41,'(E21.14)') P_UP_RECT(IP,KP)    
    ENDDO
    ENDDO
    WRITE(41,*) ')'
    WRITE(41,*) ';'
    WRITE(41,*) '    }'
    WRITE(51,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(51,*) '    {'
    WRITE(51,*) '        type            processor;'
    WRITE(51,*) '        value           nonuniform List<vector> '
    N=CNI*CNL*2
    WRITE(51,'(I0)') N
    WRITE(51,*) '('
    DO JP=1,CNL
    DO IP=1,CNI
    UU1=U_FRONT_RECT(1,IP,JP)
    UU2=U_FRONT_RECT(2,IP,JP)
    UU3=U_FRONT_RECT(3,IP,JP)
    WRITE(51,230) '( ',UU1,' ',UU2,' ',UU3,' )'  
    ENDDO
    ENDDO
    DO KP=1,CNL
    DO IP=1,CNI
    UU1=U_UP_RECT(1,IP,KP)
    UU2=U_UP_RECT(2,IP,KP)
    UU3=U_UP_RECT(3,IP,KP)
    WRITE(51,230) '( ',UU1,' ',UU2,' ',UU3,' )'  
    ENDDO
    ENDDO
    WRITE(51,*) ')'
    WRITE(51,*) ';'
    WRITE(51,*) '    }'
    WRITE(61,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(61,*) '    {'
    WRITE(61,*) '        type            processor;'
    WRITE(61,*) "        value           uniform 0;"
    WRITE(61,*) '    }'
    
        !DOWN FACE--PROCESSOR_MI1_TO_MI2  !CASE RECTANGLE 3-3-1
    PRO2=PNI*PNJ+PNL*PNI*(MJ2-1)+PNI*(MJ1-2)+MI-1
    WRITE(21,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(21,*) '    {'
    WRITE(21,*) '        type            processor;'
    WRITE(21,*) '        inGroups        List<word> 1(processor);'
    N=CNI*CNL
    WRITE(21,232) "        nFaces          ",N,";"
    N=F_INTERNAL_NUM+CNL*CNL+CNI*CNL*2
    WRITE(21,232) "        startFace       ",N,";"
    WRITE(21,*) '        matchTolerance  0.0001;'
    WRITE(21,232) "        myProcNo        ",PRO1,";"
    WRITE(21,232) "        neighbProcNo    ",PRO2,";"
    WRITE(21,*) '    }'
    WRITE(41,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(41,*) '    {'
    WRITE(41,*) '        type            processor;'
    WRITE(41,*) '        value           nonuniform List<scalar> '
    N=CNI*CNL
    WRITE(41,'(I0)') N
    WRITE(41,*) '('
    DO KP=1,CNL
    DO IP=1,CNI
    WRITE(41,'(E21.14)') P_DOWN_RECT(IP,KP)    
    ENDDO
    ENDDO
    WRITE(41,*) ')'
    WRITE(41,*) ';'
    WRITE(41,*) '    }'
    WRITE(51,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(51,*) '    {'
    WRITE(51,*) '        type            processor;'
    WRITE(51,*) '        value           nonuniform List<vector> '
    N=CNI*CNL
    WRITE(51,'(I0)') N
    WRITE(51,*) '('
    DO KP=1,CNL
    DO IP=1,CNI
    UU1=U_DOWN_RECT(1,IP,KP)
    UU2=U_DOWN_RECT(2,IP,KP)
    UU3=U_DOWN_RECT(3,IP,KP)
    WRITE(51,230) '( ',UU1,' ',UU2,' ',UU3,' )'  
    ENDDO
    ENDDO
    WRITE(51,*) ')'
    WRITE(51,*) ';'
    WRITE(51,*) '    }'
    WRITE(61,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(61,*) '    {'
    WRITE(61,*) '        type            processor;'
    WRITE(61,*) "        value           uniform 0;"
    WRITE(61,*) '    }'
    
        !RIGHT FACE--PROCESSOR_MI1_TO_MI2  !CASE RECTANGLE 3-3-1
    PRO2=PNI*PNJ+PNL*PNI*(MJ2-1)+PNI*(MJ1-1)+MI-2
    WRITE(21,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(21,*) '    {'
    WRITE(21,*) '        type            processor;'
    WRITE(21,*) '        inGroups        List<word> 1(processor);'
    N=CNL*CNL
    WRITE(21,232) "        nFaces          ",N,";"
    N=F_INTERNAL_NUM+CNL*CNL+CNI*CNL*2+CNI*CNL
    WRITE(21,232) "        startFace       ",N,";"
    WRITE(21,*) '        matchTolerance  0.0001;'
    WRITE(21,232) "        myProcNo        ",PRO1,";"
    WRITE(21,232) "        neighbProcNo    ",PRO2,";"
    WRITE(21,*) '    }'
    WRITE(41,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(41,*) '    {'
    WRITE(41,*) '        type            processor;'
    WRITE(41,*) '        value           nonuniform List<scalar> '
    N=CNL*CNL
    WRITE(41,'(I0)') N
    WRITE(41,*) '('
    DO KP=1,CNL
    DO JP=1,CNL
    WRITE(41,'(E21.14)') P_RIGHT_RECT(JP,KP)    
    ENDDO
    ENDDO
    WRITE(41,*) ')'
    WRITE(41,*) ';'
    WRITE(41,*) '    }'
    WRITE(51,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(51,*) '    {'
    WRITE(51,*) '        type            processor;'
    WRITE(51,*) '        value           nonuniform List<vector> '
    N=CNL*CNL
    WRITE(51,'(I0)') N
    WRITE(51,*) '('
    DO KP=1,CNL
    DO JP=1,CNL
    UU1=U_RIGHT_RECT(1,JP,KP)
    UU2=U_RIGHT_RECT(2,JP,KP)
    UU3=U_RIGHT_RECT(3,JP,KP)
    WRITE(51,230) '( ',UU1,' ',UU2,' ',UU3,' )'  
    ENDDO
    ENDDO
    WRITE(51,*) ')'
    WRITE(51,*) ';'
    WRITE(51,*) '    }'
    WRITE(61,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(61,*) '    {'
    WRITE(61,*) '        type            processor;'
    WRITE(61,*) "        value           uniform 0;"
    WRITE(61,*) '    }'
    
        !BACK FACE--PROCESSOR_MI1_TO_MI2  !CASE RECTANGLE 3-3-1
    PRO2=PNI*PNJ+PNL*PNI*MJ2+PNI*(MJ1-1)+MI-1
    WRITE(21,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(21,*) '    {'
    WRITE(21,*) '        type            processor;'
    WRITE(21,*) '        inGroups        List<word> 1(processor);'
    N=CNI*CNL
    WRITE(21,232) "        nFaces          ",N,";"
    N=F_INTERNAL_NUM+CNL*CNL+CNI*CNL*2+CNI*CNL+CNL*CNL
    WRITE(21,232) "        startFace       ",N,";"
    WRITE(21,*) '        matchTolerance  0.0001;'
    WRITE(21,232) "        myProcNo        ",PRO1,";"
    WRITE(21,232) "        neighbProcNo    ",PRO2,";"
    WRITE(21,*) '    }'
    WRITE(41,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(41,*) '    {'
    WRITE(41,*) '        type            processor;'
    WRITE(41,*) '        value           nonuniform List<scalar> '
    N=CNI*CNL
    WRITE(41,'(I0)') N
    WRITE(41,*) '('
    DO JP=1,CNL
    DO IP=1,CNI
    WRITE(41,'(E21.14)') P_BACK_RECT(IP,JP)    
    ENDDO
    ENDDO
    WRITE(41,*) ')'
    WRITE(41,*) ';'
    WRITE(41,*) '    }'
    WRITE(51,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(51,*) '    {'
    WRITE(51,*) '        type            processor;'
    WRITE(51,*) '        value           nonuniform List<vector> '
    N=CNI*CNL
    WRITE(51,'(I0)') N
    WRITE(51,*) '('
    DO JP=1,CNL
    DO IP=1,CNI
    UU1=U_BACK_RECT(1,IP,JP)
    UU2=U_BACK_RECT(2,IP,JP)
    UU3=U_BACK_RECT(3,IP,JP)
    WRITE(51,230) '( ',UU1,' ',UU2,' ',UU3,' )'  
    ENDDO
    ENDDO
    WRITE(51,*) ')'
    WRITE(51,*) ';'
    WRITE(51,*) '    }'
    WRITE(61,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(61,*) '    {'
    WRITE(61,*) '        type            processor;'
    WRITE(61,*) "        value           uniform 0;"
    WRITE(61,*) '    }'
    
    
    
    ELSEIF (MJ2.EQ.PNL) THEN
    !CASE RECTANGLE 3-3-3
        
        !LEFT FACE--HULL  !CASE RECTANGLE 3-3-3 
    IP=CNI
    I=IP+CNI*(MI-1)
    DO KP=1,CNL
    K=KP+CNL*(MJ2-1)
    DO JP=1,CNL
    J=JP+CNL*(MJ1-1)
    F11=(CNL+1)*(CNI+1)*(KP-1)+(CNI+1)*(JP-1)+IP-1
    F12=F11+1
    F13=F12+CNI+1
    F14=F13-1
    F21=F11+(CNL+1)*(CNI+1)
    F22=F21+1
    F23=F22+CNI+1
    F24=F23-1
    CALL CAL_FIJ_RECT(F11,F12,F13,F14,F21,F22,F23,F24)
    WRITE(2,235) '4(',F12,' ',F13,' ',F23,' ',F22,')'
    FAD1=F_INTERNAL_NUM_ALL+NL*NL+NJ*NK+NI*NK+NJ*NK+NL*(K-1)+J
    WRITE(3,231) FAD1
    OWN=CELL_ID_RECT(IP,JP,KP)
    WRITE(4,231) OWN
    ENDDO   
    ENDDO
    
        !UP FACE--PROCESSOR_MI1_TO_MI2  !CASE RECTANGLE 3-3-3
    JP=CNL
    J=JP+CNL*(MJ1-1)
    DO KP=1,CNL
    K=KP+CNL*(MJ2-1)
    DO IP=1,CNI
    I=IP+CNI*(MI-1)
    F11=(CNL+1)*(CNI+1)*(KP-1)+(CNI+1)*(JP-1)+IP-1
    F12=F11+1
    F13=F12+CNI+1
    F14=F13-1
    F21=F11+(CNL+1)*(CNI+1)
    F22=F21+1
    F23=F22+CNI+1
    F24=F23-1
    CALL CAL_FIJ_RECT(F11,F12,F13,F14,F21,F22,F23,F24)
    WRITE(2,235) '4(',F14,' ',F24,' ',F23,' ',F13,')'
    FAD1=F_UP_RECT(IP,KP)
    WRITE(3,231) FAD1
    OWN=CELL_ID_RECT(IP,JP,KP)
    WRITE(4,231) OWN
    ENDDO   
    ENDDO
    
        !BACK FACE--PROCESSOR_MI1_TO_MI2  !CASE RECTANGLE 3-3-3
    KP=CNL
    K=KP+CNL*(MJ2-1)
    DO JP=CNL,1,-1
    J=JP+CNL*(MJ1-1)
    DO IP=1,CNI
    I=IP+CNI*(MI-1)
    F11=(CNL+1)*(CNI+1)*(KP-1)+(CNI+1)*(JP-1)+IP-1
    F12=F11+1
    F13=F12+CNI+1
    F14=F13-1
    F21=F11+(CNL+1)*(CNI+1)
    F22=F21+1
    F23=F22+CNI+1
    F24=F23-1
    CALL CAL_FIJ_RECT(F11,F12,F13,F14,F21,F22,F23,F24)
    !WRITE(2,235) '4(',F21,' ',F22,' ',F23,' ',F24,')'
    WRITE(2,235) '4(',F24,' ',F21,' ',F22,' ',F23,')'
    FAD1=F_BACK_RECT(IP,JP)
    WRITE(3,231) FAD1
    OWN=CELL_ID_RECT(IP,JP,KP)
    WRITE(4,231) OWN
    ENDDO   
    ENDDO
    
        !FRONT FACE--PROCESSOR_MI1_TO_MI2  !CASE RECTANGLE 3-3-3
    KP=1
    K=KP+CNL*(MJ2-1)
    DO JP=1,CNL
    J=JP+CNL*(MJ1-1)
    DO IP=1,CNI
    I=IP+CNI*(MI-1)
    F11=(CNL+1)*(CNI+1)*(KP-1)+(CNI+1)*(JP-1)+IP-1
    F12=F11+1
    F13=F12+CNI+1
    F14=F13-1
    F21=F11+(CNL+1)*(CNI+1)
    F22=F21+1
    F23=F22+CNI+1
    F24=F23-1
    CALL CAL_FIJ_RECT(F11,F12,F13,F14,F21,F22,F23,F24)
    WRITE(2,235) '4(',F11,' ',F14,' ',F13,' ',F12,')'
    FAD1=F_FRONT_RECT(IP,JP)
    WRITE(3,231) FAD1
    OWN=CELL_ID_RECT(IP,JP,KP)
    WRITE(4,231) OWN
    ENDDO   
    ENDDO
    
        !DOWN FACE--PROCESSOR_MI1_TO_MI2  !CASE RECTANGLE 3-3-3
    JP=1
    J=JP+CNL*(MJ1-1)
    DO KP=1,CNL
    K=KP+CNL*(MJ2-1)
    DO IP=1,CNI
    I=IP+CNI*(MI-1)
    F11=(CNL+1)*(CNI+1)*(KP-1)+(CNI+1)*(JP-1)+IP-1
    F12=F11+1
    F13=F12+CNI+1
    F14=F13-1
    F21=F11+(CNL+1)*(CNI+1)
    F22=F21+1
    F23=F22+CNI+1
    F24=F23-1
    CALL CAL_FIJ_RECT(F11,F12,F13,F14,F21,F22,F23,F24)
    WRITE(2,235) '4(',F11,' ',F12,' ',F22,' ',F21,')'
    FAD1=F_DOWN_RECT(IP,KP)
    WRITE(3,231) FAD1
    OWN=CELL_ID_RECT(IP,JP,KP)
    WRITE(4,231) OWN
    ENDDO   
    ENDDO
    
        !RIGHT FACE--PROCESSOR_MI1_TO_MI2  !CASE RECTANGLE 3-3-3 
    IP=1
    I=IP+CNI*(MI-1)
    DO KP=1,CNL
    K=KP+CNL*(MJ2-1)
    DO JP=1,CNL
    J=JP+CNL*(MJ1-1)
    F11=(CNL+1)*(CNI+1)*(KP-1)+(CNI+1)*(JP-1)+IP-1
    F12=F11+1
    F13=F12+CNI+1
    F14=F13-1
    F21=F11+(CNL+1)*(CNI+1)
    F22=F21+1
    F23=F22+CNI+1
    F24=F23-1
    CALL CAL_FIJ_RECT(F11,F12,F13,F14,F21,F22,F23,F24)
    WRITE(2,235) '4(',F11,' ',F21,' ',F24,' ',F14,')'
    FAD1=F_RIGHT_RECT(JP,KP)
    WRITE(3,231) FAD1
    OWN=CELL_ID_RECT(IP,JP,KP)
    WRITE(4,231) OWN
    ENDDO   
    ENDDO
    
    N=4
    WRITE(31,231) N+4
    WRITE(31,*) '('  
    WRITE(31,*) '0'
    WRITE(31,*) '1'
    WRITE(31,*) '2' 
    WRITE(31,*) '3' 
    DO I=1,N
    WRITE(31,*) '-1' 
    ENDDO
    WRITE(31,*) ')'
    
    WRITE(21,231) N+4
    WRITE(21,*) "("    
    WRITE(21,*) "    INLET1"
    WRITE(21,*) "    {"    
    WRITE(21,*) "        type            patch;"
    N=0
    WRITE(21,232) "        nFaces          ",N,";"
    N=F_INTERNAL_NUM
    WRITE(21,232) "        startFace       ",N,";"
    WRITE(21,*) "    }"
    WRITE(21,*) "    INLET2"
    WRITE(21,*) "    {"    
    WRITE(21,*) "        type            patch;"
    N=0
    WRITE(21,232) "        nFaces          ",N,";"
    N=F_INTERNAL_NUM
    WRITE(21,232) "        startFace       ",N,";"
    WRITE(21,*) "    }"
    WRITE(21,*) "    OUTLET"
    WRITE(21,*) "    {"
    WRITE(21,*) "        type            patch;"
    N=0
    WRITE(21,232) "        nFaces          ",N,";"
    N=F_INTERNAL_NUM
    WRITE(21,232) "        startFace       ",N,";"
    WRITE(21,*) "    }"
    WRITE(21,*) "    HULL"
    WRITE(21,*) "    {"
    WRITE(21,*) "        type            wall;"
    WRITE(21,*) "        inGroups        1(wall);"
    N=CNL*CNL
    WRITE(21,232) "        nFaces          ",N,";"
    N=F_INTERNAL_NUM
    WRITE(21,232) "        startFace       ",N,";"
    WRITE(21,*) "    }"
    
        !UP/BACK FACE--PROCESSOR_MI1_TO_MI2  !CASE RECTANGLE 3-3-3
    PRO1=PNI*PNJ+PNL*PNI*(MJ2-1)+PNI*(MJ1-1)+MI-1
    MJ=1
    PRO2=PNI*(MJ-1)+MI-1
    WRITE(21,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(21,*) '    {'
    WRITE(21,*) '        type            processor;'
    WRITE(21,*) '        inGroups        List<word> 1(processor);'
    N=CNI*CNL*2
    WRITE(21,232) "        nFaces          ",N,";"
    N=F_INTERNAL_NUM+CNL*CNL
    WRITE(21,232) "        startFace       ",N,";"
    WRITE(21,*) '        matchTolerance  0.0001;'
    WRITE(21,232) "        myProcNo        ",PRO1,";"
    WRITE(21,232) "        neighbProcNo    ",PRO2,";"
    WRITE(21,*) '    }'
    WRITE(41,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(41,*) '    {'
    WRITE(41,*) '        type            processor;'
    WRITE(41,*) '        value           nonuniform List<scalar> '
    N=CNI*CNL*2
    WRITE(41,'(I0)') N
    WRITE(41,*) '('
    DO KP=1,CNL
    DO IP=1,CNI
    WRITE(41,'(E21.14)') P_UP_RECT(IP,KP)    
    ENDDO
    ENDDO
    DO JP=1,CNL
    DO IP=1,CNI
    WRITE(41,'(E21.14)') P_BACK_RECT(IP,JP)    
    ENDDO
    ENDDO
    WRITE(41,*) ')'
    WRITE(41,*) ';'
    WRITE(41,*) '    }'
    WRITE(51,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(51,*) '    {'
    WRITE(51,*) '        type            processor;'
    WRITE(51,*) '        value           nonuniform List<vector> '
    N=CNI*CNL*2
    WRITE(51,'(I0)') N
    WRITE(51,*) '('
    DO KP=1,CNL
    DO IP=1,CNI
    UU1=U_UP_RECT(1,IP,KP)
    UU2=U_UP_RECT(2,IP,KP)
    UU3=U_UP_RECT(3,IP,KP)
    WRITE(51,230) '( ',UU1,' ',UU2,' ',UU3,' )'  
    ENDDO
    ENDDO
    DO JP=1,CNL
    DO IP=1,CNI
    UU1=U_BACK_RECT(1,IP,JP)
    UU2=U_BACK_RECT(2,IP,JP)
    UU3=U_BACK_RECT(3,IP,JP)
    WRITE(51,230) '( ',UU1,' ',UU2,' ',UU3,' )'  
    ENDDO
    ENDDO
    WRITE(51,*) ')'
    WRITE(51,*) ';'
    WRITE(51,*) '    }'
    WRITE(61,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(61,*) '    {'
    WRITE(61,*) '        type            processor;'
    WRITE(61,*) "        value           uniform 0;"
    WRITE(61,*) '    }'
    
        !FRONT FACE--PROCESSOR_MI1_TO_MI2  !CASE RECTANGLE 3-3-3
    PRO2=PNI*PNJ+PNL*PNI*(MJ2-2)+PNI*(MJ1-1)+MI-1
    WRITE(21,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(21,*) '    {'
    WRITE(21,*) '        type            processor;'
    WRITE(21,*) '        inGroups        List<word> 1(processor);'
    N=CNI*CNL
    WRITE(21,232) "        nFaces          ",N,";"
    N=F_INTERNAL_NUM+CNL*CNL+CNI*CNL*2
    WRITE(21,232) "        startFace       ",N,";"
    WRITE(21,*) '        matchTolerance  0.0001;'
    WRITE(21,232) "        myProcNo        ",PRO1,";"
    WRITE(21,232) "        neighbProcNo    ",PRO2,";"
    WRITE(21,*) '    }'
    WRITE(41,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(41,*) '    {'
    WRITE(41,*) '        type            processor;'
    WRITE(41,*) '        value           nonuniform List<scalar> '
    N=CNI*CNL
    WRITE(41,'(I0)') N
    WRITE(41,*) '('
    DO JP=1,CNL
    DO IP=1,CNI
    WRITE(41,'(E21.14)') P_FRONT_RECT(IP,JP)    
    ENDDO
    ENDDO
    WRITE(41,*) ')'
    WRITE(41,*) ';'
    WRITE(41,*) '    }'
    WRITE(51,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(51,*) '    {'
    WRITE(51,*) '        type            processor;'
    WRITE(51,*) '        value           nonuniform List<vector> '
    N=CNI*CNL
    WRITE(51,'(I0)') N
    WRITE(51,*) '('
    DO JP=1,CNL
    DO IP=1,CNI
    UU1=U_FRONT_RECT(1,IP,JP)
    UU2=U_FRONT_RECT(2,IP,JP)
    UU3=U_FRONT_RECT(3,IP,JP)
    WRITE(51,230) '( ',UU1,' ',UU2,' ',UU3,' )'  
    ENDDO
    ENDDO
    WRITE(51,*) ')'
    WRITE(51,*) ';'
    WRITE(51,*) '    }'
    WRITE(61,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(61,*) '    {'
    WRITE(61,*) '        type            processor;'
    WRITE(61,*) "        value           uniform 0;"
    WRITE(61,*) '    }'
    
        !DOWN FACE--PROCESSOR_MI1_TO_MI2  !CASE RECTANGLE 3-3-3
    PRO2=PNI*PNJ+PNL*PNI*(MJ2-1)+PNI*(MJ1-2)+MI-1
    WRITE(21,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(21,*) '    {'
    WRITE(21,*) '        type            processor;'
    WRITE(21,*) '        inGroups        List<word> 1(processor);'
    N=CNI*CNL
    WRITE(21,232) "        nFaces          ",N,";"
    N=F_INTERNAL_NUM+CNL*CNL+CNI*CNL*2+CNI*CNL
    WRITE(21,232) "        startFace       ",N,";"
    WRITE(21,*) '        matchTolerance  0.0001;'
    WRITE(21,232) "        myProcNo        ",PRO1,";"
    WRITE(21,232) "        neighbProcNo    ",PRO2,";"
    WRITE(21,*) '    }'
    WRITE(41,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(41,*) '    {'
    WRITE(41,*) '        type            processor;'
    WRITE(41,*) '        value           nonuniform List<scalar> '
    N=CNI*CNL
    WRITE(41,'(I0)') N
    WRITE(41,*) '('
    DO KP=1,CNL
    DO IP=1,CNI
    WRITE(41,'(E21.14)') P_DOWN_RECT(IP,KP)    
    ENDDO
    ENDDO
    WRITE(41,*) ')'
    WRITE(41,*) ';'
    WRITE(41,*) '    }'
    WRITE(51,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(51,*) '    {'
    WRITE(51,*) '        type            processor;'
    WRITE(51,*) '        value           nonuniform List<vector> '
    N=CNI*CNL
    WRITE(51,'(I0)') N
    WRITE(51,*) '('
    DO KP=1,CNL
    DO IP=1,CNI
    UU1=U_DOWN_RECT(1,IP,KP)
    UU2=U_DOWN_RECT(2,IP,KP)
    UU3=U_DOWN_RECT(3,IP,KP)
    WRITE(51,230) '( ',UU1,' ',UU2,' ',UU3,' )'  
    ENDDO
    ENDDO
    WRITE(51,*) ')'
    WRITE(51,*) ';'
    WRITE(51,*) '    }'
    WRITE(61,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(61,*) '    {'
    WRITE(61,*) '        type            processor;'
    WRITE(61,*) "        value           uniform 0;"
    WRITE(61,*) '    }'
    
        !RIGHT FACE--PROCESSOR_MI1_TO_MI2  !CASE RECTANGLE 3-3-3 
    PRO2=PNI*PNJ+PNL*PNI*(MJ2-1)+PNI*(MJ1-1)+MI-2
    WRITE(21,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(21,*) '    {'
    WRITE(21,*) '        type            processor;'
    WRITE(21,*) '        inGroups        List<word> 1(processor);'
    N=CNL*CNL
    WRITE(21,232) "        nFaces          ",N,";"
    N=F_INTERNAL_NUM+CNL*CNL+CNI*CNL*2+CNI*CNL+CNI*CNL
    WRITE(21,232) "        startFace       ",N,";"
    WRITE(21,*) '        matchTolerance  0.0001;'
    WRITE(21,232) "        myProcNo        ",PRO1,";"
    WRITE(21,232) "        neighbProcNo    ",PRO2,";"
    WRITE(21,*) '    }'
    WRITE(41,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(41,*) '    {'
    WRITE(41,*) '        type            processor;'
    WRITE(41,*) '        value           nonuniform List<scalar> '
    N=CNL*CNL
    WRITE(41,'(I0)') N
    WRITE(41,*) '('
    DO KP=1,CNL
    DO JP=1,CNL
    WRITE(41,'(E21.14)') P_RIGHT_RECT(JP,KP)    
    ENDDO
    ENDDO
    WRITE(41,*) ')'
    WRITE(41,*) ';'
    WRITE(41,*) '    }'
    WRITE(51,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(51,*) '    {'
    WRITE(51,*) '        type            processor;'
    WRITE(51,*) '        value           nonuniform List<vector> '
    N=CNL*CNL
    WRITE(51,'(I0)') N
    WRITE(51,*) '('
    DO KP=1,CNL
    DO JP=1,CNL
    UU1=U_RIGHT_RECT(1,JP,KP)
    UU2=U_RIGHT_RECT(2,JP,KP)
    UU3=U_RIGHT_RECT(3,JP,KP)
    WRITE(51,230) '( ',UU1,' ',UU2,' ',UU3,' )'  
    ENDDO
    ENDDO
    WRITE(51,*) ')'
    WRITE(51,*) ';'
    WRITE(51,*) '    }'
    WRITE(61,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(61,*) '    {'
    WRITE(61,*) '        type            processor;'
    WRITE(61,*) "        value           uniform 0;"
    WRITE(61,*) '    }'
    
    
    
    
    
    
    
    ELSE
    !CASE RECTANGLE 3-3-2
        
        !LEFT FACE--HULL  !CASE RECTANGLE 3-3-2 
    IP=CNI
    I=IP+CNI*(MI-1)
    DO KP=1,CNL
    K=KP+CNL*(MJ2-1)
    DO JP=1,CNL
    J=JP+CNL*(MJ1-1)
    F11=(CNL+1)*(CNI+1)*(KP-1)+(CNI+1)*(JP-1)+IP-1
    F12=F11+1
    F13=F12+CNI+1
    F14=F13-1
    F21=F11+(CNL+1)*(CNI+1)
    F22=F21+1
    F23=F22+CNI+1
    F24=F23-1
    CALL CAL_FIJ_RECT(F11,F12,F13,F14,F21,F22,F23,F24)
    WRITE(2,235) '4(',F12,' ',F13,' ',F23,' ',F22,')'
    FAD1=F_INTERNAL_NUM_ALL+NL*NL+NJ*NK+NI*NK+NJ*NK+NL*(K-1)+J
    WRITE(3,231) FAD1
    OWN=CELL_ID_RECT(IP,JP,KP)
    WRITE(4,231) OWN
    ENDDO   
    ENDDO
    
        !UP FACE--PROCESSOR_MI1_TO_MI2  !CASE RECTANGLE 3-3-2
    JP=CNL
    J=JP+CNL*(MJ1-1)
    DO KP=1,CNL
    K=KP+CNL*(MJ2-1)
    DO IP=1,CNI
    I=IP+CNI*(MI-1)
    F11=(CNL+1)*(CNI+1)*(KP-1)+(CNI+1)*(JP-1)+IP-1
    F12=F11+1
    F13=F12+CNI+1
    F14=F13-1
    F21=F11+(CNL+1)*(CNI+1)
    F22=F21+1
    F23=F22+CNI+1
    F24=F23-1
    CALL CAL_FIJ_RECT(F11,F12,F13,F14,F21,F22,F23,F24)
    WRITE(2,235) '4(',F14,' ',F24,' ',F23,' ',F13,')'
    FAD1=F_UP_RECT(IP,KP)
    WRITE(3,231) FAD1
    OWN=CELL_ID_RECT(IP,JP,KP)
    WRITE(4,231) OWN
    ENDDO   
    ENDDO
    
        !FRONT FACE--PROCESSOR_MI1_TO_MI2  !CASE RECTANGLE 3-3-2
    KP=1
    K=KP+CNL*(MJ2-1)
    DO JP=1,CNL
    J=JP+CNL*(MJ1-1)
    DO IP=1,CNI
    I=IP+CNI*(MI-1)
    F11=(CNL+1)*(CNI+1)*(KP-1)+(CNI+1)*(JP-1)+IP-1
    F12=F11+1
    F13=F12+CNI+1
    F14=F13-1
    F21=F11+(CNL+1)*(CNI+1)
    F22=F21+1
    F23=F22+CNI+1
    F24=F23-1
    CALL CAL_FIJ_RECT(F11,F12,F13,F14,F21,F22,F23,F24)
    WRITE(2,235) '4(',F11,' ',F14,' ',F13,' ',F12,')'
    FAD1=F_FRONT_RECT(IP,JP)
    WRITE(3,231) FAD1
    OWN=CELL_ID_RECT(IP,JP,KP)
    WRITE(4,231) OWN
    ENDDO   
    ENDDO
    
        !DOWN FACE--PROCESSOR_MI1_TO_MI2  !CASE RECTANGLE 3-3-2
    JP=1
    J=JP+CNL*(MJ1-1)
    DO KP=1,CNL
    K=KP+CNL*(MJ2-1)
    DO IP=1,CNI
    I=IP+CNI*(MI-1)
    F11=(CNL+1)*(CNI+1)*(KP-1)+(CNI+1)*(JP-1)+IP-1
    F12=F11+1
    F13=F12+CNI+1
    F14=F13-1
    F21=F11+(CNL+1)*(CNI+1)
    F22=F21+1
    F23=F22+CNI+1
    F24=F23-1
    CALL CAL_FIJ_RECT(F11,F12,F13,F14,F21,F22,F23,F24)
    WRITE(2,235) '4(',F11,' ',F12,' ',F22,' ',F21,')'
    FAD1=F_DOWN_RECT(IP,KP)
    WRITE(3,231) FAD1
    OWN=CELL_ID_RECT(IP,JP,KP)
    WRITE(4,231) OWN
    ENDDO   
    ENDDO
    
        !RIGHT FACE--PROCESSOR_MI1_TO_MI2  !CASE RECTANGLE 3-3-2 
    IP=1
    I=IP+CNI*(MI-1)
    DO KP=1,CNL
    K=KP+CNL*(MJ2-1)
    DO JP=1,CNL
    J=JP+CNL*(MJ1-1)
    F11=(CNL+1)*(CNI+1)*(KP-1)+(CNI+1)*(JP-1)+IP-1
    F12=F11+1
    F13=F12+CNI+1
    F14=F13-1
    F21=F11+(CNL+1)*(CNI+1)
    F22=F21+1
    F23=F22+CNI+1
    F24=F23-1
    CALL CAL_FIJ_RECT(F11,F12,F13,F14,F21,F22,F23,F24)
    WRITE(2,235) '4(',F11,' ',F21,' ',F24,' ',F14,')'
    FAD1=F_RIGHT_RECT(JP,KP)
    WRITE(3,231) FAD1
    OWN=CELL_ID_RECT(IP,JP,KP)
    WRITE(4,231) OWN
    ENDDO   
    ENDDO
    
        !BACK FACE--PROCESSOR_MI1_TO_MI2  !CASE RECTANGLE 3-3-2
    KP=CNL
    K=KP+CNL*(MJ2-1)
    DO JP=1,CNL
    J=JP+CNL*(MJ1-1)
    DO IP=1,CNI
    I=IP+CNI*(MI-1)
    F11=(CNL+1)*(CNI+1)*(KP-1)+(CNI+1)*(JP-1)+IP-1
    F12=F11+1
    F13=F12+CNI+1
    F14=F13-1
    F21=F11+(CNL+1)*(CNI+1)
    F22=F21+1
    F23=F22+CNI+1
    F24=F23-1
    CALL CAL_FIJ_RECT(F11,F12,F13,F14,F21,F22,F23,F24)
    WRITE(2,235) '4(',F21,' ',F22,' ',F23,' ',F24,')'
    FAD1=F_BACK_RECT(IP,JP)
    WRITE(3,231) FAD1
    OWN=CELL_ID_RECT(IP,JP,KP)
    WRITE(4,231) OWN
    ENDDO   
    ENDDO
    
    N=5
    WRITE(31,231) N+4
    WRITE(31,*) '('  
    WRITE(31,*) '0'
    WRITE(31,*) '1'
    WRITE(31,*) '2' 
    WRITE(31,*) '3' 
    DO I=1,N
    WRITE(31,*) '-1' 
    ENDDO
    WRITE(31,*) ')'
    
    WRITE(21,231) N+4
    WRITE(21,*) "("    
    WRITE(21,*) "    INLET1"
    WRITE(21,*) "    {"    
    WRITE(21,*) "        type            patch;"
    N=0
    WRITE(21,232) "        nFaces          ",N,";"
    N=F_INTERNAL_NUM
    WRITE(21,232) "        startFace       ",N,";"
    WRITE(21,*) "    }"
    WRITE(21,*) "    INLET2"
    WRITE(21,*) "    {"    
    WRITE(21,*) "        type            patch;"
    N=0
    WRITE(21,232) "        nFaces          ",N,";"
    N=F_INTERNAL_NUM
    WRITE(21,232) "        startFace       ",N,";"
    WRITE(21,*) "    }"
    WRITE(21,*) "    OUTLET"
    WRITE(21,*) "    {"
    WRITE(21,*) "        type            patch;"
    N=0
    WRITE(21,232) "        nFaces          ",N,";"
    N=F_INTERNAL_NUM
    WRITE(21,232) "        startFace       ",N,";"
    WRITE(21,*) "    }"
    WRITE(21,*) "    HULL"
    WRITE(21,*) "    {"
    WRITE(21,*) "        type            wall;"
    WRITE(21,*) "        inGroups        1(wall);"
    N=CNL*CNL
    WRITE(21,232) "        nFaces          ",N,";"
    N=F_INTERNAL_NUM
    WRITE(21,232) "        startFace       ",N,";"
    WRITE(21,*) "    }"
    
        !UP FACE--PROCESSOR_MI1_TO_MI2  !CASE RECTANGLE 3-3-2
    PRO1=PNI*PNJ+PNL*PNI*(MJ2-1)+PNI*(MJ1-1)+MI-1
    MJ=1
    PRO2=PNI*(MJ-1)+MI-1
    WRITE(21,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(21,*) '    {'
    WRITE(21,*) '        type            processor;'
    WRITE(21,*) '        inGroups        List<word> 1(processor);'
    N=CNI*CNL
    WRITE(21,232) "        nFaces          ",N,";"
    N=F_INTERNAL_NUM+CNL*CNL
    WRITE(21,232) "        startFace       ",N,";"
    WRITE(21,*) '        matchTolerance  0.0001;'
    WRITE(21,232) "        myProcNo        ",PRO1,";"
    WRITE(21,232) "        neighbProcNo    ",PRO2,";"
    WRITE(21,*) '    }'
    WRITE(41,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(41,*) '    {'
    WRITE(41,*) '        type            processor;'
    WRITE(41,*) '        value           nonuniform List<scalar> '
    N=CNI*CNL
    WRITE(41,'(I0)') N
    WRITE(41,*) '('
    DO KP=1,CNL
    DO IP=1,CNI
    WRITE(41,'(E21.14)') P_UP_RECT(IP,KP)    
    ENDDO
    ENDDO
    WRITE(41,*) ')'
    WRITE(41,*) ';'
    WRITE(41,*) '    }'
    WRITE(51,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(51,*) '    {'
    WRITE(51,*) '        type            processor;'
    WRITE(51,*) '        value           nonuniform List<vector> '
    N=CNI*CNL
    WRITE(51,'(I0)') N
    WRITE(51,*) '('
    DO KP=1,CNL
    DO IP=1,CNI
    UU1=U_UP_RECT(1,IP,KP)
    UU2=U_UP_RECT(2,IP,KP)
    UU3=U_UP_RECT(3,IP,KP)
    WRITE(51,230) '( ',UU1,' ',UU2,' ',UU3,' )'  
    ENDDO
    ENDDO
    WRITE(51,*) ')'
    WRITE(51,*) ';'
    WRITE(51,*) '    }'
    WRITE(61,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(61,*) '    {'
    WRITE(61,*) '        type            processor;'
    WRITE(61,*) "        value           uniform 0;"
    WRITE(61,*) '    }'
    
        !FRONT FACE--PROCESSOR_MI1_TO_MI2  !CASE RECTANGLE 3-3-2
    PRO2=PNI*PNJ+PNL*PNI*(MJ2-2)+PNI*(MJ1-1)+MI-1
    WRITE(21,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(21,*) '    {'
    WRITE(21,*) '        type            processor;'
    WRITE(21,*) '        inGroups        List<word> 1(processor);'
    N=CNI*CNL
    WRITE(21,232) "        nFaces          ",N,";"
    N=F_INTERNAL_NUM+CNL*CNL+CNI*CNL
    WRITE(21,232) "        startFace       ",N,";"
    WRITE(21,*) '        matchTolerance  0.0001;'
    WRITE(21,232) "        myProcNo        ",PRO1,";"
    WRITE(21,232) "        neighbProcNo    ",PRO2,";"
    WRITE(21,*) '    }'
    WRITE(41,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(41,*) '    {'
    WRITE(41,*) '        type            processor;'
    WRITE(41,*) '        value           nonuniform List<scalar> '
    N=CNI*CNL
    WRITE(41,'(I0)') N
    WRITE(41,*) '('
    DO JP=1,CNL
    DO IP=1,CNI
    WRITE(41,'(E21.14)') P_FRONT_RECT(IP,JP)    
    ENDDO
    ENDDO
    WRITE(41,*) ')'
    WRITE(41,*) ';'
    WRITE(41,*) '    }'
    WRITE(51,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(51,*) '    {'
    WRITE(51,*) '        type            processor;'
    WRITE(51,*) '        value           nonuniform List<vector> '
    N=CNI*CNL
    WRITE(51,'(I0)') N
    WRITE(51,*) '('
    DO JP=1,CNL
    DO IP=1,CNI
    UU1=U_FRONT_RECT(1,IP,JP)
    UU2=U_FRONT_RECT(2,IP,JP)
    UU3=U_FRONT_RECT(3,IP,JP)
    WRITE(51,230) '( ',UU1,' ',UU2,' ',UU3,' )'  
    ENDDO
    ENDDO
    WRITE(51,*) ')'
    WRITE(51,*) ';'
    WRITE(51,*) '    }'
    WRITE(61,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(61,*) '    {'
    WRITE(61,*) '        type            processor;'
    WRITE(61,*) "        value           uniform 0;"
    WRITE(61,*) '    }'
    
        !DOWN FACE--PROCESSOR_MI1_TO_MI2  !CASE RECTANGLE 3-3-2
    PRO2=PNI*PNJ+PNL*PNI*(MJ2-1)+PNI*(MJ1-2)+MI-1
    WRITE(21,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(21,*) '    {'
    WRITE(21,*) '        type            processor;'
    WRITE(21,*) '        inGroups        List<word> 1(processor);'
    N=CNI*CNL
    WRITE(21,232) "        nFaces          ",N,";"
    N=F_INTERNAL_NUM+CNL*CNL+CNI*CNL+CNI*CNL
    WRITE(21,232) "        startFace       ",N,";"
    WRITE(21,*) '        matchTolerance  0.0001;'
    WRITE(21,232) "        myProcNo        ",PRO1,";"
    WRITE(21,232) "        neighbProcNo    ",PRO2,";"
    WRITE(21,*) '    }'
    WRITE(41,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(41,*) '    {'
    WRITE(41,*) '        type            processor;'
    WRITE(41,*) '        value           nonuniform List<scalar> '
    N=CNI*CNL
    WRITE(41,'(I0)') N
    WRITE(41,*) '('
    DO KP=1,CNL
    DO IP=1,CNI
    WRITE(41,'(E21.14)') P_DOWN_RECT(IP,KP)    
    ENDDO
    ENDDO
    WRITE(41,*) ')'
    WRITE(41,*) ';'
    WRITE(41,*) '    }'
    WRITE(51,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(51,*) '    {'
    WRITE(51,*) '        type            processor;'
    WRITE(51,*) '        value           nonuniform List<vector> '
    N=CNI*CNL
    WRITE(51,'(I0)') N
    WRITE(51,*) '('
    DO KP=1,CNL
    DO IP=1,CNI
    UU1=U_DOWN_RECT(1,IP,KP)
    UU2=U_DOWN_RECT(2,IP,KP)
    UU3=U_DOWN_RECT(3,IP,KP)
    WRITE(51,230) '( ',UU1,' ',UU2,' ',UU3,' )'  
    ENDDO
    ENDDO
    WRITE(51,*) ')'
    WRITE(51,*) ';'
    WRITE(51,*) '    }'
    WRITE(61,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(61,*) '    {'
    WRITE(61,*) '        type            processor;'
    WRITE(61,*) "        value           uniform 0;"
    WRITE(61,*) '    }'
    
        !RIGHT FACE--PROCESSOR_MI1_TO_MI2  !CASE RECTANGLE 3-3-2 
    PRO2=PNI*PNJ+PNL*PNI*(MJ2-1)+PNI*(MJ1-1)+MI-2
    WRITE(21,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(21,*) '    {'
    WRITE(21,*) '        type            processor;'
    WRITE(21,*) '        inGroups        List<word> 1(processor);'
    N=CNL*CNL
    WRITE(21,232) "        nFaces          ",N,";"
    N=F_INTERNAL_NUM+CNL*CNL+CNI*CNL+CNI*CNL+CNI*CNL
    WRITE(21,232) "        startFace       ",N,";"
    WRITE(21,*) '        matchTolerance  0.0001;'
    WRITE(21,232) "        myProcNo        ",PRO1,";"
    WRITE(21,232) "        neighbProcNo    ",PRO2,";"
    WRITE(21,*) '    }'
    WRITE(41,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(41,*) '    {'
    WRITE(41,*) '        type            processor;'
    WRITE(41,*) '        value           nonuniform List<scalar> '
    N=CNL*CNL
    WRITE(41,'(I0)') N
    WRITE(41,*) '('
    DO KP=1,CNL
    DO JP=1,CNL
    WRITE(41,'(E21.14)') P_RIGHT_RECT(JP,KP)    
    ENDDO
    ENDDO
    WRITE(41,*) ')'
    WRITE(41,*) ';'
    WRITE(41,*) '    }'
    WRITE(51,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(51,*) '    {'
    WRITE(51,*) '        type            processor;'
    WRITE(51,*) '        value           nonuniform List<vector> '
    N=CNL*CNL
    WRITE(51,'(I0)') N
    WRITE(51,*) '('
    DO KP=1,CNL
    DO JP=1,CNL
    UU1=U_RIGHT_RECT(1,JP,KP)
    UU2=U_RIGHT_RECT(2,JP,KP)
    UU3=U_RIGHT_RECT(3,JP,KP)
    WRITE(51,230) '( ',UU1,' ',UU2,' ',UU3,' )'  
    ENDDO
    ENDDO
    WRITE(51,*) ')'
    WRITE(51,*) ';'
    WRITE(51,*) '    }'
    WRITE(61,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(61,*) '    {'
    WRITE(61,*) '        type            processor;'
    WRITE(61,*) "        value           uniform 0;"
    WRITE(61,*) '    }' 
    
        !BACK FACE--PROCESSOR_MI1_TO_MI2  !CASE RECTANGLE 3-3-2
    PRO2=PNI*PNJ+PNL*PNI*MJ2+PNI*(MJ1-1)+MI-1
    WRITE(21,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(21,*) '    {'
    WRITE(21,*) '        type            processor;'
    WRITE(21,*) '        inGroups        List<word> 1(processor);'
    N=CNI*CNL
    WRITE(21,232) "        nFaces          ",N,";"
    N=F_INTERNAL_NUM+CNL*CNL+CNI*CNL+CNI*CNL+CNI*CNL+CNL*CNL
    WRITE(21,232) "        startFace       ",N,";"
    WRITE(21,*) '        matchTolerance  0.0001;'
    WRITE(21,232) "        myProcNo        ",PRO1,";"
    WRITE(21,232) "        neighbProcNo    ",PRO2,";"
    WRITE(21,*) '    }'
    WRITE(41,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(41,*) '    {'
    WRITE(41,*) '        type            processor;'
    WRITE(41,*) '        value           nonuniform List<scalar> '
    N=CNI*CNL
    WRITE(41,'(I0)') N
    WRITE(41,*) '('
    DO JP=1,CNL
    DO IP=1,CNI
    WRITE(41,'(E21.14)') P_BACK_RECT(IP,JP)    
    ENDDO
    ENDDO
    WRITE(41,*) ')'
    WRITE(41,*) ';'
    WRITE(41,*) '    }'
    WRITE(51,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(51,*) '    {'
    WRITE(51,*) '        type            processor;'
    WRITE(51,*) '        value           nonuniform List<vector> '
    N=CNI*CNL
    WRITE(51,'(I0)') N
    WRITE(51,*) '('
    DO JP=1,CNL
    DO IP=1,CNI
    UU1=U_BACK_RECT(1,IP,JP)
    UU2=U_BACK_RECT(2,IP,JP)
    UU3=U_BACK_RECT(3,IP,JP)
    WRITE(51,230) '( ',UU1,' ',UU2,' ',UU3,' )'  
    ENDDO
    ENDDO
    WRITE(51,*) ')'
    WRITE(51,*) ';'
    WRITE(51,*) '    }'
    WRITE(61,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(61,*) '    {'
    WRITE(61,*) '        type            processor;'
    WRITE(61,*) "        value           uniform 0;"
    WRITE(61,*) '    }'
    
    
    
    
    
    
    
    ENDIF 
    ELSE
    IF (MJ2.EQ.1) THEN
    !CASE RECTANGLE 3-2-1
        
        !LEFT FACE--HULL  !CASE RECTANGLE 3-2-1 
    IP=CNI
    I=IP+CNI*(MI-1)
    DO KP=1,CNL
    K=KP+CNL*(MJ2-1)
    DO JP=1,CNL
    J=JP+CNL*(MJ1-1)
    F11=(CNL+1)*(CNI+1)*(KP-1)+(CNI+1)*(JP-1)+IP-1
    F12=F11+1
    F13=F12+CNI+1
    F14=F13-1
    F21=F11+(CNL+1)*(CNI+1)
    F22=F21+1
    F23=F22+CNI+1
    F24=F23-1
    CALL CAL_FIJ_RECT(F11,F12,F13,F14,F21,F22,F23,F24)
    WRITE(2,235) '4(',F12,' ',F13,' ',F23,' ',F22,')'
    FAD1=F_INTERNAL_NUM_ALL+NL*NL+NJ*NK+NI*NK+NJ*NK+NL*(K-1)+J
    WRITE(3,231) FAD1
    OWN=CELL_ID_RECT(IP,JP,KP)
    WRITE(4,231) OWN
    ENDDO   
    ENDDO
    
        !FRONT FACE--PROCESSOR_MI1_TO_MI2  !CASE RECTANGLE 3-2-1
    KP=1
    K=KP+CNL*(MJ2-1)
    DO JP=1,CNL
    J=JP+CNL*(MJ1-1)
    DO IP=1,CNI
    I=IP+CNI*(MI-1)
    F11=(CNL+1)*(CNI+1)*(KP-1)+(CNI+1)*(JP-1)+IP-1
    F12=F11+1
    F13=F12+CNI+1
    F14=F13-1
    F21=F11+(CNL+1)*(CNI+1)
    F22=F21+1
    F23=F22+CNI+1
    F24=F23-1
    CALL CAL_FIJ_RECT(F11,F12,F13,F14,F21,F22,F23,F24)
    WRITE(2,235) '4(',F11,' ',F14,' ',F13,' ',F12,')'
    FAD1=F_FRONT_RECT(IP,JP)
    WRITE(3,231) FAD1
    OWN=CELL_ID_RECT(IP,JP,KP)
    WRITE(4,231) OWN
    ENDDO   
    ENDDO
    
        !DOWN FACE--PROCESSOR_MI1_TO_MI2  !CASE RECTANGLE 3-2-1
    JP=1
    J=JP+CNL*(MJ1-1)
    DO KP=1,CNL
    K=KP+CNL*(MJ2-1)
    DO IP=1,CNI
    I=IP+CNI*(MI-1)
    F11=(CNL+1)*(CNI+1)*(KP-1)+(CNI+1)*(JP-1)+IP-1
    F12=F11+1
    F13=F12+CNI+1
    F14=F13-1
    F21=F11+(CNL+1)*(CNI+1)
    F22=F21+1
    F23=F22+CNI+1
    F24=F23-1
    CALL CAL_FIJ_RECT(F11,F12,F13,F14,F21,F22,F23,F24)
    WRITE(2,235) '4(',F11,' ',F12,' ',F22,' ',F21,')'
    FAD1=F_DOWN_RECT(IP,KP)
    WRITE(3,231) FAD1
    OWN=CELL_ID_RECT(IP,JP,KP)
    WRITE(4,231) OWN
    ENDDO   
    ENDDO
    
        !RIGHT FACE--PROCESSOR_MI1_TO_MI2  !CASE RECTANGLE 3-2-1 
    IP=1
    I=IP+CNI*(MI-1)
    DO KP=1,CNL
    K=KP+CNL*(MJ2-1)
    DO JP=1,CNL
    J=JP+CNL*(MJ1-1)
    F11=(CNL+1)*(CNI+1)*(KP-1)+(CNI+1)*(JP-1)+IP-1
    F12=F11+1
    F13=F12+CNI+1
    F14=F13-1
    F21=F11+(CNL+1)*(CNI+1)
    F22=F21+1
    F23=F22+CNI+1
    F24=F23-1
    CALL CAL_FIJ_RECT(F11,F12,F13,F14,F21,F22,F23,F24)
    WRITE(2,235) '4(',F11,' ',F21,' ',F24,' ',F14,')'
    FAD1=F_RIGHT_RECT(JP,KP)
    WRITE(3,231) FAD1
    OWN=CELL_ID_RECT(IP,JP,KP)
    WRITE(4,231) OWN
    ENDDO   
    ENDDO
    
        !UP FACE--PROCESSOR_MI1_TO_MI2  !CASE RECTANGLE 3-2-1
    JP=CNL
    J=JP+CNL*(MJ1-1)
    DO KP=1,CNL
    K=KP+CNL*(MJ2-1)
    DO IP=1,CNI
    I=IP+CNI*(MI-1)
    F11=(CNL+1)*(CNI+1)*(KP-1)+(CNI+1)*(JP-1)+IP-1
    F12=F11+1
    F13=F12+CNI+1
    F14=F13-1
    F21=F11+(CNL+1)*(CNI+1)
    F22=F21+1
    F23=F22+CNI+1
    F24=F23-1
    CALL CAL_FIJ_RECT(F11,F12,F13,F14,F21,F22,F23,F24)
    WRITE(2,235) '4(',F14,' ',F24,' ',F23,' ',F13,')'
    FAD1=F_UP_RECT(IP,KP)
    WRITE(3,231) FAD1
    OWN=CELL_ID_RECT(IP,JP,KP)
    WRITE(4,231) OWN
    ENDDO   
    ENDDO
    
        !BACK FACE--PROCESSOR_MI1_TO_MI2  !CASE RECTANGLE 3-2-1
    KP=CNL
    K=KP+CNL*(MJ2-1)
    DO JP=1,CNL
    J=JP+CNL*(MJ1-1)
    DO IP=1,CNI
    I=IP+CNI*(MI-1)
    F11=(CNL+1)*(CNI+1)*(KP-1)+(CNI+1)*(JP-1)+IP-1
    F12=F11+1
    F13=F12+CNI+1
    F14=F13-1
    F21=F11+(CNL+1)*(CNI+1)
    F22=F21+1
    F23=F22+CNI+1
    F24=F23-1
    CALL CAL_FIJ_RECT(F11,F12,F13,F14,F21,F22,F23,F24)
    WRITE(2,235) '4(',F21,' ',F22,' ',F23,' ',F24,')'
    FAD1=F_BACK_RECT(IP,JP)
    WRITE(3,231) FAD1
    OWN=CELL_ID_RECT(IP,JP,KP)
    WRITE(4,231) OWN
    ENDDO   
    ENDDO
    
    N=5
    WRITE(31,231) N+4
    WRITE(31,*) '('  
    WRITE(31,*) '0'
    WRITE(31,*) '1'
    WRITE(31,*) '2' 
    WRITE(31,*) '3' 
    DO I=1,N
    WRITE(31,*) '-1' 
    ENDDO
    WRITE(31,*) ')'
    
    WRITE(21,231) N+4
    WRITE(21,*) "("    
    WRITE(21,*) "    INLET1"
    WRITE(21,*) "    {"    
    WRITE(21,*) "        type            patch;"
    N=0
    WRITE(21,232) "        nFaces          ",N,";"
    N=F_INTERNAL_NUM
    WRITE(21,232) "        startFace       ",N,";"
    WRITE(21,*) "    }"
    WRITE(21,*) "    INLET2"
    WRITE(21,*) "    {"    
    WRITE(21,*) "        type            patch;"
    N=0
    WRITE(21,232) "        nFaces          ",N,";"
    N=F_INTERNAL_NUM
    WRITE(21,232) "        startFace       ",N,";"
    WRITE(21,*) "    }"
    WRITE(21,*) "    OUTLET"
    WRITE(21,*) "    {"
    WRITE(21,*) "        type            patch;"
    N=0
    WRITE(21,232) "        nFaces          ",N,";"
    N=F_INTERNAL_NUM
    WRITE(21,232) "        startFace       ",N,";"
    WRITE(21,*) "    }"
    WRITE(21,*) "    HULL"
    WRITE(21,*) "    {"
    WRITE(21,*) "        type            wall;"
    WRITE(21,*) "        inGroups        1(wall);"
    N=CNL*CNL
    WRITE(21,232) "        nFaces          ",N,";"
    N=F_INTERNAL_NUM
    WRITE(21,232) "        startFace       ",N,";"
    WRITE(21,*) "    }"
    
        !FRONT FACE--PROCESSOR_MI1_TO_MI2  !CASE RECTANGLE 3-2-1
    PRO1=PNI*PNJ+PNL*PNI*(MJ2-1)+PNI*(MJ1-1)+MI-1
    MJ=1
    PRO2=PNI*(MJ-1)+MI-1
    WRITE(21,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(21,*) '    {'
    WRITE(21,*) '        type            processor;'
    WRITE(21,*) '        inGroups        List<word> 1(processor);'
    N=CNI*CNL
    WRITE(21,232) "        nFaces          ",N,";"
    N=F_INTERNAL_NUM+CNL*CNL
    WRITE(21,232) "        startFace       ",N,";"
    WRITE(21,*) '        matchTolerance  0.0001;'
    WRITE(21,232) "        myProcNo        ",PRO1,";"
    WRITE(21,232) "        neighbProcNo    ",PRO2,";"
    WRITE(21,*) '    }'
    WRITE(41,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(41,*) '    {'
    WRITE(41,*) '        type            processor;'
    WRITE(41,*) '        value           nonuniform List<scalar> '
    N=CNI*CNL
    WRITE(41,'(I0)') N
    WRITE(41,*) '('
    DO JP=1,CNL
    DO IP=1,CNI
    WRITE(41,'(E21.14)') P_FRONT_RECT(IP,JP)    
    ENDDO
    ENDDO
    WRITE(41,*) ')'
    WRITE(41,*) ';'
    WRITE(41,*) '    }'
    WRITE(51,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(51,*) '    {'
    WRITE(51,*) '        type            processor;'
    WRITE(51,*) '        value           nonuniform List<vector> '
    N=CNI*CNL
    WRITE(51,'(I0)') N
    WRITE(51,*) '('
    DO JP=1,CNL
    DO IP=1,CNI
    UU1=U_FRONT_RECT(1,IP,JP)
    UU2=U_FRONT_RECT(2,IP,JP)
    UU3=U_FRONT_RECT(3,IP,JP)
    WRITE(51,230) '( ',UU1,' ',UU2,' ',UU3,' )'  
    ENDDO
    ENDDO
    WRITE(51,*) ')'
    WRITE(51,*) ';'
    WRITE(51,*) '    }'
    WRITE(61,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(61,*) '    {'
    WRITE(61,*) '        type            processor;'
    WRITE(61,*) "        value           uniform 0;"
    WRITE(61,*) '    }'
    
        !DOWN FACE--PROCESSOR_MI1_TO_MI2  !CASE RECTANGLE 3-2-1
    PRO2=PNI*PNJ+PNL*PNI*(MJ2-1)+PNI*(MJ1-2)+MI-1
    WRITE(21,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(21,*) '    {'
    WRITE(21,*) '        type            processor;'
    WRITE(21,*) '        inGroups        List<word> 1(processor);'
    N=CNI*CNL
    WRITE(21,232) "        nFaces          ",N,";"
    N=F_INTERNAL_NUM+CNL*CNL+CNI*CNL
    WRITE(21,232) "        startFace       ",N,";"
    WRITE(21,*) '        matchTolerance  0.0001;'
    WRITE(21,232) "        myProcNo        ",PRO1,";"
    WRITE(21,232) "        neighbProcNo    ",PRO2,";"
    WRITE(21,*) '    }'
    WRITE(41,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(41,*) '    {'
    WRITE(41,*) '        type            processor;'
    WRITE(41,*) '        value           nonuniform List<scalar> '
    N=CNI*CNL
    WRITE(41,'(I0)') N
    WRITE(41,*) '('
    DO KP=1,CNL
    DO IP=1,CNI
    WRITE(41,'(E21.14)') P_DOWN_RECT(IP,KP)    
    ENDDO
    ENDDO
    WRITE(41,*) ')'
    WRITE(41,*) ';'
    WRITE(41,*) '    }'
    WRITE(51,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(51,*) '    {'
    WRITE(51,*) '        type            processor;'
    WRITE(51,*) '        value           nonuniform List<vector> '
    N=CNI*CNL
    WRITE(51,'(I0)') N
    WRITE(51,*) '('
    DO KP=1,CNL
    DO IP=1,CNI
    UU1=U_DOWN_RECT(1,IP,KP)
    UU2=U_DOWN_RECT(2,IP,KP)
    UU3=U_DOWN_RECT(3,IP,KP)
    WRITE(51,230) '( ',UU1,' ',UU2,' ',UU3,' )'  
    ENDDO
    ENDDO
    WRITE(51,*) ')'
    WRITE(51,*) ';'
    WRITE(51,*) '    }'
    WRITE(61,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(61,*) '    {'
    WRITE(61,*) '        type            processor;'
    WRITE(61,*) "        value           uniform 0;"
    WRITE(61,*) '    }'
    
        !RIGHT FACE--PROCESSOR_MI1_TO_MI2  !CASE RECTANGLE 3-2-1
    PRO2=PNI*PNJ+PNL*PNI*(MJ2-1)+PNI*(MJ1-1)+MI-2
    WRITE(21,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(21,*) '    {'
    WRITE(21,*) '        type            processor;'
    WRITE(21,*) '        inGroups        List<word> 1(processor);'
    N=CNL*CNL
    WRITE(21,232) "        nFaces          ",N,";"
    N=F_INTERNAL_NUM+CNL*CNL+CNI*CNL+CNI*CNL
    WRITE(21,232) "        startFace       ",N,";"
    WRITE(21,*) '        matchTolerance  0.0001;'
    WRITE(21,232) "        myProcNo        ",PRO1,";"
    WRITE(21,232) "        neighbProcNo    ",PRO2,";"
    WRITE(21,*) '    }'
    WRITE(41,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(41,*) '    {'
    WRITE(41,*) '        type            processor;'
    WRITE(41,*) '        value           nonuniform List<scalar> '
    N=CNL*CNL
    WRITE(41,'(I0)') N
    WRITE(41,*) '('
    DO KP=1,CNL
    DO JP=1,CNL
    WRITE(41,'(E21.14)') P_RIGHT_RECT(JP,KP)    
    ENDDO
    ENDDO
    WRITE(41,*) ')'
    WRITE(41,*) ';'
    WRITE(41,*) '    }'
    WRITE(51,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(51,*) '    {'
    WRITE(51,*) '        type            processor;'
    WRITE(51,*) '        value           nonuniform List<vector> '
    N=CNL*CNL
    WRITE(51,'(I0)') N
    WRITE(51,*) '('
    DO KP=1,CNL
    DO JP=1,CNL
    UU1=U_RIGHT_RECT(1,JP,KP)
    UU2=U_RIGHT_RECT(2,JP,KP)
    UU3=U_RIGHT_RECT(3,JP,KP)
    WRITE(51,230) '( ',UU1,' ',UU2,' ',UU3,' )'  
    ENDDO
    ENDDO
    WRITE(51,*) ')'
    WRITE(51,*) ';'
    WRITE(51,*) '    }'
    WRITE(61,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(61,*) '    {'
    WRITE(61,*) '        type            processor;'
    WRITE(61,*) "        value           uniform 0;"
    WRITE(61,*) '    }'
    
        !UP FACE--PROCESSOR_MI1_TO_MI2  !CASE RECTANGLE 3-2-1
    PRO2=PNI*PNJ+PNL*PNI*(MJ2-1)+PNI*MJ1+MI-1
    WRITE(21,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(21,*) '    {'
    WRITE(21,*) '        type            processor;'
    WRITE(21,*) '        inGroups        List<word> 1(processor);'
    N=CNI*CNL
    WRITE(21,232) "        nFaces          ",N,";"
    N=F_INTERNAL_NUM+CNL*CNL+CNI*CNL+CNI*CNL+CNL*CNL
    WRITE(21,232) "        startFace       ",N,";"
    WRITE(21,*) '        matchTolerance  0.0001;'
    WRITE(21,232) "        myProcNo        ",PRO1,";"
    WRITE(21,232) "        neighbProcNo    ",PRO2,";"
    WRITE(21,*) '    }'
    WRITE(41,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(41,*) '    {'
    WRITE(41,*) '        type            processor;'
    WRITE(41,*) '        value           nonuniform List<scalar> '
    N=CNI*CNL
    WRITE(41,'(I0)') N
    WRITE(41,*) '('
    DO KP=1,CNL
    DO IP=1,CNI
    WRITE(41,'(E21.14)') P_UP_RECT(IP,KP)    
    ENDDO
    ENDDO
    WRITE(41,*) ')'
    WRITE(41,*) ';'
    WRITE(41,*) '    }'
    WRITE(51,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(51,*) '    {'
    WRITE(51,*) '        type            processor;'
    WRITE(51,*) '        value           nonuniform List<vector> '
    N=CNI*CNL
    WRITE(51,'(I0)') N
    WRITE(51,*) '('
    DO KP=1,CNL
    DO IP=1,CNI
    UU1=U_UP_RECT(1,IP,KP)
    UU2=U_UP_RECT(2,IP,KP)
    UU3=U_UP_RECT(3,IP,KP)
    WRITE(51,230) '( ',UU1,' ',UU2,' ',UU3,' )'  
    ENDDO
    ENDDO
    WRITE(51,*) ')'
    WRITE(51,*) ';'
    WRITE(51,*) '    }'
    WRITE(61,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(61,*) '    {'
    WRITE(61,*) '        type            processor;'
    WRITE(61,*) "        value           uniform 0;"
    WRITE(61,*) '    }'
    
        !BACK FACE--PROCESSOR_MI1_TO_MI2  !CASE RECTANGLE 3-2-1
    PRO2=PNI*PNJ+PNL*PNI*MJ2+PNI*(MJ1-1)+MI-1
    WRITE(21,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(21,*) '    {'
    WRITE(21,*) '        type            processor;'
    WRITE(21,*) '        inGroups        List<word> 1(processor);'
    N=CNI*CNL
    WRITE(21,232) "        nFaces          ",N,";"
    N=F_INTERNAL_NUM+CNL*CNL+CNI*CNL+CNI*CNL+CNL*CNL+CNI*CNL
    WRITE(21,232) "        startFace       ",N,";"
    WRITE(21,*) '        matchTolerance  0.0001;'
    WRITE(21,232) "        myProcNo        ",PRO1,";"
    WRITE(21,232) "        neighbProcNo    ",PRO2,";"
    WRITE(21,*) '    }'
    WRITE(41,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(41,*) '    {'
    WRITE(41,*) '        type            processor;'
    WRITE(41,*) '        value           nonuniform List<scalar> '
    N=CNI*CNL
    WRITE(41,'(I0)') N
    WRITE(41,*) '('
    DO JP=1,CNL
    DO IP=1,CNI
    WRITE(41,'(E21.14)') P_BACK_RECT(IP,JP)    
    ENDDO
    ENDDO
    WRITE(41,*) ')'
    WRITE(41,*) ';'
    WRITE(41,*) '    }'
    WRITE(51,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(51,*) '    {'
    WRITE(51,*) '        type            processor;'
    WRITE(51,*) '        value           nonuniform List<vector> '
    N=CNI*CNL
    WRITE(51,'(I0)') N
    WRITE(51,*) '('
    DO JP=1,CNL
    DO IP=1,CNI
    UU1=U_BACK_RECT(1,IP,JP)
    UU2=U_BACK_RECT(2,IP,JP)
    UU3=U_BACK_RECT(3,IP,JP)
    WRITE(51,230) '( ',UU1,' ',UU2,' ',UU3,' )'  
    ENDDO
    ENDDO
    WRITE(51,*) ')'
    WRITE(51,*) ';'
    WRITE(51,*) '    }'
    WRITE(61,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(61,*) '    {'
    WRITE(61,*) '        type            processor;'
    WRITE(61,*) "        value           uniform 0;"
    WRITE(61,*) '    }'
    
    
    
    
    
    
    
    
    
    ELSEIF (MJ2.EQ.PNL) THEN
    !CASE RECTANGLE 3-2-3
        
        !LEFT FACE--HULL  !CASE RECTANGLE 3-2-3 
    IP=CNI
    I=IP+CNI*(MI-1)
    DO KP=1,CNL
    K=KP+CNL*(MJ2-1)
    DO JP=1,CNL
    J=JP+CNL*(MJ1-1)
    F11=(CNL+1)*(CNI+1)*(KP-1)+(CNI+1)*(JP-1)+IP-1
    F12=F11+1
    F13=F12+CNI+1
    F14=F13-1
    F21=F11+(CNL+1)*(CNI+1)
    F22=F21+1
    F23=F22+CNI+1
    F24=F23-1
    CALL CAL_FIJ_RECT(F11,F12,F13,F14,F21,F22,F23,F24)
    WRITE(2,235) '4(',F12,' ',F13,' ',F23,' ',F22,')'
    FAD1=F_INTERNAL_NUM_ALL+NL*NL+NJ*NK+NI*NK+NJ*NK+NL*(K-1)+J
    WRITE(3,231) FAD1
    OWN=CELL_ID_RECT(IP,JP,KP)
    WRITE(4,231) OWN
    ENDDO   
    ENDDO
    
        !BACK FACE--PROCESSOR_MI1_TO_MI2  !CASE RECTANGLE 3-2-3
    KP=CNL
    K=KP+CNL*(MJ2-1)
    DO JP=CNL,1,-1
    J=JP+CNL*(MJ1-1)
    DO IP=1,CNI
    I=IP+CNI*(MI-1)
    F11=(CNL+1)*(CNI+1)*(KP-1)+(CNI+1)*(JP-1)+IP-1
    F12=F11+1
    F13=F12+CNI+1
    F14=F13-1
    F21=F11+(CNL+1)*(CNI+1)
    F22=F21+1
    F23=F22+CNI+1
    F24=F23-1
    CALL CAL_FIJ_RECT(F11,F12,F13,F14,F21,F22,F23,F24)
    !WRITE(2,235) '4(',F21,' ',F22,' ',F23,' ',F24,')'
    WRITE(2,235) '4(',F24,' ',F21,' ',F22,' ',F23,')'
    FAD1=F_BACK_RECT(IP,JP)
    WRITE(3,231) FAD1
    OWN=CELL_ID_RECT(IP,JP,KP)
    WRITE(4,231) OWN
    ENDDO   
    ENDDO
    
        !FRONT FACE--PROCESSOR_MI1_TO_MI2  !CASE RECTANGLE 3-2-3
    KP=1
    K=KP+CNL*(MJ2-1)
    DO JP=1,CNL
    J=JP+CNL*(MJ1-1)
    DO IP=1,CNI
    I=IP+CNI*(MI-1)
    F11=(CNL+1)*(CNI+1)*(KP-1)+(CNI+1)*(JP-1)+IP-1
    F12=F11+1
    F13=F12+CNI+1
    F14=F13-1
    F21=F11+(CNL+1)*(CNI+1)
    F22=F21+1
    F23=F22+CNI+1
    F24=F23-1
    CALL CAL_FIJ_RECT(F11,F12,F13,F14,F21,F22,F23,F24)
    WRITE(2,235) '4(',F11,' ',F14,' ',F13,' ',F12,')'
    FAD1=F_FRONT_RECT(IP,JP)
    WRITE(3,231) FAD1
    OWN=CELL_ID_RECT(IP,JP,KP)
    WRITE(4,231) OWN
    ENDDO   
    ENDDO
    
        !DOWN FACE--PROCESSOR_MI1_TO_MI2  !CASE RECTANGLE 3-2-3
    JP=1
    J=JP+CNL*(MJ1-1)
    DO KP=1,CNL
    K=KP+CNL*(MJ2-1)
    DO IP=1,CNI
    I=IP+CNI*(MI-1)
    F11=(CNL+1)*(CNI+1)*(KP-1)+(CNI+1)*(JP-1)+IP-1
    F12=F11+1
    F13=F12+CNI+1
    F14=F13-1
    F21=F11+(CNL+1)*(CNI+1)
    F22=F21+1
    F23=F22+CNI+1
    F24=F23-1
    CALL CAL_FIJ_RECT(F11,F12,F13,F14,F21,F22,F23,F24)
    WRITE(2,235) '4(',F11,' ',F12,' ',F22,' ',F21,')'
    FAD1=F_DOWN_RECT(IP,KP)
    WRITE(3,231) FAD1
    OWN=CELL_ID_RECT(IP,JP,KP)
    WRITE(4,231) OWN
    ENDDO   
    ENDDO
    
        !RIGHT FACE--PROCESSOR_MI1_TO_MI2  !CASE RECTANGLE 3-2-3 
    IP=1
    I=IP+CNI*(MI-1)
    DO KP=1,CNL
    K=KP+CNL*(MJ2-1)
    DO JP=1,CNL
    J=JP+CNL*(MJ1-1)
    F11=(CNL+1)*(CNI+1)*(KP-1)+(CNI+1)*(JP-1)+IP-1
    F12=F11+1
    F13=F12+CNI+1
    F14=F13-1
    F21=F11+(CNL+1)*(CNI+1)
    F22=F21+1
    F23=F22+CNI+1
    F24=F23-1
    CALL CAL_FIJ_RECT(F11,F12,F13,F14,F21,F22,F23,F24)
    WRITE(2,235) '4(',F11,' ',F21,' ',F24,' ',F14,')'
    FAD1=F_RIGHT_RECT(JP,KP)
    WRITE(3,231) FAD1
    OWN=CELL_ID_RECT(IP,JP,KP)
    WRITE(4,231) OWN
    ENDDO   
    ENDDO
    
        !UP FACE--PROCESSOR_MI1_TO_MI2  !CASE RECTANGLE 3-2-3
    JP=CNL
    J=JP+CNL*(MJ1-1)
    DO KP=1,CNL
    K=KP+CNL*(MJ2-1)
    DO IP=1,CNI
    I=IP+CNI*(MI-1)
    F11=(CNL+1)*(CNI+1)*(KP-1)+(CNI+1)*(JP-1)+IP-1
    F12=F11+1
    F13=F12+CNI+1
    F14=F13-1
    F21=F11+(CNL+1)*(CNI+1)
    F22=F21+1
    F23=F22+CNI+1
    F24=F23-1
    CALL CAL_FIJ_RECT(F11,F12,F13,F14,F21,F22,F23,F24)
    WRITE(2,235) '4(',F14,' ',F24,' ',F23,' ',F13,')'
    FAD1=F_UP_RECT(IP,KP)
    WRITE(3,231) FAD1
    OWN=CELL_ID_RECT(IP,JP,KP)
    WRITE(4,231) OWN
    ENDDO   
    ENDDO
    
    N=5
    WRITE(31,231) N+4
    WRITE(31,*) '('  
    WRITE(31,*) '0'
    WRITE(31,*) '1'
    WRITE(31,*) '2' 
    WRITE(31,*) '3' 
    DO I=1,N
    WRITE(31,*) '-1' 
    ENDDO
    WRITE(31,*) ')'
    
    WRITE(21,231) N+4
    WRITE(21,*) "("    
    WRITE(21,*) "    INLET1"
    WRITE(21,*) "    {"    
    WRITE(21,*) "        type            patch;"
    N=0
    WRITE(21,232) "        nFaces          ",N,";"
    N=F_INTERNAL_NUM
    WRITE(21,232) "        startFace       ",N,";"
    WRITE(21,*) "    }"
    WRITE(21,*) "    INLET2"
    WRITE(21,*) "    {"    
    WRITE(21,*) "        type            patch;"
    N=0
    WRITE(21,232) "        nFaces          ",N,";"
    N=F_INTERNAL_NUM
    WRITE(21,232) "        startFace       ",N,";"
    WRITE(21,*) "    }"
    WRITE(21,*) "    OUTLET"
    WRITE(21,*) "    {"
    WRITE(21,*) "        type            patch;"
    N=0
    WRITE(21,232) "        nFaces          ",N,";"
    N=F_INTERNAL_NUM
    WRITE(21,232) "        startFace       ",N,";"
    WRITE(21,*) "    }"
    WRITE(21,*) "    HULL"
    WRITE(21,*) "    {"
    WRITE(21,*) "        type            wall;"
    WRITE(21,*) "        inGroups        1(wall);"
    N=CNL*CNL
    WRITE(21,232) "        nFaces          ",N,";"
    N=F_INTERNAL_NUM
    WRITE(21,232) "        startFace       ",N,";"
    WRITE(21,*) "    }"
    
        !BACK FACE--PROCESSOR_MI1_TO_MI2  !CASE RECTANGLE 3-2-3
    PRO1=PNI*PNJ+PNL*PNI*(MJ2-1)+PNI*(MJ1-1)+MI-1
    MJ=1
    PRO2=PNI*(MJ-1)+MI-1
    WRITE(21,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(21,*) '    {'
    WRITE(21,*) '        type            processor;'
    WRITE(21,*) '        inGroups        List<word> 1(processor);'
    N=CNI*CNL
    WRITE(21,232) "        nFaces          ",N,";"
    N=F_INTERNAL_NUM+CNL*CNL
    WRITE(21,232) "        startFace       ",N,";"
    WRITE(21,*) '        matchTolerance  0.0001;'
    WRITE(21,232) "        myProcNo        ",PRO1,";"
    WRITE(21,232) "        neighbProcNo    ",PRO2,";"
    WRITE(21,*) '    }'
    WRITE(41,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(41,*) '    {'
    WRITE(41,*) '        type            processor;'
    WRITE(41,*) '        value           nonuniform List<scalar> '
    N=CNI*CNL
    WRITE(41,'(I0)') N
    WRITE(41,*) '('
    DO JP=1,CNL
    DO IP=1,CNI
    WRITE(41,'(E21.14)') P_BACK_RECT(IP,JP)    
    ENDDO
    ENDDO
    WRITE(41,*) ')'
    WRITE(41,*) ';'
    WRITE(41,*) '    }'
    WRITE(51,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(51,*) '    {'
    WRITE(51,*) '        type            processor;'
    WRITE(51,*) '        value           nonuniform List<vector> '
    N=CNI*CNL
    WRITE(51,'(I0)') N
    WRITE(51,*) '('
    DO JP=1,CNL
    DO IP=1,CNI
    UU1=U_BACK_RECT(1,IP,JP)
    UU2=U_BACK_RECT(2,IP,JP)
    UU3=U_BACK_RECT(3,IP,JP)
    WRITE(51,230) '( ',UU1,' ',UU2,' ',UU3,' )'  
    ENDDO
    ENDDO
    WRITE(51,*) ')'
    WRITE(51,*) ';'
    WRITE(51,*) '    }'
    WRITE(61,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(61,*) '    {'
    WRITE(61,*) '        type            processor;'
    WRITE(61,*) "        value           uniform 0;"
    WRITE(61,*) '    }'
    
        !FRONT FACE--PROCESSOR_MI1_TO_MI2  !CASE RECTANGLE 3-2-3
    PRO2=PNI*PNJ+PNL*PNI*(MJ2-2)+PNI*(MJ1-1)+MI-1
    WRITE(21,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(21,*) '    {'
    WRITE(21,*) '        type            processor;'
    WRITE(21,*) '        inGroups        List<word> 1(processor);'
    N=CNI*CNL
    WRITE(21,232) "        nFaces          ",N,";"
    N=F_INTERNAL_NUM+CNL*CNL+CNI*CNL
    WRITE(21,232) "        startFace       ",N,";"
    WRITE(21,*) '        matchTolerance  0.0001;'
    WRITE(21,232) "        myProcNo        ",PRO1,";"
    WRITE(21,232) "        neighbProcNo    ",PRO2,";"
    WRITE(21,*) '    }'
    WRITE(41,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(41,*) '    {'
    WRITE(41,*) '        type            processor;'
    WRITE(41,*) '        value           nonuniform List<scalar> '
    N=CNI*CNL
    WRITE(41,'(I0)') N
    WRITE(41,*) '('
    DO JP=1,CNL
    DO IP=1,CNI
    WRITE(41,'(E21.14)') P_FRONT_RECT(IP,JP)    
    ENDDO
    ENDDO
    WRITE(41,*) ')'
    WRITE(41,*) ';'
    WRITE(41,*) '    }'
    WRITE(51,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(51,*) '    {'
    WRITE(51,*) '        type            processor;'
    WRITE(51,*) '        value           nonuniform List<vector> '
    N=CNI*CNL
    WRITE(51,'(I0)') N
    WRITE(51,*) '('
    DO JP=1,CNL
    DO IP=1,CNI
    UU1=U_FRONT_RECT(1,IP,JP)
    UU2=U_FRONT_RECT(2,IP,JP)
    UU3=U_FRONT_RECT(3,IP,JP)
    WRITE(51,230) '( ',UU1,' ',UU2,' ',UU3,' )'  
    ENDDO
    ENDDO
    WRITE(51,*) ')'
    WRITE(51,*) ';'
    WRITE(51,*) '    }'
    WRITE(61,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(61,*) '    {'
    WRITE(61,*) '        type            processor;'
    WRITE(61,*) "        value           uniform 0;"
    WRITE(61,*) '    }'
    
        !DOWN FACE--PROCESSOR_MI1_TO_MI2  !CASE RECTANGLE 3-2-3
    PRO2=PNI*PNJ+PNL*PNI*(MJ2-1)+PNI*(MJ1-2)+MI-1
    WRITE(21,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(21,*) '    {'
    WRITE(21,*) '        type            processor;'
    WRITE(21,*) '        inGroups        List<word> 1(processor);'
    N=CNI*CNL
    WRITE(21,232) "        nFaces          ",N,";"
    N=F_INTERNAL_NUM+CNL*CNL+CNI*CNL+CNI*CNL
    WRITE(21,232) "        startFace       ",N,";"
    WRITE(21,*) '        matchTolerance  0.0001;'
    WRITE(21,232) "        myProcNo        ",PRO1,";"
    WRITE(21,232) "        neighbProcNo    ",PRO2,";"
    WRITE(21,*) '    }'
    WRITE(41,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(41,*) '    {'
    WRITE(41,*) '        type            processor;'
    WRITE(41,*) '        value           nonuniform List<scalar> '
    N=CNI*CNL
    WRITE(41,'(I0)') N
    WRITE(41,*) '('
    DO KP=1,CNL
    DO IP=1,CNI
    WRITE(41,'(E21.14)') P_DOWN_RECT(IP,KP)    
    ENDDO
    ENDDO
    WRITE(41,*) ')'
    WRITE(41,*) ';'
    WRITE(41,*) '    }'
    WRITE(51,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(51,*) '    {'
    WRITE(51,*) '        type            processor;'
    WRITE(51,*) '        value           nonuniform List<vector> '
    N=CNI*CNL
    WRITE(51,'(I0)') N
    WRITE(51,*) '('
    DO KP=1,CNL
    DO IP=1,CNI
    UU1=U_DOWN_RECT(1,IP,KP)
    UU2=U_DOWN_RECT(2,IP,KP)
    UU3=U_DOWN_RECT(3,IP,KP)
    WRITE(51,230) '( ',UU1,' ',UU2,' ',UU3,' )'  
    ENDDO
    ENDDO
    WRITE(51,*) ')'
    WRITE(51,*) ';'
    WRITE(51,*) '    }'
    WRITE(61,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(61,*) '    {'
    WRITE(61,*) '        type            processor;'
    WRITE(61,*) "        value           uniform 0;"
    WRITE(61,*) '    }'
    
        !RIGHT FACE--PROCESSOR_MI1_TO_MI2  !CASE RECTANGLE 3-2-3
    PRO2=PNI*PNJ+PNL*PNI*(MJ2-1)+PNI*(MJ1-1)+MI-2
    WRITE(21,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(21,*) '    {'
    WRITE(21,*) '        type            processor;'
    WRITE(21,*) '        inGroups        List<word> 1(processor);'
    N=CNL*CNL
    WRITE(21,232) "        nFaces          ",N,";"
    N=F_INTERNAL_NUM+CNL*CNL+CNI*CNL+CNI*CNL+CNI*CNL
    WRITE(21,232) "        startFace       ",N,";"
    WRITE(21,*) '        matchTolerance  0.0001;'
    WRITE(21,232) "        myProcNo        ",PRO1,";"
    WRITE(21,232) "        neighbProcNo    ",PRO2,";"
    WRITE(21,*) '    }'
    WRITE(41,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(41,*) '    {'
    WRITE(41,*) '        type            processor;'
    WRITE(41,*) '        value           nonuniform List<scalar> '
    N=CNL*CNL
    WRITE(41,'(I0)') N
    WRITE(41,*) '('
    DO KP=1,CNL
    DO JP=1,CNL
    WRITE(41,'(E21.14)') P_RIGHT_RECT(JP,KP)    
    ENDDO
    ENDDO
    WRITE(41,*) ')'
    WRITE(41,*) ';'
    WRITE(41,*) '    }'
    WRITE(51,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(51,*) '    {'
    WRITE(51,*) '        type            processor;'
    WRITE(51,*) '        value           nonuniform List<vector> '
    N=CNL*CNL
    WRITE(51,'(I0)') N
    WRITE(51,*) '('
    DO KP=1,CNL
    DO JP=1,CNL
    UU1=U_RIGHT_RECT(1,JP,KP)
    UU2=U_RIGHT_RECT(2,JP,KP)
    UU3=U_RIGHT_RECT(3,JP,KP)
    WRITE(51,230) '( ',UU1,' ',UU2,' ',UU3,' )'  
    ENDDO
    ENDDO
    WRITE(51,*) ')'
    WRITE(51,*) ';'
    WRITE(51,*) '    }'
    WRITE(61,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(61,*) '    {'
    WRITE(61,*) '        type            processor;'
    WRITE(61,*) "        value           uniform 0;"
    WRITE(61,*) '    }'
    
        !UP FACE--PROCESSOR_MI1_TO_MI2  !CASE RECTANGLE 3-2-3
    PRO2=PNI*PNJ+PNL*PNI*(MJ2-1)+PNI*MJ1+MI-1
    WRITE(21,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(21,*) '    {'
    WRITE(21,*) '        type            processor;'
    WRITE(21,*) '        inGroups        List<word> 1(processor);'
    N=CNI*CNL
    WRITE(21,232) "        nFaces          ",N,";"
    N=F_INTERNAL_NUM+CNL*CNL+CNI*CNL+CNI*CNL+CNI*CNL+CNL*CNL
    WRITE(21,232) "        startFace       ",N,";"
    WRITE(21,*) '        matchTolerance  0.0001;'
    WRITE(21,232) "        myProcNo        ",PRO1,";"
    WRITE(21,232) "        neighbProcNo    ",PRO2,";"
    WRITE(21,*) '    }'
    WRITE(41,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(41,*) '    {'
    WRITE(41,*) '        type            processor;'
    WRITE(41,*) '        value           nonuniform List<scalar> '
    N=CNI*CNL
    WRITE(41,'(I0)') N
    WRITE(41,*) '('
    DO KP=1,CNL
    DO IP=1,CNI
    WRITE(41,'(E21.14)') P_UP_RECT(IP,KP)    
    ENDDO
    ENDDO
    WRITE(41,*) ')'
    WRITE(41,*) ';'
    WRITE(41,*) '    }'
    WRITE(51,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(51,*) '    {'
    WRITE(51,*) '        type            processor;'
    WRITE(51,*) '        value           nonuniform List<vector> '
    N=CNI*CNL
    WRITE(51,'(I0)') N
    WRITE(51,*) '('
    DO KP=1,CNL
    DO IP=1,CNI
    UU1=U_UP_RECT(1,IP,KP)
    UU2=U_UP_RECT(2,IP,KP)
    UU3=U_UP_RECT(3,IP,KP)
    WRITE(51,230) '( ',UU1,' ',UU2,' ',UU3,' )'  
    ENDDO
    ENDDO
    WRITE(51,*) ')'
    WRITE(51,*) ';'
    WRITE(51,*) '    }'
    WRITE(61,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(61,*) '    {'
    WRITE(61,*) '        type            processor;'
    WRITE(61,*) "        value           uniform 0;"
    WRITE(61,*) '    }'
    
    
    
    
    
    
    
    ELSE
    !CASE RECTANGLE 3-2-2
        
        !LEFT FACE--HULL  !CASE RECTANGLE 3-2-2 
    IP=CNI
    I=IP+CNI*(MI-1)
    DO KP=1,CNL
    K=KP+CNL*(MJ2-1)
    DO JP=1,CNL
    J=JP+CNL*(MJ1-1)
    F11=(CNL+1)*(CNI+1)*(KP-1)+(CNI+1)*(JP-1)+IP-1
    F12=F11+1
    F13=F12+CNI+1
    F14=F13-1
    F21=F11+(CNL+1)*(CNI+1)
    F22=F21+1
    F23=F22+CNI+1
    F24=F23-1
    CALL CAL_FIJ_RECT(F11,F12,F13,F14,F21,F22,F23,F24)
    WRITE(2,235) '4(',F12,' ',F13,' ',F23,' ',F22,')'
    FAD1=F_INTERNAL_NUM_ALL+NL*NL+NJ*NK+NI*NK+NJ*NK+NL*(K-1)+J
    WRITE(3,231) FAD1
    OWN=CELL_ID_RECT(IP,JP,KP)
    WRITE(4,231) OWN
    ENDDO   
    ENDDO
    
        !FRONT FACE--PROCESSOR_MI1_TO_MI2  !CASE RECTANGLE 3-2-2
    KP=1
    K=KP+CNL*(MJ2-1)
    DO JP=1,CNL
    J=JP+CNL*(MJ1-1)
    DO IP=1,CNI
    I=IP+CNI*(MI-1)
    F11=(CNL+1)*(CNI+1)*(KP-1)+(CNI+1)*(JP-1)+IP-1
    F12=F11+1
    F13=F12+CNI+1
    F14=F13-1
    F21=F11+(CNL+1)*(CNI+1)
    F22=F21+1
    F23=F22+CNI+1
    F24=F23-1
    CALL CAL_FIJ_RECT(F11,F12,F13,F14,F21,F22,F23,F24)
    WRITE(2,235) '4(',F11,' ',F14,' ',F13,' ',F12,')'
    FAD1=F_FRONT_RECT(IP,JP)
    WRITE(3,231) FAD1
    OWN=CELL_ID_RECT(IP,JP,KP)
    WRITE(4,231) OWN
    ENDDO   
    ENDDO
    
        !DOWN FACE--PROCESSOR_MI1_TO_MI2  !CASE RECTANGLE 3-2-2
    JP=1
    J=JP+CNL*(MJ1-1)
    DO KP=1,CNL
    K=KP+CNL*(MJ2-1)
    DO IP=1,CNI
    I=IP+CNI*(MI-1)
    F11=(CNL+1)*(CNI+1)*(KP-1)+(CNI+1)*(JP-1)+IP-1
    F12=F11+1
    F13=F12+CNI+1
    F14=F13-1
    F21=F11+(CNL+1)*(CNI+1)
    F22=F21+1
    F23=F22+CNI+1
    F24=F23-1
    CALL CAL_FIJ_RECT(F11,F12,F13,F14,F21,F22,F23,F24)
    WRITE(2,235) '4(',F11,' ',F12,' ',F22,' ',F21,')'
    FAD1=F_DOWN_RECT(IP,KP)
    WRITE(3,231) FAD1
    OWN=CELL_ID_RECT(IP,JP,KP)
    WRITE(4,231) OWN
    ENDDO   
    ENDDO
    
        !RIGHT FACE--PROCESSOR_MI1_TO_MI2  !CASE RECTANGLE 3-2-2 
    IP=1
    I=IP+CNI*(MI-1)
    DO KP=1,CNL
    K=KP+CNL*(MJ2-1)
    DO JP=1,CNL
    J=JP+CNL*(MJ1-1)
    F11=(CNL+1)*(CNI+1)*(KP-1)+(CNI+1)*(JP-1)+IP-1
    F12=F11+1
    F13=F12+CNI+1
    F14=F13-1
    F21=F11+(CNL+1)*(CNI+1)
    F22=F21+1
    F23=F22+CNI+1
    F24=F23-1
    CALL CAL_FIJ_RECT(F11,F12,F13,F14,F21,F22,F23,F24)
    WRITE(2,235) '4(',F11,' ',F21,' ',F24,' ',F14,')'
    FAD1=F_RIGHT_RECT(JP,KP)
    WRITE(3,231) FAD1
    OWN=CELL_ID_RECT(IP,JP,KP)
    WRITE(4,231) OWN
    ENDDO   
    ENDDO
    
        !UP FACE--PROCESSOR_MI1_TO_MI2  !CASE RECTANGLE 3-2-2
    JP=CNL
    J=JP+CNL*(MJ1-1)
    DO KP=1,CNL
    K=KP+CNL*(MJ2-1)
    DO IP=1,CNI
    I=IP+CNI*(MI-1)
    F11=(CNL+1)*(CNI+1)*(KP-1)+(CNI+1)*(JP-1)+IP-1
    F12=F11+1
    F13=F12+CNI+1
    F14=F13-1
    F21=F11+(CNL+1)*(CNI+1)
    F22=F21+1
    F23=F22+CNI+1
    F24=F23-1
    CALL CAL_FIJ_RECT(F11,F12,F13,F14,F21,F22,F23,F24)
    WRITE(2,235) '4(',F14,' ',F24,' ',F23,' ',F13,')'
    FAD1=F_UP_RECT(IP,KP)
    WRITE(3,231) FAD1
    OWN=CELL_ID_RECT(IP,JP,KP)
    WRITE(4,231) OWN
    ENDDO   
    ENDDO
    
        !BACK FACE--PROCESSOR_MI1_TO_MI2  !CASE RECTANGLE 3-2-2
    KP=CNL
    K=KP+CNL*(MJ2-1)
    DO JP=1,CNL
    J=JP+CNL*(MJ1-1)
    DO IP=1,CNI
    I=IP+CNI*(MI-1)
    F11=(CNL+1)*(CNI+1)*(KP-1)+(CNI+1)*(JP-1)+IP-1
    F12=F11+1
    F13=F12+CNI+1
    F14=F13-1
    F21=F11+(CNL+1)*(CNI+1)
    F22=F21+1
    F23=F22+CNI+1
    F24=F23-1
    CALL CAL_FIJ_RECT(F11,F12,F13,F14,F21,F22,F23,F24)
    WRITE(2,235) '4(',F21,' ',F22,' ',F23,' ',F24,')'
    FAD1=F_BACK_RECT(IP,JP)
    WRITE(3,231) FAD1
    OWN=CELL_ID_RECT(IP,JP,KP)
    WRITE(4,231) OWN
    ENDDO   
    ENDDO
    
    N=5
    WRITE(31,231) N+4
    WRITE(31,*) '('  
    WRITE(31,*) '0'
    WRITE(31,*) '1'
    WRITE(31,*) '2' 
    WRITE(31,*) '3' 
    DO I=1,N
    WRITE(31,*) '-1' 
    ENDDO
    WRITE(31,*) ')'
    
    WRITE(21,231) N+4
    WRITE(21,*) "("    
    WRITE(21,*) "    INLET1"
    WRITE(21,*) "    {"    
    WRITE(21,*) "        type            patch;"
    N=0
    WRITE(21,232) "        nFaces          ",N,";"
    N=F_INTERNAL_NUM
    WRITE(21,232) "        startFace       ",N,";"
    WRITE(21,*) "    }"
    WRITE(21,*) "    INLET2"
    WRITE(21,*) "    {"    
    WRITE(21,*) "        type            patch;"
    N=0
    WRITE(21,232) "        nFaces          ",N,";"
    N=F_INTERNAL_NUM
    WRITE(21,232) "        startFace       ",N,";"
    WRITE(21,*) "    }"
    WRITE(21,*) "    OUTLET"
    WRITE(21,*) "    {"
    WRITE(21,*) "        type            patch;"
    N=0
    WRITE(21,232) "        nFaces          ",N,";"
    N=F_INTERNAL_NUM
    WRITE(21,232) "        startFace       ",N,";"
    WRITE(21,*) "    }"
    WRITE(21,*) "    HULL"
    WRITE(21,*) "    {"
    WRITE(21,*) "        type            wall;"
    WRITE(21,*) "        inGroups        1(wall);"
    N=CNL*CNL
    WRITE(21,232) "        nFaces          ",N,";"
    N=F_INTERNAL_NUM
    WRITE(21,232) "        startFace       ",N,";"
    WRITE(21,*) "    }"
    
        !FRONT FACE--PROCESSOR_MI1_TO_MI2  !CASE RECTANGLE 3-2-2
    PRO1=PNI*PNJ+PNL*PNI*(MJ2-1)+PNI*(MJ1-1)+MI-1
    PRO2=PNI*PNJ+PNL*PNI*(MJ2-2)+PNI*(MJ1-1)+MI-1
    WRITE(21,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(21,*) '    {'
    WRITE(21,*) '        type            processor;'
    WRITE(21,*) '        inGroups        List<word> 1(processor);'
    N=CNI*CNL
    WRITE(21,232) "        nFaces          ",N,";"
    N=F_INTERNAL_NUM+CNL*CNL
    WRITE(21,232) "        startFace       ",N,";"
    WRITE(21,*) '        matchTolerance  0.0001;'
    WRITE(21,232) "        myProcNo        ",PRO1,";"
    WRITE(21,232) "        neighbProcNo    ",PRO2,";"
    WRITE(21,*) '    }'
    WRITE(41,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(41,*) '    {'
    WRITE(41,*) '        type            processor;'
    WRITE(41,*) '        value           nonuniform List<scalar> '
    N=CNI*CNL
    WRITE(41,'(I0)') N
    WRITE(41,*) '('
    DO JP=1,CNL
    DO IP=1,CNI
    WRITE(41,'(E21.14)') P_FRONT_RECT(IP,JP)    
    ENDDO
    ENDDO
    WRITE(41,*) ')'
    WRITE(41,*) ';'
    WRITE(41,*) '    }'
    WRITE(51,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(51,*) '    {'
    WRITE(51,*) '        type            processor;'
    WRITE(51,*) '        value           nonuniform List<vector> '
    N=CNI*CNL
    WRITE(51,'(I0)') N
    WRITE(51,*) '('
    DO JP=1,CNL
    DO IP=1,CNI
    UU1=U_FRONT_RECT(1,IP,JP)
    UU2=U_FRONT_RECT(2,IP,JP)
    UU3=U_FRONT_RECT(3,IP,JP)
    WRITE(51,230) '( ',UU1,' ',UU2,' ',UU3,' )'  
    ENDDO
    ENDDO
    WRITE(51,*) ')'
    WRITE(51,*) ';'
    WRITE(51,*) '    }'
    WRITE(61,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(61,*) '    {'
    WRITE(61,*) '        type            processor;'
    WRITE(61,*) "        value           uniform 0;"
    WRITE(61,*) '    }'
    
        !DOWN FACE--PROCESSOR_MI1_TO_MI2  !CASE RECTANGLE 3-2-2
    PRO2=PNI*PNJ+PNL*PNI*(MJ2-1)+PNI*(MJ1-2)+MI-1
    WRITE(21,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(21,*) '    {'
    WRITE(21,*) '        type            processor;'
    WRITE(21,*) '        inGroups        List<word> 1(processor);'
    N=CNI*CNL
    WRITE(21,232) "        nFaces          ",N,";"
    N=F_INTERNAL_NUM+CNL*CNL+CNI*CNL
    WRITE(21,232) "        startFace       ",N,";"
    WRITE(21,*) '        matchTolerance  0.0001;'
    WRITE(21,232) "        myProcNo        ",PRO1,";"
    WRITE(21,232) "        neighbProcNo    ",PRO2,";"
    WRITE(21,*) '    }'
    WRITE(41,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(41,*) '    {'
    WRITE(41,*) '        type            processor;'
    WRITE(41,*) '        value           nonuniform List<scalar> '
    N=CNI*CNL
    WRITE(41,'(I0)') N
    WRITE(41,*) '('
    DO KP=1,CNL
    DO IP=1,CNI
    WRITE(41,'(E21.14)') P_DOWN_RECT(IP,KP)    
    ENDDO
    ENDDO
    WRITE(41,*) ')'
    WRITE(41,*) ';'
    WRITE(41,*) '    }'
    WRITE(51,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(51,*) '    {'
    WRITE(51,*) '        type            processor;'
    WRITE(51,*) '        value           nonuniform List<vector> '
    N=CNI*CNL
    WRITE(51,'(I0)') N
    WRITE(51,*) '('
    DO KP=1,CNL
    DO IP=1,CNI
    UU1=U_DOWN_RECT(1,IP,KP)
    UU2=U_DOWN_RECT(2,IP,KP)
    UU3=U_DOWN_RECT(3,IP,KP)
    WRITE(51,230) '( ',UU1,' ',UU2,' ',UU3,' )'  
    ENDDO
    ENDDO
    WRITE(51,*) ')'
    WRITE(51,*) ';'
    WRITE(51,*) '    }'
    WRITE(61,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(61,*) '    {'
    WRITE(61,*) '        type            processor;'
    WRITE(61,*) "        value           uniform 0;"
    WRITE(61,*) '    }'
    
        !RIGHT FACE--PROCESSOR_MI1_TO_MI2  !CASE RECTANGLE 3-2-2
    PRO2=PNI*PNJ+PNL*PNI*(MJ2-1)+PNI*(MJ1-1)+MI-2
    WRITE(21,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(21,*) '    {'
    WRITE(21,*) '        type            processor;'
    WRITE(21,*) '        inGroups        List<word> 1(processor);'
    N=CNL*CNL
    WRITE(21,232) "        nFaces          ",N,";"
    N=F_INTERNAL_NUM+CNL*CNL+CNI*CNL+CNI*CNL
    WRITE(21,232) "        startFace       ",N,";"
    WRITE(21,*) '        matchTolerance  0.0001;'
    WRITE(21,232) "        myProcNo        ",PRO1,";"
    WRITE(21,232) "        neighbProcNo    ",PRO2,";"
    WRITE(21,*) '    }'
    WRITE(41,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(41,*) '    {'
    WRITE(41,*) '        type            processor;'
    WRITE(41,*) '        value           nonuniform List<scalar> '
    N=CNL*CNL
    WRITE(41,'(I0)') N
    WRITE(41,*) '('
    DO KP=1,CNL
    DO JP=1,CNL
    WRITE(41,'(E21.14)') P_RIGHT_RECT(JP,KP)    
    ENDDO
    ENDDO
    WRITE(41,*) ')'
    WRITE(41,*) ';'
    WRITE(41,*) '    }'
    WRITE(51,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(51,*) '    {'
    WRITE(51,*) '        type            processor;'
    WRITE(51,*) '        value           nonuniform List<vector> '
    N=CNL*CNL
    WRITE(51,'(I0)') N
    WRITE(51,*) '('
    DO KP=1,CNL
    DO JP=1,CNL
    UU1=U_RIGHT_RECT(1,JP,KP)
    UU2=U_RIGHT_RECT(2,JP,KP)
    UU3=U_RIGHT_RECT(3,JP,KP)
    WRITE(51,230) '( ',UU1,' ',UU2,' ',UU3,' )'  
    ENDDO
    ENDDO
    WRITE(51,*) ')'
    WRITE(51,*) ';'
    WRITE(51,*) '    }'
    WRITE(61,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(61,*) '    {'
    WRITE(61,*) '        type            processor;'
    WRITE(61,*) "        value           uniform 0;"
    WRITE(61,*) '    }' 
    
        !UP FACE--PROCESSOR_MI1_TO_MI2  !CASE RECTANGLE 3-2-2
    PRO2=PNI*PNJ+PNL*PNI*(MJ2-1)+PNI*MJ1+MI-1
    WRITE(21,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(21,*) '    {'
    WRITE(21,*) '        type            processor;'
    WRITE(21,*) '        inGroups        List<word> 1(processor);'
    N=CNI*CNL
    WRITE(21,232) "        nFaces          ",N,";"
    N=F_INTERNAL_NUM+CNL*CNL+CNI*CNL+CNI*CNL+CNL*CNL
    WRITE(21,232) "        startFace       ",N,";"
    WRITE(21,*) '        matchTolerance  0.0001;'
    WRITE(21,232) "        myProcNo        ",PRO1,";"
    WRITE(21,232) "        neighbProcNo    ",PRO2,";"
    WRITE(21,*) '    }'
    WRITE(41,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(41,*) '    {'
    WRITE(41,*) '        type            processor;'
    WRITE(41,*) '        value           nonuniform List<scalar> '
    N=CNI*CNL
    WRITE(41,'(I0)') N
    WRITE(41,*) '('
    DO KP=1,CNL
    DO IP=1,CNI
    WRITE(41,'(E21.14)') P_UP_RECT(IP,KP)    
    ENDDO
    ENDDO
    WRITE(41,*) ')'
    WRITE(41,*) ';'
    WRITE(41,*) '    }'
    WRITE(51,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(51,*) '    {'
    WRITE(51,*) '        type            processor;'
    WRITE(51,*) '        value           nonuniform List<vector> '
    N=CNI*CNL
    WRITE(51,'(I0)') N
    WRITE(51,*) '('
    DO KP=1,CNL
    DO IP=1,CNI
    UU1=U_UP_RECT(1,IP,KP)
    UU2=U_UP_RECT(2,IP,KP)
    UU3=U_UP_RECT(3,IP,KP)
    WRITE(51,230) '( ',UU1,' ',UU2,' ',UU3,' )'  
    ENDDO
    ENDDO
    WRITE(51,*) ')'
    WRITE(51,*) ';'
    WRITE(51,*) '    }'
    WRITE(61,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(61,*) '    {'
    WRITE(61,*) '        type            processor;'
    WRITE(61,*) "        value           uniform 0;"
    WRITE(61,*) '    }'
    
        !BACK FACE--PROCESSOR_MI1_TO_MI2  !CASE RECTANGLE 3-2-2
    PRO2=PNI*PNJ+PNL*PNI*MJ2+PNI*(MJ1-1)+MI-1
    WRITE(21,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(21,*) '    {'
    WRITE(21,*) '        type            processor;'
    WRITE(21,*) '        inGroups        List<word> 1(processor);'
    N=CNI*CNL
    WRITE(21,232) "        nFaces          ",N,";"
    N=F_INTERNAL_NUM+CNL*CNL+CNI*CNL+CNI*CNL+CNL*CNL+CNI*CNL
    WRITE(21,232) "        startFace       ",N,";"
    WRITE(21,*) '        matchTolerance  0.0001;'
    WRITE(21,232) "        myProcNo        ",PRO1,";"
    WRITE(21,232) "        neighbProcNo    ",PRO2,";"
    WRITE(21,*) '    }'
    WRITE(41,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(41,*) '    {'
    WRITE(41,*) '        type            processor;'
    WRITE(41,*) '        value           nonuniform List<scalar> '
    N=CNI*CNL
    WRITE(41,'(I0)') N
    WRITE(41,*) '('
    DO JP=1,CNL
    DO IP=1,CNI
    WRITE(41,'(E21.14)') P_BACK_RECT(IP,JP)    
    ENDDO
    ENDDO
    WRITE(41,*) ')'
    WRITE(41,*) ';'
    WRITE(41,*) '    }'
    WRITE(51,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(51,*) '    {'
    WRITE(51,*) '        type            processor;'
    WRITE(51,*) '        value           nonuniform List<vector> '
    N=CNI*CNL
    WRITE(51,'(I0)') N
    WRITE(51,*) '('
    DO JP=1,CNL
    DO IP=1,CNI
    UU1=U_BACK_RECT(1,IP,JP)
    UU2=U_BACK_RECT(2,IP,JP)
    UU3=U_BACK_RECT(3,IP,JP)
    WRITE(51,230) '( ',UU1,' ',UU2,' ',UU3,' )'  
    ENDDO
    ENDDO
    WRITE(51,*) ')'
    WRITE(51,*) ';'
    WRITE(51,*) '    }'
    WRITE(61,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(61,*) '    {'
    WRITE(61,*) '        type            processor;'
    WRITE(61,*) "        value           uniform 0;"
    WRITE(61,*) '    }'
    
    
    
    
    
    
    
    
    
    ENDIF 
    ENDIF              
    ENDIF
            
    ELSE
!****************************************NI-MIDDLE****************************************
    IF (MJ.EQ.1) THEN
    !CASE CYLINDER 2-1
    
        !RIGHT FACE--PROCESSOR_MI1_TO_MI2  !CASE CYLINDER 2-1
    IP=1
    I=IP+CNI*(MI-1)
    DO KP=1,NK
    K=KP
    DO JP=1,CNJ
    J=JP+CNJ*(MJ-1)
    F11=(CNJ+1)*(CNI+1)*(KP-1)+(CNI+1)*(JP-1)+IP-1
    F12=F11+1
    F13=F12+CNI+1
    F14=F13-1
    F21=F11+(CNJ+1)*(CNI+1)
    IF (KP.EQ.NK) F21=(CNI+1)*(JP-1)+IP-1
    F22=F21+1
    F23=F22+CNI+1
    F24=F23-1
    CALL CAL_FIJ_CY(F11,F12,F13,F14,F21,F22,F23,F24)
    WRITE(2,235) '4(',F11,' ',F21,' ',F24,' ',F14,')'
    FAD1=F_RIGHT_CY(JP,KP)
    WRITE(3,231) FAD1
    OWN=CELL_ID_CY(IP,JP,KP)
    WRITE(4,231) OWN
    ENDDO   
    ENDDO
        
        !LEFT FACE--PROCESSOR_MI1_TO_MI2  !CASE CYLINDER 2-1
    IP=CNI
    I=IP+CNI*(MI-1)
    DO KP=1,NK
    K=KP
    DO JP=1,CNJ
    J=JP+CNJ*(MJ-1)
    F11=(CNJ+1)*(CNI+1)*(KP-1)+(CNI+1)*(JP-1)+IP-1
    F12=F11+1
    F13=F12+CNI+1
    F14=F13-1
    F21=F11+(CNJ+1)*(CNI+1)
    IF (KP.EQ.NK) F21=(CNI+1)*(JP-1)+IP-1
    F22=F21+1
    F23=F22+CNI+1
    F24=F23-1
    CALL CAL_FIJ_CY(F11,F12,F13,F14,F21,F22,F23,F24)
    WRITE(2,235) '4(',F12,' ',F13,' ',F23,' ',F22,')'
    FAD1=F_LEFT_CY(JP,KP)
    WRITE(3,231) FAD1
    OWN=CELL_ID_CY(IP,JP,KP)
    WRITE(4,231) OWN
    ENDDO   
    ENDDO
    
        !UP FACE--PROCESSOR_MI1_TO_MI2  !CASE CYLINDER 2-1
    JP=CNJ
    J=JP+CNJ*(MJ-1)
    DO KP=1,NK
    K=KP
    DO IP=1,CNI
    I=IP+CNI*(MI-1)
    F11=(CNJ+1)*(CNI+1)*(KP-1)+(CNI+1)*(JP-1)+IP-1
    F12=F11+1
    F13=F12+CNI+1
    F14=F13-1
    F21=F11+(CNJ+1)*(CNI+1)
    IF (KP.EQ.NK) F21=(CNI+1)*(JP-1)+IP-1
    F22=F21+1
    F23=F22+CNI+1
    F24=F23-1
    CALL CAL_FIJ_CY(F11,F12,F13,F14,F21,F22,F23,F24)
    WRITE(2,235) '4(',F14,' ',F24,' ',F23,' ',F13,')'
    FAD1=F_UP_CY(IP,KP)
    WRITE(3,231) FAD1
    OWN=CELL_ID_CY(IP,JP,KP)
    WRITE(4,231) OWN
    ENDDO   
    ENDDO
    
        !DOWN FACE--PROCESSOR_MI1_TO_MI2(RECTANGLE 1-1)  !CASE CYLINDER 2-1
    JP=1
    J=JP+CNJ*(MJ-1)
    DO KP=1,CNL
    K=KP
    DO IP=1,CNI
    I=IP+CNI*(MI-1)
    F11=(CNJ+1)*(CNI+1)*(KP-1)+(CNI+1)*(JP-1)+IP-1
    F12=F11+1
    F13=F12+CNI+1
    F14=F13-1
    F21=F11+(CNJ+1)*(CNI+1)
    IF (KP.EQ.NK) F21=(CNI+1)*(JP-1)+IP-1
    F22=F21+1
    F23=F22+CNI+1
    F24=F23-1
    CALL CAL_FIJ_CY(F11,F12,F13,F14,F21,F22,F23,F24)
    WRITE(2,235) '4(',F11,' ',F12,' ',F22,' ',F21,')'
    FAD1=F_DOWN_CY(IP,KP)
    WRITE(3,231) FAD1
    OWN=CELL_ID_CY(IP,JP,KP)
    WRITE(4,231) OWN
    ENDDO   
    ENDDO
    DO KP=NK-CNL+1,NK
    K=KP
    DO IP=1,CNI
    I=IP+CNI*(MI-1)
    F11=(CNJ+1)*(CNI+1)*(KP-1)+(CNI+1)*(JP-1)+IP-1
    F12=F11+1
    F13=F12+CNI+1
    F14=F13-1
    F21=F11+(CNJ+1)*(CNI+1)
    IF (KP.EQ.NK) F21=(CNI+1)*(JP-1)+IP-1
    F22=F21+1
    F23=F22+CNI+1
    F24=F23-1
    CALL CAL_FIJ_CY(F11,F12,F13,F14,F21,F22,F23,F24)
    WRITE(2,235) '4(',F11,' ',F12,' ',F22,' ',F21,')'
    FAD1=F_DOWN_CY(IP,KP)
    WRITE(3,231) FAD1
    OWN=CELL_ID_CY(IP,JP,KP)
    WRITE(4,231) OWN
    ENDDO   
    ENDDO
    
        !DOWN FACE--PROCESSOR_MI1_TO_MI2(RECTANGLE 2-1)  !CASE CYLINDER 2-1
    JP=1
    J=JP+CNJ*(MJ-1)
    DO NN=2,PNL-1
    DO KP=CNL*(NN-1)+1,CNL*NN
    K=KP
    DO IP=1,CNI
    I=IP+CNI*(MI-1)
    F11=(CNJ+1)*(CNI+1)*(KP-1)+(CNI+1)*(JP-1)+IP-1
    F12=F11+1
    F13=F12+CNI+1
    F14=F13-1
    F21=F11+(CNJ+1)*(CNI+1)
    IF (KP.EQ.NK) F21=(CNI+1)*(JP-1)+IP-1
    F22=F21+1
    F23=F22+CNI+1
    F24=F23-1
    CALL CAL_FIJ_CY(F11,F12,F13,F14,F21,F22,F23,F24)
    WRITE(2,235) '4(',F11,' ',F12,' ',F22,' ',F21,')'
    FAD1=F_DOWN_CY(IP,KP)
    WRITE(3,231) FAD1
    OWN=CELL_ID_CY(IP,JP,KP)
    WRITE(4,231) OWN
    ENDDO   
    ENDDO
    ENDDO
    
        !DOWN FACE--PROCESSOR_MI1_TO_MI2(RECTANGLE 3-1)  !CASE CYLINDER 2-1 
    JP=1
    J=JP+CNJ*(MJ-1)
    DO KP=NL-CNL+1,NL
    K=KP
    DO IP=1,CNI
    I=IP+CNI*(MI-1)
    F11=(CNJ+1)*(CNI+1)*(KP-1)+(CNI+1)*(JP-1)+IP-1
    F12=F11+1
    F13=F12+CNI+1
    F14=F13-1
    F21=F11+(CNJ+1)*(CNI+1)
    IF (KP.EQ.NK) F21=(CNI+1)*(JP-1)+IP-1
    F22=F21+1
    F23=F22+CNI+1
    F24=F23-1
    CALL CAL_FIJ_CY(F11,F12,F13,F14,F21,F22,F23,F24)
    WRITE(2,235) '4(',F11,' ',F12,' ',F22,' ',F21,')'
    FAD1=F_DOWN_CY(IP,KP)
    WRITE(3,231) FAD1
    OWN=CELL_ID_CY(IP,JP,KP)
    WRITE(4,231) OWN
    ENDDO   
    ENDDO
    DO KP=NL+1,NL+CNL
    K=KP
    DO IP=1,CNI
    I=IP+CNI*(MI-1)
    F11=(CNJ+1)*(CNI+1)*(KP-1)+(CNI+1)*(JP-1)+IP-1
    F12=F11+1
    F13=F12+CNI+1
    F14=F13-1
    F21=F11+(CNJ+1)*(CNI+1)
    IF (KP.EQ.NK) F21=(CNI+1)*(JP-1)+IP-1
    F22=F21+1
    F23=F22+CNI+1
    F24=F23-1
    CALL CAL_FIJ_CY(F11,F12,F13,F14,F21,F22,F23,F24)
    WRITE(2,235) '4(',F11,' ',F12,' ',F22,' ',F21,')'
    FAD1=F_DOWN_CY(IP,KP)
    WRITE(3,231) FAD1
    OWN=CELL_ID_CY(IP,JP,KP)
    WRITE(4,231) OWN
    ENDDO   
    ENDDO
        
        !DOWN FACE--PROCESSOR_MI1_TO_MI2(RECTANGLE 1-2 / 3-2)  !CASE CYLINDER 2-1  
    JP=1
    J=JP+CNJ*(MJ-1)
    DO NN=2,PNL-1
        
    DO KP=3*NL+CNL*(PNL-NN)+1,3*NL+CNL*(PNL-NN+1)
    K=KP
    DO IP=1,CNI
    I=IP+CNI*(MI-1)
    F11=(CNJ+1)*(CNI+1)*(KP-1)+(CNI+1)*(JP-1)+IP-1
    F12=F11+1
    F13=F12+CNI+1
    F14=F13-1
    F21=F11+(CNJ+1)*(CNI+1)
    IF (KP.EQ.NK) F21=(CNI+1)*(JP-1)+IP-1
    F22=F21+1
    F23=F22+CNI+1
    F24=F23-1
    CALL CAL_FIJ_CY(F11,F12,F13,F14,F21,F22,F23,F24)
    WRITE(2,235) '4(',F11,' ',F12,' ',F22,' ',F21,')'
    FAD1=F_DOWN_CY(IP,KP)
    WRITE(3,231) FAD1
    OWN=CELL_ID_CY(IP,JP,KP)
    WRITE(4,231) OWN
    ENDDO   
    ENDDO
    
    DO KP=NL+CNL*(NN-1)+1,NL+CNL*NN
    K=KP
    DO IP=1,CNI
    I=IP+CNI*(MI-1)
    F11=(CNJ+1)*(CNI+1)*(KP-1)+(CNI+1)*(JP-1)+IP-1
    F12=F11+1
    F13=F12+CNI+1
    F14=F13-1
    F21=F11+(CNJ+1)*(CNI+1)
    IF (KP.EQ.NK) F21=(CNI+1)*(JP-1)+IP-1
    F22=F21+1
    F23=F22+CNI+1
    F24=F23-1
    CALL CAL_FIJ_CY(F11,F12,F13,F14,F21,F22,F23,F24)
    WRITE(2,235) '4(',F11,' ',F12,' ',F22,' ',F21,')'
    FAD1=F_DOWN_CY(IP,KP)
    WRITE(3,231) FAD1
    OWN=CELL_ID_CY(IP,JP,KP)
    WRITE(4,231) OWN
    ENDDO   
    ENDDO
    
    ENDDO    
        
        !DOWN FACE--PROCESSOR_MI1_TO_MI2(RECTANGLE 1-3)  !CASE CYLINDER 2-1 
    JP=1
    J=JP+CNJ*(MJ-1)
    DO KP=3*NL-CNL+1,3*NL
    K=KP
    DO IP=1,CNI
    I=IP+CNI*(MI-1)
    F11=(CNJ+1)*(CNI+1)*(KP-1)+(CNI+1)*(JP-1)+IP-1
    F12=F11+1
    F13=F12+CNI+1
    F14=F13-1
    F21=F11+(CNJ+1)*(CNI+1)
    IF (KP.EQ.NK) F21=(CNI+1)*(JP-1)+IP-1
    F22=F21+1
    F23=F22+CNI+1
    F24=F23-1
    CALL CAL_FIJ_CY(F11,F12,F13,F14,F21,F22,F23,F24)
    WRITE(2,235) '4(',F11,' ',F12,' ',F22,' ',F21,')'
    FAD1=F_DOWN_CY(IP,KP)
    WRITE(3,231) FAD1
    OWN=CELL_ID_CY(IP,JP,KP)
    WRITE(4,231) OWN
    ENDDO   
    ENDDO
    DO KP=3*NL+1,3*NL+CNL
    K=KP
    DO IP=1,CNI
    I=IP+CNI*(MI-1)
    F11=(CNJ+1)*(CNI+1)*(KP-1)+(CNI+1)*(JP-1)+IP-1
    F12=F11+1
    F13=F12+CNI+1
    F14=F13-1
    F21=F11+(CNJ+1)*(CNI+1)
    IF (KP.EQ.NK) F21=(CNI+1)*(JP-1)+IP-1
    F22=F21+1
    F23=F22+CNI+1
    F24=F23-1
    CALL CAL_FIJ_CY(F11,F12,F13,F14,F21,F22,F23,F24)
    WRITE(2,235) '4(',F11,' ',F12,' ',F22,' ',F21,')'
    FAD1=F_DOWN_CY(IP,KP)
    WRITE(3,231) FAD1
    OWN=CELL_ID_CY(IP,JP,KP)
    WRITE(4,231) OWN
    ENDDO   
    ENDDO
    
        !DOWN FACE--PROCESSOR_MI1_TO_MI2(RECTANGLE 2-3)  !CASE CYLINDER 2-1  
    JP=1
    J=JP+CNJ*(MJ-1)
    DO NN=2,PNL-1
    DO KP=2*NL+CNL*(PNL-NN)+1,2*NL+CNL*(PNL-NN+1)
    K=KP
    DO IP=1,CNI
    I=IP+CNI*(MI-1)
    F11=(CNJ+1)*(CNI+1)*(KP-1)+(CNI+1)*(JP-1)+IP-1
    F12=F11+1
    F13=F12+CNI+1
    F14=F13-1
    F21=F11+(CNJ+1)*(CNI+1)
    IF (KP.EQ.NK) F21=(CNI+1)*(JP-1)+IP-1
    F22=F21+1
    F23=F22+CNI+1
    F24=F23-1
    CALL CAL_FIJ_CY(F11,F12,F13,F14,F21,F22,F23,F24)
    WRITE(2,235) '4(',F11,' ',F12,' ',F22,' ',F21,')'
    FAD1=F_DOWN_CY(IP,KP)
    WRITE(3,231) FAD1
    OWN=CELL_ID_CY(IP,JP,KP)
    WRITE(4,231) OWN
    ENDDO   
    ENDDO
    ENDDO
    
        !DOWN FACE--PROCESSOR_MI1_TO_MI2(RECTANGLE 3-3)  !CASE CYLINDER 2-1 
    JP=1
    J=JP+CNJ*(MJ-1)
    DO KP=2*NL-CNL+1,2*NL
    K=KP
    DO IP=1,CNI
    I=IP+CNI*(MI-1)
    F11=(CNJ+1)*(CNI+1)*(KP-1)+(CNI+1)*(JP-1)+IP-1
    F12=F11+1
    F13=F12+CNI+1
    F14=F13-1
    F21=F11+(CNJ+1)*(CNI+1)
    IF (KP.EQ.NK) F21=(CNI+1)*(JP-1)+IP-1
    F22=F21+1
    F23=F22+CNI+1
    F24=F23-1
    CALL CAL_FIJ_CY(F11,F12,F13,F14,F21,F22,F23,F24)
    WRITE(2,235) '4(',F11,' ',F12,' ',F22,' ',F21,')'
    FAD1=F_DOWN_CY(IP,KP)
    WRITE(3,231) FAD1
    OWN=CELL_ID_CY(IP,JP,KP)
    WRITE(4,231) OWN
    ENDDO   
    ENDDO
    DO KP=2*NL+1,2*NL+CNL
    K=KP
    DO IP=1,CNI
    I=IP+CNI*(MI-1)
    F11=(CNJ+1)*(CNI+1)*(KP-1)+(CNI+1)*(JP-1)+IP-1
    F12=F11+1
    F13=F12+CNI+1
    F14=F13-1
    F21=F11+(CNJ+1)*(CNI+1)
    IF (KP.EQ.NK) F21=(CNI+1)*(JP-1)+IP-1
    F22=F21+1
    F23=F22+CNI+1
    F24=F23-1
    CALL CAL_FIJ_CY(F11,F12,F13,F14,F21,F22,F23,F24)
    WRITE(2,235) '4(',F11,' ',F12,' ',F22,' ',F21,')'
    FAD1=F_DOWN_CY(IP,KP)
    WRITE(3,231) FAD1
    OWN=CELL_ID_CY(IP,JP,KP)
    WRITE(4,231) OWN
    ENDDO   
    ENDDO
    
    N=3+4*(PNL-1)
    WRITE(31,231) N+4
    WRITE(31,*) '('  
    WRITE(31,*) '0'
    WRITE(31,*) '1'
    WRITE(31,*) '2' 
    WRITE(31,*) '3' 
    DO I=1,N
    WRITE(31,*) '-1' 
    ENDDO
    WRITE(31,*) ')'
    
    WRITE(21,231) N+4
    WRITE(21,*) "("    
    WRITE(21,*) "    INLET1"
    WRITE(21,*) "    {"    
    WRITE(21,*) "        type            patch;"
    N=0
    WRITE(21,232) "        nFaces          ",N,";"
    N=F_INTERNAL_NUM
    WRITE(21,232) "        startFace       ",N,";"
    WRITE(21,*) "    }"
    WRITE(21,*) "    INLET2"
    WRITE(21,*) "    {"    
    WRITE(21,*) "        type            patch;"
    N=0
    WRITE(21,232) "        nFaces          ",N,";"
    N=F_INTERNAL_NUM
    WRITE(21,232) "        startFace       ",N,";"
    WRITE(21,*) "    }"
    WRITE(21,*) "    OUTLET"
    WRITE(21,*) "    {"
    WRITE(21,*) "        type            patch;"
    N=0
    WRITE(21,232) "        nFaces          ",N,";"
    N=F_INTERNAL_NUM
    WRITE(21,232) "        startFace       ",N,";"
    WRITE(21,*) "    }"
    WRITE(21,*) "    HULL"
    WRITE(21,*) "    {"
    WRITE(21,*) "        type            wall;"
    WRITE(21,*) "        inGroups        1(wall);"
    N=0
    WRITE(21,232) "        nFaces          ",N,";"
    N=F_INTERNAL_NUM
    WRITE(21,232) "        startFace       ",N,";"
    WRITE(21,*) "    }"
    
        !RIGHT FACE--PROCESSOR_MI1_TO_MI2  !CASE CYLINDER 2-1 
    PRO1=PNI*(MJ-1)+MI-1
    PRO2=PNI*(MJ-1)+MI-2
    WRITE(21,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(21,*) '    {'
    WRITE(21,*) '        type            processor;'
    WRITE(21,*) '        inGroups        List<word> 1(processor);'
    N=CNJ*NK
    WRITE(21,232) "        nFaces          ",N,";"
    N=F_INTERNAL_NUM
    WRITE(21,232) "        startFace       ",N,";"
    WRITE(21,*) '        matchTolerance  0.0001;'
    WRITE(21,232) "        myProcNo        ",PRO1,";"
    WRITE(21,232) "        neighbProcNo    ",PRO2,";"
    WRITE(21,*) '    }'
    WRITE(41,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(41,*) '    {'
    WRITE(41,*) '        type            processor;'
    WRITE(41,*) '        value           nonuniform List<scalar> '
    N=CNJ*NK
    WRITE(41,'(I0)') N
    WRITE(41,*) '('
    DO KP=1,NK
    DO JP=1,CNJ
    WRITE(41,'(E21.14)') P_RIGHT_CY(JP,KP)    
    ENDDO
    ENDDO
    WRITE(41,*) ')'
    WRITE(41,*) ';'
    WRITE(41,*) '    }'
    WRITE(51,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(51,*) '    {'
    WRITE(51,*) '        type            processor;'
    WRITE(51,*) '        value           nonuniform List<vector> '
    N=CNJ*NK
    WRITE(51,'(I0)') N
    WRITE(51,*) '('
    DO KP=1,NK
    DO JP=1,CNJ
    UU1=U_RIGHT_CY(1,JP,KP)
    UU2=U_RIGHT_CY(2,JP,KP)
    UU3=U_RIGHT_CY(3,JP,KP)
    WRITE(51,230) '( ',UU1,' ',UU2,' ',UU3,' )'  
    ENDDO
    ENDDO
    WRITE(51,*) ')'
    WRITE(51,*) ';'
    WRITE(51,*) '    }'
    WRITE(61,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(61,*) '    {'
    WRITE(61,*) '        type            processor;'
    WRITE(61,*) "        value           uniform 0;"
    WRITE(61,*) '    }'
        !LEFT FACE--PROCESSOR_MI1_TO_MI2  !CASE CYLINDER 2-1 
    PRO2=PNI*(MJ-1)+MI
    WRITE(21,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(21,*) '    {'
    WRITE(21,*) '        type            processor;'
    WRITE(21,*) '        inGroups        List<word> 1(processor);'
    N=CNJ*NK
    WRITE(21,232) "        nFaces          ",N,";"
    N=F_INTERNAL_NUM+CNJ*NK
    WRITE(21,232) "        startFace       ",N,";"
    WRITE(21,*) '        matchTolerance  0.0001;'
    WRITE(21,232) "        myProcNo        ",PRO1,";"
    WRITE(21,232) "        neighbProcNo    ",PRO2,";"
    WRITE(21,*) '    }'
    WRITE(41,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(41,*) '    {'
    WRITE(41,*) '        type            processor;'
    WRITE(41,*) '        value           nonuniform List<scalar> '
    N=CNJ*NK
    WRITE(41,'(I0)') N
    WRITE(41,*) '('
    DO KP=1,NK
    DO JP=1,CNJ
    WRITE(41,'(E21.14)') P_LEFT_CY(JP,KP)    
    ENDDO
    ENDDO
    WRITE(41,*) ')'
    WRITE(41,*) ';'
    WRITE(41,*) '    }'
    WRITE(51,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(51,*) '    {'
    WRITE(51,*) '        type            processor;'
    WRITE(51,*) '        value           nonuniform List<vector> '
    N=CNJ*NK
    WRITE(51,'(I0)') N
    WRITE(51,*) '('
    DO KP=1,NK
    DO JP=1,CNJ
    UU1=U_LEFT_CY(1,JP,KP)
    UU2=U_LEFT_CY(2,JP,KP)
    UU3=U_LEFT_CY(3,JP,KP)
    WRITE(51,230) '( ',UU1,' ',UU2,' ',UU3,' )'  
    ENDDO
    ENDDO
    WRITE(51,*) ')'
    WRITE(51,*) ';'
    WRITE(51,*) '    }'
    WRITE(61,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(61,*) '    {'
    WRITE(61,*) '        type            processor;'
    WRITE(61,*) "        value           uniform 0;"
    WRITE(61,*) '    }'
        !UP FACE--PROCESSOR_MI1_TO_MI2  !CASE CYLINDER 2-1
    PRO2=PNI*MJ+MI-1
    WRITE(21,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(21,*) '    {'
    WRITE(21,*) '        type            processor;'
    WRITE(21,*) '        inGroups        List<word> 1(processor);'
    N=CNI*NK
    WRITE(21,232) "        nFaces          ",N,";"
    N=F_INTERNAL_NUM+CNJ*NK+CNJ*NK
    WRITE(21,232) "        startFace       ",N,";"
    WRITE(21,*) '        matchTolerance  0.0001;'
    WRITE(21,232) "        myProcNo        ",PRO1,";"
    WRITE(21,232) "        neighbProcNo    ",PRO2,";"
    WRITE(21,*) '    }'
    WRITE(41,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(41,*) '    {'
    WRITE(41,*) '        type            processor;'
    WRITE(41,*) '        value           nonuniform List<scalar> '
    N=CNI*NK
    WRITE(41,'(I0)') N
    WRITE(41,*) '('
    DO KP=1,NK
    DO IP=1,CNI
    WRITE(41,'(E21.14)') P_LEFT_CY(IP,KP)    
    ENDDO
    ENDDO
    WRITE(41,*) ')'
    WRITE(41,*) ';'
    WRITE(41,*) '    }'
    WRITE(51,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(51,*) '    {'
    WRITE(51,*) '        type            processor;'
    WRITE(51,*) '        value           nonuniform List<vector> '
    N=CNI*NK
    WRITE(51,'(I0)') N
    WRITE(51,*) '('
    DO KP=1,NK
    DO IP=1,CNI
    UU1=U_UP_CY(1,IP,KP)
    UU2=U_UP_CY(2,IP,KP)
    UU3=U_UP_CY(3,IP,KP)
    WRITE(51,230) '( ',UU1,' ',UU2,' ',UU3,' )'  
    ENDDO
    ENDDO
    WRITE(51,*) ')'
    WRITE(51,*) ';'
    WRITE(51,*) '    }'
    WRITE(61,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(61,*) '    {'
    WRITE(61,*) '        type            processor;'
    WRITE(61,*) "        value           uniform 0;"
    WRITE(61,*) '    }'
        !DOWN FACE--PROCESSOR_MI1_TO_MI2(RECTANGLE 1-1)  !CASE CYLINDER 2-1
    MJ1=1
    MJ2=1
    PRO2=PNI*PNJ+PNL*PNI*(MJ2-1)+PNI*(MJ1-1)+MI-1
    WRITE(21,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(21,*) '    {'
    WRITE(21,*) '        type            processor;'
    WRITE(21,*) '        inGroups        List<word> 1(processor);'
    N=CNI*CNL*2
    WRITE(21,232) "        nFaces          ",N,";"
    N=F_INTERNAL_NUM+CNJ*NK+CNJ*NK+CNI*NK
    WRITE(21,232) "        startFace       ",N,";"
    WRITE(21,*) '        matchTolerance  0.0001;'
    WRITE(21,232) "        myProcNo        ",PRO1,";"
    WRITE(21,232) "        neighbProcNo    ",PRO2,";"
    WRITE(21,*) '    }'
    WRITE(41,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(41,*) '    {'
    WRITE(41,*) '        type            processor;'
    WRITE(41,*) '        value           nonuniform List<scalar> '
    N=CNI*CNL*2
    WRITE(41,'(I0)') N
    WRITE(41,*) '('
    DO KP=1,CNL
    DO IP=1,CNI
    WRITE(41,'(E21.14)') P_DOWN_CY(IP,KP)    
    ENDDO
    ENDDO
    DO KP=NK-CNL+1,NK
    DO IP=1,CNI
    WRITE(41,'(E21.14)') P_DOWN_CY(IP,KP)    
    ENDDO
    ENDDO
    WRITE(41,*) ')'
    WRITE(41,*) ';'
    WRITE(41,*) '    }'
    WRITE(51,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(51,*) '    {'
    WRITE(51,*) '        type            processor;'
    WRITE(51,*) '        value           nonuniform List<vector> '
    N=CNI*CNL*2
    WRITE(51,'(I0)') N
    WRITE(51,*) '('
    DO KP=1,CNL
    DO IP=1,CNI
    UU1=U_DOWN_CY(1,IP,KP)
    UU2=U_DOWN_CY(2,IP,KP)
    UU3=U_DOWN_CY(3,IP,KP)
    WRITE(51,230) '( ',UU1,' ',UU2,' ',UU3,' )'  
    ENDDO
    ENDDO
    DO KP=NK-CNL+1,NK
    DO IP=1,CNI
    UU1=U_DOWN_CY(1,IP,KP)
    UU2=U_DOWN_CY(2,IP,KP)
    UU3=U_DOWN_CY(3,IP,KP)
    WRITE(51,230) '( ',UU1,' ',UU2,' ',UU3,' )'  
    ENDDO
    ENDDO
    WRITE(51,*) ')'
    WRITE(51,*) ';'
    WRITE(51,*) '    }'
    WRITE(61,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(61,*) '    {'
    WRITE(61,*) '        type            processor;'
    WRITE(61,*) "        value           uniform 0;"
    WRITE(61,*) '    }'
    
        !DOWN FACE--PROCESSOR_MI1_TO_MI2(RECTANGLE 2-1)  !CASE CYLINDER 2-1
    DO MJ1=2,PNL-1
    MJ2=1
    PRO2=PNI*PNJ+PNL*PNI*(MJ2-1)+PNI*(MJ1-1)+MI-1
    WRITE(21,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(21,*) '    {'
    WRITE(21,*) '        type            processor;'
    WRITE(21,*) '        inGroups        List<word> 1(processor);'
    N=CNI*CNL
    WRITE(21,232) "        nFaces          ",N,";"
    N=F_INTERNAL_NUM+CNJ*NK+CNJ*NK+CNI*NK+CNI*CNL*2+(MJ1-2)*CNI*CNL
    WRITE(21,232) "        startFace       ",N,";"
    WRITE(21,*) '        matchTolerance  0.0001;'
    WRITE(21,232) "        myProcNo        ",PRO1,";"
    WRITE(21,232) "        neighbProcNo    ",PRO2,";"
    WRITE(21,*) '    }'
    WRITE(41,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(41,*) '    {'
    WRITE(41,*) '        type            processor;'
    WRITE(41,*) '        value           nonuniform List<scalar> '
    N=CNI*CNL
    WRITE(41,'(I0)') N
    WRITE(41,*) '('
    DO KP=CNL*(MJ1-1)+1,CNL*MJ1
    DO IP=1,CNI
    WRITE(41,'(E21.14)') P_DOWN_CY(IP,KP)    
    ENDDO
    ENDDO
    WRITE(41,*) ')'
    WRITE(41,*) ';'
    WRITE(41,*) '    }'
    WRITE(51,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(51,*) '    {'
    WRITE(51,*) '        type            processor;'
    WRITE(51,*) '        value           nonuniform List<vector> '
    N=CNI*CNL
    WRITE(51,'(I0)') N
    WRITE(51,*) '('
    DO KP=CNL*(MJ1-1)+1,CNL*MJ1
    DO IP=1,CNI
    UU1=U_DOWN_CY(1,IP,KP)
    UU2=U_DOWN_CY(2,IP,KP)
    UU3=U_DOWN_CY(3,IP,KP)
    WRITE(51,230) '( ',UU1,' ',UU2,' ',UU3,' )'  
    ENDDO
    ENDDO
    WRITE(51,*) ')'
    WRITE(51,*) ';'
    WRITE(51,*) '    }'
    WRITE(61,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(61,*) '    {'
    WRITE(61,*) '        type            processor;'
    WRITE(61,*) "        value           uniform 0;"
    WRITE(61,*) '    }'    
    ENDDO
    
        !DOWN FACE--PROCESSOR_MI1_TO_MI2(RECTANGLE 3-1)  !CASE CYLINDER 2-1
    MJ1=PNL
    MJ2=1
    PRO2=PNI*PNJ+PNL*PNI*(MJ2-1)+PNI*(MJ1-1)+MI-1
    WRITE(21,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(21,*) '    {'
    WRITE(21,*) '        type            processor;'
    WRITE(21,*) '        inGroups        List<word> 1(processor);'
    N=CNI*CNL*2
    WRITE(21,232) "        nFaces          ",N,";"
    N=F_INTERNAL_NUM+CNJ*NK+CNJ*NK+CNI*NK+CNI*CNL*2+(PNL-2)*CNI*CNL
    WRITE(21,232) "        startFace       ",N,";"
    WRITE(21,*) '        matchTolerance  0.0001;'
    WRITE(21,232) "        myProcNo        ",PRO1,";"
    WRITE(21,232) "        neighbProcNo    ",PRO2,";"
    WRITE(21,*) '    }'
    WRITE(41,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(41,*) '    {'
    WRITE(41,*) '        type            processor;'
    WRITE(41,*) '        value           nonuniform List<scalar> '
    N=CNI*CNL*2
    WRITE(41,'(I0)') N
    WRITE(41,*) '('
    DO KP=NL-CNL+1,NL
    DO IP=1,CNI
    WRITE(41,'(E21.14)') P_DOWN_CY(IP,KP)    
    ENDDO
    ENDDO
    DO KP=NL+1,NL+CNL
    DO IP=1,CNI
    WRITE(41,'(E21.14)') P_DOWN_CY(IP,KP)    
    ENDDO
    ENDDO
    WRITE(41,*) ')'
    WRITE(41,*) ';'
    WRITE(41,*) '    }'
    WRITE(51,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(51,*) '    {'
    WRITE(51,*) '        type            processor;'
    WRITE(51,*) '        value           nonuniform List<vector> '
    N=CNI*CNL*2
    WRITE(51,'(I0)') N
    WRITE(51,*) '('
    DO KP=NL-CNL+1,NL
    DO IP=1,CNI
    UU1=U_DOWN_CY(1,IP,KP)
    UU2=U_DOWN_CY(2,IP,KP)
    UU3=U_DOWN_CY(3,IP,KP)
    WRITE(51,230) '( ',UU1,' ',UU2,' ',UU3,' )'  
    ENDDO
    ENDDO
    DO KP=NL+1,NL+CNL
    DO IP=1,CNI
    UU1=U_DOWN_CY(1,IP,KP)
    UU2=U_DOWN_CY(2,IP,KP)
    UU3=U_DOWN_CY(3,IP,KP)
    WRITE(51,230) '( ',UU1,' ',UU2,' ',UU3,' )'  
    ENDDO
    ENDDO
    WRITE(51,*) ')'
    WRITE(51,*) ';'
    WRITE(51,*) '    }'
    WRITE(61,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(61,*) '    {'
    WRITE(61,*) '        type            processor;'
    WRITE(61,*) "        value           uniform 0;"
    WRITE(61,*) '    }'
    
        !DOWN FACE--PROCESSOR_MI1_TO_MI2(RECTANGLE 1-2 / 3-2)  !CASE CYLINDER 2-1
    DO MJ2=2,PNL-1
    MJ1=1
    PRO2=PNI*PNJ+PNL*PNI*(MJ2-1)+PNI*(MJ1-1)+MI-1
    WRITE(21,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(21,*) '    {'
    WRITE(21,*) '        type            processor;'
    WRITE(21,*) '        inGroups        List<word> 1(processor);'
    N=CNI*CNL
    WRITE(21,232) "        nFaces          ",N,";"
    N=F_INTERNAL_NUM+CNJ*NK+CNJ*NK+CNI*NK+CNI*CNL*2+(PNL-2)*CNI*CNL+CNI*CNL*2+(MJ2-2)*CNI*CNL*2
    WRITE(21,232) "        startFace       ",N,";"
    WRITE(21,*) '        matchTolerance  0.0001;'
    WRITE(21,232) "        myProcNo        ",PRO1,";"
    WRITE(21,232) "        neighbProcNo    ",PRO2,";"
    WRITE(21,*) '    }'
    WRITE(41,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(41,*) '    {'
    WRITE(41,*) '        type            processor;'
    WRITE(41,*) '        value           nonuniform List<scalar> '
    N=CNI*CNL
    WRITE(41,'(I0)') N
    WRITE(41,*) '('
    DO KP=3*NL+CNL*(PNL-MJ2)+1,3*NL+CNL*(PNL-MJ2+1)
    DO IP=1,CNI
    WRITE(41,'(E21.14)') P_DOWN_CY(IP,KP)    
    ENDDO
    ENDDO
    WRITE(41,*) ')'
    WRITE(41,*) ';'
    WRITE(41,*) '    }'
    WRITE(51,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(51,*) '    {'
    WRITE(51,*) '        type            processor;'
    WRITE(51,*) '        value           nonuniform List<vector> '
    N=CNI*CNL
    WRITE(51,'(I0)') N
    WRITE(51,*) '('
    DO KP=3*NL+CNL*(PNL-MJ2)+1,3*NL+CNL*(PNL-MJ2+1)
    DO IP=1,CNI
    UU1=U_DOWN_CY(1,IP,KP)
    UU2=U_DOWN_CY(2,IP,KP)
    UU3=U_DOWN_CY(3,IP,KP)
    WRITE(51,230) '( ',UU1,' ',UU2,' ',UU3,' )'  
    ENDDO
    ENDDO
    WRITE(51,*) ')'
    WRITE(51,*) ';'
    WRITE(51,*) '    }'
    WRITE(61,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(61,*) '    {'
    WRITE(61,*) '        type            processor;'
    WRITE(61,*) "        value           uniform 0;"
    WRITE(61,*) '    }' 
    
    MJ1=PNL
    PRO2=PNI*PNJ+PNL*PNI*(MJ2-1)+PNI*(MJ1-1)+MI-1
    WRITE(21,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(21,*) '    {'
    WRITE(21,*) '        type            processor;'
    WRITE(21,*) '        inGroups        List<word> 1(processor);'
    N=CNI*CNL
    WRITE(21,232) "        nFaces          ",N,";"
    N=F_INTERNAL_NUM+CNJ*NK+CNJ*NK+CNI*NK+CNI*CNL*2+(PNL-2)*CNI*CNL+CNI*CNL*2+(MJ2-2)*CNI*CNL*2+CNI*CNL
    WRITE(21,232) "        startFace       ",N,";"
    WRITE(21,*) '        matchTolerance  0.0001;'
    WRITE(21,232) "        myProcNo        ",PRO1,";"
    WRITE(21,232) "        neighbProcNo    ",PRO2,";"
    WRITE(21,*) '    }'
    WRITE(41,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(41,*) '    {'
    WRITE(41,*) '        type            processor;'
    WRITE(41,*) '        value           nonuniform List<scalar> '
    N=CNI*CNL
    WRITE(41,'(I0)') N
    WRITE(41,*) '('
    DO KP=NL+CNL*(MJ2-1)+1,NL+CNL*MJ2
    DO IP=1,CNI
    WRITE(41,'(E21.14)') P_DOWN_CY(IP,KP)    
    ENDDO
    ENDDO
    WRITE(41,*) ')'
    WRITE(41,*) ';'
    WRITE(41,*) '    }'
    WRITE(51,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(51,*) '    {'
    WRITE(51,*) '        type            processor;'
    WRITE(51,*) '        value           nonuniform List<vector> '
    N=CNI*CNL
    WRITE(51,'(I0)') N
    WRITE(51,*) '('
    DO KP=NL+CNL*(MJ2-1)+1,NL+CNL*MJ2
    DO IP=1,CNI
    UU1=U_DOWN_CY(1,IP,KP)
    UU2=U_DOWN_CY(2,IP,KP)
    UU3=U_DOWN_CY(3,IP,KP)
    WRITE(51,230) '( ',UU1,' ',UU2,' ',UU3,' )'  
    ENDDO
    ENDDO
    WRITE(51,*) ')'
    WRITE(51,*) ';'
    WRITE(51,*) '    }'
    WRITE(61,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(61,*) '    {'
    WRITE(61,*) '        type            processor;'
    WRITE(61,*) "        value           uniform 0;"
    WRITE(61,*) '    }'
        
    ENDDO
    
    
        !DOWN FACE--PROCESSOR_MI1_TO_MI2(RECTANGLE 1-3)  !CASE CYLINDER 2-1
    MJ1=1
    MJ2=PNL
    PRO2=PNI*PNJ+PNL*PNI*(MJ2-1)+PNI*(MJ1-1)+MI-1
    WRITE(21,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(21,*) '    {'
    WRITE(21,*) '        type            processor;'
    WRITE(21,*) '        inGroups        List<word> 1(processor);'
    N=CNI*CNL*2
    WRITE(21,232) "        nFaces          ",N,";"
    N=F_INTERNAL_NUM+CNJ*NK+CNJ*NK+CNI*NK+CNI*CNL*2+(PNL-2)*CNI*CNL+CNI*CNL*2+(PNL-2)*CNI*CNL*2
    WRITE(21,232) "        startFace       ",N,";"
    WRITE(21,*) '        matchTolerance  0.0001;'
    WRITE(21,232) "        myProcNo        ",PRO1,";"
    WRITE(21,232) "        neighbProcNo    ",PRO2,";"
    WRITE(21,*) '    }'
    WRITE(41,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(41,*) '    {'
    WRITE(41,*) '        type            processor;'
    WRITE(41,*) '        value           nonuniform List<scalar> '
    N=CNI*CNL*2
    WRITE(41,'(I0)') N
    WRITE(41,*) '('
    DO KP=3*NL-CNL+1,3*NL
    DO IP=1,CNI
    WRITE(41,'(E21.14)') P_DOWN_CY(IP,KP)    
    ENDDO
    ENDDO
    DO KP=3*NL+1,3*NL+CNL
    DO IP=1,CNI
    WRITE(41,'(E21.14)') P_DOWN_CY(IP,KP)    
    ENDDO
    ENDDO
    WRITE(41,*) ')'
    WRITE(41,*) ';'
    WRITE(41,*) '    }'
    WRITE(51,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(51,*) '    {'
    WRITE(51,*) '        type            processor;'
    WRITE(51,*) '        value           nonuniform List<vector> '
    N=CNI*CNL*2
    WRITE(51,'(I0)') N
    WRITE(51,*) '('
    DO KP=3*NL-CNL+1,3*NL
    DO IP=1,CNI
    UU1=U_DOWN_CY(1,IP,KP)
    UU2=U_DOWN_CY(2,IP,KP)
    UU3=U_DOWN_CY(3,IP,KP)
    WRITE(51,230) '( ',UU1,' ',UU2,' ',UU3,' )'  
    ENDDO
    ENDDO
    DO KP=3*NL+1,3*NL+CNL
    DO IP=1,CNI
    UU1=U_DOWN_CY(1,IP,KP)
    UU2=U_DOWN_CY(2,IP,KP)
    UU3=U_DOWN_CY(3,IP,KP)
    WRITE(51,230) '( ',UU1,' ',UU2,' ',UU3,' )'  
    ENDDO
    ENDDO
    WRITE(51,*) ')'
    WRITE(51,*) ';'
    WRITE(51,*) '    }'
    WRITE(61,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(61,*) '    {'
    WRITE(61,*) '        type            processor;'
    WRITE(61,*) "        value           uniform 0;"
    WRITE(61,*) '    }'
    
    
        !DOWN FACE--PROCESSOR_MI1_TO_MI2(RECTANGLE 2-3)  !CASE CYLINDER 2-1
    DO MJ1=2,PNL-1
    MJ2=PNL
    PRO2=PNI*PNJ+PNL*PNI*(MJ2-1)+PNI*(MJ1-1)+MI-1
    WRITE(21,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(21,*) '    {'
    WRITE(21,*) '        type            processor;'
    WRITE(21,*) '        inGroups        List<word> 1(processor);'
    N=CNI*CNL
    WRITE(21,232) "        nFaces          ",N,";"
    N=F_INTERNAL_NUM+CNJ*NK+CNJ*NK+CNI*NK+CNI*CNL*2+(PNL-2)*CNI*CNL+CNI*CNL*2+(PNL-2)*CNI*CNL*2+CNI*CNL*2+(MJ1-2)*CNI*CNL
    WRITE(21,232) "        startFace       ",N,";"
    WRITE(21,*) '        matchTolerance  0.0001;'
    WRITE(21,232) "        myProcNo        ",PRO1,";"
    WRITE(21,232) "        neighbProcNo    ",PRO2,";"
    WRITE(21,*) '    }'
    WRITE(41,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(41,*) '    {'
    WRITE(41,*) '        type            processor;'
    WRITE(41,*) '        value           nonuniform List<scalar> '
    N=CNI*CNL
    WRITE(41,'(I0)') N
    WRITE(41,*) '('
    DO KP=2*NL+CNL*(PNL-MJ1)+1,2*NL+CNL*(PNL-MJ1+1)
    DO IP=1,CNI
    WRITE(41,'(E21.14)') P_DOWN_CY(IP,KP)    
    ENDDO
    ENDDO
    WRITE(41,*) ')'
    WRITE(41,*) ';'
    WRITE(41,*) '    }'
    WRITE(51,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(51,*) '    {'
    WRITE(51,*) '        type            processor;'
    WRITE(51,*) '        value           nonuniform List<vector> '
    N=CNI*CNL
    WRITE(51,'(I0)') N
    WRITE(51,*) '('
    DO KP=2*NL+CNL*(PNL-MJ1)+1,2*NL+CNL*(PNL-MJ1+1)
    DO IP=1,CNI
    UU1=U_DOWN_CY(1,IP,KP)
    UU2=U_DOWN_CY(2,IP,KP)
    UU3=U_DOWN_CY(3,IP,KP)
    WRITE(51,230) '( ',UU1,' ',UU2,' ',UU3,' )'  
    ENDDO
    ENDDO
    WRITE(51,*) ')'
    WRITE(51,*) ';'
    WRITE(51,*) '    }'
    WRITE(61,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(61,*) '    {'
    WRITE(61,*) '        type            processor;'
    WRITE(61,*) "        value           uniform 0;"
    WRITE(61,*) '    }'    
    ENDDO
    
    
    
        !DOWN FACE--PROCESSOR_MI1_TO_MI2(RECTANGLE 3-3)  !CASE CYLINDER 2-1
    MJ1=PNL
    MJ2=PNL
    PRO2=PNI*PNJ+PNL*PNI*(MJ2-1)+PNI*(MJ1-1)+MI-1
    WRITE(21,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(21,*) '    {'
    WRITE(21,*) '        type            processor;'
    WRITE(21,*) '        inGroups        List<word> 1(processor);'
    N=CNI*CNL*2
    WRITE(21,232) "        nFaces          ",N,";"
    N=F_INTERNAL_NUM+CNJ*NK+CNJ*NK+CNI*NK+CNI*CNL*2+(PNL-2)*CNI*CNL+CNI*CNL*2+(PNL-2)*CNI*CNL*2+CNI*CNL*2+(PNL-2)*CNI*CNL
    WRITE(21,232) "        startFace       ",N,";"
    WRITE(21,*) '        matchTolerance  0.0001;'
    WRITE(21,232) "        myProcNo        ",PRO1,";"
    WRITE(21,232) "        neighbProcNo    ",PRO2,";"
    WRITE(21,*) '    }'
    WRITE(41,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(41,*) '    {'
    WRITE(41,*) '        type            processor;'
    WRITE(41,*) '        value           nonuniform List<scalar> '
    N=CNI*CNL*2
    WRITE(41,'(I0)') N
    WRITE(41,*) '('
    DO KP=2*NL-CNL+1,2*NL
    DO IP=1,CNI
    WRITE(41,'(E21.14)') P_DOWN_CY(IP,KP)    
    ENDDO
    ENDDO
    DO KP=2*NL+1,2*NL+CNL
    DO IP=1,CNI
    WRITE(41,'(E21.14)') P_DOWN_CY(IP,KP)    
    ENDDO
    ENDDO
    WRITE(41,*) ')'
    WRITE(41,*) ';'
    WRITE(41,*) '    }'
    WRITE(51,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(51,*) '    {'
    WRITE(51,*) '        type            processor;'
    WRITE(51,*) '        value           nonuniform List<vector> '
    N=CNI*CNL*2
    WRITE(51,'(I0)') N
    WRITE(51,*) '('
    DO KP=2*NL-CNL+1,2*NL
    DO IP=1,CNI
    UU1=U_DOWN_CY(1,IP,KP)
    UU2=U_DOWN_CY(2,IP,KP)
    UU3=U_DOWN_CY(3,IP,KP)
    WRITE(51,230) '( ',UU1,' ',UU2,' ',UU3,' )'  
    ENDDO
    ENDDO
    DO KP=2*NL+1,2*NL+CNL
    DO IP=1,CNI
    UU1=U_DOWN_CY(1,IP,KP)
    UU2=U_DOWN_CY(2,IP,KP)
    UU3=U_DOWN_CY(3,IP,KP)
    WRITE(51,230) '( ',UU1,' ',UU2,' ',UU3,' )'  
    ENDDO
    ENDDO
    WRITE(51,*) ')'
    WRITE(51,*) ';'
    WRITE(51,*) '    }'
    WRITE(61,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(61,*) '    {'
    WRITE(61,*) '        type            processor;'
    WRITE(61,*) "        value           uniform 0;"
    WRITE(61,*) '    }'
    
    
    
    
    
    
    
    
    
    
    ELSEIF (MJ.GT.1 .AND. MJ.LT.PNJ) THEN
    !CASE CYLINDER 2-2
    
        !DOWN FACE--PROCESSOR_MI1_TO_MI2  !CASE CYLINDER 2-2
    JP=1
    J=JP+CNJ*(MJ-1)
    DO KP=1,NK
    K=KP
    DO IP=1,CNI
    I=IP+CNI*(MI-1)
    F11=(CNJ+1)*(CNI+1)*(KP-1)+(CNI+1)*(JP-1)+IP-1
    F12=F11+1
    F13=F12+CNI+1
    F14=F13-1
    F21=F11+(CNJ+1)*(CNI+1)
    IF (KP.EQ.NK) F21=(CNI+1)*(JP-1)+IP-1
    F22=F21+1
    F23=F22+CNI+1
    F24=F23-1
    CALL CAL_FIJ_CY(F11,F12,F13,F14,F21,F22,F23,F24)
    WRITE(2,235) '4(',F11,' ',F12,' ',F22,' ',F21,')'
    FAD1=F_DOWN_CY(IP,KP)
    WRITE(3,231) FAD1
    OWN=CELL_ID_CY(IP,JP,KP)
    WRITE(4,231) OWN
    ENDDO   
    ENDDO
        
        !RIGHT FACE--PROCESSOR_MI1_TO_MI2  !CASE CYLINDER 2-2
    IP=1
    I=IP+CNI*(MI-1)
    DO KP=1,NK
    K=KP
    DO JP=1,CNJ
    J=JP+CNJ*(MJ-1)
    F11=(CNJ+1)*(CNI+1)*(KP-1)+(CNI+1)*(JP-1)+IP-1
    F12=F11+1
    F13=F12+CNI+1
    F14=F13-1
    F21=F11+(CNJ+1)*(CNI+1)
    IF (KP.EQ.NK) F21=(CNI+1)*(JP-1)+IP-1
    F22=F21+1
    F23=F22+CNI+1
    F24=F23-1
    CALL CAL_FIJ_CY(F11,F12,F13,F14,F21,F22,F23,F24)
    WRITE(2,235) '4(',F11,' ',F21,' ',F24,' ',F14,')'
    FAD1=F_RIGHT_CY(JP,KP)
    WRITE(3,231) FAD1
    OWN=CELL_ID_CY(IP,JP,KP)
    WRITE(4,231) OWN
    ENDDO   
    ENDDO
        
        !LEFT FACE--PROCESSOR_MI1_TO_MI2  !CASE CYLINDER 2-2
    IP=CNI
    I=IP+CNI*(MI-1)
    DO KP=1,NK
    K=KP
    DO JP=1,CNJ
    J=JP+CNJ*(MJ-1)
    F11=(CNJ+1)*(CNI+1)*(KP-1)+(CNI+1)*(JP-1)+IP-1
    F12=F11+1
    F13=F12+CNI+1
    F14=F13-1
    F21=F11+(CNJ+1)*(CNI+1)
    IF (KP.EQ.NK) F21=(CNI+1)*(JP-1)+IP-1
    F22=F21+1
    F23=F22+CNI+1
    F24=F23-1
    CALL CAL_FIJ_CY(F11,F12,F13,F14,F21,F22,F23,F24)
    WRITE(2,235) '4(',F12,' ',F13,' ',F23,' ',F22,')'
    FAD1=F_LEFT_CY(JP,KP)
    WRITE(3,231) FAD1
    OWN=CELL_ID_CY(IP,JP,KP)
    WRITE(4,231) OWN
    ENDDO   
    ENDDO
    
        !UP FACE--PROCESSOR_MI1_TO_MI2  !CASE CYLINDER 2-2
    JP=CNJ
    J=JP+CNJ*(MJ-1)
    DO KP=1,NK
    K=KP
    DO IP=1,CNI
    I=IP+CNI*(MI-1)
    F11=(CNJ+1)*(CNI+1)*(KP-1)+(CNI+1)*(JP-1)+IP-1
    F12=F11+1
    F13=F12+CNI+1
    F14=F13-1
    F21=F11+(CNJ+1)*(CNI+1)
    IF (KP.EQ.NK) F21=(CNI+1)*(JP-1)+IP-1
    F22=F21+1
    F23=F22+CNI+1
    F24=F23-1
    CALL CAL_FIJ_CY(F11,F12,F13,F14,F21,F22,F23,F24)
    WRITE(2,235) '4(',F14,' ',F24,' ',F23,' ',F13,')'
    FAD1=F_UP_CY(IP,KP)
    WRITE(3,231) FAD1
    OWN=CELL_ID_CY(IP,JP,KP)
    WRITE(4,231) OWN
    ENDDO   
    ENDDO
    
    N=4
    WRITE(31,231) N+4
    WRITE(31,*) '('  
    WRITE(31,*) '0'
    WRITE(31,*) '1'
    WRITE(31,*) '2' 
    WRITE(31,*) '3' 
    DO I=1,N
    WRITE(31,*) '-1' 
    ENDDO
    WRITE(31,*) ')'
    
    WRITE(21,231) N+4
    WRITE(21,*) "("    
    WRITE(21,*) "    INLET1"
    WRITE(21,*) "    {"    
    WRITE(21,*) "        type            patch;"
    N=0
    WRITE(21,232) "        nFaces          ",N,";"
    N=F_INTERNAL_NUM
    WRITE(21,232) "        startFace       ",N,";"
    WRITE(21,*) "    }"
    WRITE(21,*) "    INLET2"
    WRITE(21,*) "    {"    
    WRITE(21,*) "        type            patch;"
    N=0
    WRITE(21,232) "        nFaces          ",N,";"
    N=F_INTERNAL_NUM
    WRITE(21,232) "        startFace       ",N,";"
    WRITE(21,*) "    }"
    WRITE(21,*) "    OUTLET"
    WRITE(21,*) "    {"
    WRITE(21,*) "        type            patch;"
    N=0
    WRITE(21,232) "        nFaces          ",N,";"
    N=F_INTERNAL_NUM
    WRITE(21,232) "        startFace       ",N,";"
    WRITE(21,*) "    }"
    WRITE(21,*) "    HULL"
    WRITE(21,*) "    {"
    WRITE(21,*) "        type            wall;"
    WRITE(21,*) "        inGroups        1(wall);"
    N=0
    WRITE(21,232) "        nFaces          ",N,";"
    N=F_INTERNAL_NUM
    WRITE(21,232) "        startFace       ",N,";"
    WRITE(21,*) "    }"
    
        !DOWN FACE--PROCESSOR_MI1_TO_MI2 !CASE CYLINDER 2-2
    PRO1=PNI*(MJ-1)+MI-1
    PRO2=PNI*(MJ-2)+MI-1
    WRITE(21,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(21,*) '    {'
    WRITE(21,*) '        type            processor;'
    WRITE(21,*) '        inGroups        List<word> 1(processor);'
    N=CNI*NK
    WRITE(21,232) "        nFaces          ",N,";"
    N=F_INTERNAL_NUM
    WRITE(21,232) "        startFace       ",N,";"
    WRITE(21,*) '        matchTolerance  0.0001;'
    WRITE(21,232) "        myProcNo        ",PRO1,";"
    WRITE(21,232) "        neighbProcNo    ",PRO2,";"
    WRITE(21,*) '    }'
    WRITE(41,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(41,*) '    {'
    WRITE(41,*) '        type            processor;'
    WRITE(41,*) '        value           nonuniform List<scalar> '
    N=CNI*NK
    WRITE(41,'(I0)') N
    WRITE(41,*) '('
    DO KP=1,NK
    DO IP=1,CNI
    WRITE(41,'(E21.14)') P_DOWN_CY(IP,KP)    
    ENDDO
    ENDDO
    WRITE(41,*) ')'
    WRITE(41,*) ';'
    WRITE(41,*) '    }'
    WRITE(51,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(51,*) '    {'
    WRITE(51,*) '        type            processor;'
    WRITE(51,*) '        value           nonuniform List<vector> '
    N=CNI*NK
    WRITE(51,'(I0)') N
    WRITE(51,*) '('
    DO KP=1,NK
    DO IP=1,CNI
    UU1=U_DOWN_CY(1,IP,KP)
    UU2=U_DOWN_CY(2,IP,KP)
    UU3=U_DOWN_CY(3,IP,KP)
    WRITE(51,230) '( ',UU1,' ',UU2,' ',UU3,' )'  
    ENDDO
    ENDDO
    WRITE(51,*) ')'
    WRITE(51,*) ';'
    WRITE(51,*) '    }'
    WRITE(61,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(61,*) '    {'
    WRITE(61,*) '        type            processor;'
    WRITE(61,*) "        value           uniform 0;"
    WRITE(61,*) '    }'
    
        !RIGHT FACE--PROCESSOR_MI1_TO_MI2  !CASE CYLINDER 2-2
    PRO2=PNI*(MJ-1)+MI-2
    WRITE(21,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(21,*) '    {'
    WRITE(21,*) '        type            processor;'
    WRITE(21,*) '        inGroups        List<word> 1(processor);'
    N=CNJ*NK
    WRITE(21,232) "        nFaces          ",N,";"
    N=F_INTERNAL_NUM+CNI*NK
    WRITE(21,232) "        startFace       ",N,";"
    WRITE(21,*) '        matchTolerance  0.0001;'
    WRITE(21,232) "        myProcNo        ",PRO1,";"
    WRITE(21,232) "        neighbProcNo    ",PRO2,";"
    WRITE(21,*) '    }'
    WRITE(41,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(41,*) '    {'
    WRITE(41,*) '        type            processor;'
    WRITE(41,*) '        value           nonuniform List<scalar> '
    N=CNJ*NK
    WRITE(41,'(I0)') N
    WRITE(41,*) '('
    DO KP=1,NK
    DO JP=1,CNJ
    WRITE(41,'(E21.14)') P_RIGHT_CY(JP,KP)    
    ENDDO
    ENDDO
    WRITE(41,*) ')'
    WRITE(41,*) ';'
    WRITE(41,*) '    }'
    WRITE(51,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(51,*) '    {'
    WRITE(51,*) '        type            processor;'
    WRITE(51,*) '        value           nonuniform List<vector> '
    N=CNJ*NK
    WRITE(51,'(I0)') N
    WRITE(51,*) '('
    DO KP=1,NK
    DO JP=1,CNJ
    UU1=U_RIGHT_CY(1,JP,KP)
    UU2=U_RIGHT_CY(2,JP,KP)
    UU3=U_RIGHT_CY(3,JP,KP)
    WRITE(51,230) '( ',UU1,' ',UU2,' ',UU3,' )'  
    ENDDO
    ENDDO
    WRITE(51,*) ')'
    WRITE(51,*) ';'
    WRITE(51,*) '    }'
    WRITE(61,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(61,*) '    {'
    WRITE(61,*) '        type            processor;'
    WRITE(61,*) "        value           uniform 0;"
    WRITE(61,*) '    }'
    
        !LEFT FACE--PROCESSOR_MI1_TO_MI2  !CASE CYLINDER 2-2
    PRO2=PNI*(MJ-1)+MI
    WRITE(21,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(21,*) '    {'
    WRITE(21,*) '        type            processor;'
    WRITE(21,*) '        inGroups        List<word> 1(processor);'
    N=CNJ*NK
    WRITE(21,232) "        nFaces          ",N,";"
    N=F_INTERNAL_NUM+CNI*NK+CNJ*NK
    WRITE(21,232) "        startFace       ",N,";"
    WRITE(21,*) '        matchTolerance  0.0001;'
    WRITE(21,232) "        myProcNo        ",PRO1,";"
    WRITE(21,232) "        neighbProcNo    ",PRO2,";"
    WRITE(21,*) '    }'
    WRITE(41,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(41,*) '    {'
    WRITE(41,*) '        type            processor;'
    WRITE(41,*) '        value           nonuniform List<scalar> '
    N=CNJ*NK
    WRITE(41,'(I0)') N
    WRITE(41,*) '('
    DO KP=1,NK
    DO JP=1,CNJ
    WRITE(41,'(E21.14)') P_LEFT_CY(JP,KP)    
    ENDDO
    ENDDO
    WRITE(41,*) ')'
    WRITE(41,*) ';'
    WRITE(41,*) '    }'
    WRITE(51,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(51,*) '    {'
    WRITE(51,*) '        type            processor;'
    WRITE(51,*) '        value           nonuniform List<vector> '
    N=CNJ*NK
    WRITE(51,'(I0)') N
    WRITE(51,*) '('
    DO KP=1,NK
    DO JP=1,CNJ
    UU1=U_LEFT_CY(1,JP,KP)
    UU2=U_LEFT_CY(2,JP,KP)
    UU3=U_LEFT_CY(3,JP,KP)
    WRITE(51,230) '( ',UU1,' ',UU2,' ',UU3,' )'  
    ENDDO
    ENDDO
    WRITE(51,*) ')'
    WRITE(51,*) ';'
    WRITE(51,*) '    }'
    WRITE(61,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(61,*) '    {'
    WRITE(61,*) '        type            processor;'
    WRITE(61,*) "        value           uniform 0;"
    WRITE(61,*) '    }'
    
        !UP FACE--PROCESSOR_MI1_TO_MI2  !CASE CYLINDER 2-2 
    PRO2=PNI*MJ+MI-1
    WRITE(21,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(21,*) '    {'
    WRITE(21,*) '        type            processor;'
    WRITE(21,*) '        inGroups        List<word> 1(processor);'
    N=CNI*NK
    WRITE(21,232) "        nFaces          ",N,";"
    N=F_INTERNAL_NUM+CNI*NK+CNJ*NK+CNJ*NK
    WRITE(21,232) "        startFace       ",N,";"
    WRITE(21,*) '        matchTolerance  0.0001;'
    WRITE(21,232) "        myProcNo        ",PRO1,";"
    WRITE(21,232) "        neighbProcNo    ",PRO2,";"
    WRITE(21,*) '    }'
    WRITE(41,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(41,*) '    {'
    WRITE(41,*) '        type            processor;'
    WRITE(41,*) '        value           nonuniform List<scalar> '
    N=CNI*NK
    WRITE(41,'(I0)') N
    WRITE(41,*) '('
    DO KP=1,NK
    DO IP=1,CNI
    WRITE(41,'(E21.14)') P_LEFT_CY(IP,KP)    
    ENDDO
    ENDDO
    WRITE(41,*) ')'
    WRITE(41,*) ';'
    WRITE(41,*) '    }'
    WRITE(51,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(51,*) '    {'
    WRITE(51,*) '        type            processor;'
    WRITE(51,*) '        value           nonuniform List<vector> '
    N=CNI*NK
    WRITE(51,'(I0)') N
    WRITE(51,*) '('
    DO KP=1,NK
    DO IP=1,CNI
    UU1=U_UP_CY(1,IP,KP)
    UU2=U_UP_CY(2,IP,KP)
    UU3=U_UP_CY(3,IP,KP)
    WRITE(51,230) '( ',UU1,' ',UU2,' ',UU3,' )'  
    ENDDO
    ENDDO
    WRITE(51,*) ')'
    WRITE(51,*) ';'
    WRITE(51,*) '    }'
    WRITE(61,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(61,*) '    {'
    WRITE(61,*) '        type            processor;'
    WRITE(61,*) "        value           uniform 0;"
    WRITE(61,*) '    }'
    
    
    
    
    
    
    
    
    
    
    
    ELSEIF (MJ.EQ.PNJ) THEN
    !CASE CYLINDER 2-3
    
        !UP FACE--OUTLET  !CASE CYLINDER 2-3
    !��������ɾ������ʱ����
    !JP=CNJ
    !J=JP+CNJ*(MJ-1)
    !DO KP=1,NK
    !K=KP
    !DO IP=1,CNI
    !I=IP+CNI*(MI-1)
    !F11=(CNJ+1)*(CNI+1)*(KP-1)+(CNI+1)*(JP-1)+IP-1
    !F12=F11+1
    !F13=F12+CNI+1
    !F14=F13-1
    !F21=F11+(CNJ+1)*(CNI+1)
    !IF (KP.EQ.NK) F21=(CNI+1)*(JP-1)+IP-1
    !F22=F21+1
    !F23=F22+CNI+1
    !F24=F23-1
    !CALL CAL_FIJ_CY(F11,F12,F13,F14,F21,F22,F23,F24)
    !WRITE(2,235) '4(',F14,' ',F24,' ',F23,' ',F13,')'
    !FAD1=F_INTERNAL_NUM_ALL+NL*NL+NJ*NK+NI*(K-1)+I
    !WRITE(3,231) FAD1
    !OWN=CELL_ID_CY(IP,JP,KP)
    !WRITE(4,231) OWN
    !ENDDO   
    !ENDDO
    JP=CNJ
    J=JP+CNJ*(MJ-1)
    DO IP=1,CNI
    I=IP+CNI*(MI-1)
    DO KP=1,NL
    K=KP
    F11=(CNJ+1)*(CNI+1)*(KP-1)+(CNI+1)*(JP-1)+IP-1
    F12=F11+1
    F13=F12+CNI+1
    F14=F13-1
    F21=F11+(CNJ+1)*(CNI+1)
    IF (KP.EQ.NK) F21=(CNI+1)*(JP-1)+IP-1
    F22=F21+1
    F23=F22+CNI+1
    F24=F23-1
    CALL CAL_FIJ_CY(F11,F12,F13,F14,F21,F22,F23,F24)
    WRITE(2,235) '4(',F14,' ',F24,' ',F23,' ',F13,')'
    FAD1=F_INTERNAL_NUM_ALL+NL*NL+NJ*NK+NL*(I-1)+K
    WRITE(3,231) FAD1
    OWN=CELL_ID_CY(IP,JP,KP)
    WRITE(4,231) OWN
    ENDDO   
    ENDDO
    DO IP=1,CNI
    I=IP+CNI*(MI-1)
    DO KP=1+NL,2*NL
    K=KP
    F11=(CNJ+1)*(CNI+1)*(KP-1)+(CNI+1)*(JP-1)+IP-1
    F12=F11+1
    F13=F12+CNI+1
    F14=F13-1
    F21=F11+(CNJ+1)*(CNI+1)
    IF (KP.EQ.NK) F21=(CNI+1)*(JP-1)+IP-1
    F22=F21+1
    F23=F22+CNI+1
    F24=F23-1
    CALL CAL_FIJ_CY(F11,F12,F13,F14,F21,F22,F23,F24)
    WRITE(2,235) '4(',F14,' ',F24,' ',F23,' ',F13,')'
    FAD1=F_INTERNAL_NUM_ALL+NL*NL+NJ*NK+NI*NL+NL*(I-1)+(K-NL)
    WRITE(3,231) FAD1
    OWN=CELL_ID_CY(IP,JP,KP)
    WRITE(4,231) OWN
    ENDDO   
    ENDDO
    DO IP=1,CNI
    I=IP+CNI*(MI-1)
    DO KP=1+2*NL,3*NL
    K=KP
    F11=(CNJ+1)*(CNI+1)*(KP-1)+(CNI+1)*(JP-1)+IP-1
    F12=F11+1
    F13=F12+CNI+1
    F14=F13-1
    F21=F11+(CNJ+1)*(CNI+1)
    IF (KP.EQ.NK) F21=(CNI+1)*(JP-1)+IP-1
    F22=F21+1
    F23=F22+CNI+1
    F24=F23-1
    CALL CAL_FIJ_CY(F11,F12,F13,F14,F21,F22,F23,F24)
    WRITE(2,235) '4(',F14,' ',F24,' ',F23,' ',F13,')'
    FAD1=F_INTERNAL_NUM_ALL+NL*NL+NJ*NK+NI*NL*2+NL*(I-1)+(K-2*NL)
    WRITE(3,231) FAD1
    OWN=CELL_ID_CY(IP,JP,KP)
    WRITE(4,231) OWN
    ENDDO   
    ENDDO
    DO IP=1,CNI
    I=IP+CNI*(MI-1)
    DO KP=1+3*NL,NK
    K=KP
    F11=(CNJ+1)*(CNI+1)*(KP-1)+(CNI+1)*(JP-1)+IP-1
    F12=F11+1
    F13=F12+CNI+1
    F14=F13-1
    F21=F11+(CNJ+1)*(CNI+1)
    IF (KP.EQ.NK) F21=(CNI+1)*(JP-1)+IP-1
    F22=F21+1
    F23=F22+CNI+1
    F24=F23-1
    CALL CAL_FIJ_CY(F11,F12,F13,F14,F21,F22,F23,F24)
    WRITE(2,235) '4(',F14,' ',F24,' ',F23,' ',F13,')'
    FAD1=F_INTERNAL_NUM_ALL+NL*NL+NJ*NK+NI*NL*3+NL*(I-1)+(K-3*NL)
    WRITE(3,231) FAD1
    OWN=CELL_ID_CY(IP,JP,KP)
    WRITE(4,231) OWN
    ENDDO   
    ENDDO
    
        !DOWN FACE--PROCESSOR_MI1_TO_MI2  !CASE CYLINDER 2-3
    JP=1
    J=JP+CNJ*(MJ-1)
    DO KP=1,NK
    K=KP
    DO IP=1,CNI
    I=IP+CNI*(MI-1)
    F11=(CNJ+1)*(CNI+1)*(KP-1)+(CNI+1)*(JP-1)+IP-1
    F12=F11+1
    F13=F12+CNI+1
    F14=F13-1
    F21=F11+(CNJ+1)*(CNI+1)
    IF (KP.EQ.NK) F21=(CNI+1)*(JP-1)+IP-1
    F22=F21+1
    F23=F22+CNI+1
    F24=F23-1
    CALL CAL_FIJ_CY(F11,F12,F13,F14,F21,F22,F23,F24)
    WRITE(2,235) '4(',F11,' ',F12,' ',F22,' ',F21,')'
    FAD1=F_DOWN_CY(IP,KP)
    WRITE(3,231) FAD1
    OWN=CELL_ID_CY(IP,JP,KP)
    WRITE(4,231) OWN
    ENDDO   
    ENDDO
        
        !RIGHT FACE--PROCESSOR_MI1_TO_MI2  !CASE CYLINDER 2-3
    IP=1
    I=IP+CNI*(MI-1)
    DO KP=1,NK
    K=KP
    DO JP=1,CNJ
    J=JP+CNJ*(MJ-1)
    F11=(CNJ+1)*(CNI+1)*(KP-1)+(CNI+1)*(JP-1)+IP-1
    F12=F11+1
    F13=F12+CNI+1
    F14=F13-1
    F21=F11+(CNJ+1)*(CNI+1)
    IF (KP.EQ.NK) F21=(CNI+1)*(JP-1)+IP-1
    F22=F21+1
    F23=F22+CNI+1
    F24=F23-1
    CALL CAL_FIJ_CY(F11,F12,F13,F14,F21,F22,F23,F24)
    WRITE(2,235) '4(',F11,' ',F21,' ',F24,' ',F14,')'
    FAD1=F_RIGHT_CY(JP,KP)
    WRITE(3,231) FAD1
    OWN=CELL_ID_CY(IP,JP,KP)
    WRITE(4,231) OWN
    ENDDO   
    ENDDO
        
        !LEFT FACE--PROCESSOR_MI1_TO_MI2  !CASE CYLINDER 2-3
    IP=CNI
    I=IP+CNI*(MI-1)
    DO KP=1,NK
    K=KP
    DO JP=1,CNJ
    J=JP+CNJ*(MJ-1)
    F11=(CNJ+1)*(CNI+1)*(KP-1)+(CNI+1)*(JP-1)+IP-1
    F12=F11+1
    F13=F12+CNI+1
    F14=F13-1
    F21=F11+(CNJ+1)*(CNI+1)
    IF (KP.EQ.NK) F21=(CNI+1)*(JP-1)+IP-1
    F22=F21+1
    F23=F22+CNI+1
    F24=F23-1
    CALL CAL_FIJ_CY(F11,F12,F13,F14,F21,F22,F23,F24)
    WRITE(2,235) '4(',F12,' ',F13,' ',F23,' ',F22,')'
    FAD1=F_LEFT_CY(JP,KP)
    WRITE(3,231) FAD1
    OWN=CELL_ID_CY(IP,JP,KP)
    WRITE(4,231) OWN
    ENDDO   
    ENDDO
    
    N=3
    WRITE(31,231) N+4
    WRITE(31,*) '('  
    WRITE(31,*) '0'
    WRITE(31,*) '1'
    WRITE(31,*) '2' 
    WRITE(31,*) '3' 
    DO I=1,N
    WRITE(31,*) '-1' 
    ENDDO
    WRITE(31,*) ')'
    
    WRITE(21,231) N+4
    WRITE(21,*) "("    
    WRITE(21,*) "    INLET1"
    WRITE(21,*) "    {"    
    WRITE(21,*) "        type            patch;"
    N=0
    WRITE(21,232) "        nFaces          ",N,";"
    N=F_INTERNAL_NUM
    WRITE(21,232) "        startFace       ",N,";"
    WRITE(21,*) "    }"
    WRITE(21,*) "    INLET2"
    WRITE(21,*) "    {"    
    WRITE(21,*) "        type            patch;"
    N=0
    WRITE(21,232) "        nFaces          ",N,";"
    N=F_INTERNAL_NUM
    WRITE(21,232) "        startFace       ",N,";"
    WRITE(21,*) "    }"
    WRITE(21,*) "    OUTLET"
    WRITE(21,*) "    {"
    WRITE(21,*) "        type            patch;"
    N=CNI*NK
    WRITE(21,232) "        nFaces          ",N,";"
    N=F_INTERNAL_NUM
    WRITE(21,232) "        startFace       ",N,";"
    WRITE(21,*) "    }"
    WRITE(21,*) "    HULL"
    WRITE(21,*) "    {"
    WRITE(21,*) "        type            wall;"
    WRITE(21,*) "        inGroups        1(wall);"
    N=0
    WRITE(21,232) "        nFaces          ",N,";"
    N=F_INTERNAL_NUM+CNI*NK
    WRITE(21,232) "        startFace       ",N,";"
    WRITE(21,*) "    }"
    
        !DOWN FACE--PROCESSOR_MI1_TO_MI2 !CASE CYLINDER 2-3
    PRO1=PNI*(MJ-1)+MI-1
    PRO2=PNI*(MJ-2)+MI-1
    WRITE(21,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(21,*) '    {'
    WRITE(21,*) '        type            processor;'
    WRITE(21,*) '        inGroups        List<word> 1(processor);'
    N=CNI*NK
    WRITE(21,232) "        nFaces          ",N,";"
    N=F_INTERNAL_NUM+CNI*NK
    WRITE(21,232) "        startFace       ",N,";"
    WRITE(21,*) '        matchTolerance  0.0001;'
    WRITE(21,232) "        myProcNo        ",PRO1,";"
    WRITE(21,232) "        neighbProcNo    ",PRO2,";"
    WRITE(21,*) '    }'
    WRITE(41,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(41,*) '    {'
    WRITE(41,*) '        type            processor;'
    WRITE(41,*) '        value           nonuniform List<scalar> '
    N=CNI*NK
    WRITE(41,'(I0)') N
    WRITE(41,*) '('
    DO KP=1,NK
    DO IP=1,CNI
    WRITE(41,'(E21.14)') P_DOWN_CY(IP,KP)    
    ENDDO
    ENDDO
    WRITE(41,*) ')'
    WRITE(41,*) ';'
    WRITE(41,*) '    }'
    WRITE(51,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(51,*) '    {'
    WRITE(51,*) '        type            processor;'
    WRITE(51,*) '        value           nonuniform List<vector> '
    N=CNI*NK
    WRITE(51,'(I0)') N
    WRITE(51,*) '('
    DO KP=1,NK
    DO IP=1,CNI
    UU1=U_DOWN_CY(1,IP,KP)
    UU2=U_DOWN_CY(2,IP,KP)
    UU3=U_DOWN_CY(3,IP,KP)
    WRITE(51,230) '( ',UU1,' ',UU2,' ',UU3,' )'  
    ENDDO
    ENDDO
    WRITE(51,*) ')'
    WRITE(51,*) ';'
    WRITE(51,*) '    }'
    WRITE(61,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(61,*) '    {'
    WRITE(61,*) '        type            processor;'
    WRITE(61,*) "        value           uniform 0;"
    WRITE(61,*) '    }'
    
        !RIGHT FACE--PROCESSOR_MI1_TO_MI2  !CASE CYLINDER 2-3
    PRO2=PNI*(MJ-1)+MI-2
    WRITE(21,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(21,*) '    {'
    WRITE(21,*) '        type            processor;'
    WRITE(21,*) '        inGroups        List<word> 1(processor);'
    N=CNJ*NK
    WRITE(21,232) "        nFaces          ",N,";"
    N=F_INTERNAL_NUM+CNI*NK+CNI*NK
    WRITE(21,232) "        startFace       ",N,";"
    WRITE(21,*) '        matchTolerance  0.0001;'
    WRITE(21,232) "        myProcNo        ",PRO1,";"
    WRITE(21,232) "        neighbProcNo    ",PRO2,";"
    WRITE(21,*) '    }'
    WRITE(41,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(41,*) '    {'
    WRITE(41,*) '        type            processor;'
    WRITE(41,*) '        value           nonuniform List<scalar> '
    N=CNJ*NK
    WRITE(41,'(I0)') N
    WRITE(41,*) '('
    DO KP=1,NK
    DO JP=1,CNJ
    WRITE(41,'(E21.14)') P_RIGHT_CY(JP,KP)    
    ENDDO
    ENDDO
    WRITE(41,*) ')'
    WRITE(41,*) ';'
    WRITE(41,*) '    }'
    WRITE(51,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(51,*) '    {'
    WRITE(51,*) '        type            processor;'
    WRITE(51,*) '        value           nonuniform List<vector> '
    N=CNJ*NK
    WRITE(51,'(I0)') N
    WRITE(51,*) '('
    DO KP=1,NK
    DO JP=1,CNJ
    UU1=U_RIGHT_CY(1,JP,KP)
    UU2=U_RIGHT_CY(2,JP,KP)
    UU3=U_RIGHT_CY(3,JP,KP)
    WRITE(51,230) '( ',UU1,' ',UU2,' ',UU3,' )'  
    ENDDO
    ENDDO
    WRITE(51,*) ')'
    WRITE(51,*) ';'
    WRITE(51,*) '    }'
    WRITE(61,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(61,*) '    {'
    WRITE(61,*) '        type            processor;'
    WRITE(61,*) "        value           uniform 0;"
    WRITE(61,*) '    }'
    
        !LEFT FACE--PROCESSOR_MI1_TO_MI2  !CASE CYLINDER 2-3
    PRO2=PNI*(MJ-1)+MI
    WRITE(21,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(21,*) '    {'
    WRITE(21,*) '        type            processor;'
    WRITE(21,*) '        inGroups        List<word> 1(processor);'
    N=CNJ*NK
    WRITE(21,232) "        nFaces          ",N,";"
    N=F_INTERNAL_NUM+CNI*NK+CNI*NK+CNJ*NK
    WRITE(21,232) "        startFace       ",N,";"
    WRITE(21,*) '        matchTolerance  0.0001;'
    WRITE(21,232) "        myProcNo        ",PRO1,";"
    WRITE(21,232) "        neighbProcNo    ",PRO2,";"
    WRITE(21,*) '    }'
    WRITE(41,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(41,*) '    {'
    WRITE(41,*) '        type            processor;'
    WRITE(41,*) '        value           nonuniform List<scalar> '
    N=CNJ*NK
    WRITE(41,'(I0)') N
    WRITE(41,*) '('
    DO KP=1,NK
    DO JP=1,CNJ
    WRITE(41,'(E21.14)') P_LEFT_CY(JP,KP)    
    ENDDO
    ENDDO
    WRITE(41,*) ')'
    WRITE(41,*) ';'
    WRITE(41,*) '    }'
    WRITE(51,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(51,*) '    {'
    WRITE(51,*) '        type            processor;'
    WRITE(51,*) '        value           nonuniform List<vector> '
    N=CNJ*NK
    WRITE(51,'(I0)') N
    WRITE(51,*) '('
    DO KP=1,NK
    DO JP=1,CNJ
    UU1=U_LEFT_CY(1,JP,KP)
    UU2=U_LEFT_CY(2,JP,KP)
    UU3=U_LEFT_CY(3,JP,KP)
    WRITE(51,230) '( ',UU1,' ',UU2,' ',UU3,' )'  
    ENDDO
    ENDDO
    WRITE(51,*) ')'
    WRITE(51,*) ';'
    WRITE(51,*) '    }'
    WRITE(61,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(61,*) '    {'
    WRITE(61,*) '        type            processor;'
    WRITE(61,*) "        value           uniform 0;"
    WRITE(61,*) '    }'
    
    
    
    
    
    
    
    
    ELSE
    MJ1=( (COORD(1)-PNJ)-MOD(COORD(1)-PNJ,PNL) )/PNL+1
    MJ2=MJ-PNJ-(MJ1-1)*PNL
    IF (MJ1.EQ.1) THEN
    IF (MJ2.EQ.1) THEN
    !CASE RECTANGLE 2-1-1
    
        !FRONT FACE--PROCESSOR_MI1_TO_MI2  !CASE RECTANGLE 2-1-1
    KP=1
    K=KP+CNL*(MJ2-1)
    DO JP=1,CNL
    J=JP+CNL*(MJ1-1)
    DO IP=1,CNI
    I=IP+CNI*(MI-1)
    F11=(CNL+1)*(CNI+1)*(KP-1)+(CNI+1)*(JP-1)+IP-1
    F12=F11+1
    F13=F12+CNI+1
    F14=F13-1
    F21=F11+(CNL+1)*(CNI+1)
    F22=F21+1
    F23=F22+CNI+1
    F24=F23-1
    CALL CAL_FIJ_RECT(F11,F12,F13,F14,F21,F22,F23,F24)
    WRITE(2,235) '4(',F11,' ',F14,' ',F13,' ',F12,')'
    FAD1=F_FRONT_RECT(IP,JP)
    WRITE(3,231) FAD1
    OWN=CELL_ID_RECT(IP,JP,KP)
    WRITE(4,231) OWN
    ENDDO   
    ENDDO
    
        !DOWN FACE--PROCESSOR_MI1_TO_MI2  !CASE RECTANGLE 2-1-1
    JP=1
    J=JP+CNL*(MJ1-1)
    DO KP=CNL,1,-1
    K=KP+CNL*(MJ2-1)
    DO IP=1,CNI
    I=IP+CNI*(MI-1)
    F11=(CNL+1)*(CNI+1)*(KP-1)+(CNI+1)*(JP-1)+IP-1
    F12=F11+1
    F13=F12+CNI+1
    F14=F13-1
    F21=F11+(CNL+1)*(CNI+1)
    F22=F21+1
    F23=F22+CNI+1
    F24=F23-1
    CALL CAL_FIJ_RECT(F11,F12,F13,F14,F21,F22,F23,F24)
    !WRITE(2,235) '4(',F11,' ',F12,' ',F22,' ',F21,')'
    WRITE(2,235) '4(',F21,' ',F11,' ',F12,' ',F22,')'
    FAD1=F_DOWN_RECT(IP,KP)
    WRITE(3,231) FAD1
    OWN=CELL_ID_RECT(IP,JP,KP)
    WRITE(4,231) OWN
    ENDDO   
    ENDDO
    
        !RIGHT FACE--PROCESSOR_MI1_TO_MI2  !CASE RECTANGLE 2-1-1 
    IP=1
    I=IP+CNI*(MI-1)
    DO KP=1,CNL
    K=KP+CNL*(MJ2-1)
    DO JP=1,CNL
    J=JP+CNL*(MJ1-1)
    F11=(CNL+1)*(CNI+1)*(KP-1)+(CNI+1)*(JP-1)+IP-1
    F12=F11+1
    F13=F12+CNI+1
    F14=F13-1
    F21=F11+(CNL+1)*(CNI+1)
    F22=F21+1
    F23=F22+CNI+1
    F24=F23-1
    CALL CAL_FIJ_RECT(F11,F12,F13,F14,F21,F22,F23,F24)
    WRITE(2,235) '4(',F11,' ',F21,' ',F24,' ',F14,')'
    FAD1=F_RIGHT_RECT(JP,KP)
    WRITE(3,231) FAD1
    OWN=CELL_ID_RECT(IP,JP,KP)
    WRITE(4,231) OWN
    ENDDO   
    ENDDO
    
        !LEFT FACE--PROCESSOR_MI1_TO_MI2  !CASE RECTANGLE 2-1-1 
    IP=CNI
    I=IP+CNI*(MI-1)
    DO KP=1,CNL
    K=KP+CNL*(MJ2-1)
    DO JP=1,CNL
    J=JP+CNL*(MJ1-1)
    F11=(CNL+1)*(CNI+1)*(KP-1)+(CNI+1)*(JP-1)+IP-1
    F12=F11+1
    F13=F12+CNI+1
    F14=F13-1
    F21=F11+(CNL+1)*(CNI+1)
    F22=F21+1
    F23=F22+CNI+1
    F24=F23-1
    CALL CAL_FIJ_RECT(F11,F12,F13,F14,F21,F22,F23,F24)
    WRITE(2,235) '4(',F12,' ',F13,' ',F23,' ',F22,')'
    FAD1=F_LEFT_RECT(JP,KP)
    WRITE(3,231) FAD1
    OWN=CELL_ID_RECT(IP,JP,KP)
    WRITE(4,231) OWN
    ENDDO   
    ENDDO
    
        !UP FACE--PROCESSOR_MI1_TO_MI2  !CASE RECTANGLE 2-1-1
    JP=CNL
    J=JP+CNL*(MJ1-1)
    DO KP=1,CNL
    K=KP+CNL*(MJ2-1)
    DO IP=1,CNI
    I=IP+CNI*(MI-1)
    F11=(CNL+1)*(CNI+1)*(KP-1)+(CNI+1)*(JP-1)+IP-1
    F12=F11+1
    F13=F12+CNI+1
    F14=F13-1
    F21=F11+(CNL+1)*(CNI+1)
    F22=F21+1
    F23=F22+CNI+1
    F24=F23-1
    CALL CAL_FIJ_RECT(F11,F12,F13,F14,F21,F22,F23,F24)
    WRITE(2,235) '4(',F14,' ',F24,' ',F23,' ',F13,')'
    FAD1=F_UP_RECT(IP,KP)
    WRITE(3,231) FAD1
    OWN=CELL_ID_RECT(IP,JP,KP)
    WRITE(4,231) OWN
    ENDDO   
    ENDDO
    
        !BACK FACE--PROCESSOR_MI1_TO_MI2  !CASE RECTANGLE 2-1-1
    KP=CNL
    K=KP+CNL*(MJ2-1)
    DO JP=1,CNL
    J=JP+CNL*(MJ1-1)
    DO IP=1,CNI
    I=IP+CNI*(MI-1)
    F11=(CNL+1)*(CNI+1)*(KP-1)+(CNI+1)*(JP-1)+IP-1
    F12=F11+1
    F13=F12+CNI+1
    F14=F13-1
    F21=F11+(CNL+1)*(CNI+1)
    F22=F21+1
    F23=F22+CNI+1
    F24=F23-1
    CALL CAL_FIJ_RECT(F11,F12,F13,F14,F21,F22,F23,F24)
    WRITE(2,235) '4(',F21,' ',F22,' ',F23,' ',F24,')'
    FAD1=F_BACK_RECT(IP,JP)
    WRITE(3,231) FAD1
    OWN=CELL_ID_RECT(IP,JP,KP)
    WRITE(4,231) OWN
    ENDDO   
    ENDDO    
        
    N=5
    WRITE(31,231) N+4
    WRITE(31,*) '('  
    WRITE(31,*) '0'
    WRITE(31,*) '1'
    WRITE(31,*) '2' 
    WRITE(31,*) '3' 
    DO I=1,N
    WRITE(31,*) '-1' 
    ENDDO
    WRITE(31,*) ')'
    
    WRITE(21,231) N+4
    WRITE(21,*) "("    
    WRITE(21,*) "    INLET1"
    WRITE(21,*) "    {"    
    WRITE(21,*) "        type            patch;"
    N=0
    WRITE(21,232) "        nFaces          ",N,";"
    N=F_INTERNAL_NUM
    WRITE(21,232) "        startFace       ",N,";"
    WRITE(21,*) "    }"
    WRITE(21,*) "    INLET2"
    WRITE(21,*) "    {"    
    WRITE(21,*) "        type            patch;"
    N=0
    WRITE(21,232) "        nFaces          ",N,";"
    N=F_INTERNAL_NUM
    WRITE(21,232) "        startFace       ",N,";"
    WRITE(21,*) "    }"
    WRITE(21,*) "    OUTLET"
    WRITE(21,*) "    {"
    WRITE(21,*) "        type            patch;"
    N=0
    WRITE(21,232) "        nFaces          ",N,";"
    N=F_INTERNAL_NUM
    WRITE(21,232) "        startFace       ",N,";"
    WRITE(21,*) "    }"
    WRITE(21,*) "    HULL"
    WRITE(21,*) "    {"
    WRITE(21,*) "        type            wall;"
    WRITE(21,*) "        inGroups        1(wall);"
    N=0
    WRITE(21,232) "        nFaces          ",N,";"
    N=F_INTERNAL_NUM
    WRITE(21,232) "        startFace       ",N,";"
    WRITE(21,*) "    }"    
        
        !FRONT/DOWN FACE--PROCESSOR_MI1_TO_MI2  !CASE RECTANGLE 2-1-1
    PRO1=PNI*PNJ+PNL*PNI*(MJ2-1)+PNI*(MJ1-1)+MI-1
    MJ=1
    PRO2=PNI*(MJ-1)+MI-1
    WRITE(21,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(21,*) '    {'
    WRITE(21,*) '        type            processor;'
    WRITE(21,*) '        inGroups        List<word> 1(processor);'
    N=CNI*CNL*2
    WRITE(21,232) "        nFaces          ",N,";"
    N=F_INTERNAL_NUM
    WRITE(21,232) "        startFace       ",N,";"
    WRITE(21,*) '        matchTolerance  0.0001;'
    WRITE(21,232) "        myProcNo        ",PRO1,";"
    WRITE(21,232) "        neighbProcNo    ",PRO2,";"
    WRITE(21,*) '    }'
    WRITE(41,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(41,*) '    {'
    WRITE(41,*) '        type            processor;'
    WRITE(41,*) '        value           nonuniform List<scalar> '
    N=CNI*CNL*2
    WRITE(41,'(I0)') N
    WRITE(41,*) '('
    DO JP=1,CNL
    DO IP=1,CNI
    WRITE(41,'(E21.14)') P_FRONT_RECT(IP,JP)    
    ENDDO
    ENDDO
    DO KP=1,CNL
    DO IP=1,CNI
    WRITE(41,'(E21.14)') P_DOWN_RECT(IP,KP)    
    ENDDO
    ENDDO
    WRITE(41,*) ')'
    WRITE(41,*) ';'
    WRITE(41,*) '    }'
    WRITE(51,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(51,*) '    {'
    WRITE(51,*) '        type            processor;'
    WRITE(51,*) '        value           nonuniform List<vector> '
    N=CNI*CNL*2
    WRITE(51,'(I0)') N
    WRITE(51,*) '('
    DO JP=1,CNL
    DO IP=1,CNI
    UU1=U_FRONT_RECT(1,IP,JP)
    UU2=U_FRONT_RECT(2,IP,JP)
    UU3=U_FRONT_RECT(3,IP,JP)
    WRITE(51,230) '( ',UU1,' ',UU2,' ',UU3,' )'  
    ENDDO
    ENDDO
    DO KP=1,CNL
    DO IP=1,CNI
    UU1=U_DOWN_RECT(1,IP,KP)
    UU2=U_DOWN_RECT(2,IP,KP)
    UU3=U_DOWN_RECT(3,IP,KP)
    WRITE(51,230) '( ',UU1,' ',UU2,' ',UU3,' )'  
    ENDDO
    ENDDO
    WRITE(51,*) ')'
    WRITE(51,*) ';'
    WRITE(51,*) '    }'
    WRITE(61,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(61,*) '    {'
    WRITE(61,*) '        type            processor;'
    WRITE(61,*) "        value           uniform 0;"
    WRITE(61,*) '    }'
    
        !RIGHT FACE--PROCESSOR_MI1_TO_MI2  !CASE RECTANGLE 2-1-1
    PRO2=PNI*PNJ+PNL*PNI*(MJ2-1)+PNI*(MJ1-1)+MI-2
    WRITE(21,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(21,*) '    {'
    WRITE(21,*) '        type            processor;'
    WRITE(21,*) '        inGroups        List<word> 1(processor);'
    N=CNL*CNL
    WRITE(21,232) "        nFaces          ",N,";"
    N=F_INTERNAL_NUM+CNI*CNL*2
    WRITE(21,232) "        startFace       ",N,";"
    WRITE(21,*) '        matchTolerance  0.0001;'
    WRITE(21,232) "        myProcNo        ",PRO1,";"
    WRITE(21,232) "        neighbProcNo    ",PRO2,";"
    WRITE(21,*) '    }'
    WRITE(41,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(41,*) '    {'
    WRITE(41,*) '        type            processor;'
    WRITE(41,*) '        value           nonuniform List<scalar> '
    N=CNL*CNL
    WRITE(41,'(I0)') N
    WRITE(41,*) '('
    DO KP=1,CNL
    DO JP=1,CNL
    WRITE(41,'(E21.14)') P_RIGHT_RECT(JP,KP)    
    ENDDO
    ENDDO
    WRITE(41,*) ')'
    WRITE(41,*) ';'
    WRITE(41,*) '    }'
    WRITE(51,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(51,*) '    {'
    WRITE(51,*) '        type            processor;'
    WRITE(51,*) '        value           nonuniform List<vector> '
    N=CNL*CNL
    WRITE(51,'(I0)') N
    WRITE(51,*) '('
    DO KP=1,CNL
    DO JP=1,CNL
    UU1=U_RIGHT_RECT(1,JP,KP)
    UU2=U_RIGHT_RECT(2,JP,KP)
    UU3=U_RIGHT_RECT(3,JP,KP)
    WRITE(51,230) '( ',UU1,' ',UU2,' ',UU3,' )'  
    ENDDO
    ENDDO
    WRITE(51,*) ')'
    WRITE(51,*) ';'
    WRITE(51,*) '    }'
    WRITE(61,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(61,*) '    {'
    WRITE(61,*) '        type            processor;'
    WRITE(61,*) "        value           uniform 0;"
    WRITE(61,*) '    }'
    
        !LEFT FACE--PROCESSOR_MI1_TO_MI2  !CASE RECTANGLE 2-1-1 
    PRO2=PNI*PNJ+PNL*PNI*(MJ2-1)+PNI*(MJ1-1)+MI
    WRITE(21,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(21,*) '    {'
    WRITE(21,*) '        type            processor;'
    WRITE(21,*) '        inGroups        List<word> 1(processor);'
    N=CNL*CNL
    WRITE(21,232) "        nFaces          ",N,";"
    N=F_INTERNAL_NUM+CNI*CNL*2+CNL*CNL
    WRITE(21,232) "        startFace       ",N,";"
    WRITE(21,*) '        matchTolerance  0.0001;'
    WRITE(21,232) "        myProcNo        ",PRO1,";"
    WRITE(21,232) "        neighbProcNo    ",PRO2,";"
    WRITE(21,*) '    }'
    WRITE(41,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(41,*) '    {'
    WRITE(41,*) '        type            processor;'
    WRITE(41,*) '        value           nonuniform List<scalar> '
    N=CNL*CNL
    WRITE(41,'(I0)') N
    WRITE(41,*) '('
    DO KP=1,CNL
    DO JP=1,CNL
    WRITE(41,'(E21.14)') P_LEFT_RECT(JP,KP)    
    ENDDO
    ENDDO
    WRITE(41,*) ')'
    WRITE(41,*) ';'
    WRITE(41,*) '    }'
    WRITE(51,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(51,*) '    {'
    WRITE(51,*) '        type            processor;'
    WRITE(51,*) '        value           nonuniform List<vector> '
    N=CNL*CNL
    WRITE(51,'(I0)') N
    WRITE(51,*) '('
    DO KP=1,CNL
    DO JP=1,CNL
    UU1=U_LEFT_RECT(1,JP,KP)
    UU2=U_LEFT_RECT(2,JP,KP)
    UU3=U_LEFT_RECT(3,JP,KP)
    WRITE(51,230) '( ',UU1,' ',UU2,' ',UU3,' )'  
    ENDDO
    ENDDO
    WRITE(51,*) ')'
    WRITE(51,*) ';'
    WRITE(51,*) '    }'
    WRITE(61,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(61,*) '    {'
    WRITE(61,*) '        type            processor;'
    WRITE(61,*) "        value           uniform 0;"
    WRITE(61,*) '    }'
    
        !UP FACE--PROCESSOR_MI1_TO_MI2  !CASE RECTANGLE 2-1-1
    PRO2=PNI*PNJ+PNL*PNI*(MJ2-1)+PNI*MJ1+MI-1
    WRITE(21,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(21,*) '    {'
    WRITE(21,*) '        type            processor;'
    WRITE(21,*) '        inGroups        List<word> 1(processor);'
    N=CNI*CNL
    WRITE(21,232) "        nFaces          ",N,";"
    N=F_INTERNAL_NUM+CNI*CNL*2+CNL*CNL+CNL*CNL
    WRITE(21,232) "        startFace       ",N,";"
    WRITE(21,*) '        matchTolerance  0.0001;'
    WRITE(21,232) "        myProcNo        ",PRO1,";"
    WRITE(21,232) "        neighbProcNo    ",PRO2,";"
    WRITE(21,*) '    }'
    WRITE(41,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(41,*) '    {'
    WRITE(41,*) '        type            processor;'
    WRITE(41,*) '        value           nonuniform List<scalar> '
    N=CNI*CNL
    WRITE(41,'(I0)') N
    WRITE(41,*) '('
    DO KP=1,CNL
    DO IP=1,CNI
    WRITE(41,'(E21.14)') P_UP_RECT(IP,KP)    
    ENDDO
    ENDDO
    WRITE(41,*) ')'
    WRITE(41,*) ';'
    WRITE(41,*) '    }'
    WRITE(51,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(51,*) '    {'
    WRITE(51,*) '        type            processor;'
    WRITE(51,*) '        value           nonuniform List<vector> '
    N=CNI*CNL
    WRITE(51,'(I0)') N
    WRITE(51,*) '('
    DO KP=1,CNL
    DO IP=1,CNI
    UU1=U_UP_RECT(1,IP,KP)
    UU2=U_UP_RECT(2,IP,KP)
    UU3=U_UP_RECT(3,IP,KP)
    WRITE(51,230) '( ',UU1,' ',UU2,' ',UU3,' )'  
    ENDDO
    ENDDO
    WRITE(51,*) ')'
    WRITE(51,*) ';'
    WRITE(51,*) '    }'
    WRITE(61,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(61,*) '    {'
    WRITE(61,*) '        type            processor;'
    WRITE(61,*) "        value           uniform 0;"
    WRITE(61,*) '    }'
    
        !BACK FACE--PROCESSOR_MI1_TO_MI2  !CASE RECTANGLE 2-1-1
    PRO2=PNI*PNJ+PNL*PNI*MJ2+PNI*(MJ1-1)+MI-1
    WRITE(21,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(21,*) '    {'
    WRITE(21,*) '        type            processor;'
    WRITE(21,*) '        inGroups        List<word> 1(processor);'
    N=CNI*CNL
    WRITE(21,232) "        nFaces          ",N,";"
    N=F_INTERNAL_NUM+CNI*CNL*2+CNL*CNL+CNL*CNL+CNI*CNL
    WRITE(21,232) "        startFace       ",N,";"
    WRITE(21,*) '        matchTolerance  0.0001;'
    WRITE(21,232) "        myProcNo        ",PRO1,";"
    WRITE(21,232) "        neighbProcNo    ",PRO2,";"
    WRITE(21,*) '    }'
    WRITE(41,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(41,*) '    {'
    WRITE(41,*) '        type            processor;'
    WRITE(41,*) '        value           nonuniform List<scalar> '
    N=CNI*CNL
    WRITE(41,'(I0)') N
    WRITE(41,*) '('
    DO JP=1,CNL
    DO IP=1,CNI
    WRITE(41,'(E21.14)') P_BACK_RECT(IP,JP)    
    ENDDO
    ENDDO
    WRITE(41,*) ')'
    WRITE(41,*) ';'
    WRITE(41,*) '    }'
    WRITE(51,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(51,*) '    {'
    WRITE(51,*) '        type            processor;'
    WRITE(51,*) '        value           nonuniform List<vector> '
    N=CNI*CNL
    WRITE(51,'(I0)') N
    WRITE(51,*) '('
    DO JP=1,CNL
    DO IP=1,CNI
    UU1=U_BACK_RECT(1,IP,JP)
    UU2=U_BACK_RECT(2,IP,JP)
    UU3=U_BACK_RECT(3,IP,JP)
    WRITE(51,230) '( ',UU1,' ',UU2,' ',UU3,' )'  
    ENDDO
    ENDDO
    WRITE(51,*) ')'
    WRITE(51,*) ';'
    WRITE(51,*) '    }'
    WRITE(61,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(61,*) '    {'
    WRITE(61,*) '        type            processor;'
    WRITE(61,*) "        value           uniform 0;"
    WRITE(61,*) '    }'    
        
        
        
        
        
        
        
        
        
        
    ELSEIF (MJ2.EQ.PNL) THEN
    !CASE RECTANGLE 2-1-3
    
        !BACK FACE--PROCESSOR_MI1_TO_MI2  !CASE RECTANGLE 2-1-3
    KP=CNL
    K=KP+CNL*(MJ2-1)
    DO JP=CNL,1,-1
    J=JP+CNL*(MJ1-1)
    DO IP=1,CNI
    I=IP+CNI*(MI-1)
    F11=(CNL+1)*(CNI+1)*(KP-1)+(CNI+1)*(JP-1)+IP-1
    F12=F11+1
    F13=F12+CNI+1
    F14=F13-1
    F21=F11+(CNL+1)*(CNI+1)
    F22=F21+1
    F23=F22+CNI+1
    F24=F23-1
    CALL CAL_FIJ_RECT(F11,F12,F13,F14,F21,F22,F23,F24)
    !WRITE(2,235) '4(',F21,' ',F22,' ',F23,' ',F24,')'
    WRITE(2,235) '4(',F24,' ',F21,' ',F22,' ',F23,')'
    FAD1=F_BACK_RECT(IP,JP)
    WRITE(3,231) FAD1
    OWN=CELL_ID_RECT(IP,JP,KP)
    WRITE(4,231) OWN
    ENDDO   
    ENDDO    
        
        !DOWN FACE--PROCESSOR_MI1_TO_MI2  !CASE RECTANGLE 2-1-3
    JP=1
    J=JP+CNL*(MJ1-1)
    DO KP=CNL,1,-1
    K=KP+CNL*(MJ2-1)
    DO IP=1,CNI
    I=IP+CNI*(MI-1)
    F11=(CNL+1)*(CNI+1)*(KP-1)+(CNI+1)*(JP-1)+IP-1
    F12=F11+1
    F13=F12+CNI+1
    F14=F13-1
    F21=F11+(CNL+1)*(CNI+1)
    F22=F21+1
    F23=F22+CNI+1
    F24=F23-1
    CALL CAL_FIJ_RECT(F11,F12,F13,F14,F21,F22,F23,F24)
    !WRITE(2,235) '4(',F11,' ',F12,' ',F22,' ',F21,')'
    WRITE(2,235) '4(',F21,' ',F11,' ',F12,' ',F22,')'
    FAD1=F_DOWN_RECT(IP,KP)
    WRITE(3,231) FAD1
    OWN=CELL_ID_RECT(IP,JP,KP)
    WRITE(4,231) OWN
    ENDDO   
    ENDDO    
        
        !FRONT FACE--PROCESSOR_MI1_TO_MI2  !CASE RECTANGLE 2-1-3
    KP=1
    K=KP+CNL*(MJ2-1)
    DO JP=1,CNL
    J=JP+CNL*(MJ1-1)
    DO IP=1,CNI
    I=IP+CNI*(MI-1)
    F11=(CNL+1)*(CNI+1)*(KP-1)+(CNI+1)*(JP-1)+IP-1
    F12=F11+1
    F13=F12+CNI+1
    F14=F13-1
    F21=F11+(CNL+1)*(CNI+1)
    F22=F21+1
    F23=F22+CNI+1
    F24=F23-1
    CALL CAL_FIJ_RECT(F11,F12,F13,F14,F21,F22,F23,F24)
    WRITE(2,235) '4(',F11,' ',F14,' ',F13,' ',F12,')'
    FAD1=F_FRONT_RECT(IP,JP)
    WRITE(3,231) FAD1
    OWN=CELL_ID_RECT(IP,JP,KP)
    WRITE(4,231) OWN
    ENDDO   
    ENDDO
    
        !RIGHT FACE--PROCESSOR_MI1_TO_MI2  !CASE RECTANGLE 2-1-3 
    IP=1
    I=IP+CNI*(MI-1)
    DO KP=1,CNL
    K=KP+CNL*(MJ2-1)
    DO JP=1,CNL
    J=JP+CNL*(MJ1-1)
    F11=(CNL+1)*(CNI+1)*(KP-1)+(CNI+1)*(JP-1)+IP-1
    F12=F11+1
    F13=F12+CNI+1
    F14=F13-1
    F21=F11+(CNL+1)*(CNI+1)
    F22=F21+1
    F23=F22+CNI+1
    F24=F23-1
    CALL CAL_FIJ_RECT(F11,F12,F13,F14,F21,F22,F23,F24)
    WRITE(2,235) '4(',F11,' ',F21,' ',F24,' ',F14,')'
    FAD1=F_RIGHT_RECT(JP,KP)
    WRITE(3,231) FAD1
    OWN=CELL_ID_RECT(IP,JP,KP)
    WRITE(4,231) OWN
    ENDDO   
    ENDDO
    
        !LEFT FACE--PROCESSOR_MI1_TO_MI2  !CASE RECTANGLE 2-1-3 
    IP=CNI
    I=IP+CNI*(MI-1)
    DO KP=1,CNL
    K=KP+CNL*(MJ2-1)
    DO JP=1,CNL
    J=JP+CNL*(MJ1-1)
    F11=(CNL+1)*(CNI+1)*(KP-1)+(CNI+1)*(JP-1)+IP-1
    F12=F11+1
    F13=F12+CNI+1
    F14=F13-1
    F21=F11+(CNL+1)*(CNI+1)
    F22=F21+1
    F23=F22+CNI+1
    F24=F23-1
    CALL CAL_FIJ_RECT(F11,F12,F13,F14,F21,F22,F23,F24)
    WRITE(2,235) '4(',F12,' ',F13,' ',F23,' ',F22,')'
    FAD1=F_LEFT_RECT(JP,KP)
    WRITE(3,231) FAD1
    OWN=CELL_ID_RECT(IP,JP,KP)
    WRITE(4,231) OWN
    ENDDO   
    ENDDO
    
        !UP FACE--PROCESSOR_MI1_TO_MI2  !CASE RECTANGLE 2-1-3
    JP=CNL
    J=JP+CNL*(MJ1-1)
    DO KP=1,CNL
    K=KP+CNL*(MJ2-1)
    DO IP=1,CNI
    I=IP+CNI*(MI-1)
    F11=(CNL+1)*(CNI+1)*(KP-1)+(CNI+1)*(JP-1)+IP-1
    F12=F11+1
    F13=F12+CNI+1
    F14=F13-1
    F21=F11+(CNL+1)*(CNI+1)
    F22=F21+1
    F23=F22+CNI+1
    F24=F23-1
    CALL CAL_FIJ_RECT(F11,F12,F13,F14,F21,F22,F23,F24)
    WRITE(2,235) '4(',F14,' ',F24,' ',F23,' ',F13,')'
    FAD1=F_UP_RECT(IP,KP)
    WRITE(3,231) FAD1
    OWN=CELL_ID_RECT(IP,JP,KP)
    WRITE(4,231) OWN
    ENDDO   
    ENDDO    
        
    N=5
    WRITE(31,231) N+4
    WRITE(31,*) '('  
    WRITE(31,*) '0'
    WRITE(31,*) '1'
    WRITE(31,*) '2' 
    WRITE(31,*) '3' 
    DO I=1,N
    WRITE(31,*) '-1' 
    ENDDO
    WRITE(31,*) ')'
    
    WRITE(21,231) N+4
    WRITE(21,*) "("    
    WRITE(21,*) "    INLET1"
    WRITE(21,*) "    {"    
    WRITE(21,*) "        type            patch;"
    N=0
    WRITE(21,232) "        nFaces          ",N,";"
    N=F_INTERNAL_NUM
    WRITE(21,232) "        startFace       ",N,";"
    WRITE(21,*) "    }"
    WRITE(21,*) "    INLET2"
    WRITE(21,*) "    {"    
    WRITE(21,*) "        type            patch;"
    N=0
    WRITE(21,232) "        nFaces          ",N,";"
    N=F_INTERNAL_NUM
    WRITE(21,232) "        startFace       ",N,";"
    WRITE(21,*) "    }"
    WRITE(21,*) "    OUTLET"
    WRITE(21,*) "    {"
    WRITE(21,*) "        type            patch;"
    N=0
    WRITE(21,232) "        nFaces          ",N,";"
    N=F_INTERNAL_NUM
    WRITE(21,232) "        startFace       ",N,";"
    WRITE(21,*) "    }"
    WRITE(21,*) "    HULL"
    WRITE(21,*) "    {"
    WRITE(21,*) "        type            wall;"
    WRITE(21,*) "        inGroups        1(wall);"
    N=0
    WRITE(21,232) "        nFaces          ",N,";"
    N=F_INTERNAL_NUM
    WRITE(21,232) "        startFace       ",N,";"
    WRITE(21,*) "    }"    
        
        !BACK/DOWN FACE--PROCESSOR_MI1_TO_MI2  !CASE RECTANGLE 2-1-3
    PRO1=PNI*PNJ+PNL*PNI*(MJ2-1)+PNI*(MJ1-1)+MI-1
    MJ=1
    PRO2=PNI*(MJ-1)+MI-1
    WRITE(21,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(21,*) '    {'
    WRITE(21,*) '        type            processor;'
    WRITE(21,*) '        inGroups        List<word> 1(processor);'
    N=CNI*CNL*2
    WRITE(21,232) "        nFaces          ",N,";"
    N=F_INTERNAL_NUM
    WRITE(21,232) "        startFace       ",N,";"
    WRITE(21,*) '        matchTolerance  0.0001;'
    WRITE(21,232) "        myProcNo        ",PRO1,";"
    WRITE(21,232) "        neighbProcNo    ",PRO2,";"
    WRITE(21,*) '    }'
    WRITE(41,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(41,*) '    {'
    WRITE(41,*) '        type            processor;'
    WRITE(41,*) '        value           nonuniform List<scalar> '
    N=CNI*CNL*2
    WRITE(41,'(I0)') N
    WRITE(41,*) '('
    DO JP=1,CNL
    DO IP=1,CNI
    WRITE(41,'(E21.14)') P_BACK_RECT(IP,JP)    
    ENDDO
    ENDDO
    DO KP=1,CNL
    DO IP=1,CNI
    WRITE(41,'(E21.14)') P_DOWN_RECT(IP,KP)    
    ENDDO
    ENDDO
    WRITE(41,*) ')'
    WRITE(41,*) ';'
    WRITE(41,*) '    }'
    WRITE(51,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(51,*) '    {'
    WRITE(51,*) '        type            processor;'
    WRITE(51,*) '        value           nonuniform List<vector> '
    N=CNI*CNL*2
    WRITE(51,'(I0)') N
    WRITE(51,*) '('
    DO JP=1,CNL
    DO IP=1,CNI
    UU1=U_BACK_RECT(1,IP,JP)
    UU2=U_BACK_RECT(2,IP,JP)
    UU3=U_BACK_RECT(3,IP,JP)
    WRITE(51,230) '( ',UU1,' ',UU2,' ',UU3,' )'  
    ENDDO
    ENDDO
    DO KP=1,CNL
    DO IP=1,CNI
    UU1=U_DOWN_RECT(1,IP,KP)
    UU2=U_DOWN_RECT(2,IP,KP)
    UU3=U_DOWN_RECT(3,IP,KP)
    WRITE(51,230) '( ',UU1,' ',UU2,' ',UU3,' )'  
    ENDDO
    ENDDO
    WRITE(51,*) ')'
    WRITE(51,*) ';'
    WRITE(51,*) '    }'
    WRITE(61,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(61,*) '    {'
    WRITE(61,*) '        type            processor;'
    WRITE(61,*) "        value           uniform 0;"
    WRITE(61,*) '    }'
    
        !FRONT FACE--PROCESSOR_MI1_TO_MI2  !CASE RECTANGLE 2-1-3
    PRO2=PNI*PNJ+PNL*PNI*(MJ2-2)+PNI*(MJ1-1)+MI-1
    WRITE(21,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(21,*) '    {'
    WRITE(21,*) '        type            processor;'
    WRITE(21,*) '        inGroups        List<word> 1(processor);'
    N=CNI*CNL
    WRITE(21,232) "        nFaces          ",N,";"
    N=F_INTERNAL_NUM+CNI*CNL*2
    WRITE(21,232) "        startFace       ",N,";"
    WRITE(21,*) '        matchTolerance  0.0001;'
    WRITE(21,232) "        myProcNo        ",PRO1,";"
    WRITE(21,232) "        neighbProcNo    ",PRO2,";"
    WRITE(21,*) '    }'
    WRITE(41,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(41,*) '    {'
    WRITE(41,*) '        type            processor;'
    WRITE(41,*) '        value           nonuniform List<scalar> '
    N=CNI*CNL
    WRITE(41,'(I0)') N
    WRITE(41,*) '('
    DO JP=1,CNL
    DO IP=1,CNI
    WRITE(41,'(E21.14)') P_FRONT_RECT(IP,JP)    
    ENDDO
    ENDDO
    WRITE(41,*) ')'
    WRITE(41,*) ';'
    WRITE(41,*) '    }'
    WRITE(51,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(51,*) '    {'
    WRITE(51,*) '        type            processor;'
    WRITE(51,*) '        value           nonuniform List<vector> '
    N=CNI*CNL
    WRITE(51,'(I0)') N
    WRITE(51,*) '('
    DO JP=1,CNL
    DO IP=1,CNI
    UU1=U_FRONT_RECT(1,IP,JP)
    UU2=U_FRONT_RECT(2,IP,JP)
    UU3=U_FRONT_RECT(3,IP,JP)
    WRITE(51,230) '( ',UU1,' ',UU2,' ',UU3,' )'  
    ENDDO
    ENDDO
    WRITE(51,*) ')'
    WRITE(51,*) ';'
    WRITE(51,*) '    }'
    WRITE(61,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(61,*) '    {'
    WRITE(61,*) '        type            processor;'
    WRITE(61,*) "        value           uniform 0;"
    WRITE(61,*) '    }'
    
        !RIGHT FACE--PROCESSOR_MI1_TO_MI2  !CASE RECTANGLE 2-1-3
    PRO2=PNI*PNJ+PNL*PNI*(MJ2-1)+PNI*(MJ1-1)+MI-2
    WRITE(21,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(21,*) '    {'
    WRITE(21,*) '        type            processor;'
    WRITE(21,*) '        inGroups        List<word> 1(processor);'
    N=CNL*CNL
    WRITE(21,232) "        nFaces          ",N,";"
    N=F_INTERNAL_NUM+CNI*CNL*2+CNI*CNL
    WRITE(21,232) "        startFace       ",N,";"
    WRITE(21,*) '        matchTolerance  0.0001;'
    WRITE(21,232) "        myProcNo        ",PRO1,";"
    WRITE(21,232) "        neighbProcNo    ",PRO2,";"
    WRITE(21,*) '    }'
    WRITE(41,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(41,*) '    {'
    WRITE(41,*) '        type            processor;'
    WRITE(41,*) '        value           nonuniform List<scalar> '
    N=CNL*CNL
    WRITE(41,'(I0)') N
    WRITE(41,*) '('
    DO KP=1,CNL
    DO JP=1,CNL
    WRITE(41,'(E21.14)') P_RIGHT_RECT(JP,KP)    
    ENDDO
    ENDDO
    WRITE(41,*) ')'
    WRITE(41,*) ';'
    WRITE(41,*) '    }'
    WRITE(51,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(51,*) '    {'
    WRITE(51,*) '        type            processor;'
    WRITE(51,*) '        value           nonuniform List<vector> '
    N=CNL*CNL
    WRITE(51,'(I0)') N
    WRITE(51,*) '('
    DO KP=1,CNL
    DO JP=1,CNL
    UU1=U_RIGHT_RECT(1,JP,KP)
    UU2=U_RIGHT_RECT(2,JP,KP)
    UU3=U_RIGHT_RECT(3,JP,KP)
    WRITE(51,230) '( ',UU1,' ',UU2,' ',UU3,' )'  
    ENDDO
    ENDDO
    WRITE(51,*) ')'
    WRITE(51,*) ';'
    WRITE(51,*) '    }'
    WRITE(61,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(61,*) '    {'
    WRITE(61,*) '        type            processor;'
    WRITE(61,*) "        value           uniform 0;"
    WRITE(61,*) '    }' 
    
        !LEFT FACE--PROCESSOR_MI1_TO_MI2  !CASE RECTANGLE 2-1-3
    PRO2=PNI*PNJ+PNL*PNI*(MJ2-1)+PNI*(MJ1-1)+MI
    WRITE(21,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(21,*) '    {'
    WRITE(21,*) '        type            processor;'
    WRITE(21,*) '        inGroups        List<word> 1(processor);'
    N=CNL*CNL
    WRITE(21,232) "        nFaces          ",N,";"
    N=F_INTERNAL_NUM+CNI*CNL*2+CNI*CNL+CNL*CNL
    WRITE(21,232) "        startFace       ",N,";"
    WRITE(21,*) '        matchTolerance  0.0001;'
    WRITE(21,232) "        myProcNo        ",PRO1,";"
    WRITE(21,232) "        neighbProcNo    ",PRO2,";"
    WRITE(21,*) '    }'
    WRITE(41,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(41,*) '    {'
    WRITE(41,*) '        type            processor;'
    WRITE(41,*) '        value           nonuniform List<scalar> '
    N=CNL*CNL
    WRITE(41,'(I0)') N
    WRITE(41,*) '('
    DO KP=1,CNL
    DO JP=1,CNL
    WRITE(41,'(E21.14)') P_LEFT_RECT(JP,KP)    
    ENDDO
    ENDDO
    WRITE(41,*) ')'
    WRITE(41,*) ';'
    WRITE(41,*) '    }'
    WRITE(51,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(51,*) '    {'
    WRITE(51,*) '        type            processor;'
    WRITE(51,*) '        value           nonuniform List<vector> '
    N=CNL*CNL
    WRITE(51,'(I0)') N
    WRITE(51,*) '('
    DO KP=1,CNL
    DO JP=1,CNL
    UU1=U_LEFT_RECT(1,JP,KP)
    UU2=U_LEFT_RECT(2,JP,KP)
    UU3=U_LEFT_RECT(3,JP,KP)
    WRITE(51,230) '( ',UU1,' ',UU2,' ',UU3,' )'  
    ENDDO
    ENDDO
    WRITE(51,*) ')'
    WRITE(51,*) ';'
    WRITE(51,*) '    }'
    WRITE(61,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(61,*) '    {'
    WRITE(61,*) '        type            processor;'
    WRITE(61,*) "        value           uniform 0;"
    WRITE(61,*) '    }'
    
        !UP FACE--PROCESSOR_MI1_TO_MI2  !CASE RECTANGLE 2-1-3
    PRO2=PNI*PNJ+PNL*PNI*(MJ2-1)+PNI*MJ1+MI-1
    WRITE(21,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(21,*) '    {'
    WRITE(21,*) '        type            processor;'
    WRITE(21,*) '        inGroups        List<word> 1(processor);'
    N=CNI*CNL
    WRITE(21,232) "        nFaces          ",N,";"
    N=F_INTERNAL_NUM+CNI*CNL*2+CNI*CNL+CNL*CNL+CNL*CNL
    WRITE(21,232) "        startFace       ",N,";"
    WRITE(21,*) '        matchTolerance  0.0001;'
    WRITE(21,232) "        myProcNo        ",PRO1,";"
    WRITE(21,232) "        neighbProcNo    ",PRO2,";"
    WRITE(21,*) '    }'
    WRITE(41,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(41,*) '    {'
    WRITE(41,*) '        type            processor;'
    WRITE(41,*) '        value           nonuniform List<scalar> '
    N=CNI*CNL
    WRITE(41,'(I0)') N
    WRITE(41,*) '('
    DO KP=1,CNL
    DO IP=1,CNI
    WRITE(41,'(E21.14)') P_UP_RECT(IP,KP)    
    ENDDO
    ENDDO
    WRITE(41,*) ')'
    WRITE(41,*) ';'
    WRITE(41,*) '    }'
    WRITE(51,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(51,*) '    {'
    WRITE(51,*) '        type            processor;'
    WRITE(51,*) '        value           nonuniform List<vector> '
    N=CNI*CNL
    WRITE(51,'(I0)') N
    WRITE(51,*) '('
    DO KP=1,CNL
    DO IP=1,CNI
    UU1=U_UP_RECT(1,IP,KP)
    UU2=U_UP_RECT(2,IP,KP)
    UU3=U_UP_RECT(3,IP,KP)
    WRITE(51,230) '( ',UU1,' ',UU2,' ',UU3,' )'  
    ENDDO
    ENDDO
    WRITE(51,*) ')'
    WRITE(51,*) ';'
    WRITE(51,*) '    }'
    WRITE(61,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(61,*) '    {'
    WRITE(61,*) '        type            processor;'
    WRITE(61,*) "        value           uniform 0;"
    WRITE(61,*) '    }'    
        
        
        
        
        
        
    ELSE
    !CASE RECTANGLE 2-1-2
    
        !DOWN FACE--PROCESSOR_MI1_TO_MI2  !CASE RECTANGLE 2-1-2
    JP=1
    J=JP+CNL*(MJ1-1)
    DO KP=CNL,1,-1
    K=KP+CNL*(MJ2-1)
    DO IP=1,CNI
    I=IP+CNI*(MI-1)
    F11=(CNL+1)*(CNI+1)*(KP-1)+(CNI+1)*(JP-1)+IP-1
    F12=F11+1
    F13=F12+CNI+1
    F14=F13-1
    F21=F11+(CNL+1)*(CNI+1)
    F22=F21+1
    F23=F22+CNI+1
    F24=F23-1
    CALL CAL_FIJ_RECT(F11,F12,F13,F14,F21,F22,F23,F24)
    !WRITE(2,235) '4(',F11,' ',F12,' ',F22,' ',F21,')'
    WRITE(2,235) '4(',F21,' ',F11,' ',F12,' ',F22,')'
    FAD1=F_DOWN_RECT(IP,KP)
    WRITE(3,231) FAD1
    OWN=CELL_ID_RECT(IP,JP,KP)
    WRITE(4,231) OWN
    ENDDO   
    ENDDO    
        
        !FRONT FACE--PROCESSOR_MI1_TO_MI2  !CASE RECTANGLE 2-1-2
    KP=1
    K=KP+CNL*(MJ2-1)
    DO JP=1,CNL
    J=JP+CNL*(MJ1-1)
    DO IP=1,CNI
    I=IP+CNI*(MI-1)
    F11=(CNL+1)*(CNI+1)*(KP-1)+(CNI+1)*(JP-1)+IP-1
    F12=F11+1
    F13=F12+CNI+1
    F14=F13-1
    F21=F11+(CNL+1)*(CNI+1)
    F22=F21+1
    F23=F22+CNI+1
    F24=F23-1
    CALL CAL_FIJ_RECT(F11,F12,F13,F14,F21,F22,F23,F24)
    WRITE(2,235) '4(',F11,' ',F14,' ',F13,' ',F12,')'
    FAD1=F_FRONT_RECT(IP,JP)
    WRITE(3,231) FAD1
    OWN=CELL_ID_RECT(IP,JP,KP)
    WRITE(4,231) OWN
    ENDDO   
    ENDDO
    
        !RIGHT FACE--PROCESSOR_MI1_TO_MI2  !CASE RECTANGLE 2-1-2 
    IP=1
    I=IP+CNI*(MI-1)
    DO KP=1,CNL
    K=KP+CNL*(MJ2-1)
    DO JP=1,CNL
    J=JP+CNL*(MJ1-1)
    F11=(CNL+1)*(CNI+1)*(KP-1)+(CNI+1)*(JP-1)+IP-1
    F12=F11+1
    F13=F12+CNI+1
    F14=F13-1
    F21=F11+(CNL+1)*(CNI+1)
    F22=F21+1
    F23=F22+CNI+1
    F24=F23-1
    CALL CAL_FIJ_RECT(F11,F12,F13,F14,F21,F22,F23,F24)
    WRITE(2,235) '4(',F11,' ',F21,' ',F24,' ',F14,')'
    FAD1=F_RIGHT_RECT(JP,KP)
    WRITE(3,231) FAD1
    OWN=CELL_ID_RECT(IP,JP,KP)
    WRITE(4,231) OWN
    ENDDO   
    ENDDO
    
        !LEFT FACE--PROCESSOR_MI1_TO_MI2  !CASE RECTANGLE 2-1-2 
    IP=CNI
    I=IP+CNI*(MI-1)
    DO KP=1,CNL
    K=KP+CNL*(MJ2-1)
    DO JP=1,CNL
    J=JP+CNL*(MJ1-1)
    F11=(CNL+1)*(CNI+1)*(KP-1)+(CNI+1)*(JP-1)+IP-1
    F12=F11+1
    F13=F12+CNI+1
    F14=F13-1
    F21=F11+(CNL+1)*(CNI+1)
    F22=F21+1
    F23=F22+CNI+1
    F24=F23-1
    CALL CAL_FIJ_RECT(F11,F12,F13,F14,F21,F22,F23,F24)
    WRITE(2,235) '4(',F12,' ',F13,' ',F23,' ',F22,')'
    FAD1=F_LEFT_RECT(JP,KP)
    WRITE(3,231) FAD1
    OWN=CELL_ID_RECT(IP,JP,KP)
    WRITE(4,231) OWN
    ENDDO   
    ENDDO
    
        !UP FACE--PROCESSOR_MI1_TO_MI2  !CASE RECTANGLE 2-1-2
    JP=CNL
    J=JP+CNL*(MJ1-1)
    DO KP=1,CNL
    K=KP+CNL*(MJ2-1)
    DO IP=1,CNI
    I=IP+CNI*(MI-1)
    F11=(CNL+1)*(CNI+1)*(KP-1)+(CNI+1)*(JP-1)+IP-1
    F12=F11+1
    F13=F12+CNI+1
    F14=F13-1
    F21=F11+(CNL+1)*(CNI+1)
    F22=F21+1
    F23=F22+CNI+1
    F24=F23-1
    CALL CAL_FIJ_RECT(F11,F12,F13,F14,F21,F22,F23,F24)
    WRITE(2,235) '4(',F14,' ',F24,' ',F23,' ',F13,')'
    FAD1=F_UP_RECT(IP,KP)
    WRITE(3,231) FAD1
    OWN=CELL_ID_RECT(IP,JP,KP)
    WRITE(4,231) OWN
    ENDDO   
    ENDDO
    
        !BACK FACE--PROCESSOR_MI1_TO_MI2  !CASE RECTANGLE 2-1-2
    KP=CNL
    K=KP+CNL*(MJ2-1)
    DO JP=1,CNL
    J=JP+CNL*(MJ1-1)
    DO IP=1,CNI
    I=IP+CNI*(MI-1)
    F11=(CNL+1)*(CNI+1)*(KP-1)+(CNI+1)*(JP-1)+IP-1
    F12=F11+1
    F13=F12+CNI+1
    F14=F13-1
    F21=F11+(CNL+1)*(CNI+1)
    F22=F21+1
    F23=F22+CNI+1
    F24=F23-1
    CALL CAL_FIJ_RECT(F11,F12,F13,F14,F21,F22,F23,F24)
    WRITE(2,235) '4(',F21,' ',F22,' ',F23,' ',F24,')'
    FAD1=F_BACK_RECT(IP,JP)
    WRITE(3,231) FAD1
    OWN=CELL_ID_RECT(IP,JP,KP)
    WRITE(4,231) OWN
    ENDDO   
    ENDDO    
        
    N=6
    WRITE(31,231) N+4
    WRITE(31,*) '('  
    WRITE(31,*) '0'
    WRITE(31,*) '1'
    WRITE(31,*) '2' 
    WRITE(31,*) '3' 
    DO I=1,N
    WRITE(31,*) '-1' 
    ENDDO
    WRITE(31,*) ')'
    
    WRITE(21,231) N+4
    WRITE(21,*) "("    
    WRITE(21,*) "    INLET1"
    WRITE(21,*) "    {"    
    WRITE(21,*) "        type            patch;"
    N=0
    WRITE(21,232) "        nFaces          ",N,";"
    N=F_INTERNAL_NUM
    WRITE(21,232) "        startFace       ",N,";"
    WRITE(21,*) "    }"
    WRITE(21,*) "    INLET2"
    WRITE(21,*) "    {"    
    WRITE(21,*) "        type            patch;"
    N=0
    WRITE(21,232) "        nFaces          ",N,";"
    N=F_INTERNAL_NUM
    WRITE(21,232) "        startFace       ",N,";"
    WRITE(21,*) "    }"
    WRITE(21,*) "    OUTLET"
    WRITE(21,*) "    {"
    WRITE(21,*) "        type            patch;"
    N=0
    WRITE(21,232) "        nFaces          ",N,";"
    N=F_INTERNAL_NUM
    WRITE(21,232) "        startFace       ",N,";"
    WRITE(21,*) "    }"
    WRITE(21,*) "    HULL"
    WRITE(21,*) "    {"
    WRITE(21,*) "        type            wall;"
    WRITE(21,*) "        inGroups        1(wall);"
    N=0
    WRITE(21,232) "        nFaces          ",N,";"
    N=F_INTERNAL_NUM
    WRITE(21,232) "        startFace       ",N,";"
    WRITE(21,*) "    }"    
        
        !DOWN FACE--PROCESSOR_MI1_TO_MI2  !CASE RECTANGLE 2-1-2
    PRO1=PNI*PNJ+PNL*PNI*(MJ2-1)+PNI*(MJ1-1)+MI-1
    MJ=1
    PRO2=PNI*(MJ-1)+MI-1
    WRITE(21,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(21,*) '    {'
    WRITE(21,*) '        type            processor;'
    WRITE(21,*) '        inGroups        List<word> 1(processor);'
    N=CNI*CNL
    WRITE(21,232) "        nFaces          ",N,";"
    N=F_INTERNAL_NUM
    WRITE(21,232) "        startFace       ",N,";"
    WRITE(21,*) '        matchTolerance  0.0001;'
    WRITE(21,232) "        myProcNo        ",PRO1,";"
    WRITE(21,232) "        neighbProcNo    ",PRO2,";"
    WRITE(21,*) '    }'
    WRITE(41,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(41,*) '    {'
    WRITE(41,*) '        type            processor;'
    WRITE(41,*) '        value           nonuniform List<scalar> '
    N=CNI*CNL
    WRITE(41,'(I0)') N
    WRITE(41,*) '('
    DO KP=1,CNL
    DO IP=1,CNI
    WRITE(41,'(E21.14)') P_DOWN_RECT(IP,KP)    
    ENDDO
    ENDDO
    WRITE(41,*) ')'
    WRITE(41,*) ';'
    WRITE(41,*) '    }'
    WRITE(51,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(51,*) '    {'
    WRITE(51,*) '        type            processor;'
    WRITE(51,*) '        value           nonuniform List<vector> '
    N=CNI*CNL
    WRITE(51,'(I0)') N
    WRITE(51,*) '('
    DO KP=1,CNL
    DO IP=1,CNI
    UU1=U_DOWN_RECT(1,IP,KP)
    UU2=U_DOWN_RECT(2,IP,KP)
    UU3=U_DOWN_RECT(3,IP,KP)
    WRITE(51,230) '( ',UU1,' ',UU2,' ',UU3,' )'  
    ENDDO
    ENDDO
    WRITE(51,*) ')'
    WRITE(51,*) ';'
    WRITE(51,*) '    }'
    WRITE(61,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(61,*) '    {'
    WRITE(61,*) '        type            processor;'
    WRITE(61,*) "        value           uniform 0;"
    WRITE(61,*) '    }'
    
        !FRONT FACE--PROCESSOR_MI1_TO_MI2  !CASE RECTANGLE 2-1-2
    PRO2=PNI*PNJ+PNL*PNI*(MJ2-2)+PNI*(MJ1-1)+MI-1
    WRITE(21,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(21,*) '    {'
    WRITE(21,*) '        type            processor;'
    WRITE(21,*) '        inGroups        List<word> 1(processor);'
    N=CNI*CNL
    WRITE(21,232) "        nFaces          ",N,";"
    N=F_INTERNAL_NUM+CNI*CNL
    WRITE(21,232) "        startFace       ",N,";"
    WRITE(21,*) '        matchTolerance  0.0001;'
    WRITE(21,232) "        myProcNo        ",PRO1,";"
    WRITE(21,232) "        neighbProcNo    ",PRO2,";"
    WRITE(21,*) '    }'
    WRITE(41,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(41,*) '    {'
    WRITE(41,*) '        type            processor;'
    WRITE(41,*) '        value           nonuniform List<scalar> '
    N=CNI*CNL
    WRITE(41,'(I0)') N
    WRITE(41,*) '('
    DO JP=1,CNL
    DO IP=1,CNI
    WRITE(41,'(E21.14)') P_FRONT_RECT(IP,JP)    
    ENDDO
    ENDDO
    WRITE(41,*) ')'
    WRITE(41,*) ';'
    WRITE(41,*) '    }'
    WRITE(51,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(51,*) '    {'
    WRITE(51,*) '        type            processor;'
    WRITE(51,*) '        value           nonuniform List<vector> '
    N=CNI*CNL
    WRITE(51,'(I0)') N
    WRITE(51,*) '('
    DO JP=1,CNL
    DO IP=1,CNI
    UU1=U_FRONT_RECT(1,IP,JP)
    UU2=U_FRONT_RECT(2,IP,JP)
    UU3=U_FRONT_RECT(3,IP,JP)
    WRITE(51,230) '( ',UU1,' ',UU2,' ',UU3,' )'  
    ENDDO
    ENDDO
    WRITE(51,*) ')'
    WRITE(51,*) ';'
    WRITE(51,*) '    }'
    WRITE(61,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(61,*) '    {'
    WRITE(61,*) '        type            processor;'
    WRITE(61,*) "        value           uniform 0;"
    WRITE(61,*) '    }'
    
        !RIGHT FACE--PROCESSOR_MI1_TO_MI2  !CASE RECTANGLE 2-1-2
    PRO2=PNI*PNJ+PNL*PNI*(MJ2-1)+PNI*(MJ1-1)+MI-2
    WRITE(21,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(21,*) '    {'
    WRITE(21,*) '        type            processor;'
    WRITE(21,*) '        inGroups        List<word> 1(processor);'
    N=CNL*CNL
    WRITE(21,232) "        nFaces          ",N,";"
    N=F_INTERNAL_NUM+CNI*CNL+CNI*CNL
    WRITE(21,232) "        startFace       ",N,";"
    WRITE(21,*) '        matchTolerance  0.0001;'
    WRITE(21,232) "        myProcNo        ",PRO1,";"
    WRITE(21,232) "        neighbProcNo    ",PRO2,";"
    WRITE(21,*) '    }'
    WRITE(41,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(41,*) '    {'
    WRITE(41,*) '        type            processor;'
    WRITE(41,*) '        value           nonuniform List<scalar> '
    N=CNL*CNL
    WRITE(41,'(I0)') N
    WRITE(41,*) '('
    DO KP=1,CNL
    DO JP=1,CNL
    WRITE(41,'(E21.14)') P_RIGHT_RECT(JP,KP)    
    ENDDO
    ENDDO
    WRITE(41,*) ')'
    WRITE(41,*) ';'
    WRITE(41,*) '    }'
    WRITE(51,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(51,*) '    {'
    WRITE(51,*) '        type            processor;'
    WRITE(51,*) '        value           nonuniform List<vector> '
    N=CNL*CNL
    WRITE(51,'(I0)') N
    WRITE(51,*) '('
    DO KP=1,CNL
    DO JP=1,CNL
    UU1=U_RIGHT_RECT(1,JP,KP)
    UU2=U_RIGHT_RECT(2,JP,KP)
    UU3=U_RIGHT_RECT(3,JP,KP)
    WRITE(51,230) '( ',UU1,' ',UU2,' ',UU3,' )'  
    ENDDO
    ENDDO
    WRITE(51,*) ')'
    WRITE(51,*) ';'
    WRITE(51,*) '    }'
    WRITE(61,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(61,*) '    {'
    WRITE(61,*) '        type            processor;'
    WRITE(61,*) "        value           uniform 0;"
    WRITE(61,*) '    }'
    
        !LEFT FACE--PROCESSOR_MI1_TO_MI2  !CASE RECTANGLE 2-1-2
    PRO2=PNI*PNJ+PNL*PNI*(MJ2-1)+PNI*(MJ1-1)+MI
    WRITE(21,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(21,*) '    {'
    WRITE(21,*) '        type            processor;'
    WRITE(21,*) '        inGroups        List<word> 1(processor);'
    N=CNL*CNL
    WRITE(21,232) "        nFaces          ",N,";"
    N=F_INTERNAL_NUM+CNI*CNL+CNI*CNL+CNL*CNL
    WRITE(21,232) "        startFace       ",N,";"
    WRITE(21,*) '        matchTolerance  0.0001;'
    WRITE(21,232) "        myProcNo        ",PRO1,";"
    WRITE(21,232) "        neighbProcNo    ",PRO2,";"
    WRITE(21,*) '    }'
    WRITE(41,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(41,*) '    {'
    WRITE(41,*) '        type            processor;'
    WRITE(41,*) '        value           nonuniform List<scalar> '
    N=CNL*CNL
    WRITE(41,'(I0)') N
    WRITE(41,*) '('
    DO KP=1,CNL
    DO JP=1,CNL
    WRITE(41,'(E21.14)') P_LEFT_RECT(JP,KP)    
    ENDDO
    ENDDO
    WRITE(41,*) ')'
    WRITE(41,*) ';'
    WRITE(41,*) '    }'
    WRITE(51,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(51,*) '    {'
    WRITE(51,*) '        type            processor;'
    WRITE(51,*) '        value           nonuniform List<vector> '
    N=CNL*CNL
    WRITE(51,'(I0)') N
    WRITE(51,*) '('
    DO KP=1,CNL
    DO JP=1,CNL
    UU1=U_LEFT_RECT(1,JP,KP)
    UU2=U_LEFT_RECT(2,JP,KP)
    UU3=U_LEFT_RECT(3,JP,KP)
    WRITE(51,230) '( ',UU1,' ',UU2,' ',UU3,' )'  
    ENDDO
    ENDDO
    WRITE(51,*) ')'
    WRITE(51,*) ';'
    WRITE(51,*) '    }'
    WRITE(61,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(61,*) '    {'
    WRITE(61,*) '        type            processor;'
    WRITE(61,*) "        value           uniform 0;"
    WRITE(61,*) '    }' 
    
        !UP FACE--PROCESSOR_MI1_TO_MI2  !CASE RECTANGLE 2-1-2
    PRO2=PNI*PNJ+PNL*PNI*(MJ2-1)+PNI*MJ1+MI-1
    WRITE(21,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(21,*) '    {'
    WRITE(21,*) '        type            processor;'
    WRITE(21,*) '        inGroups        List<word> 1(processor);'
    N=CNI*CNL
    WRITE(21,232) "        nFaces          ",N,";"
    N=F_INTERNAL_NUM+CNI*CNL+CNI*CNL+CNL*CNL+CNL*CNL
    WRITE(21,232) "        startFace       ",N,";"
    WRITE(21,*) '        matchTolerance  0.0001;'
    WRITE(21,232) "        myProcNo        ",PRO1,";"
    WRITE(21,232) "        neighbProcNo    ",PRO2,";"
    WRITE(21,*) '    }'
    WRITE(41,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(41,*) '    {'
    WRITE(41,*) '        type            processor;'
    WRITE(41,*) '        value           nonuniform List<scalar> '
    N=CNI*CNL
    WRITE(41,'(I0)') N
    WRITE(41,*) '('
    DO KP=1,CNL
    DO IP=1,CNI
    WRITE(41,'(E21.14)') P_UP_RECT(IP,KP)    
    ENDDO
    ENDDO
    WRITE(41,*) ')'
    WRITE(41,*) ';'
    WRITE(41,*) '    }'
    WRITE(51,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(51,*) '    {'
    WRITE(51,*) '        type            processor;'
    WRITE(51,*) '        value           nonuniform List<vector> '
    N=CNI*CNL
    WRITE(51,'(I0)') N
    WRITE(51,*) '('
    DO KP=1,CNL
    DO IP=1,CNI
    UU1=U_UP_RECT(1,IP,KP)
    UU2=U_UP_RECT(2,IP,KP)
    UU3=U_UP_RECT(3,IP,KP)
    WRITE(51,230) '( ',UU1,' ',UU2,' ',UU3,' )'  
    ENDDO
    ENDDO
    WRITE(51,*) ')'
    WRITE(51,*) ';'
    WRITE(51,*) '    }'
    WRITE(61,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(61,*) '    {'
    WRITE(61,*) '        type            processor;'
    WRITE(61,*) "        value           uniform 0;"
    WRITE(61,*) '    }'
    
        !BACK FACE--PROCESSOR_MI1_TO_MI2  !CASE RECTANGLE 2-1-2
    PRO2=PNI*PNJ+PNL*PNI*MJ2+PNI*(MJ1-1)+MI-1
    WRITE(21,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(21,*) '    {'
    WRITE(21,*) '        type            processor;'
    WRITE(21,*) '        inGroups        List<word> 1(processor);'
    N=CNI*CNL
    WRITE(21,232) "        nFaces          ",N,";"
    N=F_INTERNAL_NUM+CNI*CNL+CNI*CNL+CNL*CNL+CNL*CNL+CNI*CNL
    WRITE(21,232) "        startFace       ",N,";"
    WRITE(21,*) '        matchTolerance  0.0001;'
    WRITE(21,232) "        myProcNo        ",PRO1,";"
    WRITE(21,232) "        neighbProcNo    ",PRO2,";"
    WRITE(21,*) '    }'
    WRITE(41,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(41,*) '    {'
    WRITE(41,*) '        type            processor;'
    WRITE(41,*) '        value           nonuniform List<scalar> '
    N=CNI*CNL
    WRITE(41,'(I0)') N
    WRITE(41,*) '('
    DO JP=1,CNL
    DO IP=1,CNI
    WRITE(41,'(E21.14)') P_BACK_RECT(IP,JP)    
    ENDDO
    ENDDO
    WRITE(41,*) ')'
    WRITE(41,*) ';'
    WRITE(41,*) '    }'
    WRITE(51,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(51,*) '    {'
    WRITE(51,*) '        type            processor;'
    WRITE(51,*) '        value           nonuniform List<vector> '
    N=CNI*CNL
    WRITE(51,'(I0)') N
    WRITE(51,*) '('
    DO JP=1,CNL
    DO IP=1,CNI
    UU1=U_BACK_RECT(1,IP,JP)
    UU2=U_BACK_RECT(2,IP,JP)
    UU3=U_BACK_RECT(3,IP,JP)
    WRITE(51,230) '( ',UU1,' ',UU2,' ',UU3,' )'  
    ENDDO
    ENDDO
    WRITE(51,*) ')'
    WRITE(51,*) ';'
    WRITE(51,*) '    }'
    WRITE(61,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(61,*) '    {'
    WRITE(61,*) '        type            processor;'
    WRITE(61,*) "        value           uniform 0;"
    WRITE(61,*) '    }'    
        
        
        
        
        
    ENDIF    
    ELSEIF (MJ1.EQ.PNL) THEN
    IF (MJ2.EQ.1) THEN
    !CASE RECTANGLE 2-3-1
    
        !FRONT FACE--PROCESSOR_MI1_TO_MI2  !CASE RECTANGLE 2-3-1
    KP=1
    K=KP+CNL*(MJ2-1)
    DO JP=1,CNL
    J=JP+CNL*(MJ1-1)
    DO IP=1,CNI
    I=IP+CNI*(MI-1)
    F11=(CNL+1)*(CNI+1)*(KP-1)+(CNI+1)*(JP-1)+IP-1
    F12=F11+1
    F13=F12+CNI+1
    F14=F13-1
    F21=F11+(CNL+1)*(CNI+1)
    F22=F21+1
    F23=F22+CNI+1
    F24=F23-1
    CALL CAL_FIJ_RECT(F11,F12,F13,F14,F21,F22,F23,F24)
    WRITE(2,235) '4(',F11,' ',F14,' ',F13,' ',F12,')'
    FAD1=F_FRONT_RECT(IP,JP)
    WRITE(3,231) FAD1
    OWN=CELL_ID_RECT(IP,JP,KP)
    WRITE(4,231) OWN
    ENDDO   
    ENDDO    
        
        !UP FACE--PROCESSOR_MI1_TO_MI2  !CASE RECTANGLE 2-3-1
    JP=CNL
    J=JP+CNL*(MJ1-1)
    DO KP=1,CNL
    K=KP+CNL*(MJ2-1)
    DO IP=1,CNI
    I=IP+CNI*(MI-1)
    F11=(CNL+1)*(CNI+1)*(KP-1)+(CNI+1)*(JP-1)+IP-1
    F12=F11+1
    F13=F12+CNI+1
    F14=F13-1
    F21=F11+(CNL+1)*(CNI+1)
    F22=F21+1
    F23=F22+CNI+1
    F24=F23-1
    CALL CAL_FIJ_RECT(F11,F12,F13,F14,F21,F22,F23,F24)
    WRITE(2,235) '4(',F14,' ',F24,' ',F23,' ',F13,')'
    FAD1=F_UP_RECT(IP,KP)
    WRITE(3,231) FAD1
    OWN=CELL_ID_RECT(IP,JP,KP)
    WRITE(4,231) OWN
    ENDDO   
    ENDDO    
        
        !DOWN FACE--PROCESSOR_MI1_TO_MI2  !CASE RECTANGLE 2-3-1
    JP=1
    J=JP+CNL*(MJ1-1)
    DO KP=1,CNL
    K=KP+CNL*(MJ2-1)
    DO IP=1,CNI
    I=IP+CNI*(MI-1)
    F11=(CNL+1)*(CNI+1)*(KP-1)+(CNI+1)*(JP-1)+IP-1
    F12=F11+1
    F13=F12+CNI+1
    F14=F13-1
    F21=F11+(CNL+1)*(CNI+1)
    F22=F21+1
    F23=F22+CNI+1
    F24=F23-1
    CALL CAL_FIJ_RECT(F11,F12,F13,F14,F21,F22,F23,F24)
    WRITE(2,235) '4(',F11,' ',F12,' ',F22,' ',F21,')'
    FAD1=F_DOWN_RECT(IP,KP)
    WRITE(3,231) FAD1
    OWN=CELL_ID_RECT(IP,JP,KP)
    WRITE(4,231) OWN
    ENDDO   
    ENDDO
    
        !RIGHT FACE--PROCESSOR_MI1_TO_MI2  !CASE RECTANGLE 2-3-1 
    IP=1
    I=IP+CNI*(MI-1)
    DO KP=1,CNL
    K=KP+CNL*(MJ2-1)
    DO JP=1,CNL
    J=JP+CNL*(MJ1-1)
    F11=(CNL+1)*(CNI+1)*(KP-1)+(CNI+1)*(JP-1)+IP-1
    F12=F11+1
    F13=F12+CNI+1
    F14=F13-1
    F21=F11+(CNL+1)*(CNI+1)
    F22=F21+1
    F23=F22+CNI+1
    F24=F23-1
    CALL CAL_FIJ_RECT(F11,F12,F13,F14,F21,F22,F23,F24)
    WRITE(2,235) '4(',F11,' ',F21,' ',F24,' ',F14,')'
    FAD1=F_RIGHT_RECT(JP,KP)
    WRITE(3,231) FAD1
    OWN=CELL_ID_RECT(IP,JP,KP)
    WRITE(4,231) OWN
    ENDDO   
    ENDDO
    
        !LEFT FACE--PROCESSOR_MI1_TO_MI2  !CASE RECTANGLE 2-3-1 
    IP=CNI
    I=IP+CNI*(MI-1)
    DO KP=1,CNL
    K=KP+CNL*(MJ2-1)
    DO JP=1,CNL
    J=JP+CNL*(MJ1-1)
    F11=(CNL+1)*(CNI+1)*(KP-1)+(CNI+1)*(JP-1)+IP-1
    F12=F11+1
    F13=F12+CNI+1
    F14=F13-1
    F21=F11+(CNL+1)*(CNI+1)
    F22=F21+1
    F23=F22+CNI+1
    F24=F23-1
    CALL CAL_FIJ_RECT(F11,F12,F13,F14,F21,F22,F23,F24)
    WRITE(2,235) '4(',F12,' ',F13,' ',F23,' ',F22,')'
    FAD1=F_LEFT_RECT(JP,KP)
    WRITE(3,231) FAD1
    OWN=CELL_ID_RECT(IP,JP,KP)
    WRITE(4,231) OWN
    ENDDO   
    ENDDO
    
        !BACK FACE--PROCESSOR_MI1_TO_MI2  !CASE RECTANGLE 2-3-1
    KP=CNL
    K=KP+CNL*(MJ2-1)
    DO JP=1,CNL
    J=JP+CNL*(MJ1-1)
    DO IP=1,CNI
    I=IP+CNI*(MI-1)
    F11=(CNL+1)*(CNI+1)*(KP-1)+(CNI+1)*(JP-1)+IP-1
    F12=F11+1
    F13=F12+CNI+1
    F14=F13-1
    F21=F11+(CNL+1)*(CNI+1)
    F22=F21+1
    F23=F22+CNI+1
    F24=F23-1
    CALL CAL_FIJ_RECT(F11,F12,F13,F14,F21,F22,F23,F24)
    WRITE(2,235) '4(',F21,' ',F22,' ',F23,' ',F24,')'
    FAD1=F_BACK_RECT(IP,JP)
    WRITE(3,231) FAD1
    OWN=CELL_ID_RECT(IP,JP,KP)
    WRITE(4,231) OWN
    ENDDO   
    ENDDO    
        
    N=5
    WRITE(31,231) N+4
    WRITE(31,*) '('  
    WRITE(31,*) '0'
    WRITE(31,*) '1'
    WRITE(31,*) '2' 
    WRITE(31,*) '3' 
    DO I=1,N
    WRITE(31,*) '-1' 
    ENDDO
    WRITE(31,*) ')'
    
    WRITE(21,231) N+4
    WRITE(21,*) "("    
    WRITE(21,*) "    INLET1"
    WRITE(21,*) "    {"    
    WRITE(21,*) "        type            patch;"
    N=0
    WRITE(21,232) "        nFaces          ",N,";"
    N=F_INTERNAL_NUM
    WRITE(21,232) "        startFace       ",N,";"
    WRITE(21,*) "    }"
    WRITE(21,*) "    INLET2"
    WRITE(21,*) "    {"    
    WRITE(21,*) "        type            patch;"
    N=0
    WRITE(21,232) "        nFaces          ",N,";"
    N=F_INTERNAL_NUM
    WRITE(21,232) "        startFace       ",N,";"
    WRITE(21,*) "    }"
    WRITE(21,*) "    OUTLET"
    WRITE(21,*) "    {"
    WRITE(21,*) "        type            patch;"
    N=0
    WRITE(21,232) "        nFaces          ",N,";"
    N=F_INTERNAL_NUM
    WRITE(21,232) "        startFace       ",N,";"
    WRITE(21,*) "    }"
    WRITE(21,*) "    HULL"
    WRITE(21,*) "    {"
    WRITE(21,*) "        type            wall;"
    WRITE(21,*) "        inGroups        1(wall);"
    N=0
    WRITE(21,232) "        nFaces          ",N,";"
    N=F_INTERNAL_NUM
    WRITE(21,232) "        startFace       ",N,";"
    WRITE(21,*) "    }"    
        
        !FRONT/UP FACE--PROCESSOR_MI1_TO_MI2  !CASE RECTANGLE 2-3-1
    PRO1=PNI*PNJ+PNL*PNI*(MJ2-1)+PNI*(MJ1-1)+MI-1
    MJ=1
    PRO2=PNI*(MJ-1)+MI-1
    WRITE(21,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(21,*) '    {'
    WRITE(21,*) '        type            processor;'
    WRITE(21,*) '        inGroups        List<word> 1(processor);'
    N=CNI*CNL*2
    WRITE(21,232) "        nFaces          ",N,";"
    N=F_INTERNAL_NUM
    WRITE(21,232) "        startFace       ",N,";"
    WRITE(21,*) '        matchTolerance  0.0001;'
    WRITE(21,232) "        myProcNo        ",PRO1,";"
    WRITE(21,232) "        neighbProcNo    ",PRO2,";"
    WRITE(21,*) '    }'
    WRITE(41,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(41,*) '    {'
    WRITE(41,*) '        type            processor;'
    WRITE(41,*) '        value           nonuniform List<scalar> '
    N=CNI*CNL*2
    WRITE(41,'(I0)') N
    WRITE(41,*) '('
    DO JP=1,CNL
    DO IP=1,CNI
    WRITE(41,'(E21.14)') P_FRONT_RECT(IP,JP)    
    ENDDO
    ENDDO
    DO KP=1,CNL
    DO IP=1,CNI
    WRITE(41,'(E21.14)') P_UP_RECT(IP,KP)    
    ENDDO
    ENDDO
    WRITE(41,*) ')'
    WRITE(41,*) ';'
    WRITE(41,*) '    }'
    WRITE(51,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(51,*) '    {'
    WRITE(51,*) '        type            processor;'
    WRITE(51,*) '        value           nonuniform List<vector> '
    N=CNI*CNL*2
    WRITE(51,'(I0)') N
    WRITE(51,*) '('
    DO JP=1,CNL
    DO IP=1,CNI
    UU1=U_FRONT_RECT(1,IP,JP)
    UU2=U_FRONT_RECT(2,IP,JP)
    UU3=U_FRONT_RECT(3,IP,JP)
    WRITE(51,230) '( ',UU1,' ',UU2,' ',UU3,' )'  
    ENDDO
    ENDDO
    DO KP=1,CNL
    DO IP=1,CNI
    UU1=U_UP_RECT(1,IP,KP)
    UU2=U_UP_RECT(2,IP,KP)
    UU3=U_UP_RECT(3,IP,KP)
    WRITE(51,230) '( ',UU1,' ',UU2,' ',UU3,' )'  
    ENDDO
    ENDDO
    WRITE(51,*) ')'
    WRITE(51,*) ';'
    WRITE(51,*) '    }'
    WRITE(61,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(61,*) '    {'
    WRITE(61,*) '        type            processor;'
    WRITE(61,*) "        value           uniform 0;"
    WRITE(61,*) '    }'
    
        !DOWN FACE--PROCESSOR_MI1_TO_MI2  !CASE RECTANGLE 2-3-1
    PRO2=PNI*PNJ+PNL*PNI*(MJ2-1)+PNI*(MJ1-2)+MI-1
    WRITE(21,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(21,*) '    {'
    WRITE(21,*) '        type            processor;'
    WRITE(21,*) '        inGroups        List<word> 1(processor);'
    N=CNI*CNL
    WRITE(21,232) "        nFaces          ",N,";"
    N=F_INTERNAL_NUM+CNI*CNL*2
    WRITE(21,232) "        startFace       ",N,";"
    WRITE(21,*) '        matchTolerance  0.0001;'
    WRITE(21,232) "        myProcNo        ",PRO1,";"
    WRITE(21,232) "        neighbProcNo    ",PRO2,";"
    WRITE(21,*) '    }'
    WRITE(41,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(41,*) '    {'
    WRITE(41,*) '        type            processor;'
    WRITE(41,*) '        value           nonuniform List<scalar> '
    N=CNI*CNL
    WRITE(41,'(I0)') N
    WRITE(41,*) '('
    DO KP=1,CNL
    DO IP=1,CNI
    WRITE(41,'(E21.14)') P_DOWN_RECT(IP,KP)    
    ENDDO
    ENDDO
    WRITE(41,*) ')'
    WRITE(41,*) ';'
    WRITE(41,*) '    }'
    WRITE(51,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(51,*) '    {'
    WRITE(51,*) '        type            processor;'
    WRITE(51,*) '        value           nonuniform List<vector> '
    N=CNI*CNL
    WRITE(51,'(I0)') N
    WRITE(51,*) '('
    DO KP=1,CNL
    DO IP=1,CNI
    UU1=U_DOWN_RECT(1,IP,KP)
    UU2=U_DOWN_RECT(2,IP,KP)
    UU3=U_DOWN_RECT(3,IP,KP)
    WRITE(51,230) '( ',UU1,' ',UU2,' ',UU3,' )'  
    ENDDO
    ENDDO
    WRITE(51,*) ')'
    WRITE(51,*) ';'
    WRITE(51,*) '    }'
    WRITE(61,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(61,*) '    {'
    WRITE(61,*) '        type            processor;'
    WRITE(61,*) "        value           uniform 0;"
    WRITE(61,*) '    }'
    
        !RIGHT FACE--PROCESSOR_MI1_TO_MI2  !CASE RECTANGLE 2-3-1
    PRO2=PNI*PNJ+PNL*PNI*(MJ2-1)+PNI*(MJ1-1)+MI-2
    WRITE(21,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(21,*) '    {'
    WRITE(21,*) '        type            processor;'
    WRITE(21,*) '        inGroups        List<word> 1(processor);'
    N=CNL*CNL
    WRITE(21,232) "        nFaces          ",N,";"
    N=F_INTERNAL_NUM+CNI*CNL*2+CNI*CNL
    WRITE(21,232) "        startFace       ",N,";"
    WRITE(21,*) '        matchTolerance  0.0001;'
    WRITE(21,232) "        myProcNo        ",PRO1,";"
    WRITE(21,232) "        neighbProcNo    ",PRO2,";"
    WRITE(21,*) '    }'
    WRITE(41,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(41,*) '    {'
    WRITE(41,*) '        type            processor;'
    WRITE(41,*) '        value           nonuniform List<scalar> '
    N=CNL*CNL
    WRITE(41,'(I0)') N
    WRITE(41,*) '('
    DO KP=1,CNL
    DO JP=1,CNL
    WRITE(41,'(E21.14)') P_RIGHT_RECT(JP,KP)    
    ENDDO
    ENDDO
    WRITE(41,*) ')'
    WRITE(41,*) ';'
    WRITE(41,*) '    }'
    WRITE(51,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(51,*) '    {'
    WRITE(51,*) '        type            processor;'
    WRITE(51,*) '        value           nonuniform List<vector> '
    N=CNL*CNL
    WRITE(51,'(I0)') N
    WRITE(51,*) '('
    DO KP=1,CNL
    DO JP=1,CNL
    UU1=U_RIGHT_RECT(1,JP,KP)
    UU2=U_RIGHT_RECT(2,JP,KP)
    UU3=U_RIGHT_RECT(3,JP,KP)
    WRITE(51,230) '( ',UU1,' ',UU2,' ',UU3,' )'  
    ENDDO
    ENDDO
    WRITE(51,*) ')'
    WRITE(51,*) ';'
    WRITE(51,*) '    }'
    WRITE(61,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(61,*) '    {'
    WRITE(61,*) '        type            processor;'
    WRITE(61,*) "        value           uniform 0;"
    WRITE(61,*) '    }'
    
        !LEFT FACE--PROCESSOR_MI1_TO_MI2  !CASE RECTANGLE 2-3-1
    PRO2=PNI*PNJ+PNL*PNI*(MJ2-1)+PNI*(MJ1-1)+MI
    WRITE(21,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(21,*) '    {'
    WRITE(21,*) '        type            processor;'
    WRITE(21,*) '        inGroups        List<word> 1(processor);'
    N=CNL*CNL
    WRITE(21,232) "        nFaces          ",N,";"
    N=F_INTERNAL_NUM+CNI*CNL*2+CNI*CNL+CNL*CNL
    WRITE(21,232) "        startFace       ",N,";"
    WRITE(21,*) '        matchTolerance  0.0001;'
    WRITE(21,232) "        myProcNo        ",PRO1,";"
    WRITE(21,232) "        neighbProcNo    ",PRO2,";"
    WRITE(21,*) '    }'
    WRITE(41,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(41,*) '    {'
    WRITE(41,*) '        type            processor;'
    WRITE(41,*) '        value           nonuniform List<scalar> '
    N=CNL*CNL
    WRITE(41,'(I0)') N
    WRITE(41,*) '('
    DO KP=1,CNL
    DO JP=1,CNL
    WRITE(41,'(E21.14)') P_LEFT_RECT(JP,KP)    
    ENDDO
    ENDDO
    WRITE(41,*) ')'
    WRITE(41,*) ';'
    WRITE(41,*) '    }'
    WRITE(51,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(51,*) '    {'
    WRITE(51,*) '        type            processor;'
    WRITE(51,*) '        value           nonuniform List<vector> '
    N=CNL*CNL
    WRITE(51,'(I0)') N
    WRITE(51,*) '('
    DO KP=1,CNL
    DO JP=1,CNL
    UU1=U_LEFT_RECT(1,JP,KP)
    UU2=U_LEFT_RECT(2,JP,KP)
    UU3=U_LEFT_RECT(3,JP,KP)
    WRITE(51,230) '( ',UU1,' ',UU2,' ',UU3,' )'  
    ENDDO
    ENDDO
    WRITE(51,*) ')'
    WRITE(51,*) ';'
    WRITE(51,*) '    }'
    WRITE(61,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(61,*) '    {'
    WRITE(61,*) '        type            processor;'
    WRITE(61,*) "        value           uniform 0;"
    WRITE(61,*) '    }'
    
        !BACK FACE--PROCESSOR_MI1_TO_MI2  !CASE RECTANGLE 2-3-1
    PRO2=PNI*PNJ+PNL*PNI*MJ2+PNI*(MJ1-1)+MI-1
    WRITE(21,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(21,*) '    {'
    WRITE(21,*) '        type            processor;'
    WRITE(21,*) '        inGroups        List<word> 1(processor);'
    N=CNI*CNL
    WRITE(21,232) "        nFaces          ",N,";"
    N=F_INTERNAL_NUM+CNI*CNL*2+CNI*CNL+CNL*CNL+CNL*CNL
    WRITE(21,232) "        startFace       ",N,";"
    WRITE(21,*) '        matchTolerance  0.0001;'
    WRITE(21,232) "        myProcNo        ",PRO1,";"
    WRITE(21,232) "        neighbProcNo    ",PRO2,";"
    WRITE(21,*) '    }'
    WRITE(41,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(41,*) '    {'
    WRITE(41,*) '        type            processor;'
    WRITE(41,*) '        value           nonuniform List<scalar> '
    N=CNI*CNL
    WRITE(41,'(I0)') N
    WRITE(41,*) '('
    DO JP=1,CNL
    DO IP=1,CNI
    WRITE(41,'(E21.14)') P_BACK_RECT(IP,JP)    
    ENDDO
    ENDDO
    WRITE(41,*) ')'
    WRITE(41,*) ';'
    WRITE(41,*) '    }'
    WRITE(51,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(51,*) '    {'
    WRITE(51,*) '        type            processor;'
    WRITE(51,*) '        value           nonuniform List<vector> '
    N=CNI*CNL
    WRITE(51,'(I0)') N
    WRITE(51,*) '('
    DO JP=1,CNL
    DO IP=1,CNI
    UU1=U_BACK_RECT(1,IP,JP)
    UU2=U_BACK_RECT(2,IP,JP)
    UU3=U_BACK_RECT(3,IP,JP)
    WRITE(51,230) '( ',UU1,' ',UU2,' ',UU3,' )'  
    ENDDO
    ENDDO
    WRITE(51,*) ')'
    WRITE(51,*) ';'
    WRITE(51,*) '    }'
    WRITE(61,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(61,*) '    {'
    WRITE(61,*) '        type            processor;'
    WRITE(61,*) "        value           uniform 0;"
    WRITE(61,*) '    }'    
        
        
        
        
        
        
        
        
    ELSEIF (MJ2.EQ.PNL) THEN
    !CASE RECTANGLE 2-3-3
    
        !UP FACE--PROCESSOR_MI1_TO_MI2  !CASE RECTANGLE 2-3-3
    JP=CNL
    J=JP+CNL*(MJ1-1)
    DO KP=1,CNL
    K=KP+CNL*(MJ2-1)
    DO IP=1,CNI
    I=IP+CNI*(MI-1)
    F11=(CNL+1)*(CNI+1)*(KP-1)+(CNI+1)*(JP-1)+IP-1
    F12=F11+1
    F13=F12+CNI+1
    F14=F13-1
    F21=F11+(CNL+1)*(CNI+1)
    F22=F21+1
    F23=F22+CNI+1
    F24=F23-1
    CALL CAL_FIJ_RECT(F11,F12,F13,F14,F21,F22,F23,F24)
    WRITE(2,235) '4(',F14,' ',F24,' ',F23,' ',F13,')'
    FAD1=F_UP_RECT(IP,KP)
    WRITE(3,231) FAD1
    OWN=CELL_ID_RECT(IP,JP,KP)
    WRITE(4,231) OWN
    ENDDO   
    ENDDO
    
        !BACK FACE--PROCESSOR_MI1_TO_MI2  !CASE RECTANGLE 2-3-3
    KP=CNL
    K=KP+CNL*(MJ2-1)
    DO JP=CNL,1,-1
    J=JP+CNL*(MJ1-1)
    DO IP=1,CNI
    I=IP+CNI*(MI-1)
    F11=(CNL+1)*(CNI+1)*(KP-1)+(CNI+1)*(JP-1)+IP-1
    F12=F11+1
    F13=F12+CNI+1
    F14=F13-1
    F21=F11+(CNL+1)*(CNI+1)
    F22=F21+1
    F23=F22+CNI+1
    F24=F23-1
    CALL CAL_FIJ_RECT(F11,F12,F13,F14,F21,F22,F23,F24)
    !WRITE(2,235) '4(',F21,' ',F22,' ',F23,' ',F24,')'
    WRITE(2,235) '4(',F24,' ',F21,' ',F22,' ',F23,')'
    FAD1=F_BACK_RECT(IP,JP)
    WRITE(3,231) FAD1
    OWN=CELL_ID_RECT(IP,JP,KP)
    WRITE(4,231) OWN
    ENDDO   
    ENDDO    
        
        !FRONT FACE--PROCESSOR_MI1_TO_MI2  !CASE RECTANGLE 2-3-3
    KP=1
    K=KP+CNL*(MJ2-1)
    DO JP=1,CNL
    J=JP+CNL*(MJ1-1)
    DO IP=1,CNI
    I=IP+CNI*(MI-1)
    F11=(CNL+1)*(CNI+1)*(KP-1)+(CNI+1)*(JP-1)+IP-1
    F12=F11+1
    F13=F12+CNI+1
    F14=F13-1
    F21=F11+(CNL+1)*(CNI+1)
    F22=F21+1
    F23=F22+CNI+1
    F24=F23-1
    CALL CAL_FIJ_RECT(F11,F12,F13,F14,F21,F22,F23,F24)
    WRITE(2,235) '4(',F11,' ',F14,' ',F13,' ',F12,')'
    FAD1=F_FRONT_RECT(IP,JP)
    WRITE(3,231) FAD1
    OWN=CELL_ID_RECT(IP,JP,KP)
    WRITE(4,231) OWN
    ENDDO   
    ENDDO
    
        !DOWN FACE--PROCESSOR_MI1_TO_MI2  !CASE RECTANGLE 2-3-3
    JP=1
    J=JP+CNL*(MJ1-1)
    DO KP=1,CNL
    K=KP+CNL*(MJ2-1)
    DO IP=1,CNI
    I=IP+CNI*(MI-1)
    F11=(CNL+1)*(CNI+1)*(KP-1)+(CNI+1)*(JP-1)+IP-1
    F12=F11+1
    F13=F12+CNI+1
    F14=F13-1
    F21=F11+(CNL+1)*(CNI+1)
    F22=F21+1
    F23=F22+CNI+1
    F24=F23-1
    CALL CAL_FIJ_RECT(F11,F12,F13,F14,F21,F22,F23,F24)
    WRITE(2,235) '4(',F11,' ',F12,' ',F22,' ',F21,')'
    FAD1=F_DOWN_RECT(IP,KP)
    WRITE(3,231) FAD1
    OWN=CELL_ID_RECT(IP,JP,KP)
    WRITE(4,231) OWN
    ENDDO   
    ENDDO
    
        !RIGHT FACE--PROCESSOR_MI1_TO_MI2  !CASE RECTANGLE 2-3-3 
    IP=1
    I=IP+CNI*(MI-1)
    DO KP=1,CNL
    K=KP+CNL*(MJ2-1)
    DO JP=1,CNL
    J=JP+CNL*(MJ1-1)
    F11=(CNL+1)*(CNI+1)*(KP-1)+(CNI+1)*(JP-1)+IP-1
    F12=F11+1
    F13=F12+CNI+1
    F14=F13-1
    F21=F11+(CNL+1)*(CNI+1)
    F22=F21+1
    F23=F22+CNI+1
    F24=F23-1
    CALL CAL_FIJ_RECT(F11,F12,F13,F14,F21,F22,F23,F24)
    WRITE(2,235) '4(',F11,' ',F21,' ',F24,' ',F14,')'
    FAD1=F_RIGHT_RECT(JP,KP)
    WRITE(3,231) FAD1
    OWN=CELL_ID_RECT(IP,JP,KP)
    WRITE(4,231) OWN
    ENDDO   
    ENDDO
    
        !LEFT FACE--PROCESSOR_MI1_TO_MI2  !CASE RECTANGLE 2-3-3 
    IP=CNI
    I=IP+CNI*(MI-1)
    DO KP=1,CNL
    K=KP+CNL*(MJ2-1)
    DO JP=1,CNL
    J=JP+CNL*(MJ1-1)
    F11=(CNL+1)*(CNI+1)*(KP-1)+(CNI+1)*(JP-1)+IP-1
    F12=F11+1
    F13=F12+CNI+1
    F14=F13-1
    F21=F11+(CNL+1)*(CNI+1)
    F22=F21+1
    F23=F22+CNI+1
    F24=F23-1
    CALL CAL_FIJ_RECT(F11,F12,F13,F14,F21,F22,F23,F24)
    WRITE(2,235) '4(',F12,' ',F13,' ',F23,' ',F22,')'
    FAD1=F_LEFT_RECT(JP,KP)
    WRITE(3,231) FAD1
    OWN=CELL_ID_RECT(IP,JP,KP)
    WRITE(4,231) OWN
    ENDDO   
    ENDDO    
        
    N=5
    WRITE(31,231) N+4
    WRITE(31,*) '('  
    WRITE(31,*) '0'
    WRITE(31,*) '1'
    WRITE(31,*) '2' 
    WRITE(31,*) '3' 
    DO I=1,N
    WRITE(31,*) '-1' 
    ENDDO
    WRITE(31,*) ')'
    
    WRITE(21,231) N+4
    WRITE(21,*) "("    
    WRITE(21,*) "    INLET1"
    WRITE(21,*) "    {"    
    WRITE(21,*) "        type            patch;"
    N=0
    WRITE(21,232) "        nFaces          ",N,";"
    N=F_INTERNAL_NUM
    WRITE(21,232) "        startFace       ",N,";"
    WRITE(21,*) "    }"
    WRITE(21,*) "    INLET2"
    WRITE(21,*) "    {"    
    WRITE(21,*) "        type            patch;"
    N=0
    WRITE(21,232) "        nFaces          ",N,";"
    N=F_INTERNAL_NUM
    WRITE(21,232) "        startFace       ",N,";"
    WRITE(21,*) "    }"
    WRITE(21,*) "    OUTLET"
    WRITE(21,*) "    {"
    WRITE(21,*) "        type            patch;"
    N=0
    WRITE(21,232) "        nFaces          ",N,";"
    N=F_INTERNAL_NUM
    WRITE(21,232) "        startFace       ",N,";"
    WRITE(21,*) "    }"
    WRITE(21,*) "    HULL"
    WRITE(21,*) "    {"
    WRITE(21,*) "        type            wall;"
    WRITE(21,*) "        inGroups        1(wall);"
    N=0
    WRITE(21,232) "        nFaces          ",N,";"
    N=F_INTERNAL_NUM
    WRITE(21,232) "        startFace       ",N,";"
    WRITE(21,*) "    }"    
        
        !UP/BACK FACE--PROCESSOR_MI1_TO_MI2  !CASE RECTANGLE 2-3-3
    PRO1=PNI*PNJ+PNL*PNI*(MJ2-1)+PNI*(MJ1-1)+MI-1
    MJ=1
    PRO2=PNI*(MJ-1)+MI-1
    WRITE(21,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(21,*) '    {'
    WRITE(21,*) '        type            processor;'
    WRITE(21,*) '        inGroups        List<word> 1(processor);'
    N=CNI*CNL*2
    WRITE(21,232) "        nFaces          ",N,";"
    N=F_INTERNAL_NUM
    WRITE(21,232) "        startFace       ",N,";"
    WRITE(21,*) '        matchTolerance  0.0001;'
    WRITE(21,232) "        myProcNo        ",PRO1,";"
    WRITE(21,232) "        neighbProcNo    ",PRO2,";"
    WRITE(21,*) '    }'
    WRITE(41,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(41,*) '    {'
    WRITE(41,*) '        type            processor;'
    WRITE(41,*) '        value           nonuniform List<scalar> '
    N=CNI*CNL*2
    WRITE(41,'(I0)') N
    WRITE(41,*) '('
    DO KP=1,CNL
    DO IP=1,CNI
    WRITE(41,'(E21.14)') P_UP_RECT(IP,KP)    
    ENDDO
    ENDDO
    DO JP=1,CNL
    DO IP=1,CNI
    WRITE(41,'(E21.14)') P_BACK_RECT(IP,JP)    
    ENDDO
    ENDDO
    WRITE(41,*) ')'
    WRITE(41,*) ';'
    WRITE(41,*) '    }'
    WRITE(51,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(51,*) '    {'
    WRITE(51,*) '        type            processor;'
    WRITE(51,*) '        value           nonuniform List<vector> '
    N=CNI*CNL*2
    WRITE(51,'(I0)') N
    WRITE(51,*) '('
    DO KP=1,CNL
    DO IP=1,CNI
    UU1=U_UP_RECT(1,IP,KP)
    UU2=U_UP_RECT(2,IP,KP)
    UU3=U_UP_RECT(3,IP,KP)
    WRITE(51,230) '( ',UU1,' ',UU2,' ',UU3,' )'  
    ENDDO
    ENDDO
    DO JP=1,CNL
    DO IP=1,CNI
    UU1=U_BACK_RECT(1,IP,JP)
    UU2=U_BACK_RECT(2,IP,JP)
    UU3=U_BACK_RECT(3,IP,JP)
    WRITE(51,230) '( ',UU1,' ',UU2,' ',UU3,' )'  
    ENDDO
    ENDDO
    WRITE(51,*) ')'
    WRITE(51,*) ';'
    WRITE(51,*) '    }'
    WRITE(61,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(61,*) '    {'
    WRITE(61,*) '        type            processor;'
    WRITE(61,*) "        value           uniform 0;"
    WRITE(61,*) '    }'
    
        !FRONT FACE--PROCESSOR_MI1_TO_MI2  !CASE RECTANGLE 2-3-3
    PRO2=PNI*PNJ+PNL*PNI*(MJ2-2)+PNI*(MJ1-1)+MI-1
    WRITE(21,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(21,*) '    {'
    WRITE(21,*) '        type            processor;'
    WRITE(21,*) '        inGroups        List<word> 1(processor);'
    N=CNI*CNL
    WRITE(21,232) "        nFaces          ",N,";"
    N=F_INTERNAL_NUM+CNI*CNL*2
    WRITE(21,232) "        startFace       ",N,";"
    WRITE(21,*) '        matchTolerance  0.0001;'
    WRITE(21,232) "        myProcNo        ",PRO1,";"
    WRITE(21,232) "        neighbProcNo    ",PRO2,";"
    WRITE(21,*) '    }'
    WRITE(41,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(41,*) '    {'
    WRITE(41,*) '        type            processor;'
    WRITE(41,*) '        value           nonuniform List<scalar> '
    N=CNI*CNL
    WRITE(41,'(I0)') N
    WRITE(41,*) '('
    DO JP=1,CNL
    DO IP=1,CNI
    WRITE(41,'(E21.14)') P_FRONT_RECT(IP,JP)    
    ENDDO
    ENDDO
    WRITE(41,*) ')'
    WRITE(41,*) ';'
    WRITE(41,*) '    }'
    WRITE(51,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(51,*) '    {'
    WRITE(51,*) '        type            processor;'
    WRITE(51,*) '        value           nonuniform List<vector> '
    N=CNI*CNL
    WRITE(51,'(I0)') N
    WRITE(51,*) '('
    DO JP=1,CNL
    DO IP=1,CNI
    UU1=U_FRONT_RECT(1,IP,JP)
    UU2=U_FRONT_RECT(2,IP,JP)
    UU3=U_FRONT_RECT(3,IP,JP)
    WRITE(51,230) '( ',UU1,' ',UU2,' ',UU3,' )'  
    ENDDO
    ENDDO
    WRITE(51,*) ')'
    WRITE(51,*) ';'
    WRITE(51,*) '    }'
    WRITE(61,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(61,*) '    {'
    WRITE(61,*) '        type            processor;'
    WRITE(61,*) "        value           uniform 0;"
    WRITE(61,*) '    }'
    
        !DOWN FACE--PROCESSOR_MI1_TO_MI2  !CASE RECTANGLE 2-3-3
    PRO2=PNI*PNJ+PNL*PNI*(MJ2-1)+PNI*(MJ1-2)+MI-1
    WRITE(21,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(21,*) '    {'
    WRITE(21,*) '        type            processor;'
    WRITE(21,*) '        inGroups        List<word> 1(processor);'
    N=CNI*CNL
    WRITE(21,232) "        nFaces          ",N,";"
    N=F_INTERNAL_NUM+CNI*CNL*2+CNI*CNL
    WRITE(21,232) "        startFace       ",N,";"
    WRITE(21,*) '        matchTolerance  0.0001;'
    WRITE(21,232) "        myProcNo        ",PRO1,";"
    WRITE(21,232) "        neighbProcNo    ",PRO2,";"
    WRITE(21,*) '    }'
    WRITE(41,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(41,*) '    {'
    WRITE(41,*) '        type            processor;'
    WRITE(41,*) '        value           nonuniform List<scalar> '
    N=CNI*CNL
    WRITE(41,'(I0)') N
    WRITE(41,*) '('
    DO KP=1,CNL
    DO IP=1,CNI
    WRITE(41,'(E21.14)') P_DOWN_RECT(IP,KP)    
    ENDDO
    ENDDO
    WRITE(41,*) ')'
    WRITE(41,*) ';'
    WRITE(41,*) '    }'
    WRITE(51,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(51,*) '    {'
    WRITE(51,*) '        type            processor;'
    WRITE(51,*) '        value           nonuniform List<vector> '
    N=CNI*CNL
    WRITE(51,'(I0)') N
    WRITE(51,*) '('
    DO KP=1,CNL
    DO IP=1,CNI
    UU1=U_DOWN_RECT(1,IP,KP)
    UU2=U_DOWN_RECT(2,IP,KP)
    UU3=U_DOWN_RECT(3,IP,KP)
    WRITE(51,230) '( ',UU1,' ',UU2,' ',UU3,' )'  
    ENDDO
    ENDDO
    WRITE(51,*) ')'
    WRITE(51,*) ';'
    WRITE(51,*) '    }'
    WRITE(61,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(61,*) '    {'
    WRITE(61,*) '        type            processor;'
    WRITE(61,*) "        value           uniform 0;"
    WRITE(61,*) '    }'
    
        !RIGHT FACE--PROCESSOR_MI1_TO_MI2  !CASE RECTANGLE 2-3-3
    PRO2=PNI*PNJ+PNL*PNI*(MJ2-1)+PNI*(MJ1-1)+MI-2
    WRITE(21,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(21,*) '    {'
    WRITE(21,*) '        type            processor;'
    WRITE(21,*) '        inGroups        List<word> 1(processor);'
    N=CNL*CNL
    WRITE(21,232) "        nFaces          ",N,";"
    N=F_INTERNAL_NUM+CNI*CNL*2+CNI*CNL+CNI*CNL
    WRITE(21,232) "        startFace       ",N,";"
    WRITE(21,*) '        matchTolerance  0.0001;'
    WRITE(21,232) "        myProcNo        ",PRO1,";"
    WRITE(21,232) "        neighbProcNo    ",PRO2,";"
    WRITE(21,*) '    }'
    WRITE(41,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(41,*) '    {'
    WRITE(41,*) '        type            processor;'
    WRITE(41,*) '        value           nonuniform List<scalar> '
    N=CNL*CNL
    WRITE(41,'(I0)') N
    WRITE(41,*) '('
    DO KP=1,CNL
    DO JP=1,CNL
    WRITE(41,'(E21.14)') P_RIGHT_RECT(JP,KP)    
    ENDDO
    ENDDO
    WRITE(41,*) ')'
    WRITE(41,*) ';'
    WRITE(41,*) '    }'
    WRITE(51,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(51,*) '    {'
    WRITE(51,*) '        type            processor;'
    WRITE(51,*) '        value           nonuniform List<vector> '
    N=CNL*CNL
    WRITE(51,'(I0)') N
    WRITE(51,*) '('
    DO KP=1,CNL
    DO JP=1,CNL
    UU1=U_RIGHT_RECT(1,JP,KP)
    UU2=U_RIGHT_RECT(2,JP,KP)
    UU3=U_RIGHT_RECT(3,JP,KP)
    WRITE(51,230) '( ',UU1,' ',UU2,' ',UU3,' )'  
    ENDDO
    ENDDO
    WRITE(51,*) ')'
    WRITE(51,*) ';'
    WRITE(51,*) '    }'
    WRITE(61,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(61,*) '    {'
    WRITE(61,*) '        type            processor;'
    WRITE(61,*) "        value           uniform 0;"
    WRITE(61,*) '    }' 
    
        !LEFT FACE--PROCESSOR_MI1_TO_MI2  !CASE RECTANGLE 2-3-3
    PRO2=PNI*PNJ+PNL*PNI*(MJ2-1)+PNI*(MJ1-1)+MI
    WRITE(21,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(21,*) '    {'
    WRITE(21,*) '        type            processor;'
    WRITE(21,*) '        inGroups        List<word> 1(processor);'
    N=CNL*CNL
    WRITE(21,232) "        nFaces          ",N,";"
    N=F_INTERNAL_NUM+CNI*CNL*2+CNI*CNL+CNI*CNL+CNL*CNL
    WRITE(21,232) "        startFace       ",N,";"
    WRITE(21,*) '        matchTolerance  0.0001;'
    WRITE(21,232) "        myProcNo        ",PRO1,";"
    WRITE(21,232) "        neighbProcNo    ",PRO2,";"
    WRITE(21,*) '    }'
    WRITE(41,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(41,*) '    {'
    WRITE(41,*) '        type            processor;'
    WRITE(41,*) '        value           nonuniform List<scalar> '
    N=CNL*CNL
    WRITE(41,'(I0)') N
    WRITE(41,*) '('
    DO KP=1,CNL
    DO JP=1,CNL
    WRITE(41,'(E21.14)') P_LEFT_RECT(JP,KP)    
    ENDDO
    ENDDO
    WRITE(41,*) ')'
    WRITE(41,*) ';'
    WRITE(41,*) '    }'
    WRITE(51,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(51,*) '    {'
    WRITE(51,*) '        type            processor;'
    WRITE(51,*) '        value           nonuniform List<vector> '
    N=CNL*CNL
    WRITE(51,'(I0)') N
    WRITE(51,*) '('
    DO KP=1,CNL
    DO JP=1,CNL
    UU1=U_LEFT_RECT(1,JP,KP)
    UU2=U_LEFT_RECT(2,JP,KP)
    UU3=U_LEFT_RECT(3,JP,KP)
    WRITE(51,230) '( ',UU1,' ',UU2,' ',UU3,' )'  
    ENDDO
    ENDDO
    WRITE(51,*) ')'
    WRITE(51,*) ';'
    WRITE(51,*) '    }'
    WRITE(61,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(61,*) '    {'
    WRITE(61,*) '        type            processor;'
    WRITE(61,*) "        value           uniform 0;"
    WRITE(61,*) '    }'      
        
        
    ELSE
    !CASE RECTANGLE 2-3-2
    
        !UP FACE--PROCESSOR_MI1_TO_MI2  !CASE RECTANGLE 2-3-2
    JP=CNL
    J=JP+CNL*(MJ1-1)
    DO KP=1,CNL
    K=KP+CNL*(MJ2-1)
    DO IP=1,CNI
    I=IP+CNI*(MI-1)
    F11=(CNL+1)*(CNI+1)*(KP-1)+(CNI+1)*(JP-1)+IP-1
    F12=F11+1
    F13=F12+CNI+1
    F14=F13-1
    F21=F11+(CNL+1)*(CNI+1)
    F22=F21+1
    F23=F22+CNI+1
    F24=F23-1
    CALL CAL_FIJ_RECT(F11,F12,F13,F14,F21,F22,F23,F24)
    WRITE(2,235) '4(',F14,' ',F24,' ',F23,' ',F13,')'
    FAD1=F_UP_RECT(IP,KP)
    WRITE(3,231) FAD1
    OWN=CELL_ID_RECT(IP,JP,KP)
    WRITE(4,231) OWN
    ENDDO   
    ENDDO    
        
        !FRONT FACE--PROCESSOR_MI1_TO_MI2  !CASE RECTANGLE 2-3-2
    KP=1
    K=KP+CNL*(MJ2-1)
    DO JP=1,CNL
    J=JP+CNL*(MJ1-1)
    DO IP=1,CNI
    I=IP+CNI*(MI-1)
    F11=(CNL+1)*(CNI+1)*(KP-1)+(CNI+1)*(JP-1)+IP-1
    F12=F11+1
    F13=F12+CNI+1
    F14=F13-1
    F21=F11+(CNL+1)*(CNI+1)
    F22=F21+1
    F23=F22+CNI+1
    F24=F23-1
    CALL CAL_FIJ_RECT(F11,F12,F13,F14,F21,F22,F23,F24)
    WRITE(2,235) '4(',F11,' ',F14,' ',F13,' ',F12,')'
    FAD1=F_FRONT_RECT(IP,JP)
    WRITE(3,231) FAD1
    OWN=CELL_ID_RECT(IP,JP,KP)
    WRITE(4,231) OWN
    ENDDO   
    ENDDO
    
        !DOWN FACE--PROCESSOR_MI1_TO_MI2  !CASE RECTANGLE 2-3-2
    JP=1
    J=JP+CNL*(MJ1-1)
    DO KP=1,CNL
    K=KP+CNL*(MJ2-1)
    DO IP=1,CNI
    I=IP+CNI*(MI-1)
    F11=(CNL+1)*(CNI+1)*(KP-1)+(CNI+1)*(JP-1)+IP-1
    F12=F11+1
    F13=F12+CNI+1
    F14=F13-1
    F21=F11+(CNL+1)*(CNI+1)
    F22=F21+1
    F23=F22+CNI+1
    F24=F23-1
    CALL CAL_FIJ_RECT(F11,F12,F13,F14,F21,F22,F23,F24)
    WRITE(2,235) '4(',F11,' ',F12,' ',F22,' ',F21,')'
    FAD1=F_DOWN_RECT(IP,KP)
    WRITE(3,231) FAD1
    OWN=CELL_ID_RECT(IP,JP,KP)
    WRITE(4,231) OWN
    ENDDO   
    ENDDO
    
        !RIGHT FACE--PROCESSOR_MI1_TO_MI2  !CASE RECTANGLE 2-3-2 
    IP=1
    I=IP+CNI*(MI-1)
    DO KP=1,CNL
    K=KP+CNL*(MJ2-1)
    DO JP=1,CNL
    J=JP+CNL*(MJ1-1)
    F11=(CNL+1)*(CNI+1)*(KP-1)+(CNI+1)*(JP-1)+IP-1
    F12=F11+1
    F13=F12+CNI+1
    F14=F13-1
    F21=F11+(CNL+1)*(CNI+1)
    F22=F21+1
    F23=F22+CNI+1
    F24=F23-1
    CALL CAL_FIJ_RECT(F11,F12,F13,F14,F21,F22,F23,F24)
    WRITE(2,235) '4(',F11,' ',F21,' ',F24,' ',F14,')'
    FAD1=F_RIGHT_RECT(JP,KP)
    WRITE(3,231) FAD1
    OWN=CELL_ID_RECT(IP,JP,KP)
    WRITE(4,231) OWN
    ENDDO   
    ENDDO
    
        !LEFT FACE--PROCESSOR_MI1_TO_MI2  !CASE RECTANGLE 2-3-2 
    IP=CNI
    I=IP+CNI*(MI-1)
    DO KP=1,CNL
    K=KP+CNL*(MJ2-1)
    DO JP=1,CNL
    J=JP+CNL*(MJ1-1)
    F11=(CNL+1)*(CNI+1)*(KP-1)+(CNI+1)*(JP-1)+IP-1
    F12=F11+1
    F13=F12+CNI+1
    F14=F13-1
    F21=F11+(CNL+1)*(CNI+1)
    F22=F21+1
    F23=F22+CNI+1
    F24=F23-1
    CALL CAL_FIJ_RECT(F11,F12,F13,F14,F21,F22,F23,F24)
    WRITE(2,235) '4(',F12,' ',F13,' ',F23,' ',F22,')'
    FAD1=F_LEFT_RECT(JP,KP)
    WRITE(3,231) FAD1
    OWN=CELL_ID_RECT(IP,JP,KP)
    WRITE(4,231) OWN
    ENDDO   
    ENDDO
    
        !BACK FACE--PROCESSOR_MI1_TO_MI2  !CASE RECTANGLE 2-3-2
    KP=CNL
    K=KP+CNL*(MJ2-1)
    DO JP=1,CNL
    J=JP+CNL*(MJ1-1)
    DO IP=1,CNI
    I=IP+CNI*(MI-1)
    F11=(CNL+1)*(CNI+1)*(KP-1)+(CNI+1)*(JP-1)+IP-1
    F12=F11+1
    F13=F12+CNI+1
    F14=F13-1
    F21=F11+(CNL+1)*(CNI+1)
    F22=F21+1
    F23=F22+CNI+1
    F24=F23-1
    CALL CAL_FIJ_RECT(F11,F12,F13,F14,F21,F22,F23,F24)
    WRITE(2,235) '4(',F21,' ',F22,' ',F23,' ',F24,')'
    FAD1=F_BACK_RECT(IP,JP)
    WRITE(3,231) FAD1
    OWN=CELL_ID_RECT(IP,JP,KP)
    WRITE(4,231) OWN
    ENDDO   
    ENDDO    
        
    N=6
    WRITE(31,231) N+4
    WRITE(31,*) '('  
    WRITE(31,*) '0'
    WRITE(31,*) '1'
    WRITE(31,*) '2' 
    WRITE(31,*) '3' 
    DO I=1,N
    WRITE(31,*) '-1' 
    ENDDO
    WRITE(31,*) ')'
    
    WRITE(21,231) N+4
    WRITE(21,*) "("    
    WRITE(21,*) "    INLET1"
    WRITE(21,*) "    {"    
    WRITE(21,*) "        type            patch;"
    N=0
    WRITE(21,232) "        nFaces          ",N,";"
    N=F_INTERNAL_NUM
    WRITE(21,232) "        startFace       ",N,";"
    WRITE(21,*) "    }"
    WRITE(21,*) "    INLET2"
    WRITE(21,*) "    {"    
    WRITE(21,*) "        type            patch;"
    N=0
    WRITE(21,232) "        nFaces          ",N,";"
    N=F_INTERNAL_NUM
    WRITE(21,232) "        startFace       ",N,";"
    WRITE(21,*) "    }"
    WRITE(21,*) "    OUTLET"
    WRITE(21,*) "    {"
    WRITE(21,*) "        type            patch;"
    N=0
    WRITE(21,232) "        nFaces          ",N,";"
    N=F_INTERNAL_NUM
    WRITE(21,232) "        startFace       ",N,";"
    WRITE(21,*) "    }"
    WRITE(21,*) "    HULL"
    WRITE(21,*) "    {"
    WRITE(21,*) "        type            wall;"
    WRITE(21,*) "        inGroups        1(wall);"
    N=0
    WRITE(21,232) "        nFaces          ",N,";"
    N=F_INTERNAL_NUM
    WRITE(21,232) "        startFace       ",N,";"
    WRITE(21,*) "    }"    
        
        !UP FACE--PROCESSOR_MI1_TO_MI2  !CASE RECTANGLE 2-3-2
    PRO1=PNI*PNJ+PNL*PNI*(MJ2-1)+PNI*(MJ1-1)+MI-1
    MJ=1
    PRO2=PNI*(MJ-1)+MI-1
    WRITE(21,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(21,*) '    {'
    WRITE(21,*) '        type            processor;'
    WRITE(21,*) '        inGroups        List<word> 1(processor);'
    N=CNI*CNL
    WRITE(21,232) "        nFaces          ",N,";"
    N=F_INTERNAL_NUM
    WRITE(21,232) "        startFace       ",N,";"
    WRITE(21,*) '        matchTolerance  0.0001;'
    WRITE(21,232) "        myProcNo        ",PRO1,";"
    WRITE(21,232) "        neighbProcNo    ",PRO2,";"
    WRITE(21,*) '    }'
    WRITE(41,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(41,*) '    {'
    WRITE(41,*) '        type            processor;'
    WRITE(41,*) '        value           nonuniform List<scalar> '
    N=CNI*CNL
    WRITE(41,'(I0)') N
    WRITE(41,*) '('
    DO KP=1,CNL
    DO IP=1,CNI
    WRITE(41,'(E21.14)') P_UP_RECT(IP,KP)    
    ENDDO
    ENDDO
    WRITE(41,*) ')'
    WRITE(41,*) ';'
    WRITE(41,*) '    }'
    WRITE(51,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(51,*) '    {'
    WRITE(51,*) '        type            processor;'
    WRITE(51,*) '        value           nonuniform List<vector> '
    N=CNI*CNL
    WRITE(51,'(I0)') N
    WRITE(51,*) '('
    DO KP=1,CNL
    DO IP=1,CNI
    UU1=U_UP_RECT(1,IP,KP)
    UU2=U_UP_RECT(2,IP,KP)
    UU3=U_UP_RECT(3,IP,KP)
    WRITE(51,230) '( ',UU1,' ',UU2,' ',UU3,' )'  
    ENDDO
    ENDDO
    WRITE(51,*) ')'
    WRITE(51,*) ';'
    WRITE(51,*) '    }'
    WRITE(61,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(61,*) '    {'
    WRITE(61,*) '        type            processor;'
    WRITE(61,*) "        value           uniform 0;"
    WRITE(61,*) '    }'
    
        !FRONT FACE--PROCESSOR_MI1_TO_MI2  !CASE RECTANGLE 2-3-2
    PRO2=PNI*PNJ+PNL*PNI*(MJ2-2)+PNI*(MJ1-1)+MI-1
    WRITE(21,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(21,*) '    {'
    WRITE(21,*) '        type            processor;'
    WRITE(21,*) '        inGroups        List<word> 1(processor);'
    N=CNI*CNL
    WRITE(21,232) "        nFaces          ",N,";"
    N=F_INTERNAL_NUM+CNI*CNL
    WRITE(21,232) "        startFace       ",N,";"
    WRITE(21,*) '        matchTolerance  0.0001;'
    WRITE(21,232) "        myProcNo        ",PRO1,";"
    WRITE(21,232) "        neighbProcNo    ",PRO2,";"
    WRITE(21,*) '    }'
    WRITE(41,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(41,*) '    {'
    WRITE(41,*) '        type            processor;'
    WRITE(41,*) '        value           nonuniform List<scalar> '
    N=CNI*CNL
    WRITE(41,'(I0)') N
    WRITE(41,*) '('
    DO JP=1,CNL
    DO IP=1,CNI
    WRITE(41,'(E21.14)') P_FRONT_RECT(IP,JP)    
    ENDDO
    ENDDO
    WRITE(41,*) ')'
    WRITE(41,*) ';'
    WRITE(41,*) '    }'
    WRITE(51,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(51,*) '    {'
    WRITE(51,*) '        type            processor;'
    WRITE(51,*) '        value           nonuniform List<vector> '
    N=CNI*CNL
    WRITE(51,'(I0)') N
    WRITE(51,*) '('
    DO JP=1,CNL
    DO IP=1,CNI
    UU1=U_FRONT_RECT(1,IP,JP)
    UU2=U_FRONT_RECT(2,IP,JP)
    UU3=U_FRONT_RECT(3,IP,JP)
    WRITE(51,230) '( ',UU1,' ',UU2,' ',UU3,' )'  
    ENDDO
    ENDDO
    WRITE(51,*) ')'
    WRITE(51,*) ';'
    WRITE(51,*) '    }'
    WRITE(61,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(61,*) '    {'
    WRITE(61,*) '        type            processor;'
    WRITE(61,*) "        value           uniform 0;"
    WRITE(61,*) '    }'
    
        !DOWN FACE--PROCESSOR_MI1_TO_MI2  !CASE RECTANGLE 2-3-2
    PRO2=PNI*PNJ+PNL*PNI*(MJ2-1)+PNI*(MJ1-2)+MI-1
    WRITE(21,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(21,*) '    {'
    WRITE(21,*) '        type            processor;'
    WRITE(21,*) '        inGroups        List<word> 1(processor);'
    N=CNI*CNL
    WRITE(21,232) "        nFaces          ",N,";"
    N=F_INTERNAL_NUM+CNI*CNL+CNI*CNL
    WRITE(21,232) "        startFace       ",N,";"
    WRITE(21,*) '        matchTolerance  0.0001;'
    WRITE(21,232) "        myProcNo        ",PRO1,";"
    WRITE(21,232) "        neighbProcNo    ",PRO2,";"
    WRITE(21,*) '    }'
    WRITE(41,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(41,*) '    {'
    WRITE(41,*) '        type            processor;'
    WRITE(41,*) '        value           nonuniform List<scalar> '
    N=CNI*CNL
    WRITE(41,'(I0)') N
    WRITE(41,*) '('
    DO KP=1,CNL
    DO IP=1,CNI
    WRITE(41,'(E21.14)') P_DOWN_RECT(IP,KP)    
    ENDDO
    ENDDO
    WRITE(41,*) ')'
    WRITE(41,*) ';'
    WRITE(41,*) '    }'
    WRITE(51,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(51,*) '    {'
    WRITE(51,*) '        type            processor;'
    WRITE(51,*) '        value           nonuniform List<vector> '
    N=CNI*CNL
    WRITE(51,'(I0)') N
    WRITE(51,*) '('
    DO KP=1,CNL
    DO IP=1,CNI
    UU1=U_DOWN_RECT(1,IP,KP)
    UU2=U_DOWN_RECT(2,IP,KP)
    UU3=U_DOWN_RECT(3,IP,KP)
    WRITE(51,230) '( ',UU1,' ',UU2,' ',UU3,' )'  
    ENDDO
    ENDDO
    WRITE(51,*) ')'
    WRITE(51,*) ';'
    WRITE(51,*) '    }'
    WRITE(61,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(61,*) '    {'
    WRITE(61,*) '        type            processor;'
    WRITE(61,*) "        value           uniform 0;"
    WRITE(61,*) '    }'
    
        !RIGHT FACE--PROCESSOR_MI1_TO_MI2  !CASE RECTANGLE 2-3-2
    PRO2=PNI*PNJ+PNL*PNI*(MJ2-1)+PNI*(MJ1-1)+MI-2
    WRITE(21,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(21,*) '    {'
    WRITE(21,*) '        type            processor;'
    WRITE(21,*) '        inGroups        List<word> 1(processor);'
    N=CNL*CNL
    WRITE(21,232) "        nFaces          ",N,";"
    N=F_INTERNAL_NUM+CNI*CNL+CNI*CNL+CNI*CNL
    WRITE(21,232) "        startFace       ",N,";"
    WRITE(21,*) '        matchTolerance  0.0001;'
    WRITE(21,232) "        myProcNo        ",PRO1,";"
    WRITE(21,232) "        neighbProcNo    ",PRO2,";"
    WRITE(21,*) '    }'
    WRITE(41,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(41,*) '    {'
    WRITE(41,*) '        type            processor;'
    WRITE(41,*) '        value           nonuniform List<scalar> '
    N=CNL*CNL
    WRITE(41,'(I0)') N
    WRITE(41,*) '('
    DO KP=1,CNL
    DO JP=1,CNL
    WRITE(41,'(E21.14)') P_RIGHT_RECT(JP,KP)    
    ENDDO
    ENDDO
    WRITE(41,*) ')'
    WRITE(41,*) ';'
    WRITE(41,*) '    }'
    WRITE(51,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(51,*) '    {'
    WRITE(51,*) '        type            processor;'
    WRITE(51,*) '        value           nonuniform List<vector> '
    N=CNL*CNL
    WRITE(51,'(I0)') N
    WRITE(51,*) '('
    DO KP=1,CNL
    DO JP=1,CNL
    UU1=U_RIGHT_RECT(1,JP,KP)
    UU2=U_RIGHT_RECT(2,JP,KP)
    UU3=U_RIGHT_RECT(3,JP,KP)
    WRITE(51,230) '( ',UU1,' ',UU2,' ',UU3,' )'  
    ENDDO
    ENDDO
    WRITE(51,*) ')'
    WRITE(51,*) ';'
    WRITE(51,*) '    }'
    WRITE(61,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(61,*) '    {'
    WRITE(61,*) '        type            processor;'
    WRITE(61,*) "        value           uniform 0;"
    WRITE(61,*) '    }'
    
        !LEFT FACE--PROCESSOR_MI1_TO_MI2  !CASE RECTANGLE 2-3-2
    PRO2=PNI*PNJ+PNL*PNI*(MJ2-1)+PNI*(MJ1-1)+MI
    WRITE(21,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(21,*) '    {'
    WRITE(21,*) '        type            processor;'
    WRITE(21,*) '        inGroups        List<word> 1(processor);'
    N=CNL*CNL
    WRITE(21,232) "        nFaces          ",N,";"
    N=F_INTERNAL_NUM+CNI*CNL+CNI*CNL+CNI*CNL+CNL*CNL
    WRITE(21,232) "        startFace       ",N,";"
    WRITE(21,*) '        matchTolerance  0.0001;'
    WRITE(21,232) "        myProcNo        ",PRO1,";"
    WRITE(21,232) "        neighbProcNo    ",PRO2,";"
    WRITE(21,*) '    }'
    WRITE(41,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(41,*) '    {'
    WRITE(41,*) '        type            processor;'
    WRITE(41,*) '        value           nonuniform List<scalar> '
    N=CNL*CNL
    WRITE(41,'(I0)') N
    WRITE(41,*) '('
    DO KP=1,CNL
    DO JP=1,CNL
    WRITE(41,'(E21.14)') P_LEFT_RECT(JP,KP)    
    ENDDO
    ENDDO
    WRITE(41,*) ')'
    WRITE(41,*) ';'
    WRITE(41,*) '    }'
    WRITE(51,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(51,*) '    {'
    WRITE(51,*) '        type            processor;'
    WRITE(51,*) '        value           nonuniform List<vector> '
    N=CNL*CNL
    WRITE(51,'(I0)') N
    WRITE(51,*) '('
    DO KP=1,CNL
    DO JP=1,CNL
    UU1=U_LEFT_RECT(1,JP,KP)
    UU2=U_LEFT_RECT(2,JP,KP)
    UU3=U_LEFT_RECT(3,JP,KP)
    WRITE(51,230) '( ',UU1,' ',UU2,' ',UU3,' )'  
    ENDDO
    ENDDO
    WRITE(51,*) ')'
    WRITE(51,*) ';'
    WRITE(51,*) '    }'
    WRITE(61,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(61,*) '    {'
    WRITE(61,*) '        type            processor;'
    WRITE(61,*) "        value           uniform 0;"
    WRITE(61,*) '    }'
    
        !BACK FACE--PROCESSOR_MI1_TO_MI2  !CASE RECTANGLE 2-3-2
    PRO2=PNI*PNJ+PNL*PNI*MJ2+PNI*(MJ1-1)+MI-1
    WRITE(21,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(21,*) '    {'
    WRITE(21,*) '        type            processor;'
    WRITE(21,*) '        inGroups        List<word> 1(processor);'
    N=CNI*CNL
    WRITE(21,232) "        nFaces          ",N,";"
    N=F_INTERNAL_NUM+CNI*CNL+CNI*CNL+CNI*CNL+CNL*CNL+CNL*CNL
    WRITE(21,232) "        startFace       ",N,";"
    WRITE(21,*) '        matchTolerance  0.0001;'
    WRITE(21,232) "        myProcNo        ",PRO1,";"
    WRITE(21,232) "        neighbProcNo    ",PRO2,";"
    WRITE(21,*) '    }'
    WRITE(41,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(41,*) '    {'
    WRITE(41,*) '        type            processor;'
    WRITE(41,*) '        value           nonuniform List<scalar> '
    N=CNI*CNL
    WRITE(41,'(I0)') N
    WRITE(41,*) '('
    DO JP=1,CNL
    DO IP=1,CNI
    WRITE(41,'(E21.14)') P_BACK_RECT(IP,JP)    
    ENDDO
    ENDDO
    WRITE(41,*) ')'
    WRITE(41,*) ';'
    WRITE(41,*) '    }'
    WRITE(51,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(51,*) '    {'
    WRITE(51,*) '        type            processor;'
    WRITE(51,*) '        value           nonuniform List<vector> '
    N=CNI*CNL
    WRITE(51,'(I0)') N
    WRITE(51,*) '('
    DO JP=1,CNL
    DO IP=1,CNI
    UU1=U_BACK_RECT(1,IP,JP)
    UU2=U_BACK_RECT(2,IP,JP)
    UU3=U_BACK_RECT(3,IP,JP)
    WRITE(51,230) '( ',UU1,' ',UU2,' ',UU3,' )'  
    ENDDO
    ENDDO
    WRITE(51,*) ')'
    WRITE(51,*) ';'
    WRITE(51,*) '    }'
    WRITE(61,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(61,*) '    {'
    WRITE(61,*) '        type            processor;'
    WRITE(61,*) "        value           uniform 0;"
    WRITE(61,*) '    }'       
        
        
        
    ENDIF 
    ELSE
    IF (MJ2.EQ.1) THEN
    !CASE RECTANGLE 2-2-1
    
        !FRONT FACE--PROCESSOR_MI1_TO_MI2  !CASE RECTANGLE 2-2-1
    KP=1
    K=KP+CNL*(MJ2-1)
    DO JP=1,CNL
    J=JP+CNL*(MJ1-1)
    DO IP=1,CNI
    I=IP+CNI*(MI-1)
    F11=(CNL+1)*(CNI+1)*(KP-1)+(CNI+1)*(JP-1)+IP-1
    F12=F11+1
    F13=F12+CNI+1
    F14=F13-1
    F21=F11+(CNL+1)*(CNI+1)
    F22=F21+1
    F23=F22+CNI+1
    F24=F23-1
    CALL CAL_FIJ_RECT(F11,F12,F13,F14,F21,F22,F23,F24)
    WRITE(2,235) '4(',F11,' ',F14,' ',F13,' ',F12,')'
    FAD1=F_FRONT_RECT(IP,JP)
    WRITE(3,231) FAD1
    OWN=CELL_ID_RECT(IP,JP,KP)
    WRITE(4,231) OWN
    ENDDO   
    ENDDO
    
        !DOWN FACE--PROCESSOR_MI1_TO_MI2  !CASE RECTANGLE 2-2-1
    JP=1
    J=JP+CNL*(MJ1-1)
    DO KP=1,CNL
    K=KP+CNL*(MJ2-1)
    DO IP=1,CNI
    I=IP+CNI*(MI-1)
    F11=(CNL+1)*(CNI+1)*(KP-1)+(CNI+1)*(JP-1)+IP-1
    F12=F11+1
    F13=F12+CNI+1
    F14=F13-1
    F21=F11+(CNL+1)*(CNI+1)
    F22=F21+1
    F23=F22+CNI+1
    F24=F23-1
    CALL CAL_FIJ_RECT(F11,F12,F13,F14,F21,F22,F23,F24)
    WRITE(2,235) '4(',F11,' ',F12,' ',F22,' ',F21,')'
    FAD1=F_DOWN_RECT(IP,KP)
    WRITE(3,231) FAD1
    OWN=CELL_ID_RECT(IP,JP,KP)
    WRITE(4,231) OWN
    ENDDO   
    ENDDO
    
        !RIGHT FACE--PROCESSOR_MI1_TO_MI2  !CASE RECTANGLE 2-2-1 
    IP=1
    I=IP+CNI*(MI-1)
    DO KP=1,CNL
    K=KP+CNL*(MJ2-1)
    DO JP=1,CNL
    J=JP+CNL*(MJ1-1)
    F11=(CNL+1)*(CNI+1)*(KP-1)+(CNI+1)*(JP-1)+IP-1
    F12=F11+1
    F13=F12+CNI+1
    F14=F13-1
    F21=F11+(CNL+1)*(CNI+1)
    F22=F21+1
    F23=F22+CNI+1
    F24=F23-1
    CALL CAL_FIJ_RECT(F11,F12,F13,F14,F21,F22,F23,F24)
    WRITE(2,235) '4(',F11,' ',F21,' ',F24,' ',F14,')'
    FAD1=F_RIGHT_RECT(JP,KP)
    WRITE(3,231) FAD1
    OWN=CELL_ID_RECT(IP,JP,KP)
    WRITE(4,231) OWN
    ENDDO   
    ENDDO
    
        !LEFT FACE--PROCESSOR_MI1_TO_MI2  !CASE RECTANGLE 2-2-1 
    IP=CNI
    I=IP+CNI*(MI-1)
    DO KP=1,CNL
    K=KP+CNL*(MJ2-1)
    DO JP=1,CNL
    J=JP+CNL*(MJ1-1)
    F11=(CNL+1)*(CNI+1)*(KP-1)+(CNI+1)*(JP-1)+IP-1
    F12=F11+1
    F13=F12+CNI+1
    F14=F13-1
    F21=F11+(CNL+1)*(CNI+1)
    F22=F21+1
    F23=F22+CNI+1
    F24=F23-1
    CALL CAL_FIJ_RECT(F11,F12,F13,F14,F21,F22,F23,F24)
    WRITE(2,235) '4(',F12,' ',F13,' ',F23,' ',F22,')'
    FAD1=F_LEFT_RECT(JP,KP)
    WRITE(3,231) FAD1
    OWN=CELL_ID_RECT(IP,JP,KP)
    WRITE(4,231) OWN
    ENDDO   
    ENDDO
    
        !UP FACE--PROCESSOR_MI1_TO_MI2  !CASE RECTANGLE 2-2-1
    JP=CNL
    J=JP+CNL*(MJ1-1)
    DO KP=1,CNL
    K=KP+CNL*(MJ2-1)
    DO IP=1,CNI
    I=IP+CNI*(MI-1)
    F11=(CNL+1)*(CNI+1)*(KP-1)+(CNI+1)*(JP-1)+IP-1
    F12=F11+1
    F13=F12+CNI+1
    F14=F13-1
    F21=F11+(CNL+1)*(CNI+1)
    F22=F21+1
    F23=F22+CNI+1
    F24=F23-1
    CALL CAL_FIJ_RECT(F11,F12,F13,F14,F21,F22,F23,F24)
    WRITE(2,235) '4(',F14,' ',F24,' ',F23,' ',F13,')'
    FAD1=F_UP_RECT(IP,KP)
    WRITE(3,231) FAD1
    OWN=CELL_ID_RECT(IP,JP,KP)
    WRITE(4,231) OWN
    ENDDO   
    ENDDO
    
        !BACK FACE--PROCESSOR_MI1_TO_MI2  !CASE RECTANGLE 2-2-1
    KP=CNL
    K=KP+CNL*(MJ2-1)
    DO JP=1,CNL
    J=JP+CNL*(MJ1-1)
    DO IP=1,CNI
    I=IP+CNI*(MI-1)
    F11=(CNL+1)*(CNI+1)*(KP-1)+(CNI+1)*(JP-1)+IP-1
    F12=F11+1
    F13=F12+CNI+1
    F14=F13-1
    F21=F11+(CNL+1)*(CNI+1)
    F22=F21+1
    F23=F22+CNI+1
    F24=F23-1
    CALL CAL_FIJ_RECT(F11,F12,F13,F14,F21,F22,F23,F24)
    WRITE(2,235) '4(',F21,' ',F22,' ',F23,' ',F24,')'
    FAD1=F_BACK_RECT(IP,JP)
    WRITE(3,231) FAD1
    OWN=CELL_ID_RECT(IP,JP,KP)
    WRITE(4,231) OWN
    ENDDO   
    ENDDO    
        
    N=6
    WRITE(31,231) N+4
    WRITE(31,*) '('  
    WRITE(31,*) '0'
    WRITE(31,*) '1'
    WRITE(31,*) '2' 
    WRITE(31,*) '3' 
    DO I=1,N
    WRITE(31,*) '-1' 
    ENDDO
    WRITE(31,*) ')'
    
    WRITE(21,231) N+4
    WRITE(21,*) "("    
    WRITE(21,*) "    INLET1"
    WRITE(21,*) "    {"    
    WRITE(21,*) "        type            patch;"
    N=0
    WRITE(21,232) "        nFaces          ",N,";"
    N=F_INTERNAL_NUM
    WRITE(21,232) "        startFace       ",N,";"
    WRITE(21,*) "    }"
    WRITE(21,*) "    INLET2"
    WRITE(21,*) "    {"    
    WRITE(21,*) "        type            patch;"
    N=0
    WRITE(21,232) "        nFaces          ",N,";"
    N=F_INTERNAL_NUM
    WRITE(21,232) "        startFace       ",N,";"
    WRITE(21,*) "    }"
    WRITE(21,*) "    OUTLET"
    WRITE(21,*) "    {"
    WRITE(21,*) "        type            patch;"
    N=0
    WRITE(21,232) "        nFaces          ",N,";"
    N=F_INTERNAL_NUM
    WRITE(21,232) "        startFace       ",N,";"
    WRITE(21,*) "    }"
    WRITE(21,*) "    HULL"
    WRITE(21,*) "    {"
    WRITE(21,*) "        type            wall;"
    WRITE(21,*) "        inGroups        1(wall);"
    N=0
    WRITE(21,232) "        nFaces          ",N,";"
    N=F_INTERNAL_NUM
    WRITE(21,232) "        startFace       ",N,";"
    WRITE(21,*) "    }"    
        
        !FRONT FACE--PROCESSOR_MI1_TO_MI2  !CASE RECTANGLE 2-2-1
    PRO1=PNI*PNJ+PNL*PNI*(MJ2-1)+PNI*(MJ1-1)+MI-1
    MJ=1
    PRO2=PNI*(MJ-1)+MI-1
    WRITE(21,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(21,*) '    {'
    WRITE(21,*) '        type            processor;'
    WRITE(21,*) '        inGroups        List<word> 1(processor);'
    N=CNI*CNL
    WRITE(21,232) "        nFaces          ",N,";"
    N=F_INTERNAL_NUM
    WRITE(21,232) "        startFace       ",N,";"
    WRITE(21,*) '        matchTolerance  0.0001;'
    WRITE(21,232) "        myProcNo        ",PRO1,";"
    WRITE(21,232) "        neighbProcNo    ",PRO2,";"
    WRITE(21,*) '    }'
    WRITE(41,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(41,*) '    {'
    WRITE(41,*) '        type            processor;'
    WRITE(41,*) '        value           nonuniform List<scalar> '
    N=CNI*CNL
    WRITE(41,'(I0)') N
    WRITE(41,*) '('
    DO JP=1,CNL
    DO IP=1,CNI
    WRITE(41,'(E21.14)') P_FRONT_RECT(IP,JP)    
    ENDDO
    ENDDO
    WRITE(41,*) ')'
    WRITE(41,*) ';'
    WRITE(41,*) '    }'
    WRITE(51,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(51,*) '    {'
    WRITE(51,*) '        type            processor;'
    WRITE(51,*) '        value           nonuniform List<vector> '
    N=CNI*CNL
    WRITE(51,'(I0)') N
    WRITE(51,*) '('
    DO JP=1,CNL
    DO IP=1,CNI
    UU1=U_FRONT_RECT(1,IP,JP)
    UU2=U_FRONT_RECT(2,IP,JP)
    UU3=U_FRONT_RECT(3,IP,JP)
    WRITE(51,230) '( ',UU1,' ',UU2,' ',UU3,' )'  
    ENDDO
    ENDDO
    WRITE(51,*) ')'
    WRITE(51,*) ';'
    WRITE(51,*) '    }'
    WRITE(61,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(61,*) '    {'
    WRITE(61,*) '        type            processor;'
    WRITE(61,*) "        value           uniform 0;"
    WRITE(61,*) '    }'
    
        !DOWN FACE--PROCESSOR_MI1_TO_MI2  !CASE RECTANGLE 2-2-1
    PRO2=PNI*PNJ+PNL*PNI*(MJ2-1)+PNI*(MJ1-2)+MI-1
    WRITE(21,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(21,*) '    {'
    WRITE(21,*) '        type            processor;'
    WRITE(21,*) '        inGroups        List<word> 1(processor);'
    N=CNI*CNL
    WRITE(21,232) "        nFaces          ",N,";"
    N=F_INTERNAL_NUM+CNI*CNL
    WRITE(21,232) "        startFace       ",N,";"
    WRITE(21,*) '        matchTolerance  0.0001;'
    WRITE(21,232) "        myProcNo        ",PRO1,";"
    WRITE(21,232) "        neighbProcNo    ",PRO2,";"
    WRITE(21,*) '    }'
    WRITE(41,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(41,*) '    {'
    WRITE(41,*) '        type            processor;'
    WRITE(41,*) '        value           nonuniform List<scalar> '
    N=CNI*CNL
    WRITE(41,'(I0)') N
    WRITE(41,*) '('
    DO KP=1,CNL
    DO IP=1,CNI
    WRITE(41,'(E21.14)') P_DOWN_RECT(IP,KP)    
    ENDDO
    ENDDO
    WRITE(41,*) ')'
    WRITE(41,*) ';'
    WRITE(41,*) '    }'
    WRITE(51,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(51,*) '    {'
    WRITE(51,*) '        type            processor;'
    WRITE(51,*) '        value           nonuniform List<vector> '
    N=CNI*CNL
    WRITE(51,'(I0)') N
    WRITE(51,*) '('
    DO KP=1,CNL
    DO IP=1,CNI
    UU1=U_DOWN_RECT(1,IP,KP)
    UU2=U_DOWN_RECT(2,IP,KP)
    UU3=U_DOWN_RECT(3,IP,KP)
    WRITE(51,230) '( ',UU1,' ',UU2,' ',UU3,' )'  
    ENDDO
    ENDDO
    WRITE(51,*) ')'
    WRITE(51,*) ';'
    WRITE(51,*) '    }'
    WRITE(61,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(61,*) '    {'
    WRITE(61,*) '        type            processor;'
    WRITE(61,*) "        value           uniform 0;"
    WRITE(61,*) '    }'
    
        !RIGHT FACE--PROCESSOR_MI1_TO_MI2  !CASE RECTANGLE 2-2-1
    PRO2=PNI*PNJ+PNL*PNI*(MJ2-1)+PNI*(MJ1-1)+MI-2
    WRITE(21,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(21,*) '    {'
    WRITE(21,*) '        type            processor;'
    WRITE(21,*) '        inGroups        List<word> 1(processor);'
    N=CNL*CNL
    WRITE(21,232) "        nFaces          ",N,";"
    N=F_INTERNAL_NUM+CNI*CNL+CNI*CNL
    WRITE(21,232) "        startFace       ",N,";"
    WRITE(21,*) '        matchTolerance  0.0001;'
    WRITE(21,232) "        myProcNo        ",PRO1,";"
    WRITE(21,232) "        neighbProcNo    ",PRO2,";"
    WRITE(21,*) '    }'
    WRITE(41,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(41,*) '    {'
    WRITE(41,*) '        type            processor;'
    WRITE(41,*) '        value           nonuniform List<scalar> '
    N=CNL*CNL
    WRITE(41,'(I0)') N
    WRITE(41,*) '('
    DO KP=1,CNL
    DO JP=1,CNL
    WRITE(41,'(E21.14)') P_RIGHT_RECT(JP,KP)    
    ENDDO
    ENDDO
    WRITE(41,*) ')'
    WRITE(41,*) ';'
    WRITE(41,*) '    }'
    WRITE(51,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(51,*) '    {'
    WRITE(51,*) '        type            processor;'
    WRITE(51,*) '        value           nonuniform List<vector> '
    N=CNL*CNL
    WRITE(51,'(I0)') N
    WRITE(51,*) '('
    DO KP=1,CNL
    DO JP=1,CNL
    UU1=U_RIGHT_RECT(1,JP,KP)
    UU2=U_RIGHT_RECT(2,JP,KP)
    UU3=U_RIGHT_RECT(3,JP,KP)
    WRITE(51,230) '( ',UU1,' ',UU2,' ',UU3,' )'  
    ENDDO
    ENDDO
    WRITE(51,*) ')'
    WRITE(51,*) ';'
    WRITE(51,*) '    }'
    WRITE(61,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(61,*) '    {'
    WRITE(61,*) '        type            processor;'
    WRITE(61,*) "        value           uniform 0;"
    WRITE(61,*) '    }' 
    
        !LEFT FACE--PROCESSOR_MI1_TO_MI2  !CASE RECTANGLE 2-2-1
    PRO2=PNI*PNJ+PNL*PNI*(MJ2-1)+PNI*(MJ1-1)+MI
    WRITE(21,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(21,*) '    {'
    WRITE(21,*) '        type            processor;'
    WRITE(21,*) '        inGroups        List<word> 1(processor);'
    N=CNL*CNL
    WRITE(21,232) "        nFaces          ",N,";"
    N=F_INTERNAL_NUM+CNI*CNL+CNI*CNL+CNL*CNL
    WRITE(21,232) "        startFace       ",N,";"
    WRITE(21,*) '        matchTolerance  0.0001;'
    WRITE(21,232) "        myProcNo        ",PRO1,";"
    WRITE(21,232) "        neighbProcNo    ",PRO2,";"
    WRITE(21,*) '    }'
    WRITE(41,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(41,*) '    {'
    WRITE(41,*) '        type            processor;'
    WRITE(41,*) '        value           nonuniform List<scalar> '
    N=CNL*CNL
    WRITE(41,'(I0)') N
    WRITE(41,*) '('
    DO KP=1,CNL
    DO JP=1,CNL
    WRITE(41,'(E21.14)') P_LEFT_RECT(JP,KP)    
    ENDDO
    ENDDO
    WRITE(41,*) ')'
    WRITE(41,*) ';'
    WRITE(41,*) '    }'
    WRITE(51,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(51,*) '    {'
    WRITE(51,*) '        type            processor;'
    WRITE(51,*) '        value           nonuniform List<vector> '
    N=CNL*CNL
    WRITE(51,'(I0)') N
    WRITE(51,*) '('
    DO KP=1,CNL
    DO JP=1,CNL
    UU1=U_LEFT_RECT(1,JP,KP)
    UU2=U_LEFT_RECT(2,JP,KP)
    UU3=U_LEFT_RECT(3,JP,KP)
    WRITE(51,230) '( ',UU1,' ',UU2,' ',UU3,' )'  
    ENDDO
    ENDDO
    WRITE(51,*) ')'
    WRITE(51,*) ';'
    WRITE(51,*) '    }'
    WRITE(61,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(61,*) '    {'
    WRITE(61,*) '        type            processor;'
    WRITE(61,*) "        value           uniform 0;"
    WRITE(61,*) '    }' 
    
        !UP FACE--PROCESSOR_MI1_TO_MI2  !CASE RECTANGLE 2-2-1
    PRO2=PNI*PNJ+PNL*PNI*(MJ2-1)+PNI*MJ1+MI-1
    WRITE(21,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(21,*) '    {'
    WRITE(21,*) '        type            processor;'
    WRITE(21,*) '        inGroups        List<word> 1(processor);'
    N=CNI*CNL
    WRITE(21,232) "        nFaces          ",N,";"
    N=F_INTERNAL_NUM+CNI*CNL+CNI*CNL+CNL*CNL+CNL*CNL
    WRITE(21,232) "        startFace       ",N,";"
    WRITE(21,*) '        matchTolerance  0.0001;'
    WRITE(21,232) "        myProcNo        ",PRO1,";"
    WRITE(21,232) "        neighbProcNo    ",PRO2,";"
    WRITE(21,*) '    }'
    WRITE(41,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(41,*) '    {'
    WRITE(41,*) '        type            processor;'
    WRITE(41,*) '        value           nonuniform List<scalar> '
    N=CNI*CNL
    WRITE(41,'(I0)') N
    WRITE(41,*) '('
    DO KP=1,CNL
    DO IP=1,CNI
    WRITE(41,'(E21.14)') P_UP_RECT(IP,KP)    
    ENDDO
    ENDDO
    WRITE(41,*) ')'
    WRITE(41,*) ';'
    WRITE(41,*) '    }'
    WRITE(51,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(51,*) '    {'
    WRITE(51,*) '        type            processor;'
    WRITE(51,*) '        value           nonuniform List<vector> '
    N=CNI*CNL
    WRITE(51,'(I0)') N
    WRITE(51,*) '('
    DO KP=1,CNL
    DO IP=1,CNI
    UU1=U_UP_RECT(1,IP,KP)
    UU2=U_UP_RECT(2,IP,KP)
    UU3=U_UP_RECT(3,IP,KP)
    WRITE(51,230) '( ',UU1,' ',UU2,' ',UU3,' )'  
    ENDDO
    ENDDO
    WRITE(51,*) ')'
    WRITE(51,*) ';'
    WRITE(51,*) '    }'
    WRITE(61,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(61,*) '    {'
    WRITE(61,*) '        type            processor;'
    WRITE(61,*) "        value           uniform 0;"
    WRITE(61,*) '    }'
    
        !BACK FACE--PROCESSOR_MI1_TO_MI2  !CASE RECTANGLE 2-2-1
    PRO2=PNI*PNJ+PNL*PNI*MJ2+PNI*(MJ1-1)+MI-1
    WRITE(21,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(21,*) '    {'
    WRITE(21,*) '        type            processor;'
    WRITE(21,*) '        inGroups        List<word> 1(processor);'
    N=CNI*CNL
    WRITE(21,232) "        nFaces          ",N,";"
    N=F_INTERNAL_NUM+CNI*CNL+CNI*CNL+CNL*CNL+CNL*CNL+CNI*CNL
    WRITE(21,232) "        startFace       ",N,";"
    WRITE(21,*) '        matchTolerance  0.0001;'
    WRITE(21,232) "        myProcNo        ",PRO1,";"
    WRITE(21,232) "        neighbProcNo    ",PRO2,";"
    WRITE(21,*) '    }'
    WRITE(41,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(41,*) '    {'
    WRITE(41,*) '        type            processor;'
    WRITE(41,*) '        value           nonuniform List<scalar> '
    N=CNI*CNL
    WRITE(41,'(I0)') N
    WRITE(41,*) '('
    DO JP=1,CNL
    DO IP=1,CNI
    WRITE(41,'(E21.14)') P_BACK_RECT(IP,JP)    
    ENDDO
    ENDDO
    WRITE(41,*) ')'
    WRITE(41,*) ';'
    WRITE(41,*) '    }'
    WRITE(51,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(51,*) '    {'
    WRITE(51,*) '        type            processor;'
    WRITE(51,*) '        value           nonuniform List<vector> '
    N=CNI*CNL
    WRITE(51,'(I0)') N
    WRITE(51,*) '('
    DO JP=1,CNL
    DO IP=1,CNI
    UU1=U_BACK_RECT(1,IP,JP)
    UU2=U_BACK_RECT(2,IP,JP)
    UU3=U_BACK_RECT(3,IP,JP)
    WRITE(51,230) '( ',UU1,' ',UU2,' ',UU3,' )'  
    ENDDO
    ENDDO
    WRITE(51,*) ')'
    WRITE(51,*) ';'
    WRITE(51,*) '    }'
    WRITE(61,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(61,*) '    {'
    WRITE(61,*) '        type            processor;'
    WRITE(61,*) "        value           uniform 0;"
    WRITE(61,*) '    }'      
        
        
    ELSEIF (MJ2.EQ.PNL) THEN
    !CASE RECTANGLE 2-2-3
    
        !BACK FACE--PROCESSOR_MI1_TO_MI2  !CASE RECTANGLE 2-2-3
    KP=CNL
    K=KP+CNL*(MJ2-1)
    DO JP=CNL,1,-1
    J=JP+CNL*(MJ1-1)
    DO IP=1,CNI
    I=IP+CNI*(MI-1)
    F11=(CNL+1)*(CNI+1)*(KP-1)+(CNI+1)*(JP-1)+IP-1
    F12=F11+1
    F13=F12+CNI+1
    F14=F13-1
    F21=F11+(CNL+1)*(CNI+1)
    F22=F21+1
    F23=F22+CNI+1
    F24=F23-1
    CALL CAL_FIJ_RECT(F11,F12,F13,F14,F21,F22,F23,F24)
    !WRITE(2,235) '4(',F21,' ',F22,' ',F23,' ',F24,')'
    WRITE(2,235) '4(',F24,' ',F21,' ',F22,' ',F23,')'
    FAD1=F_BACK_RECT(IP,JP)
    WRITE(3,231) FAD1
    OWN=CELL_ID_RECT(IP,JP,KP)
    WRITE(4,231) OWN
    ENDDO   
    ENDDO    
        
        !FRONT FACE--PROCESSOR_MI1_TO_MI2  !CASE RECTANGLE 2-2-3
    KP=1
    K=KP+CNL*(MJ2-1)
    DO JP=1,CNL
    J=JP+CNL*(MJ1-1)
    DO IP=1,CNI
    I=IP+CNI*(MI-1)
    F11=(CNL+1)*(CNI+1)*(KP-1)+(CNI+1)*(JP-1)+IP-1
    F12=F11+1
    F13=F12+CNI+1
    F14=F13-1
    F21=F11+(CNL+1)*(CNI+1)
    F22=F21+1
    F23=F22+CNI+1
    F24=F23-1
    CALL CAL_FIJ_RECT(F11,F12,F13,F14,F21,F22,F23,F24)
    WRITE(2,235) '4(',F11,' ',F14,' ',F13,' ',F12,')'
    FAD1=F_FRONT_RECT(IP,JP)
    WRITE(3,231) FAD1
    OWN=CELL_ID_RECT(IP,JP,KP)
    WRITE(4,231) OWN
    ENDDO   
    ENDDO
    
        !DOWN FACE--PROCESSOR_MI1_TO_MI2  !CASE RECTANGLE 2-2-3
    JP=1
    J=JP+CNL*(MJ1-1)
    DO KP=1,CNL
    K=KP+CNL*(MJ2-1)
    DO IP=1,CNI
    I=IP+CNI*(MI-1)
    F11=(CNL+1)*(CNI+1)*(KP-1)+(CNI+1)*(JP-1)+IP-1
    F12=F11+1
    F13=F12+CNI+1
    F14=F13-1
    F21=F11+(CNL+1)*(CNI+1)
    F22=F21+1
    F23=F22+CNI+1
    F24=F23-1
    CALL CAL_FIJ_RECT(F11,F12,F13,F14,F21,F22,F23,F24)
    WRITE(2,235) '4(',F11,' ',F12,' ',F22,' ',F21,')'
    FAD1=F_DOWN_RECT(IP,KP)
    WRITE(3,231) FAD1
    OWN=CELL_ID_RECT(IP,JP,KP)
    WRITE(4,231) OWN
    ENDDO   
    ENDDO
    
        !RIGHT FACE--PROCESSOR_MI1_TO_MI2  !CASE RECTANGLE 2-2-3 
    IP=1
    I=IP+CNI*(MI-1)
    DO KP=1,CNL
    K=KP+CNL*(MJ2-1)
    DO JP=1,CNL
    J=JP+CNL*(MJ1-1)
    F11=(CNL+1)*(CNI+1)*(KP-1)+(CNI+1)*(JP-1)+IP-1
    F12=F11+1
    F13=F12+CNI+1
    F14=F13-1
    F21=F11+(CNL+1)*(CNI+1)
    F22=F21+1
    F23=F22+CNI+1
    F24=F23-1
    CALL CAL_FIJ_RECT(F11,F12,F13,F14,F21,F22,F23,F24)
    WRITE(2,235) '4(',F11,' ',F21,' ',F24,' ',F14,')'
    FAD1=F_RIGHT_RECT(JP,KP)
    WRITE(3,231) FAD1
    OWN=CELL_ID_RECT(IP,JP,KP)
    WRITE(4,231) OWN
    ENDDO   
    ENDDO
    
        !LEFT FACE--PROCESSOR_MI1_TO_MI2  !CASE RECTANGLE 2-2-3 
    IP=CNI
    I=IP+CNI*(MI-1)
    DO KP=1,CNL
    K=KP+CNL*(MJ2-1)
    DO JP=1,CNL
    J=JP+CNL*(MJ1-1)
    F11=(CNL+1)*(CNI+1)*(KP-1)+(CNI+1)*(JP-1)+IP-1
    F12=F11+1
    F13=F12+CNI+1
    F14=F13-1
    F21=F11+(CNL+1)*(CNI+1)
    F22=F21+1
    F23=F22+CNI+1
    F24=F23-1
    CALL CAL_FIJ_RECT(F11,F12,F13,F14,F21,F22,F23,F24)
    WRITE(2,235) '4(',F12,' ',F13,' ',F23,' ',F22,')'
    FAD1=F_LEFT_RECT(JP,KP)
    WRITE(3,231) FAD1
    OWN=CELL_ID_RECT(IP,JP,KP)
    WRITE(4,231) OWN
    ENDDO   
    ENDDO
    
        !UP FACE--PROCESSOR_MI1_TO_MI2  !CASE RECTANGLE 2-2-3
    JP=CNL
    J=JP+CNL*(MJ1-1)
    DO KP=1,CNL
    K=KP+CNL*(MJ2-1)
    DO IP=1,CNI
    I=IP+CNI*(MI-1)
    F11=(CNL+1)*(CNI+1)*(KP-1)+(CNI+1)*(JP-1)+IP-1
    F12=F11+1
    F13=F12+CNI+1
    F14=F13-1
    F21=F11+(CNL+1)*(CNI+1)
    F22=F21+1
    F23=F22+CNI+1
    F24=F23-1
    CALL CAL_FIJ_RECT(F11,F12,F13,F14,F21,F22,F23,F24)
    WRITE(2,235) '4(',F14,' ',F24,' ',F23,' ',F13,')'
    FAD1=F_UP_RECT(IP,KP)
    WRITE(3,231) FAD1
    OWN=CELL_ID_RECT(IP,JP,KP)
    WRITE(4,231) OWN
    ENDDO   
    ENDDO    
        
    N=6
    WRITE(31,231) N+4
    WRITE(31,*) '('  
    WRITE(31,*) '0'
    WRITE(31,*) '1'
    WRITE(31,*) '2' 
    WRITE(31,*) '3' 
    DO I=1,N
    WRITE(31,*) '-1' 
    ENDDO
    WRITE(31,*) ')'
    
    WRITE(21,231) N+4
    WRITE(21,*) "("    
    WRITE(21,*) "    INLET1"
    WRITE(21,*) "    {"    
    WRITE(21,*) "        type            patch;"
    N=0
    WRITE(21,232) "        nFaces          ",N,";"
    N=F_INTERNAL_NUM
    WRITE(21,232) "        startFace       ",N,";"
    WRITE(21,*) "    }"
    WRITE(21,*) "    INLET2"
    WRITE(21,*) "    {"    
    WRITE(21,*) "        type            patch;"
    N=0
    WRITE(21,232) "        nFaces          ",N,";"
    N=F_INTERNAL_NUM
    WRITE(21,232) "        startFace       ",N,";"
    WRITE(21,*) "    }"
    WRITE(21,*) "    OUTLET"
    WRITE(21,*) "    {"
    WRITE(21,*) "        type            patch;"
    N=0
    WRITE(21,232) "        nFaces          ",N,";"
    N=F_INTERNAL_NUM
    WRITE(21,232) "        startFace       ",N,";"
    WRITE(21,*) "    }"
    WRITE(21,*) "    HULL"
    WRITE(21,*) "    {"
    WRITE(21,*) "        type            wall;"
    WRITE(21,*) "        inGroups        1(wall);"
    N=0
    WRITE(21,232) "        nFaces          ",N,";"
    N=F_INTERNAL_NUM
    WRITE(21,232) "        startFace       ",N,";"
    WRITE(21,*) "    }"    
        
        !BACK FACE--PROCESSOR_MI1_TO_MI2  !CASE RECTANGLE 2-2-3
    PRO1=PNI*PNJ+PNL*PNI*(MJ2-1)+PNI*(MJ1-1)+MI-1
    MJ=1
    PRO2=PNI*(MJ-1)+MI-1
    WRITE(21,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(21,*) '    {'
    WRITE(21,*) '        type            processor;'
    WRITE(21,*) '        inGroups        List<word> 1(processor);'
    N=CNI*CNL
    WRITE(21,232) "        nFaces          ",N,";"
    N=F_INTERNAL_NUM
    WRITE(21,232) "        startFace       ",N,";"
    WRITE(21,*) '        matchTolerance  0.0001;'
    WRITE(21,232) "        myProcNo        ",PRO1,";"
    WRITE(21,232) "        neighbProcNo    ",PRO2,";"
    WRITE(21,*) '    }'
    WRITE(41,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(41,*) '    {'
    WRITE(41,*) '        type            processor;'
    WRITE(41,*) '        value           nonuniform List<scalar> '
    N=CNI*CNL
    WRITE(41,'(I0)') N
    WRITE(41,*) '('
    DO JP=1,CNL
    DO IP=1,CNI
    WRITE(41,'(E21.14)') P_BACK_RECT(IP,JP)    
    ENDDO
    ENDDO
    WRITE(41,*) ')'
    WRITE(41,*) ';'
    WRITE(41,*) '    }'
    WRITE(51,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(51,*) '    {'
    WRITE(51,*) '        type            processor;'
    WRITE(51,*) '        value           nonuniform List<vector> '
    N=CNI*CNL
    WRITE(51,'(I0)') N
    WRITE(51,*) '('
    DO JP=1,CNL
    DO IP=1,CNI
    UU1=U_BACK_RECT(1,IP,JP)
    UU2=U_BACK_RECT(2,IP,JP)
    UU3=U_BACK_RECT(3,IP,JP)
    WRITE(51,230) '( ',UU1,' ',UU2,' ',UU3,' )'  
    ENDDO
    ENDDO
    WRITE(51,*) ')'
    WRITE(51,*) ';'
    WRITE(51,*) '    }'
    WRITE(61,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(61,*) '    {'
    WRITE(61,*) '        type            processor;'
    WRITE(61,*) "        value           uniform 0;"
    WRITE(61,*) '    }'
    
        !FRONT FACE--PROCESSOR_MI1_TO_MI2  !CASE RECTANGLE 2-2-3
    PRO2=PNI*PNJ+PNL*PNI*(MJ2-2)+PNI*(MJ1-1)+MI-1
    WRITE(21,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(21,*) '    {'
    WRITE(21,*) '        type            processor;'
    WRITE(21,*) '        inGroups        List<word> 1(processor);'
    N=CNI*CNL
    WRITE(21,232) "        nFaces          ",N,";"
    N=F_INTERNAL_NUM+CNI*CNL
    WRITE(21,232) "        startFace       ",N,";"
    WRITE(21,*) '        matchTolerance  0.0001;'
    WRITE(21,232) "        myProcNo        ",PRO1,";"
    WRITE(21,232) "        neighbProcNo    ",PRO2,";"
    WRITE(21,*) '    }'
    WRITE(41,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(41,*) '    {'
    WRITE(41,*) '        type            processor;'
    WRITE(41,*) '        value           nonuniform List<scalar> '
    N=CNI*CNL
    WRITE(41,'(I0)') N
    WRITE(41,*) '('
    DO JP=1,CNL
    DO IP=1,CNI
    WRITE(41,'(E21.14)') P_FRONT_RECT(IP,JP)    
    ENDDO
    ENDDO
    WRITE(41,*) ')'
    WRITE(41,*) ';'
    WRITE(41,*) '    }'
    WRITE(51,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(51,*) '    {'
    WRITE(51,*) '        type            processor;'
    WRITE(51,*) '        value           nonuniform List<vector> '
    N=CNI*CNL
    WRITE(51,'(I0)') N
    WRITE(51,*) '('
    DO JP=1,CNL
    DO IP=1,CNI
    UU1=U_FRONT_RECT(1,IP,JP)
    UU2=U_FRONT_RECT(2,IP,JP)
    UU3=U_FRONT_RECT(3,IP,JP)
    WRITE(51,230) '( ',UU1,' ',UU2,' ',UU3,' )'  
    ENDDO
    ENDDO
    WRITE(51,*) ')'
    WRITE(51,*) ';'
    WRITE(51,*) '    }'
    WRITE(61,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(61,*) '    {'
    WRITE(61,*) '        type            processor;'
    WRITE(61,*) "        value           uniform 0;"
    WRITE(61,*) '    }'
    
        !DOWN FACE--PROCESSOR_MI1_TO_MI2  !CASE RECTANGLE 2-2-3
    PRO2=PNI*PNJ+PNL*PNI*(MJ2-1)+PNI*(MJ1-2)+MI-1
    WRITE(21,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(21,*) '    {'
    WRITE(21,*) '        type            processor;'
    WRITE(21,*) '        inGroups        List<word> 1(processor);'
    N=CNI*CNL
    WRITE(21,232) "        nFaces          ",N,";"
    N=F_INTERNAL_NUM+CNI*CNL+CNI*CNL
    WRITE(21,232) "        startFace       ",N,";"
    WRITE(21,*) '        matchTolerance  0.0001;'
    WRITE(21,232) "        myProcNo        ",PRO1,";"
    WRITE(21,232) "        neighbProcNo    ",PRO2,";"
    WRITE(21,*) '    }'
    WRITE(41,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(41,*) '    {'
    WRITE(41,*) '        type            processor;'
    WRITE(41,*) '        value           nonuniform List<scalar> '
    N=CNI*CNL
    WRITE(41,'(I0)') N
    WRITE(41,*) '('
    DO KP=1,CNL
    DO IP=1,CNI
    WRITE(41,'(E21.14)') P_DOWN_RECT(IP,KP)    
    ENDDO
    ENDDO
    WRITE(41,*) ')'
    WRITE(41,*) ';'
    WRITE(41,*) '    }'
    WRITE(51,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(51,*) '    {'
    WRITE(51,*) '        type            processor;'
    WRITE(51,*) '        value           nonuniform List<vector> '
    N=CNI*CNL
    WRITE(51,'(I0)') N
    WRITE(51,*) '('
    DO KP=1,CNL
    DO IP=1,CNI
    UU1=U_DOWN_RECT(1,IP,KP)
    UU2=U_DOWN_RECT(2,IP,KP)
    UU3=U_DOWN_RECT(3,IP,KP)
    WRITE(51,230) '( ',UU1,' ',UU2,' ',UU3,' )'  
    ENDDO
    ENDDO
    WRITE(51,*) ')'
    WRITE(51,*) ';'
    WRITE(51,*) '    }'
    WRITE(61,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(61,*) '    {'
    WRITE(61,*) '        type            processor;'
    WRITE(61,*) "        value           uniform 0;"
    WRITE(61,*) '    }'
    
        !RIGHT FACE--PROCESSOR_MI1_TO_MI2  !CASE RECTANGLE 2-2-3
    PRO2=PNI*PNJ+PNL*PNI*(MJ2-1)+PNI*(MJ1-1)+MI-2
    WRITE(21,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(21,*) '    {'
    WRITE(21,*) '        type            processor;'
    WRITE(21,*) '        inGroups        List<word> 1(processor);'
    N=CNL*CNL
    WRITE(21,232) "        nFaces          ",N,";"
    N=F_INTERNAL_NUM+CNI*CNL+CNI*CNL+CNI*CNL
    WRITE(21,232) "        startFace       ",N,";"
    WRITE(21,*) '        matchTolerance  0.0001;'
    WRITE(21,232) "        myProcNo        ",PRO1,";"
    WRITE(21,232) "        neighbProcNo    ",PRO2,";"
    WRITE(21,*) '    }'
    WRITE(41,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(41,*) '    {'
    WRITE(41,*) '        type            processor;'
    WRITE(41,*) '        value           nonuniform List<scalar> '
    N=CNL*CNL
    WRITE(41,'(I0)') N
    WRITE(41,*) '('
    DO KP=1,CNL
    DO JP=1,CNL
    WRITE(41,'(E21.14)') P_RIGHT_RECT(JP,KP)    
    ENDDO
    ENDDO
    WRITE(41,*) ')'
    WRITE(41,*) ';'
    WRITE(41,*) '    }'
    WRITE(51,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(51,*) '    {'
    WRITE(51,*) '        type            processor;'
    WRITE(51,*) '        value           nonuniform List<vector> '
    N=CNL*CNL
    WRITE(51,'(I0)') N
    WRITE(51,*) '('
    DO KP=1,CNL
    DO JP=1,CNL
    UU1=U_RIGHT_RECT(1,JP,KP)
    UU2=U_RIGHT_RECT(2,JP,KP)
    UU3=U_RIGHT_RECT(3,JP,KP)
    WRITE(51,230) '( ',UU1,' ',UU2,' ',UU3,' )'  
    ENDDO
    ENDDO
    WRITE(51,*) ')'
    WRITE(51,*) ';'
    WRITE(51,*) '    }'
    WRITE(61,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(61,*) '    {'
    WRITE(61,*) '        type            processor;'
    WRITE(61,*) "        value           uniform 0;"
    WRITE(61,*) '    }' 
    
        !LEFT FACE--PROCESSOR_MI1_TO_MI2  !CASE RECTANGLE 2-2-3
    PRO2=PNI*PNJ+PNL*PNI*(MJ2-1)+PNI*(MJ1-1)+MI
    WRITE(21,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(21,*) '    {'
    WRITE(21,*) '        type            processor;'
    WRITE(21,*) '        inGroups        List<word> 1(processor);'
    N=CNL*CNL
    WRITE(21,232) "        nFaces          ",N,";"
    N=F_INTERNAL_NUM+CNI*CNL+CNI*CNL+CNI*CNL+CNL*CNL
    WRITE(21,232) "        startFace       ",N,";"
    WRITE(21,*) '        matchTolerance  0.0001;'
    WRITE(21,232) "        myProcNo        ",PRO1,";"
    WRITE(21,232) "        neighbProcNo    ",PRO2,";"
    WRITE(21,*) '    }'
    WRITE(41,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(41,*) '    {'
    WRITE(41,*) '        type            processor;'
    WRITE(41,*) '        value           nonuniform List<scalar> '
    N=CNL*CNL
    WRITE(41,'(I0)') N
    WRITE(41,*) '('
    DO KP=1,CNL
    DO JP=1,CNL
    WRITE(41,'(E21.14)') P_LEFT_RECT(JP,KP)    
    ENDDO
    ENDDO
    WRITE(41,*) ')'
    WRITE(41,*) ';'
    WRITE(41,*) '    }'
    WRITE(51,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(51,*) '    {'
    WRITE(51,*) '        type            processor;'
    WRITE(51,*) '        value           nonuniform List<vector> '
    N=CNL*CNL
    WRITE(51,'(I0)') N
    WRITE(51,*) '('
    DO KP=1,CNL
    DO JP=1,CNL
    UU1=U_LEFT_RECT(1,JP,KP)
    UU2=U_LEFT_RECT(2,JP,KP)
    UU3=U_LEFT_RECT(3,JP,KP)
    WRITE(51,230) '( ',UU1,' ',UU2,' ',UU3,' )'  
    ENDDO
    ENDDO
    WRITE(51,*) ')'
    WRITE(51,*) ';'
    WRITE(51,*) '    }'
    WRITE(61,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(61,*) '    {'
    WRITE(61,*) '        type            processor;'
    WRITE(61,*) "        value           uniform 0;"
    WRITE(61,*) '    }' 
    
        !UP FACE--PROCESSOR_MI1_TO_MI2  !CASE RECTANGLE 2-2-3
    PRO2=PNI*PNJ+PNL*PNI*(MJ2-1)+PNI*MJ1+MI-1
    WRITE(21,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(21,*) '    {'
    WRITE(21,*) '        type            processor;'
    WRITE(21,*) '        inGroups        List<word> 1(processor);'
    N=CNI*CNL
    WRITE(21,232) "        nFaces          ",N,";"
    N=F_INTERNAL_NUM+CNI*CNL+CNI*CNL+CNI*CNL+CNL*CNL+CNL*CNL
    WRITE(21,232) "        startFace       ",N,";"
    WRITE(21,*) '        matchTolerance  0.0001;'
    WRITE(21,232) "        myProcNo        ",PRO1,";"
    WRITE(21,232) "        neighbProcNo    ",PRO2,";"
    WRITE(21,*) '    }'
    WRITE(41,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(41,*) '    {'
    WRITE(41,*) '        type            processor;'
    WRITE(41,*) '        value           nonuniform List<scalar> '
    N=CNI*CNL
    WRITE(41,'(I0)') N
    WRITE(41,*) '('
    DO KP=1,CNL
    DO IP=1,CNI
    WRITE(41,'(E21.14)') P_UP_RECT(IP,KP)    
    ENDDO
    ENDDO
    WRITE(41,*) ')'
    WRITE(41,*) ';'
    WRITE(41,*) '    }'
    WRITE(51,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(51,*) '    {'
    WRITE(51,*) '        type            processor;'
    WRITE(51,*) '        value           nonuniform List<vector> '
    N=CNI*CNL
    WRITE(51,'(I0)') N
    WRITE(51,*) '('
    DO KP=1,CNL
    DO IP=1,CNI
    UU1=U_UP_RECT(1,IP,KP)
    UU2=U_UP_RECT(2,IP,KP)
    UU3=U_UP_RECT(3,IP,KP)
    WRITE(51,230) '( ',UU1,' ',UU2,' ',UU3,' )'  
    ENDDO
    ENDDO
    WRITE(51,*) ')'
    WRITE(51,*) ';'
    WRITE(51,*) '    }'
    WRITE(61,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(61,*) '    {'
    WRITE(61,*) '        type            processor;'
    WRITE(61,*) "        value           uniform 0;"
    WRITE(61,*) '    }'    
        
        
        
        
    ELSE
    !CASE RECTANGLE 2-2-2
    
        !FRONT FACE--PROCESSOR_MI1_TO_MI2  !CASE RECTANGLE 2-2-2
    KP=1
    K=KP+CNL*(MJ2-1)
    DO JP=1,CNL
    J=JP+CNL*(MJ1-1)
    DO IP=1,CNI
    I=IP+CNI*(MI-1)
    F11=(CNL+1)*(CNI+1)*(KP-1)+(CNI+1)*(JP-1)+IP-1
    F12=F11+1
    F13=F12+CNI+1
    F14=F13-1
    F21=F11+(CNL+1)*(CNI+1)
    F22=F21+1
    F23=F22+CNI+1
    F24=F23-1
    CALL CAL_FIJ_RECT(F11,F12,F13,F14,F21,F22,F23,F24)
    WRITE(2,235) '4(',F11,' ',F14,' ',F13,' ',F12,')'
    FAD1=F_FRONT_RECT(IP,JP)
    WRITE(3,231) FAD1
    OWN=CELL_ID_RECT(IP,JP,KP)
    WRITE(4,231) OWN
    ENDDO   
    ENDDO
    
        !DOWN FACE--PROCESSOR_MI1_TO_MI2  !CASE RECTANGLE 2-2-2
    JP=1
    J=JP+CNL*(MJ1-1)
    DO KP=1,CNL
    K=KP+CNL*(MJ2-1)
    DO IP=1,CNI
    I=IP+CNI*(MI-1)
    F11=(CNL+1)*(CNI+1)*(KP-1)+(CNI+1)*(JP-1)+IP-1
    F12=F11+1
    F13=F12+CNI+1
    F14=F13-1
    F21=F11+(CNL+1)*(CNI+1)
    F22=F21+1
    F23=F22+CNI+1
    F24=F23-1
    CALL CAL_FIJ_RECT(F11,F12,F13,F14,F21,F22,F23,F24)
    WRITE(2,235) '4(',F11,' ',F12,' ',F22,' ',F21,')'
    FAD1=F_DOWN_RECT(IP,KP)
    WRITE(3,231) FAD1
    OWN=CELL_ID_RECT(IP,JP,KP)
    WRITE(4,231) OWN
    ENDDO   
    ENDDO
    
        !RIGHT FACE--PROCESSOR_MI1_TO_MI2  !CASE RECTANGLE 2-2-2 
    IP=1
    I=IP+CNI*(MI-1)
    DO KP=1,CNL
    K=KP+CNL*(MJ2-1)
    DO JP=1,CNL
    J=JP+CNL*(MJ1-1)
    F11=(CNL+1)*(CNI+1)*(KP-1)+(CNI+1)*(JP-1)+IP-1
    F12=F11+1
    F13=F12+CNI+1
    F14=F13-1
    F21=F11+(CNL+1)*(CNI+1)
    F22=F21+1
    F23=F22+CNI+1
    F24=F23-1
    CALL CAL_FIJ_RECT(F11,F12,F13,F14,F21,F22,F23,F24)
    WRITE(2,235) '4(',F11,' ',F21,' ',F24,' ',F14,')'
    FAD1=F_RIGHT_RECT(JP,KP)
    WRITE(3,231) FAD1
    OWN=CELL_ID_RECT(IP,JP,KP)
    WRITE(4,231) OWN
    ENDDO   
    ENDDO
    
        !LEFT FACE--PROCESSOR_MI1_TO_MI2  !CASE RECTANGLE 2-2-2 
    IP=CNI
    I=IP+CNI*(MI-1)
    DO KP=1,CNL
    K=KP+CNL*(MJ2-1)
    DO JP=1,CNL
    J=JP+CNL*(MJ1-1)
    F11=(CNL+1)*(CNI+1)*(KP-1)+(CNI+1)*(JP-1)+IP-1
    F12=F11+1
    F13=F12+CNI+1
    F14=F13-1
    F21=F11+(CNL+1)*(CNI+1)
    F22=F21+1
    F23=F22+CNI+1
    F24=F23-1
    CALL CAL_FIJ_RECT(F11,F12,F13,F14,F21,F22,F23,F24)
    WRITE(2,235) '4(',F12,' ',F13,' ',F23,' ',F22,')'
    FAD1=F_LEFT_RECT(JP,KP)
    WRITE(3,231) FAD1
    OWN=CELL_ID_RECT(IP,JP,KP)
    WRITE(4,231) OWN
    ENDDO   
    ENDDO
    
        !UP FACE--PROCESSOR_MI1_TO_MI2  !CASE RECTANGLE 2-2-2
    JP=CNL
    J=JP+CNL*(MJ1-1)
    DO KP=1,CNL
    K=KP+CNL*(MJ2-1)
    DO IP=1,CNI
    I=IP+CNI*(MI-1)
    F11=(CNL+1)*(CNI+1)*(KP-1)+(CNI+1)*(JP-1)+IP-1
    F12=F11+1
    F13=F12+CNI+1
    F14=F13-1
    F21=F11+(CNL+1)*(CNI+1)
    F22=F21+1
    F23=F22+CNI+1
    F24=F23-1
    CALL CAL_FIJ_RECT(F11,F12,F13,F14,F21,F22,F23,F24)
    WRITE(2,235) '4(',F14,' ',F24,' ',F23,' ',F13,')'
    FAD1=F_UP_RECT(IP,KP)
    WRITE(3,231) FAD1
    OWN=CELL_ID_RECT(IP,JP,KP)
    WRITE(4,231) OWN
    ENDDO   
    ENDDO
    
        !BACK FACE--PROCESSOR_MI1_TO_MI2  !CASE RECTANGLE 2-2-2
    KP=CNL
    K=KP+CNL*(MJ2-1)
    DO JP=1,CNL
    J=JP+CNL*(MJ1-1)
    DO IP=1,CNI
    I=IP+CNI*(MI-1)
    F11=(CNL+1)*(CNI+1)*(KP-1)+(CNI+1)*(JP-1)+IP-1
    F12=F11+1
    F13=F12+CNI+1
    F14=F13-1
    F21=F11+(CNL+1)*(CNI+1)
    F22=F21+1
    F23=F22+CNI+1
    F24=F23-1
    CALL CAL_FIJ_RECT(F11,F12,F13,F14,F21,F22,F23,F24)
    WRITE(2,235) '4(',F21,' ',F22,' ',F23,' ',F24,')'
    FAD1=F_BACK_RECT(IP,JP)
    WRITE(3,231) FAD1
    OWN=CELL_ID_RECT(IP,JP,KP)
    WRITE(4,231) OWN
    ENDDO   
    ENDDO    
        
    N=6
    WRITE(31,231) N+4
    WRITE(31,*) '('  
    WRITE(31,*) '0'
    WRITE(31,*) '1'
    WRITE(31,*) '2' 
    WRITE(31,*) '3' 
    DO I=1,N
    WRITE(31,*) '-1' 
    ENDDO
    WRITE(31,*) ')'
    
    WRITE(21,231) N+4
    WRITE(21,*) "("    
    WRITE(21,*) "    INLET1"
    WRITE(21,*) "    {"    
    WRITE(21,*) "        type            patch;"
    N=0
    WRITE(21,232) "        nFaces          ",N,";"
    N=F_INTERNAL_NUM
    WRITE(21,232) "        startFace       ",N,";"
    WRITE(21,*) "    }"
    WRITE(21,*) "    INLET2"
    WRITE(21,*) "    {"    
    WRITE(21,*) "        type            patch;"
    N=0
    WRITE(21,232) "        nFaces          ",N,";"
    N=F_INTERNAL_NUM
    WRITE(21,232) "        startFace       ",N,";"
    WRITE(21,*) "    }"
    WRITE(21,*) "    OUTLET"
    WRITE(21,*) "    {"
    WRITE(21,*) "        type            patch;"
    N=0
    WRITE(21,232) "        nFaces          ",N,";"
    N=F_INTERNAL_NUM
    WRITE(21,232) "        startFace       ",N,";"
    WRITE(21,*) "    }"
    WRITE(21,*) "    HULL"
    WRITE(21,*) "    {"
    WRITE(21,*) "        type            wall;"
    WRITE(21,*) "        inGroups        1(wall);"
    N=0
    WRITE(21,232) "        nFaces          ",N,";"
    N=F_INTERNAL_NUM
    WRITE(21,232) "        startFace       ",N,";"
    WRITE(21,*) "    }"    
        
        !FRONT FACE--PROCESSOR_MI1_TO_MI2  !CASE RECTANGLE 2-2-2
    PRO1=PNI*PNJ+PNL*PNI*(MJ2-1)+PNI*(MJ1-1)+MI-1
    PRO2=PNI*PNJ+PNL*PNI*(MJ2-2)+PNI*(MJ1-1)+MI-1
    WRITE(21,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(21,*) '    {'
    WRITE(21,*) '        type            processor;'
    WRITE(21,*) '        inGroups        List<word> 1(processor);'
    N=CNI*CNL
    WRITE(21,232) "        nFaces          ",N,";"
    N=F_INTERNAL_NUM
    WRITE(21,232) "        startFace       ",N,";"
    WRITE(21,*) '        matchTolerance  0.0001;'
    WRITE(21,232) "        myProcNo        ",PRO1,";"
    WRITE(21,232) "        neighbProcNo    ",PRO2,";"
    WRITE(21,*) '    }'
    WRITE(41,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(41,*) '    {'
    WRITE(41,*) '        type            processor;'
    WRITE(41,*) '        value           nonuniform List<scalar> '
    N=CNI*CNL
    WRITE(41,'(I0)') N
    WRITE(41,*) '('
    DO JP=1,CNL
    DO IP=1,CNI
    WRITE(41,'(E21.14)') P_FRONT_RECT(IP,JP)    
    ENDDO
    ENDDO
    WRITE(41,*) ')'
    WRITE(41,*) ';'
    WRITE(41,*) '    }'
    WRITE(51,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(51,*) '    {'
    WRITE(51,*) '        type            processor;'
    WRITE(51,*) '        value           nonuniform List<vector> '
    N=CNI*CNL
    WRITE(51,'(I0)') N
    WRITE(51,*) '('
    DO JP=1,CNL
    DO IP=1,CNI
    UU1=U_FRONT_RECT(1,IP,JP)
    UU2=U_FRONT_RECT(2,IP,JP)
    UU3=U_FRONT_RECT(3,IP,JP)
    WRITE(51,230) '( ',UU1,' ',UU2,' ',UU3,' )'  
    ENDDO
    ENDDO
    WRITE(51,*) ')'
    WRITE(51,*) ';'
    WRITE(51,*) '    }'
    WRITE(61,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(61,*) '    {'
    WRITE(61,*) '        type            processor;'
    WRITE(61,*) "        value           uniform 0;"
    WRITE(61,*) '    }'
    
        !DOWN FACE--PROCESSOR_MI1_TO_MI2  !CASE RECTANGLE 2-2-2
    PRO2=PNI*PNJ+PNL*PNI*(MJ2-1)+PNI*(MJ1-2)+MI-1
    WRITE(21,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(21,*) '    {'
    WRITE(21,*) '        type            processor;'
    WRITE(21,*) '        inGroups        List<word> 1(processor);'
    N=CNI*CNL
    WRITE(21,232) "        nFaces          ",N,";"
    N=F_INTERNAL_NUM+CNI*CNL
    WRITE(21,232) "        startFace       ",N,";"
    WRITE(21,*) '        matchTolerance  0.0001;'
    WRITE(21,232) "        myProcNo        ",PRO1,";"
    WRITE(21,232) "        neighbProcNo    ",PRO2,";"
    WRITE(21,*) '    }'
    WRITE(41,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(41,*) '    {'
    WRITE(41,*) '        type            processor;'
    WRITE(41,*) '        value           nonuniform List<scalar> '
    N=CNI*CNL
    WRITE(41,'(I0)') N
    WRITE(41,*) '('
    DO KP=1,CNL
    DO IP=1,CNI
    WRITE(41,'(E21.14)') P_DOWN_RECT(IP,KP)    
    ENDDO
    ENDDO
    WRITE(41,*) ')'
    WRITE(41,*) ';'
    WRITE(41,*) '    }'
    WRITE(51,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(51,*) '    {'
    WRITE(51,*) '        type            processor;'
    WRITE(51,*) '        value           nonuniform List<vector> '
    N=CNI*CNL
    WRITE(51,'(I0)') N
    WRITE(51,*) '('
    DO KP=1,CNL
    DO IP=1,CNI
    UU1=U_DOWN_RECT(1,IP,KP)
    UU2=U_DOWN_RECT(2,IP,KP)
    UU3=U_DOWN_RECT(3,IP,KP)
    WRITE(51,230) '( ',UU1,' ',UU2,' ',UU3,' )'  
    ENDDO
    ENDDO
    WRITE(51,*) ')'
    WRITE(51,*) ';'
    WRITE(51,*) '    }'
    WRITE(61,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(61,*) '    {'
    WRITE(61,*) '        type            processor;'
    WRITE(61,*) "        value           uniform 0;"
    WRITE(61,*) '    }'
    
        !RIGHT FACE--PROCESSOR_MI1_TO_MI2  !CASE RECTANGLE 2-2-2
    PRO2=PNI*PNJ+PNL*PNI*(MJ2-1)+PNI*(MJ1-1)+MI-2
    WRITE(21,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(21,*) '    {'
    WRITE(21,*) '        type            processor;'
    WRITE(21,*) '        inGroups        List<word> 1(processor);'
    N=CNL*CNL
    WRITE(21,232) "        nFaces          ",N,";"
    N=F_INTERNAL_NUM+CNI*CNL+CNI*CNL
    WRITE(21,232) "        startFace       ",N,";"
    WRITE(21,*) '        matchTolerance  0.0001;'
    WRITE(21,232) "        myProcNo        ",PRO1,";"
    WRITE(21,232) "        neighbProcNo    ",PRO2,";"
    WRITE(21,*) '    }'
    WRITE(41,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(41,*) '    {'
    WRITE(41,*) '        type            processor;'
    WRITE(41,*) '        value           nonuniform List<scalar> '
    N=CNL*CNL
    WRITE(41,'(I0)') N
    WRITE(41,*) '('
    DO KP=1,CNL
    DO JP=1,CNL
    WRITE(41,'(E21.14)') P_RIGHT_RECT(JP,KP)    
    ENDDO
    ENDDO
    WRITE(41,*) ')'
    WRITE(41,*) ';'
    WRITE(41,*) '    }'
    WRITE(51,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(51,*) '    {'
    WRITE(51,*) '        type            processor;'
    WRITE(51,*) '        value           nonuniform List<vector> '
    N=CNL*CNL
    WRITE(51,'(I0)') N
    WRITE(51,*) '('
    DO KP=1,CNL
    DO JP=1,CNL
    UU1=U_RIGHT_RECT(1,JP,KP)
    UU2=U_RIGHT_RECT(2,JP,KP)
    UU3=U_RIGHT_RECT(3,JP,KP)
    WRITE(51,230) '( ',UU1,' ',UU2,' ',UU3,' )'  
    ENDDO
    ENDDO
    WRITE(51,*) ')'
    WRITE(51,*) ';'
    WRITE(51,*) '    }'
    WRITE(61,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(61,*) '    {'
    WRITE(61,*) '        type            processor;'
    WRITE(61,*) "        value           uniform 0;"
    WRITE(61,*) '    }' 
    
        !LEFT FACE--PROCESSOR_MI1_TO_MI2  !CASE RECTANGLE 2-2-2
    PRO2=PNI*PNJ+PNL*PNI*(MJ2-1)+PNI*(MJ1-1)+MI
    WRITE(21,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(21,*) '    {'
    WRITE(21,*) '        type            processor;'
    WRITE(21,*) '        inGroups        List<word> 1(processor);'
    N=CNL*CNL
    WRITE(21,232) "        nFaces          ",N,";"
    N=F_INTERNAL_NUM+CNI*CNL+CNI*CNL+CNL*CNL
    WRITE(21,232) "        startFace       ",N,";"
    WRITE(21,*) '        matchTolerance  0.0001;'
    WRITE(21,232) "        myProcNo        ",PRO1,";"
    WRITE(21,232) "        neighbProcNo    ",PRO2,";"
    WRITE(21,*) '    }'
    WRITE(41,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(41,*) '    {'
    WRITE(41,*) '        type            processor;'
    WRITE(41,*) '        value           nonuniform List<scalar> '
    N=CNL*CNL
    WRITE(41,'(I0)') N
    WRITE(41,*) '('
    DO KP=1,CNL
    DO JP=1,CNL
    WRITE(41,'(E21.14)') P_LEFT_RECT(JP,KP)    
    ENDDO
    ENDDO
    WRITE(41,*) ')'
    WRITE(41,*) ';'
    WRITE(41,*) '    }'
    WRITE(51,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(51,*) '    {'
    WRITE(51,*) '        type            processor;'
    WRITE(51,*) '        value           nonuniform List<vector> '
    N=CNL*CNL
    WRITE(51,'(I0)') N
    WRITE(51,*) '('
    DO KP=1,CNL
    DO JP=1,CNL
    UU1=U_LEFT_RECT(1,JP,KP)
    UU2=U_LEFT_RECT(2,JP,KP)
    UU3=U_LEFT_RECT(3,JP,KP)
    WRITE(51,230) '( ',UU1,' ',UU2,' ',UU3,' )'  
    ENDDO
    ENDDO
    WRITE(51,*) ')'
    WRITE(51,*) ';'
    WRITE(51,*) '    }'
    WRITE(61,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(61,*) '    {'
    WRITE(61,*) '        type            processor;'
    WRITE(61,*) "        value           uniform 0;"
    WRITE(61,*) '    }' 
    
        !UP FACE--PROCESSOR_MI1_TO_MI2  !CASE RECTANGLE 2-2-2
    PRO2=PNI*PNJ+PNL*PNI*(MJ2-1)+PNI*MJ1+MI-1
    WRITE(21,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(21,*) '    {'
    WRITE(21,*) '        type            processor;'
    WRITE(21,*) '        inGroups        List<word> 1(processor);'
    N=CNI*CNL
    WRITE(21,232) "        nFaces          ",N,";"
    N=F_INTERNAL_NUM+CNI*CNL+CNI*CNL+CNL*CNL+CNL*CNL
    WRITE(21,232) "        startFace       ",N,";"
    WRITE(21,*) '        matchTolerance  0.0001;'
    WRITE(21,232) "        myProcNo        ",PRO1,";"
    WRITE(21,232) "        neighbProcNo    ",PRO2,";"
    WRITE(21,*) '    }'
    WRITE(41,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(41,*) '    {'
    WRITE(41,*) '        type            processor;'
    WRITE(41,*) '        value           nonuniform List<scalar> '
    N=CNI*CNL
    WRITE(41,'(I0)') N
    WRITE(41,*) '('
    DO KP=1,CNL
    DO IP=1,CNI
    WRITE(41,'(E21.14)') P_UP_RECT(IP,KP)    
    ENDDO
    ENDDO
    WRITE(41,*) ')'
    WRITE(41,*) ';'
    WRITE(41,*) '    }'
    WRITE(51,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(51,*) '    {'
    WRITE(51,*) '        type            processor;'
    WRITE(51,*) '        value           nonuniform List<vector> '
    N=CNI*CNL
    WRITE(51,'(I0)') N
    WRITE(51,*) '('
    DO KP=1,CNL
    DO IP=1,CNI
    UU1=U_UP_RECT(1,IP,KP)
    UU2=U_UP_RECT(2,IP,KP)
    UU3=U_UP_RECT(3,IP,KP)
    WRITE(51,230) '( ',UU1,' ',UU2,' ',UU3,' )'  
    ENDDO
    ENDDO
    WRITE(51,*) ')'
    WRITE(51,*) ';'
    WRITE(51,*) '    }'
    WRITE(61,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(61,*) '    {'
    WRITE(61,*) '        type            processor;'
    WRITE(61,*) "        value           uniform 0;"
    WRITE(61,*) '    }'
    
        !BACK FACE--PROCESSOR_MI1_TO_MI2  !CASE RECTANGLE 2-2-2
    PRO2=PNI*PNJ+PNL*PNI*MJ2+PNI*(MJ1-1)+MI-1
    WRITE(21,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(21,*) '    {'
    WRITE(21,*) '        type            processor;'
    WRITE(21,*) '        inGroups        List<word> 1(processor);'
    N=CNI*CNL
    WRITE(21,232) "        nFaces          ",N,";"
    N=F_INTERNAL_NUM+CNI*CNL+CNI*CNL+CNL*CNL+CNL*CNL+CNI*CNL
    WRITE(21,232) "        startFace       ",N,";"
    WRITE(21,*) '        matchTolerance  0.0001;'
    WRITE(21,232) "        myProcNo        ",PRO1,";"
    WRITE(21,232) "        neighbProcNo    ",PRO2,";"
    WRITE(21,*) '    }'
    WRITE(41,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(41,*) '    {'
    WRITE(41,*) '        type            processor;'
    WRITE(41,*) '        value           nonuniform List<scalar> '
    N=CNI*CNL
    WRITE(41,'(I0)') N
    WRITE(41,*) '('
    DO JP=1,CNL
    DO IP=1,CNI
        WRITE(41,'(E21.14)') P_BACK_RECT(IP,JP)    
    ENDDO
    ENDDO
    WRITE(41,*) ')'
    WRITE(41,*) ';'
    WRITE(41,*) '    }'
    WRITE(51,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(51,*) '    {'
    WRITE(51,*) '        type            processor;'
    WRITE(51,*) '        value           nonuniform List<vector> '
    N=CNI*CNL
    WRITE(51,'(I0)') N
    WRITE(51,*) '('
    DO JP=1,CNL
    DO IP=1,CNI
    UU1=U_BACK_RECT(1,IP,JP)
    UU2=U_BACK_RECT(2,IP,JP)
    UU3=U_BACK_RECT(3,IP,JP)
    WRITE(51,230) '( ',UU1,' ',UU2,' ',UU3,' )'  
    ENDDO
    ENDDO
    WRITE(51,*) ')'
    WRITE(51,*) ';'
    WRITE(51,*) '    }'
    WRITE(61,236) 'procBoundary',PRO1,'to',PRO2
    WRITE(61,*) '    {'
    WRITE(61,*) '        type            processor;'
    WRITE(61,*) "        value           uniform 0;"
    WRITE(61,*) '    }'        
        
        
        
        
        
    ENDIF 
    ENDIF              
    ENDIF
        
    ENDIF

    If ( iMPI_MyID.EQ.0 ) THEN
    write(*,*) 'test 05'
    ENDIF

    WRITE(2,*) ")"
    WRITE(2,*) " "
    WRITE(2,*) " "
    WRITE(2,*) "// ************************************************************************* //"
    CLOSE(2)
    WRITE(3,*) ")"
    WRITE(3,*) " "
    WRITE(3,*) " "
    WRITE(3,*) "// ************************************************************************* //"
    CLOSE(3)
    WRITE(4,*) ")"
    WRITE(4,*) " "
    WRITE(4,*) " "
    WRITE(4,*) "// ************************************************************************* //"
    CLOSE(4)
    WRITE(5,*) ")"
    WRITE(5,*) " "
    WRITE(5,*) " "
    WRITE(5,*) "// ************************************************************************* //"
    CLOSE(5)
    WRITE(21,*) ")"
    WRITE(21,*) " "
    WRITE(21,*) "// ************************************************************************* //"
    CLOSE(21)
    WRITE(31,*) " "
    WRITE(31,*) "// ************************************************************************* //"
    CLOSE(31)  
    WRITE(41,*) "}"
    WRITE(41,*) " "
    WRITE(41,*) "// ************************************************************************* //"
    CLOSE(41)
    WRITE(51,*) "}"
    WRITE(51,*) " "
    WRITE(51,*) "// ************************************************************************* //"
    CLOSE(51)
    WRITE(61,*) '}'
    WRITE(61,*) ' '
    WRITE(61,*) ' '
    WRITE(61,*) '// ************************************************************************* //'
    CLOSE(61)
    
   CALL MPI_Barrier( MPI_COMM_WORLD, iMPI_ErrorInfo )
230 FORMAT(A2,E21.14,A1,E21.14,A1,E21.14,A2)     
231 FORMAT(I0) 
232 FORMAT(A24,I0,A1) 
233 FORMAT(A26,I0,A9,I0,A9,I0,A17,I0,A2)   
234 FORMAT(A2,I0,A1,I0,A1,I0,A1) 
235 FORMAT(A2,I0,A1,I0,A1,I0,A1,I0,A1)
236 FORMAT(A12,I0,A2,I0)  
237 FORMAT(A26,E21.14,A1,E21.14,A1,E21.14,A3) 
238 FORMAT(A34,E21.14,A1,E21.14,A1,E21.14,A3)
239 FORMAT(A1,I0,A1)
!-------------------------------------------------------------------------------
!    local variable
!-------------------------------------------------------------------------------
    
!-------------------------------------------------------------------------
!    execute statement
!------------------------------------------------------------------------- 
	

!-------------------------------------------------

    END SUBROUTINE