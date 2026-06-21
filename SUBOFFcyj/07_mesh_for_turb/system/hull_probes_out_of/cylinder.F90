
    PROGRAM cylinder_main
    IMPLICIT NONE
    INCLUDE 'parameter.h'

    REAL*8 INLET_X(N2,N3),INLET_Y(N2,N3),INLET_Z(N2,N3)
    REAL*8 PAI,DX
    REAL*8 :: Y(0:NI)
    REAL*8 X
    INTEGER*8 I,J,K,I2,NUM
    CHARACTER*20 CHAR
 
    

    OPEN(1,FILE='ReadData/Y.txt',STATUS='UNKNOWN')
    DO I=0,NI
       READ(1,*) Y(I)
    ENDDO
    CLOSE(1)

    DX=(Y(NI)-Y(NI-1))/2
    WRITE(*,'(A,E22.15)') 'DX=', DX

    OPEN(1,FILE='ReadData/INLET1.plt',STATUS='OLD')   
    DO K=1,N3
    DO J=1,N2  
        READ(1,*) INLET_X(J,K),INLET_Y(J,K),INLET_Z(J,K)
    ENDDO
    ENDDO  
    CLOSE(1)


    OPEN(2,FILE='../probesInletShift.c',STATUS='UNKNOWN')
    DO K=1,N3
        DO J=1,N2
            X=INLET_X(J,K)+DX
            WRITE(2,230) '( ',X,' ',INLET_Y(J,K),' ',INLET_Z(J,K),' )'
        ENDDO
        ENDDO  
    CLOSE(2)
    

    OPEN(2,FILE='../probesInletShiftCenterLine.c',STATUS='UNKNOWN')
    DO K=1,N3
        DO J=N2/2,N2/2
            X=INLET_X(J,K)+DX
            WRITE(2,230) '( ',X,' ',INLET_Y(J,K),' ',INLET_Z(J,K),' )'
        ENDDO
        ENDDO  
    CLOSE(2)


    ! OPEN(11,FILE='../controlDict',STATUS='UNKNOWN')
    ! WRITE(11,*) "/*--------------------------------*- C++ -*----------------------------------*\"
    ! WRITE(11,*) "  =========                 |"
    ! WRITE(11,*) "  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox"
    ! WRITE(11,*) "   \\    /   O peration     | Website:  https://openfoam.org"
    ! WRITE(11,*) "    \\  /    A nd           | Version:  v2312"
    ! WRITE(11,*) "     \\/     M anipulation  |"
    ! WRITE(11,*) "\*---------------------------------------------------------------------------*/"
    ! WRITE(11,*) "FoamFile"
    ! WRITE(11,*) "{"
    ! WRITE(11,*) "    version     2.0;"
    ! WRITE(11,*) "    format      ascii;"
    ! WRITE(11,*) "    class       dictionary;"
    ! WRITE(11,*) "    object      controlDict;"
    ! WRITE(11,*) "}"
    ! WRITE(11,*) "// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //"
    ! WRITE(11,*) "application     pimpleFoam;"
    ! WRITE(11,*) "startFrom       startTime;"
    ! WRITE(11,*) "startTime       0;"
    ! WRITE(11,*) "stopAt          endTime;"
    ! WRITE(11,*) "endTime         1.0200001;"
    ! WRITE(11,*) "deltaT          2e-4;"

    ! WRITE(11,*) "writeControl    timeStep;"
    ! WRITE(11,*) "writeInterval   1000;"
    ! WRITE(11,*) "purgeWrite      0;"
    ! WRITE(11,*) "writeFormat     ascii;"
    ! WRITE(11,*) "writePrecision  14;"
    ! WRITE(11,*) "writeCompression off;"
    ! WRITE(11,*) "timeFormat      general;"
    ! WRITE(11,*) "timePrecision   6;"
    ! WRITE(11,*) "runTimeModifiable true;"
    ! WRITE(11,*) "functions"
    ! WRITE(11,*) "{"
    ! !WRITE(11,*) '    #includeFunc wallShearStress'
    ! !WRITE(11,*) '    #includeFunc yPlus'
    ! WRITE(11,*) '    #includeFunc Q'
    ! WRITE(11,*) '    #include "FOcomponents"'    
    ! WRITE(11,*) '    #include "mass_flow_cal"'    
    ! WRITE(11,*) "    probes1_single"
    ! WRITE(11,*) "    {"
    ! WRITE(11,*) "        type            probes;"
    ! WRITE(11,*) '        libs            ("libsampling.so");'
    ! WRITE(11,*) "        writeControl    timeStep;"
    ! WRITE(11,*) "        writeInterval   1;"
    ! WRITE(11,*) "        timeStart      0;"
    ! WRITE(11,*) "        timeEnd        1000000;"
    ! WRITE(11,*) "        fields"
    ! WRITE(11,*) "        ("
    ! WRITE(11,*) "            p  Ux Uy Uz divU"
    ! WRITE(11,*) "         );"
    ! WRITE(11,*) "        probeLocations"
    ! WRITE(11,*) "        ("
    ! DO K=N3/2,N3/2
    ! DO J=N2/2,N2/2
    !     X=INLET_X(J,K)+DX
    !     WRITE(11,230) '( ',X,' ',INLET_Y(J,K),' ',INLET_Z(J,K),' )'
    ! ENDDO
    ! ENDDO 
    ! WRITE(11,*) "         );"
    ! WRITE(11,*) "    }"
    
    ! WRITE(11,*) "    probes1"
    ! WRITE(11,*) "    {"
    ! WRITE(11,*) "        type            probes;"
    ! WRITE(11,*) '        libs            ("libsampling.so");'
    ! WRITE(11,*) "        writeControl    timeStep;"
    ! WRITE(11,*) "        writeInterval   1;"
    ! WRITE(11,*) "        timeStart      0;"
    ! WRITE(11,*) "        timeEnd        1000000;"
    ! WRITE(11,*) "        fields"
    ! WRITE(11,*) "        ("
    ! WRITE(11,*) "            p Ux Uy Uz"
    ! WRITE(11,*) "         );"
    ! WRITE(11,*) "        probeLocations"
    ! WRITE(11,*) "        ("     
    ! WRITE(11,*) "        #include     ""probesInletShift.c"""
    ! WRITE(11,*) "         );"
    ! WRITE(11,*) "    }"
    
    ! WRITE(11,*) "    probes2"
    ! WRITE(11,*) "    {"
    ! WRITE(11,*) "        type            probes;"
    ! WRITE(11,*) '        libs            ("libsampling.so");'
    ! WRITE(11,*) "        writeControl    timeStep;"
    ! WRITE(11,*) "        writeInterval   1;"
    ! WRITE(11,*) "        timeStart      0;"
    ! WRITE(11,*) "        timeEnd        1000000;"
    ! WRITE(11,*) "        fields"
    ! WRITE(11,*) "        ("
    ! WRITE(11,*) "            p Ux Uy Uz"
    ! WRITE(11,*) "         );"
    ! WRITE(11,*) "        probeLocations"
    ! WRITE(11,*) "        ("     
    ! WRITE(11,*) "        #include     ""probesInletShiftCenterLine.c"""
    ! WRITE(11,*) "         );"
    ! WRITE(11,*) "    }"
    ! WRITE(11,*) "}"
    ! WRITE(11,*) " "
    ! WRITE(11,*) "// ************************************************************************* //"
    ! CLOSE(11)
        
    
       
    
230 FORMAT(A2,E22.15,A1,E22.15,A1,E22.15,A2)         
231 FORMAT(I0)
232 FORMAT(A10,I3.3) 
233 FORMAT(A13,I3.3)
235 FORMAT(A14,I3.3)

    END PROGRAM


