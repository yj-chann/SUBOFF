!--------------------------------------------------------------
!     GRID GENERATION IN CHANNEL FLOW
!     X2 Clustering only
!     JICHOI
!--------------------------------------------------------------

      PROGRAM MAIN
      IMPLICIT NONE
      INTEGER N1,N2,N3
      INTEGER PN2,PN3,PRON
      INTEGER NPR,NPC
      INTEGER CN2,CN3
      PARAMETER (N1=3,N2=160,N3=160)
      PARAMETER (   NPR=8 ,NPC=8 ) ! NPP=64
      
      INTEGER I,J,K,NALL

    !   WRITE(*,*) INT(-4.3),INT(4.5),INT(4.7)

     
      NALL=N1*N2*N3
      CN2=N2/NPR
      CN3=N3/NPC
      OPEN(1,FILE='../manualDecomposition',STATUS='UNKNOWN')
      WRITE(1,*) "/*--------------------------------*- C++ -*----------------------------------*\"
      WRITE(1,*) "  =========                 |"
      WRITE(1,*) "  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox"
      WRITE(1,*) "   \\    /   O peration     | Website:  https://openfoam.org"
      WRITE(1,*) "    \\  /    A nd           | Version:  v2312"
      WRITE(1,*) "     \\/     M anipulation  |"
      WRITE(1,*) "\*---------------------------------------------------------------------------*/"
      WRITE(1,*) "FoamFile"
      WRITE(1,*) "{"
      WRITE(1,*) "    version     2.0;"
      WRITE(1,*) "    format      ascii;"
      WRITE(1,*) "    class       labelList;"
      WRITE(1,*) '    location    "constant";'
      WRITE(1,*) "    object      cellDecomposition;"
      WRITE(1,*) "}"
      WRITE(1,*) "// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //"
      WRITE(1,"(I0)") NALL
      WRITE(1,*) "("
      DO K=1,N3 !y dir
          PN3=INT((K-1)/CN3)
      DO J=1,N2 ! z dir 从(YMAX,ZMIN)->(YMAX->ZMAX) (YMAX-1,ZMIN)->(YMAX-1->ZMAX)  -> ...
          PN2=INT((J-1)/CN2)
          PRON=NPR*PN3+PN2
      DO I=1,N1
      WRITE(1,"(I0)") PRON
      ENDDO
      ENDDO
      ENDDO
      WRITE(1,*) ")"
      CLOSE(1)

        
      END
