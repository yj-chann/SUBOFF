!--------------------------------------------------------------
!     GRID GENERATION IN CHANNEL FLOW
!     X2 Clustering only
!     JICHOI
!       生成点文件和进行流场初始化
!--------------------------------------------------------------

      SUBROUTINE GET_POINTS_INI_UP(CY_X,CY_Y,CY_Z,RECT_X1,RECT_Y0,RECT_Z0,Y,COUNT1)
      INCLUDE 'head.fi'
      REAL XIELV
      REAL RECT_Y0(NL+1,NL+1),RECT_Z0(NL+1,NL+1),RECT_X1(NL+1,NL+1)
      REAL Y(0:NI)
      REAL CY_X(NK,NJ+1),CY_Y(NK,NJ+1),CY_Z(NK,NJ+1)
      REAL XP,YP,ZP
      INTEGER I,J,K
      INTEGER MI,MJ,MJ1,MJ2,PN
      INTEGER IP,JP,KP
      CHARACTER*20 CHAR1
      
      REAL UREF(2,REF_N1,REF_N2),PREF(REF_N1,REF_N2),U3
      REAL XPOINTS_CY(CNI+1,CNJ+1,NK+1),YPOINTS_CY(CNI+1,CNJ+1,NK+1),ZPOINTS_CY(CNI+1,CNJ+1,NK+1)
      REAL XCELL_CY(CNI,CNJ,NK),YCELL_CY(CNI,CNJ,NK),ZCELL_CY(CNI,CNJ,NK)
      REAL XPOINTS_RECT(CNI+1,CNL+1,CNL+1),YPOINTS_RECT(CNI+1,CNL+1,CNL+1),ZPOINTS_RECT(CNI+1,CNL+1,CNL+1)
      REAL XCELL_RECT(CNI,CNL,CNL),YCELL_RECT(CNI,CNL,CNL),ZCELL_RECT(CNI,CNL,CNL)
      REAL XLEFT_CY(CNJ,NK),YLEFT_CY(CNJ,NK),ZLEFT_CY(CNJ,NK)
      REAL XRIGHT_CY(CNJ,NK),YRIGHT_CY(CNJ,NK),ZRIGHT_CY(CNJ,NK)
      REAL XUP_CY(CNI,NK),YUP_CY(CNI,NK),ZUP_CY(CNI,NK)
      REAL XDOWN_CY(CNI,NK),YDOWN_CY(CNI,NK),ZDOWN_CY(CNI,NK)
      REAL XLEFT_RECT(CNL,CNL),YLEFT_RECT(CNL,CNL),ZLEFT_RECT(CNL,CNL)
      REAL XRIGHT_RECT(CNL,CNL),YRIGHT_RECT(CNL,CNL),ZRIGHT_RECT(CNL,CNL)
      REAL XUP_RECT(CNI,CNL),YUP_RECT(CNI,CNL),ZUP_RECT(CNI,CNL)
      REAL XDOWN_RECT(CNI,CNL),YDOWN_RECT(CNI,CNL),ZDOWN_RECT(CNI,CNL)
      REAL XFRONT_RECT(CNI,CNL),YFRONT_RECT(CNI,CNL),ZFRONT_RECT(CNI,CNL)
      REAL XBACK_RECT(CNI,CNL),YBACK_RECT(CNI,CNL),ZBACK_RECT(CNI,CNL)
      REAL XTEMP,YTEMP
      REAL P_RECT(NL,NL),U_RECT(3,NL,NL),P_CY(NJ,NK),U_CY(3,NJ,NK)
      INTEGER P_LABEL_CY0(P_NUM_CY),P_LABEL_RECT0(P_NUM_RECT)
      INTEGER P_LABEL_CY1(P_NUM_CY),P_LABEL_RECT1(P_NUM_RECT)
      INTEGER SORT_CY(P_NUM_CY),SORT_RECT(P_NUM_RECT)
      INTEGER IP1,ITEMP,JTEMP,NUM_P,JISHU,I2,J2,K2,IF_CY,I22,J22,K22
      CHARACTER*900 DIR,DIR_FOR_INLET
      CHARACTER*200 CHAR3
      INTEGER*8 COUNT1,COUNT2,COUNT3,COUNT_RATE,COUNT_MAX
  
      
      !READ P/U LOCATION
      !DIR='/es01/paratera/sce4049/zf/20221115/04_ini_2d/5000/' 
      !DIR_FOR_INLET='/es01/paratera/sce4049/zf/20221215/05_ini_2d/get_ini_p_u/'
      !DIR_FOR_INLET='/work/home/ac6w4p63bu/njy21/SUBOFF_tensity20_meshPolished/05_ini_2d/get_ini_p_u/'
      !DIR_FOR_INLET='/public1/home/scb7552/njy21/SUBOFF/tensity80_meshPolished/05_ini_2d/get_ini_p_u/'
      DIR_FOR_INLET='../05_ini_2d/get_ini_p_u/'
            
      
    !PROCESSOR NUMBER
    MI=COORD(0)+1 ! 1~4
    MJ=COORD(1)+1 ! 1~64
    IF (MJ.LE.PNJ) THEN ! 圆柱部分法向优先外缘到壁面，流向从内部到出口 PNJ=60 
        PN=PNI*(MJ-1)+MI-1
    ELSE
        MJ1=( (COORD(1)-PNJ)-MOD(COORD(1)-PNJ,PNL) )/PNL+1 ! 1~2沿Z
        MJ2=MJ-PNJ-(MJ1-1)*PNL ! 1~2 沿-Y
        PN=PNI*PNJ+PNL*PNI*(MJ2-1)+PNI*(MJ1-1)+MI-1
    ENDIF
       
      
      
   !points/pointProcAddressing
    WRITE(CHAR1,231) PN
    OPEN(2,FILE='../processor'//TRIM(ADJUSTL(CHAR1))//'/constant/polyMesh/points',STATUS='UNKNOWN')
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
    IF (MJ.LE.PNJ) THEN
    WRITE(2,231) P_NUM_CY ! 11*51*640=359,040
    ELSE
    WRITE(2,231) P_NUM_RECT ! 81*81*51=334611
    ENDIF   
    WRITE(2,*) "("
    
    IF (MJ.LE.PNJ) THEN ! CY
        DO KP=1,NK
            K=KP
        DO JP=1,CNJ+1 ! CNJ=10,CNI=50
            J=JP+(MJ-1)*CNJ
        DO IP=0,CNI
            I=NI-( IP+(MI-1)*CNI )                                            
            CALL GET_K(CY_X(K,J),XIELV)                   
            CALL GET_P(CY_X(K,J),CY_Y(K,J),CY_Z(K,J),XIELV,Y(I),XP,YP,ZP)

          XP=XP
          YP=YP
          ZP=ZP
          XPOINTS_CY(IP+1,JP,KP)=XP
          YPOINTS_CY(IP+1,JP,KP)=YP
          ZPOINTS_CY(IP+1,JP,KP)=ZP
          !WRITE(2,230) '( ',XP,' ',YP,' ',ZP,' )' 
      ENDDO
      ENDDO
      ENDDO

    !write(*,*) 'final j=:',j
    !write(*,*) 'max xpoints_cy:', maxval(XPOINTS_CY)
    !write(*,*) 'size of xpoints_cy in x-dirc:', size(XPOINTS_CY,2)-1
    !write(*,*) 'max cy_x:', maxval(CY_X)

      XPOINTS_CY(:,:,NK+1)=XPOINTS_CY(:,:,1)
      YPOINTS_CY(:,:,NK+1)=YPOINTS_CY(:,:,1)
      ZPOINTS_CY(:,:,NK+1)=ZPOINTS_CY(:,:,1)
      DO KP=1,NK
      DO JP=1,CNJ
      DO IP=1,CNI
          XCELL_CY(IP,JP,KP)=0.125*(  XPOINTS_CY(IP,JP,KP)    + XPOINTS_CY(IP,JP,KP+1) &
                                    + XPOINTS_CY(IP,JP+1,KP)  + XPOINTS_CY(IP,JP+1,KP+1) &
                                    + XPOINTS_CY(IP+1,JP,KP)  + XPOINTS_CY(IP+1,JP,KP+1) &
                                    + XPOINTS_CY(IP+1,JP+1,KP)+ XPOINTS_CY(IP+1,JP+1,KP+1) )
          YCELL_CY(IP,JP,KP)=0.125*(  YPOINTS_CY(IP,JP,KP)    + YPOINTS_CY(IP,JP,KP+1) &
                                    + YPOINTS_CY(IP,JP+1,KP)  + YPOINTS_CY(IP,JP+1,KP+1) &      
                                    + YPOINTS_CY(IP+1,JP,KP)  + YPOINTS_CY(IP+1,JP,KP+1) &      
                                    + YPOINTS_CY(IP+1,JP+1,KP)+ YPOINTS_CY(IP+1,JP+1,KP+1) )      
          ZCELL_CY(IP,JP,KP)=0.125*(  ZPOINTS_CY(IP,JP,KP)    + ZPOINTS_CY(IP,JP,KP+1) &        
                                    + ZPOINTS_CY(IP,JP+1,KP)  + ZPOINTS_CY(IP,JP+1,KP+1) &      
                                    + ZPOINTS_CY(IP+1,JP,KP)  + ZPOINTS_CY(IP+1,JP,KP+1) &      
                                    + ZPOINTS_CY(IP+1,JP+1,KP)+ ZPOINTS_CY(IP+1,JP+1,KP+1) )
      ENDDO
      ENDDO
      ENDDO
    ELSE
      DO JP=1,CNL+1 ! CNL=80 沿负Y
          J=NL+2-(JP+(MJ2-1)*CNL)
      DO KP=1,CNL+1 ! 沿正Z
          K=KP+(MJ1-1)*CNL
          IF (J.EQ.(NL/2+1) .AND. K.EQ.(NL/2+1)) THEN
              DO IP=0,CNI
              I=NI-( IP+(MI-1)*CNI )    
              XP=-Y(I)
              YP=0.0
              ZP=0.0
              XP=XP
              YP=YP
              ZP=ZP
              XPOINTS_RECT(IP+1,KP,JP)=XP
              YPOINTS_RECT(IP+1,KP,JP)=YP
              ZPOINTS_RECT(IP+1,KP,JP)=ZP
              ENDDO
          ELSE
              DO IP=0,CNI
              I=NI-( IP+(MI-1)*CNI ) 
              CALL GET_K(RECT_X1(J,K),XIELV)
              CALL GET_P(RECT_X1(J,K),RECT_Y0(J,K),RECT_Z0(J,K),XIELV,Y(I),XP,YP,ZP)
              XP=XP
              YP=YP
              ZP=ZP
              XPOINTS_RECT(IP+1,KP,JP)=XP
              YPOINTS_RECT(IP+1,KP,JP)=YP
              ZPOINTS_RECT(IP+1,KP,JP)=ZP
              ENDDO
          ENDIF
      ENDDO
      ENDDO
      
      DO KP=1,CNL
      DO JP=1,CNL
      DO IP=1,CNI
          XCELL_RECT(IP,JP,KP)=0.125*( XPOINTS_RECT(IP,JP,KP)+XPOINTS_RECT(IP,JP,KP+1) &
                                      +XPOINTS_RECT(IP,JP+1,KP)+XPOINTS_RECT(IP,JP+1,KP+1) &
                                      +XPOINTS_RECT(IP+1,JP,KP)+XPOINTS_RECT(IP+1,JP,KP+1) &
                                      +XPOINTS_RECT(IP+1,JP+1,KP)+XPOINTS_RECT(IP+1,JP+1,KP+1) )
          YCELL_RECT(IP,JP,KP)=0.125*( YPOINTS_RECT(IP,JP,KP)+YPOINTS_RECT(IP,JP,KP+1) &
                                      +YPOINTS_RECT(IP,JP+1,KP)+YPOINTS_RECT(IP,JP+1,KP+1) &
                                      +YPOINTS_RECT(IP+1,JP,KP)+YPOINTS_RECT(IP+1,JP,KP+1) &
                                      +YPOINTS_RECT(IP+1,JP+1,KP)+YPOINTS_RECT(IP+1,JP+1,KP+1) )
          ZCELL_RECT(IP,JP,KP)=0.125*( ZPOINTS_RECT(IP,JP,KP)+ZPOINTS_RECT(IP,JP,KP+1) &
                                      +ZPOINTS_RECT(IP,JP+1,KP)+ZPOINTS_RECT(IP,JP+1,KP+1) &
                                      +ZPOINTS_RECT(IP+1,JP,KP)+ZPOINTS_RECT(IP+1,JP,KP+1) &
                                      +ZPOINTS_RECT(IP+1,JP+1,KP)+ZPOINTS_RECT(IP+1,JP+1,KP+1) )
      ENDDO
      ENDDO
      ENDDO
      
    ENDIF
    
    !call mpi_bcast( XPOINTS_CY,(CNI+1)*(CNJ+1)*(NK+1), mpi_real8, 0, MPI_COMM_WORLD, iMPI_ErrorInfo )
    
    
    
    
    
    
    
    WRITE(CHAR1,231) PN
    OPEN(3,FILE='../processor'//TRIM(ADJUSTL(CHAR1))//'/constant/polyMesh/pointProcAddressing',STATUS='UNKNOWN')
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
    WRITE(3,*) "    object      pointProcAddressing;"
    WRITE(3,*) "}"
    WRITE(3,*) "// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //"
    WRITE(3,*) " "
    WRITE(3,*) " "
    IF (MJ.LE.PNJ) THEN
    WRITE(3,231) P_NUM_CY
    ELSE
    WRITE(3,231) P_NUM_RECT
    ENDIF 
    WRITE(3,*) "("
    
    
    JISHU = 0
    IF (MJ .LE. PNJ) THEN
            
      DO KP = 1, NK
        K = KP
        DO JP = 1, CNJ + 1
          J = JP + (MJ - 1) * CNJ
          DO IP = 1, CNI + 1
            I = IP + (MI - 1) * CNI
            ! I,J,K 是 CY 部分的全局结构化索引
            NUM_P = (NI + 1) * (NJ + 1) * (K - 1) + (NI + 1) * (J - 1) + I - 1
            JISHU = JISHU + 1
            P_LABEL_CY0(JISHU) = NUM_P   ! 原始全局编号 
          ENDDO
        ENDDO
      ENDDO
        
    ELSE   ! MJ .GT. PNJ
            
      ! 事实上只用到了 1-1 1-3 3-1 3-3 的case
      IF (MJ1 .EQ. 1) THEN
        IF (MJ2 .EQ. 1) THEN
          ! CASE 1-1
          !  k YMAX,ZMIN->0
          DO KP = 1, CNL + 1
            K = KP
            DO IP = 1, CNI + 1
              I = IP + (MI - 1) * CNI
              NUM_P = (NI + 1) * (NJ + 1) * (K - 1) + I - 1
              JISHU = JISHU + 1
              P_LABEL_RECT0(JISHU) = NUM_P
            ENDDO
          ENDDO
    
        DO KP = 2, CNL + 1
            K = KP + (MJ2 - 1) * CNL ! K沿负Y
            DO IP = 1, CNI + 1
              I = IP + (MI - 1) * CNI
              NUM_P = (NI + 1) * (NJ + 1) * (NK - K + 1) + I - 1
              JISHU = JISHU + 1
              P_LABEL_RECT0(JISHU) = NUM_P
            ENDDO

            DO JP = 2, CNL + 1
              J = JP + (MJ1 - 1) * CNL ! J沿正Z
              DO IP = 1, CNI + 1
                I = IP + (MI - 1) * CNI
                NUM_P = P_NUM_CY_ALL + (NI + 1) * (NL - 1) * (K - 2) + (NI + 1) * (J - 2) + I - 1
                JISHU = JISHU + 1
                P_LABEL_RECT0(JISHU) = NUM_P
              ENDDO
            ENDDO
        ENDDO
    
        ELSEIF (MJ2 .EQ. PNL) THEN
          ! CASE 1-3
          DO KP = 1, CNL
            K = KP + (MJ2 - 1) * CNL
            DO IP = 1, CNI + 1
              I = IP + (MI - 1) * CNI
              NUM_P = (NI + 1) * (NJ + 1) * (NK - K + 1) + I - 1
              JISHU = JISHU + 1
              P_LABEL_RECT0(JISHU) = NUM_P
              ! WRITE(3,231) NUM_P
            ENDDO
            DO JP = 2, CNL + 1
              J = JP + (MJ1 - 1) * CNL
              DO IP = 1, CNI + 1
                I = IP + (MI - 1) * CNI
                NUM_P = P_NUM_CY_ALL + (NI + 1) * (NL - 1) * (K - 2) + (NI + 1) * (J - 2) + I - 1
                JISHU = JISHU + 1
                P_LABEL_RECT0(JISHU) = NUM_P
                ! WRITE(3,231) NUM_P
              ENDDO
            ENDDO
          ENDDO
    
          DO KP = 3 * NL + 1, 3 * NL - CNL + 1, -1
            K = KP
            DO IP = 1, CNI + 1
              I = IP + (MI - 1) * CNI
              NUM_P = (NI + 1) * (NJ + 1) * (K - 1) + I - 1
              JISHU = JISHU + 1
              P_LABEL_RECT0(JISHU) = NUM_P
              ! WRITE(3,231) NUM_P
            ENDDO
          ENDDO
    
        ELSE
          ! CASE 1-2
          DO KP = 1, CNL + 1
            K = KP + (MJ2 - 1) * CNL
            DO IP = 1, CNI + 1
              I = IP + (MI - 1) * CNI
              NUM_P = (NI + 1) * (NJ + 1) * (NK - K + 1) + I - 1
              JISHU = JISHU + 1
              P_LABEL_RECT0(JISHU) = NUM_P
              ! WRITE(3,231) NUM_P
            ENDDO
            DO JP = 2, CNL + 1
              J = JP + (MJ1 - 1) * CNL
              DO IP = 1, CNI + 1
                I = IP + (MI - 1) * CNI
                NUM_P = P_NUM_CY_ALL + (NI + 1) * (NL - 1) * (K - 2) + (NI + 1) * (J - 2) + I - 1
                JISHU = JISHU + 1
                P_LABEL_RECT0(JISHU) = NUM_P
                ! WRITE(3,231) NUM_P
              ENDDO
            ENDDO
          ENDDO
        ENDIF
    
      ELSEIF (MJ1 .EQ. PNL) THEN
    
        IF (MJ2 .EQ. 1) THEN
          ! CASE 3-1
          DO KP = NL + 1 - CNL, NL + 1
            K = KP
            DO IP = 1, CNI + 1
              I = IP + (MI - 1) * CNI
              NUM_P = (NI + 1) * (NJ + 1) * (K - 1) + I - 1
              JISHU = JISHU + 1
              P_LABEL_RECT0(JISHU) = NUM_P
              ! WRITE(3,231) NUM_P
            ENDDO
          ENDDO
    
          DO KP = 2, CNL + 1
            K = KP + (MJ2 - 1) * CNL
            DO JP = 1, CNL
              J = JP + (MJ1 - 1) * CNL
              DO IP = 1, CNI + 1
                I = IP + (MI - 1) * CNI
                NUM_P = P_NUM_CY_ALL + (NI + 1) * (NL - 1) * (K - 2) + (NI + 1) * (J - 2) + I - 1
                JISHU = JISHU + 1
                P_LABEL_RECT0(JISHU) = NUM_P
                ! WRITE(3,231) NUM_P
              ENDDO
            ENDDO
            DO IP = 1, CNI + 1
              I = IP + (MI - 1) * CNI
              NUM_P = (NI + 1) * (NJ + 1) * (K + NL - 1) + I - 1
              JISHU = JISHU + 1
              P_LABEL_RECT0(JISHU) = NUM_P
              ! WRITE(3,231) NUM_P
            ENDDO
          ENDDO
    
        ELSEIF (MJ2 .EQ. PNL) THEN
          ! CASE 3-3
          DO KP = 1, CNL
            K = KP + (MJ2 - 1) * CNL
            DO JP = 1, CNL
              J = JP + (MJ1 - 1) * CNL
              DO IP = 1, CNI + 1
                I = IP + (MI - 1) * CNI
                NUM_P = P_NUM_CY_ALL + (NI + 1) * (NL - 1) * (K - 2) + (NI + 1) * (J - 2) + I - 1
                JISHU = JISHU + 1
                P_LABEL_RECT0(JISHU) = NUM_P
                ! WRITE(3,231) NUM_P
              ENDDO
            ENDDO
            DO IP = 1, CNI + 1
              I = IP + (MI - 1) * CNI
              NUM_P = (NI + 1) * (NJ + 1) * (K + NL - 1) + I - 1
              JISHU = JISHU + 1
              P_LABEL_RECT0(JISHU) = NUM_P
              ! WRITE(3,231) NUM_P
            ENDDO
          ENDDO
    
          DO KP = 2 * NL + 1 + CNL, 2 * NL + 1, -1
            K = KP
            DO IP = 1, CNI + 1
              I = IP + (MI - 1) * CNI
              NUM_P = (NI + 1) * (NJ + 1) * (K - 1) + I - 1
              JISHU = JISHU + 1
              P_LABEL_RECT0(JISHU) = NUM_P
              ! WRITE(3,231) NUM_P
            ENDDO
          ENDDO
    
        ELSE
          ! CASE 3-2
          DO KP = 1, CNL + 1
            K = KP + (MJ2 - 1) * CNL
            DO JP = 1, CNL
              J = JP + (MJ1 - 1) * CNL
              DO IP = 1, CNI + 1
                I = IP + (MI - 1) * CNI
                NUM_P = P_NUM_CY_ALL + (NI + 1) * (NL - 1) * (K - 2) + (NI + 1) * (J - 2) + I - 1
                JISHU = JISHU + 1
                P_LABEL_RECT0(JISHU) = NUM_P
                ! WRITE(3,231) NUM_P
              ENDDO
            ENDDO
            DO IP = 1, CNI + 1
              I = IP + (MI - 1) * CNI
              NUM_P = (NI + 1) * (NJ + 1) * (K + NL - 1) + I - 1
              JISHU = JISHU + 1
              P_LABEL_RECT0(JISHU) = NUM_P
              ! WRITE(3,231) NUM_P
            ENDDO
          ENDDO
        ENDIF
    
      ELSE
    
        IF (MJ2 .EQ. 1) THEN
          ! CASE 2-1
          DO KP = 1 + (MJ1 - 1) * CNL, 1 + MJ1 * CNL
            K = KP
            DO IP = 1, CNI + 1
              I = IP + (MI - 1) * CNI
              NUM_P = (NI + 1) * (NJ + 1) * (K - 1) + I - 1
              JISHU = JISHU + 1
              P_LABEL_RECT0(JISHU) = NUM_P
              ! WRITE(3,231) NUM_P
            ENDDO
          ENDDO
    
          DO KP = 2, CNL + 1
            K = KP + (MJ2 - 1) * CNL
            DO JP = 1, CNL + 1
              J = JP + (MJ1 - 1) * CNL
              DO IP = 1, CNI + 1
                I = IP + (MI - 1) * CNI
                NUM_P = P_NUM_CY_ALL + (NI + 1) * (NL - 1) * (K - 2) + (NI + 1) * (J - 2) + I - 1
                JISHU = JISHU + 1
                P_LABEL_RECT0(JISHU) = NUM_P
                ! WRITE(3,231) NUM_P
              ENDDO
            ENDDO
          ENDDO
    
        ELSEIF (MJ2 .EQ. PNL) THEN
          ! CASE 2-3
          DO KP = 1, CNL
            K = KP + (MJ2 - 1) * CNL
            DO JP = 1, CNL + 1
              J = JP + (MJ1 - 1) * CNL
              DO IP = 1, CNI + 1
                I = IP + (MI - 1) * CNI
                NUM_P = P_NUM_CY_ALL + (NI + 1) * (NL - 1) * (K - 2) + (NI + 1) * (J - 2) + I - 1
                JISHU = JISHU + 1
                P_LABEL_RECT0(JISHU) = NUM_P
                ! WRITE(3,231) NUM_P
              ENDDO
            ENDDO
          ENDDO
    
          DO KP = 3 * NL + 1 - (MJ1 - 1) * CNL, 3 * NL + 1 - MJ1 * CNL, -1
            K = KP
            DO IP = 1, CNI + 1
              I = IP + (MI - 1) * CNI
              NUM_P = (NI + 1) * (NJ + 1) * (K - 1) + I - 1
              JISHU = JISHU + 1
              P_LABEL_RECT0(JISHU) = NUM_P
              ! WRITE(3,231) NUM_P
            ENDDO
          ENDDO
    
        ELSE
          ! CASE 2-2
          DO KP = 1, CNL + 1
            K = KP + (MJ2 - 1) * CNL
            DO JP = 1, CNL + 1
              J = JP + (MJ1 - 1) * CNL
              DO IP = 1, CNI + 1
                I = IP + (MI - 1) * CNI
                NUM_P = P_NUM_CY_ALL + (NI + 1) * (NL - 1) * (K - 2) + (NI + 1) * (J - 2) + I - 1
                JISHU = JISHU + 1
                P_LABEL_RECT0(JISHU) = NUM_P
                IF (JISHU .EQ. 1) THEN
                  ! WRITE(*,*) 'PN=', PN, P_LABEL_RECT0(1)
                ENDIF
                ! WRITE(3,231) NUM_P
              ENDDO
            ENDDO
          ENDDO
        ENDIF
    
      ENDIF
    
    ENDIF    
 

    
    
    
    IF (MJ.LE.PNJ) THEN
        P_LABEL_CY1=P_LABEL_CY0
        CALL BUBBLE_SORT_POINTS_CY(P_LABEL_CY1,SORT_CY)
        DO I=1,P_NUM_CY
            WRITE(3,231) P_LABEL_CY1(I)    
        ENDDO
        DO I=1,P_NUM_CY
            SORT_CY2(SORT_CY(I)-1)=I-1 
        ENDDO
    ELSE
        P_LABEL_RECT1=P_LABEL_RECT0
        CALL BUBBLE_SORT_POINTS_RECT(P_LABEL_RECT1,SORT_RECT)
        DO I=1,P_NUM_RECT
            WRITE(3,231) P_LABEL_RECT1(I)  
        ENDDO 
        DO I=1,P_NUM_RECT
            SORT_RECT2(SORT_RECT(I)-1)=I-1 
        ENDDO
    ENDIF
    
    ! 写入points
    IF (MJ.LE.PNJ) THEN
        DO I=1,P_NUM_CY
            CALL GET_INDEX_POINT(P_LABEL_CY1(I),I2,J2,K2,IF_CY)
            KP=K2
            JP=J2-(MJ-1)*CNJ
            IP=I2-(MI-1)*CNI
            WRITE(2,230) '( ',XPOINTS_CY(IP,JP,KP),' ',YPOINTS_CY(IP,JP,KP),' ',ZPOINTS_CY(IP,JP,KP),' )'
        ENDDO   
    ELSE        
        DO I=1,P_NUM_RECT
            CALL GET_INDEX_POINT(P_LABEL_RECT1(I),I2,J2,K2,IF_CY)
            IF (IF_CY.EQ.1) THEN
                I22=I2
                IF (K2.GE.1 .AND. K2.LE.NL) THEN
                    J22=K2
                    K22=1
                ELSEIF (K2.GE.(NL+1) .AND. K2.LE.(2*NL)) THEN
                    J22=NL+1
                    K22=K2-NL
                ELSEIF (K2.GE.(2*NL+1) .AND. K2.LE.(3*NL)) THEN
                    J22=NL+2-(K2-2*NL)
                    K22=NL+1
                ELSE
                    J22=1
                    K22=NK+2-K2
                ENDIF        
                KP=K22-(MJ2-1)*CNL
                JP=J22-(MJ1-1)*CNL
                IP=I22-(MI-1)*CNI
                WRITE(2,230) '( ',XPOINTS_RECT(IP,JP,KP),' ',YPOINTS_RECT(IP,JP,KP),' ',ZPOINTS_RECT(IP,JP,KP),' )'
            ELSEIF (IF_CY.EQ.0) THEN
                KP=K2-(MJ2-1)*CNL
                JP=J2-(MJ1-1)*CNL
                IP=I2-(MI-1)*CNI
                WRITE(2,230) '( ',XPOINTS_RECT(IP,JP,KP),' ',YPOINTS_RECT(IP,JP,KP),' ',ZPOINTS_RECT(IP,JP,KP),' )'    
            ENDIF
        ENDDO 
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
    
    If ( iMPI_MyID.EQ.0 ) THEN
    CALL SYSTEM_CLOCK(COUNT2,COUNT_RATE,COUNT_MAX)
    WRITE(*,*) 'INI_UP_WRITE_POINTS FINISHED',(COUNT2-COUNT1)/DBLE(COUNT_RATE)  
    ENDIF
    
    IF (MJ.LE.PNJ) THEN
    DO KP=1,NK
    DO JP=1,CNJ
    DO IP=CNI,CNI
        XLEFT_CY(JP,KP)=0.25*(XPOINTS_CY(IP,JP,KP)+XPOINTS_CY(IP,JP,KP+1) &
                          +XPOINTS_CY(IP,JP+1,KP)+XPOINTS_CY(IP,JP+1,KP+1))
        YLEFT_CY(JP,KP)=0.25*(YPOINTS_CY(IP,JP,KP)+YPOINTS_CY(IP,JP,KP+1) &
                          +YPOINTS_CY(IP,JP+1,KP)+YPOINTS_CY(IP,JP+1,KP+1))
        ZLEFT_CY(JP,KP)=0.25*(ZPOINTS_CY(IP,JP,KP)+ZPOINTS_CY(IP,JP,KP+1) &
                          +ZPOINTS_CY(IP,JP+1,KP)+ZPOINTS_CY(IP,JP+1,KP+1))
    ENDDO
    DO IP=1,1
        XRIGHT_CY(JP,KP)=0.25*(XPOINTS_CY(IP,JP,KP)+XPOINTS_CY(IP,JP,KP+1) &
                          +XPOINTS_CY(IP,JP+1,KP)+XPOINTS_CY(IP,JP+1,KP+1))
        YRIGHT_CY(JP,KP)=0.25*(YPOINTS_CY(IP,JP,KP)+YPOINTS_CY(IP,JP,KP+1) &
                          +YPOINTS_CY(IP,JP+1,KP)+YPOINTS_CY(IP,JP+1,KP+1))
        ZRIGHT_CY(JP,KP)=0.25*(ZPOINTS_CY(IP,JP,KP)+ZPOINTS_CY(IP,JP,KP+1) &
                          +ZPOINTS_CY(IP,JP+1,KP)+ZPOINTS_CY(IP,JP+1,KP+1))
    ENDDO
    ENDDO
    ENDDO  
    DO KP=1,NK
    DO JP=CNJ,CNJ
    DO IP=1,CNI
        XUP_CY(IP,KP)=0.25*(XPOINTS_CY(IP,JP,KP)+XPOINTS_CY(IP,JP,KP+1) &
                          +XPOINTS_CY(IP+1,JP,KP)+XPOINTS_CY(IP+1,JP,KP+1))
        YUP_CY(IP,KP)=0.25*(YPOINTS_CY(IP,JP,KP)+YPOINTS_CY(IP,JP,KP+1) &
                          +YPOINTS_CY(IP+1,JP,KP)+YPOINTS_CY(IP+1,JP,KP+1))
        ZUP_CY(IP,KP)=0.25*(ZPOINTS_CY(IP,JP,KP)+ZPOINTS_CY(IP,JP,KP+1) &
                          +ZPOINTS_CY(IP+1,JP,KP)+ZPOINTS_CY(IP+1,JP,KP+1))
    ENDDO
    ENDDO
    DO JP=1,1
    DO IP=1,CNI
        XDOWN_CY(IP,KP)=0.25*(XPOINTS_CY(IP,JP,KP)+XPOINTS_CY(IP,JP,KP+1) &
                          +XPOINTS_CY(IP+1,JP,KP)+XPOINTS_CY(IP+1,JP,KP+1))
        YDOWN_CY(IP,KP)=0.25*(YPOINTS_CY(IP,JP,KP)+YPOINTS_CY(IP,JP,KP+1) &
                          +YPOINTS_CY(IP+1,JP,KP)+YPOINTS_CY(IP+1,JP,KP+1))
        ZDOWN_CY(IP,KP)=0.25*(ZPOINTS_CY(IP,JP,KP)+ZPOINTS_CY(IP,JP,KP+1) &
                          +ZPOINTS_CY(IP+1,JP,KP)+ZPOINTS_CY(IP+1,JP,KP+1))
    ENDDO
    ENDDO
    ENDDO
    
    ELSE
        
    DO KP=1,CNL
    DO JP=1,CNL
    DO IP=CNI,CNI
        XLEFT_RECT(JP,KP)=0.25*(XPOINTS_RECT(IP,JP,KP)+XPOINTS_RECT(IP,JP,KP+1) &
                          +XPOINTS_RECT(IP,JP+1,KP)+XPOINTS_RECT(IP,JP+1,KP+1))
        YLEFT_RECT(JP,KP)=0.25*(YPOINTS_RECT(IP,JP,KP)+YPOINTS_RECT(IP,JP,KP+1) &
                          +YPOINTS_RECT(IP,JP+1,KP)+YPOINTS_RECT(IP,JP+1,KP+1))
        ZLEFT_RECT(JP,KP)=0.25*(ZPOINTS_RECT(IP,JP,KP)+ZPOINTS_RECT(IP,JP,KP+1) &
                          +ZPOINTS_RECT(IP,JP+1,KP)+ZPOINTS_RECT(IP,JP+1,KP+1))
    ENDDO
    DO IP=1,1
        XRIGHT_RECT(JP,KP)=0.25*(XPOINTS_RECT(IP,JP,KP)+XPOINTS_RECT(IP,JP,KP+1) &
                          +XPOINTS_RECT(IP,JP+1,KP)+XPOINTS_RECT(IP,JP+1,KP+1))
        YRIGHT_RECT(JP,KP)=0.25*(YPOINTS_RECT(IP,JP,KP)+YPOINTS_RECT(IP,JP,KP+1) &
                          +YPOINTS_RECT(IP,JP+1,KP)+YPOINTS_RECT(IP,JP+1,KP+1))
        ZRIGHT_RECT(JP,KP)=0.25*(ZPOINTS_RECT(IP,JP,KP)+ZPOINTS_RECT(IP,JP,KP+1) &
                          +ZPOINTS_RECT(IP,JP+1,KP)+ZPOINTS_RECT(IP,JP+1,KP+1))
    ENDDO
    ENDDO
    ENDDO    
    DO KP=1,CNL
    DO JP=CNL,CNL
    DO IP=1,CNI
        XUP_RECT(IP,KP)=0.25*(XPOINTS_RECT(IP,JP,KP)+XPOINTS_RECT(IP,JP,KP+1) &
                          +XPOINTS_RECT(IP+1,JP,KP)+XPOINTS_RECT(IP+1,JP,KP+1))
        YUP_RECT(IP,KP)=0.25*(YPOINTS_RECT(IP,JP,KP)+YPOINTS_RECT(IP,JP,KP+1) &
                          +YPOINTS_RECT(IP+1,JP,KP)+YPOINTS_RECT(IP+1,JP,KP+1))
        ZUP_RECT(IP,KP)=0.25*(ZPOINTS_RECT(IP,JP,KP)+ZPOINTS_RECT(IP,JP,KP+1) &
                          +ZPOINTS_RECT(IP+1,JP,KP)+ZPOINTS_RECT(IP+1,JP,KP+1))
    ENDDO
    ENDDO
    DO JP=1,1
    DO IP=1,CNI
        XDOWN_RECT(IP,KP)=0.25*(XPOINTS_RECT(IP,JP,KP)+XPOINTS_RECT(IP,JP,KP+1) &
                          +XPOINTS_RECT(IP+1,JP,KP)+XPOINTS_RECT(IP+1,JP,KP+1))
        YDOWN_RECT(IP,KP)=0.25*(YPOINTS_RECT(IP,JP,KP)+YPOINTS_RECT(IP,JP,KP+1) &
                          +YPOINTS_RECT(IP+1,JP,KP)+YPOINTS_RECT(IP+1,JP,KP+1))
        ZDOWN_RECT(IP,KP)=0.25*(ZPOINTS_RECT(IP,JP,KP)+ZPOINTS_RECT(IP,JP,KP+1) &
                          +ZPOINTS_RECT(IP+1,JP,KP)+ZPOINTS_RECT(IP+1,JP,KP+1))
    ENDDO
    ENDDO
    ENDDO    
    DO KP=1,1
    DO JP=1,CNL
    DO IP=1,CNI
        XFRONT_RECT(IP,JP)=0.25*(XPOINTS_RECT(IP,JP,KP)+XPOINTS_RECT(IP+1,JP,KP) &
                          +XPOINTS_RECT(IP,JP+1,KP)+XPOINTS_RECT(IP+1,JP+1,KP))
        YFRONT_RECT(IP,JP)=0.25*(YPOINTS_RECT(IP,JP,KP)+YPOINTS_RECT(IP+1,JP,KP) &
                          +YPOINTS_RECT(IP,JP+1,KP)+YPOINTS_RECT(IP+1,JP+1,KP))
        ZFRONT_RECT(IP,JP)=0.25*(ZPOINTS_RECT(IP,JP,KP)+ZPOINTS_RECT(IP+1,JP,KP) &
                          +ZPOINTS_RECT(IP,JP+1,KP)+ZPOINTS_RECT(IP+1,JP+1,KP))
    ENDDO
    ENDDO
    ENDDO   
    DO KP=CNL,CNL
    DO JP=1,CNL
    DO IP=1,CNI
        XBACK_RECT(IP,JP)=0.25*(XPOINTS_RECT(IP,JP,KP)+XPOINTS_RECT(IP+1,JP,KP) &
                          +XPOINTS_RECT(IP,JP+1,KP)+XPOINTS_RECT(IP+1,JP+1,KP))
        YBACK_RECT(IP,JP)=0.25*(YPOINTS_RECT(IP,JP,KP)+YPOINTS_RECT(IP+1,JP,KP) &
                          +YPOINTS_RECT(IP,JP+1,KP)+YPOINTS_RECT(IP+1,JP+1,KP))
        ZBACK_RECT(IP,JP)=0.25*(ZPOINTS_RECT(IP,JP,KP)+ZPOINTS_RECT(IP+1,JP,KP) &
                          +ZPOINTS_RECT(IP,JP+1,KP)+ZPOINTS_RECT(IP+1,JP+1,KP))
    ENDDO
    ENDDO
    ENDDO 
        
        
        
    ENDIF
    

    
230 FORMAT(A2,E22.15,A1,E22.15,A1,E22.15,A2)
231 FORMAT(I0) 

! 跳过初值映射
    

      
    !INITIALIZE P AND U FIELDS
    
    !OPEN(1,FILE=TRIM(ADJUSTL(DIR))//'p',STATUS='OLD')
    !OPEN(2,FILE=TRIM(ADJUSTL(DIR))//'U',STATUS='OLD')
    !DO I=1,22
    !    READ(1,*)
    !    READ(2,*)
    !ENDDO
    !DO J=1,REF_N2
    !DO I=1,REF_N1
    !    READ(1,*) PREF(I,J)
    !    READ(2,'(A)') CHAR3
    !    IP1=INDEX(CHAR3,'(') 
	   ! CHAR3(1:IP1)=''
	   ! IP1=INDEX(CHAR3,')')
	   ! CHAR3(IP1:IP1)=''
	   ! READ(CHAR3,*) UREF(1,I,J),UREF(2,I,J),U3
    !ENDDO
    !ENDDO
    !CLOSE(1)
    !CLOSE(2)
    ! UREF(1,:,:)=UX_CONSTANT
    ! UREF(2,:,:)=0.0
    ! UREF(3,:,:)=0.0 ! 无用
    ! PREF=PREF_CONSTANT
    
    ! OPEN(1,FILE='mesh/suboff_mesh_2d_ADD.plt',STATUS='OLD')   
    !   READ(1,*)
    !   DO I=1,REF_N1+1 
    !   DO J=1,REF_N2+1
    !       READ(1,*) X0_REF0(I,J),Y0_REF0(I,J)
    !   ENDDO
    !   ENDDO  
    ! CLOSE(1)
    !   DO J=1,REF_N2
    !       X00_REF00(J)=0.5*(X0_REF0(REF_N1+1,J)+X0_REF0(REF_N1+1,J+1))
    !       Y00_REF00(J)=0.5*(Y0_REF0(REF_N1+1,J)+Y0_REF0(REF_N1+1,J+1))
    !       L2_REF(J)=SQRT((X0_REF0(REF_N1+1,J+1)-X0_REF0(REF_N1+1,J))**2+(Y0_REF0(REF_N1+1,J+1)-Y0_REF0(REF_N1+1,J))**2)
    !   ENDDO
    
    
    ! IF (iMPI_MyID.EQ.0) THEN
    ! OPEN(11,FILE='TEST2.plt',STATUS='UNKNOWN')
    ! DO KP=1,NK
    ! DO JP=1,CNJ
    ! DO IP=1,CNI
    ! XTEMP=XCELL_CY(IP,JP,KP)
    ! YTEMP=SQRT(YCELL_CY(IP,JP,KP)**2+ZCELL_CY(IP,JP,KP)**2)
    ! CALL PROJECTION(XTEMP,YTEMP,ITEMP,JTEMP)
    ! WRITE(11,'(4(E21.14,2X))') XTEMP,YTEMP,DBLE(ITEMP),DBLE(JTEMP)
    ! ENDDO
    ! ENDDO
    ! ENDDO 
    ! CLOSE(11)
    ! ENDIF
    
    
    ! IF (MJ.LE.PNJ) THEN
    ! DO KP=1,NK
    ! DO JP=1,CNJ
    ! DO IP=1,CNI
    ! XTEMP=XCELL_CY(IP,JP,KP)
    ! YTEMP=SQRT(YCELL_CY(IP,JP,KP)**2+ZCELL_CY(IP,JP,KP)**2)
    ! CALL PROJECTION(XTEMP,YTEMP,ITEMP,JTEMP)
    ! UREF_CY(1,IP,JP,KP)=UREF(1,ITEMP,JTEMP)
    ! UREF_CY(2,IP,JP,KP)=UREF(2,ITEMP,JTEMP)*YCELL_CY(IP,JP,KP)/YTEMP
    ! UREF_CY(3,IP,JP,KP)=UREF(2,ITEMP,JTEMP)*ZCELL_CY(IP,JP,KP)/YTEMP  
    ! PREF_CY(IP,JP,KP)=PREF(ITEMP,JTEMP)
    ! ENDDO
    ! ENDDO
    ! ENDDO    
    ! DO KP=1,NK
    ! DO JP=1,CNJ
    ! XTEMP=XLEFT_CY(JP,KP)
    ! YTEMP=SQRT(YLEFT_CY(JP,KP)**2+ZLEFT_CY(JP,KP)**2)
    ! CALL PROJECTION(XTEMP,YTEMP,ITEMP,JTEMP)  
    ! U_LEFT_CY(1,JP,KP)=UREF(1,ITEMP,JTEMP)
    ! U_LEFT_CY(2,JP,KP)=UREF(2,ITEMP,JTEMP)*YLEFT_CY(JP,KP)/YTEMP
    ! U_LEFT_CY(3,JP,KP)=UREF(2,ITEMP,JTEMP)*ZLEFT_CY(JP,KP)/YTEMP
    ! P_LEFT_CY(JP,KP)=PREF(ITEMP,JTEMP)
    ! XTEMP=XRIGHT_CY(JP,KP)
    ! YTEMP=SQRT(YRIGHT_CY(JP,KP)**2+ZRIGHT_CY(JP,KP)**2)
    ! CALL PROJECTION(XTEMP,YTEMP,ITEMP,JTEMP)  
    ! U_RIGHT_CY(1,JP,KP)=UREF(1,ITEMP,JTEMP)
    ! U_RIGHT_CY(2,JP,KP)=UREF(2,ITEMP,JTEMP)*YRIGHT_CY(JP,KP)/YTEMP
    ! U_RIGHT_CY(3,JP,KP)=UREF(2,ITEMP,JTEMP)*ZRIGHT_CY(JP,KP)/YTEMP
    ! P_RIGHT_CY(JP,KP)=PREF(ITEMP,JTEMP)
    ! ENDDO
    ! ENDDO 
    ! DO KP=1,NK
    ! DO IP=1,CNI
    ! XTEMP=XUP_CY(IP,KP)
    ! YTEMP=SQRT(YUP_CY(IP,KP)**2+ZUP_CY(IP,KP)**2)
    ! CALL PROJECTION(XTEMP,YTEMP,ITEMP,JTEMP)
    ! U_UP_CY(1,IP,KP)=UREF(1,ITEMP,JTEMP)
    ! U_UP_CY(2,IP,KP)=UREF(2,ITEMP,JTEMP)*YUP_CY(IP,KP)/YTEMP
    ! U_UP_CY(3,IP,KP)=UREF(2,ITEMP,JTEMP)*ZUP_CY(IP,KP)/YTEMP  
    ! P_UP_CY(IP,KP)=PREF(ITEMP,JTEMP)
    ! XTEMP=XDOWN_CY(IP,KP)
    ! YTEMP=SQRT(YDOWN_CY(IP,KP)**2+ZDOWN_CY(IP,KP)**2)
    ! CALL PROJECTION(XTEMP,YTEMP,ITEMP,JTEMP)
    ! U_DOWN_CY(1,IP,KP)=UREF(1,ITEMP,JTEMP)
    ! U_DOWN_CY(2,IP,KP)=UREF(2,ITEMP,JTEMP)*YDOWN_CY(IP,KP)/YTEMP
    ! U_DOWN_CY(3,IP,KP)=UREF(2,ITEMP,JTEMP)*ZDOWN_CY(IP,KP)/YTEMP 
    ! P_DOWN_CY(IP,KP)=PREF(ITEMP,JTEMP)
    ! ENDDO
    ! ENDDO
        
    ! ELSE
        
    ! DO KP=1,CNL
    ! DO JP=1,CNL
    ! DO IP=1,CNI
    ! XTEMP=XCELL_RECT(IP,JP,KP)
    ! YTEMP=SQRT(YCELL_RECT(IP,JP,KP)**2+ZCELL_RECT(IP,JP,KP)**2)
    ! CALL PROJECTION(XTEMP,YTEMP,ITEMP,JTEMP)
    ! UREF_RECT(1,IP,JP,KP)=UREF(1,ITEMP,JTEMP)
    ! UREF_RECT(2,IP,JP,KP)=UREF(2,ITEMP,JTEMP)*YCELL_RECT(IP,JP,KP)/YTEMP
    ! UREF_RECT(3,IP,JP,KP)=UREF(2,ITEMP,JTEMP)*ZCELL_RECT(IP,JP,KP)/YTEMP 
    ! PREF_RECT(IP,JP,KP)=PREF(ITEMP,JTEMP)
    ! ENDDO
    ! ENDDO
    ! ENDDO 
    ! DO KP=1,CNL
    ! DO JP=1,CNL
    ! XTEMP=XLEFT_RECT(JP,KP)
    ! YTEMP=SQRT(YLEFT_RECT(JP,KP)**2+ZLEFT_RECT(JP,KP)**2)
    ! CALL PROJECTION(XTEMP,YTEMP,ITEMP,JTEMP)  
    ! U_LEFT_RECT(1,JP,KP)=UREF(1,ITEMP,JTEMP)
    ! U_LEFT_RECT(2,JP,KP)=UREF(2,ITEMP,JTEMP)*YLEFT_RECT(JP,KP)/YTEMP
    ! U_LEFT_RECT(3,JP,KP)=UREF(2,ITEMP,JTEMP)*ZLEFT_RECT(JP,KP)/YTEMP
    ! P_LEFT_RECT(JP,KP)=PREF(ITEMP,JTEMP)
    ! XTEMP=XRIGHT_RECT(JP,KP)
    ! YTEMP=SQRT(YRIGHT_RECT(JP,KP)**2+ZRIGHT_RECT(JP,KP)**2)
    ! CALL PROJECTION(XTEMP,YTEMP,ITEMP,JTEMP)  
    ! U_RIGHT_RECT(1,JP,KP)=UREF(1,ITEMP,JTEMP)
    ! U_RIGHT_RECT(2,JP,KP)=UREF(2,ITEMP,JTEMP)*YRIGHT_RECT(JP,KP)/YTEMP
    ! U_RIGHT_RECT(3,JP,KP)=UREF(2,ITEMP,JTEMP)*ZRIGHT_RECT(JP,KP)/YTEMP
    ! P_RIGHT_RECT(JP,KP)=PREF(ITEMP,JTEMP)
    ! ENDDO
    ! ENDDO    
    ! DO KP=1,CNL
    ! DO IP=1,CNI
    ! XTEMP=XUP_RECT(IP,KP)
    ! YTEMP=SQRT(YUP_RECT(IP,KP)**2+ZUP_RECT(IP,KP)**2)
    ! CALL PROJECTION(XTEMP,YTEMP,ITEMP,JTEMP)
    ! U_UP_RECT(1,IP,KP)=UREF(1,ITEMP,JTEMP)
    ! U_UP_RECT(2,IP,KP)=UREF(2,ITEMP,JTEMP)*YUP_RECT(IP,KP)/YTEMP
    ! U_UP_RECT(3,IP,KP)=UREF(2,ITEMP,JTEMP)*ZUP_RECT(IP,KP)/YTEMP  
    ! P_UP_RECT(IP,KP)=PREF(ITEMP,JTEMP)
    ! XTEMP=XDOWN_RECT(IP,KP)
    ! YTEMP=SQRT(YDOWN_RECT(IP,KP)**2+ZDOWN_RECT(IP,KP)**2)
    ! CALL PROJECTION(XTEMP,YTEMP,ITEMP,JTEMP)
    ! U_DOWN_RECT(1,IP,KP)=UREF(1,ITEMP,JTEMP)
    ! U_DOWN_RECT(2,IP,KP)=UREF(2,ITEMP,JTEMP)*YDOWN_RECT(IP,KP)/YTEMP
    ! U_DOWN_RECT(3,IP,KP)=UREF(2,ITEMP,JTEMP)*ZDOWN_RECT(IP,KP)/YTEMP 
    ! P_DOWN_RECT(IP,KP)=PREF(ITEMP,JTEMP)
    ! ENDDO
    ! ENDDO    
    ! DO JP=1,CNL
    ! DO IP=1,CNI
    ! XTEMP=XFRONT_RECT(IP,JP)
    ! YTEMP=SQRT(YFRONT_RECT(IP,JP)**2+ZFRONT_RECT(IP,JP)**2)
    ! CALL PROJECTION(XTEMP,YTEMP,ITEMP,JTEMP)
    ! U_FRONT_RECT(1,IP,JP)=UREF(1,ITEMP,JTEMP)
    ! U_FRONT_RECT(2,IP,JP)=UREF(2,ITEMP,JTEMP)*YFRONT_RECT(IP,JP)/YTEMP
    ! U_FRONT_RECT(3,IP,JP)=UREF(2,ITEMP,JTEMP)*ZFRONT_RECT(IP,JP)/YTEMP  
    ! P_FRONT_RECT(IP,JP)=PREF(ITEMP,JTEMP)
    ! XTEMP=XBACK_RECT(IP,JP)
    ! YTEMP=SQRT(YBACK_RECT(IP,JP)**2+ZBACK_RECT(IP,JP)**2)
    ! CALL PROJECTION(XTEMP,YTEMP,ITEMP,JTEMP)
    ! U_BACK_RECT(1,IP,JP)=UREF(1,ITEMP,JTEMP)
    ! U_BACK_RECT(2,IP,JP)=UREF(2,ITEMP,JTEMP)*YBACK_RECT(IP,JP)/YTEMP
    ! U_BACK_RECT(3,IP,JP)=UREF(2,ITEMP,JTEMP)*ZBACK_RECT(IP,JP)/YTEMP
    ! P_BACK_RECT(IP,JP)=PREF(ITEMP,JTEMP)
    ! ENDDO
    ! ENDDO    
    ! ENDIF


! --- DIRECT CONSTANT ASSIGNMENTS ---
IF (MJ.LE.PNJ) THEN
          
  ! CYLINDER BLOCK (CY) - Cell Centers
  UREF_CY(1,:,:,:) = UX_CONSTANT
  UREF_CY(2,:,:,:) = UY_CONSTANT
  UREF_CY(3,:,:,:) = UZ_CONSTANT
  PREF_CY(:,:,:)   = PREF_CONSTANT
  
  ! CYLINDER BLOCK (CY) - Boundaries
  U_LEFT_CY(1,:,:) = UX_CONSTANT
  U_LEFT_CY(2,:,:) = UY_CONSTANT
  U_LEFT_CY(3,:,:) = UZ_CONSTANT
  P_LEFT_CY(:,:)   = PREF_CONSTANT
  
  U_RIGHT_CY(1,:,:) = UX_CONSTANT
  U_RIGHT_CY(2,:,:) = UY_CONSTANT
  U_RIGHT_CY(3,:,:) = UZ_CONSTANT
  P_RIGHT_CY(:,:)   = PREF_CONSTANT
  
  U_UP_CY(1,:,:) = UX_CONSTANT
  U_UP_CY(2,:,:) = UY_CONSTANT
  U_UP_CY(3,:,:) = UZ_CONSTANT
  P_UP_CY(:,:)   = PREF_CONSTANT
  
  U_DOWN_CY(1,:,:) = UX_CONSTANT
  U_DOWN_CY(2,:,:) = UY_CONSTANT
  U_DOWN_CY(3,:,:) = UZ_CONSTANT
  P_DOWN_CY(:,:)   = PREF_CONSTANT
  
ELSE
  
  ! RECTANGULAR BLOCK (RECT) - Cell Centers
  UREF_RECT(1,:,:,:) = UX_CONSTANT
  UREF_RECT(2,:,:,:) = UY_CONSTANT
  UREF_RECT(3,:,:,:) = UZ_CONSTANT
  PREF_RECT(:,:,:)   = PREF_CONSTANT
  
  ! RECTANGULAR BLOCK (RECT) - Boundaries
  U_LEFT_RECT(1,:,:) = UX_CONSTANT
  U_LEFT_RECT(2,:,:) = UY_CONSTANT
  U_LEFT_RECT(3,:,:) = UZ_CONSTANT
  P_LEFT_RECT(:,:)   = PREF_CONSTANT
  
  U_RIGHT_RECT(1,:,:) = UX_CONSTANT
  U_RIGHT_RECT(2,:,:) = UY_CONSTANT
  U_RIGHT_RECT(3,:,:) = UZ_CONSTANT
  P_RIGHT_RECT(:,:)   = PREF_CONSTANT
  
  U_UP_RECT(1,:,:) = UX_CONSTANT
  U_UP_RECT(2,:,:) = UY_CONSTANT
  U_UP_RECT(3,:,:) = UZ_CONSTANT
  P_UP_RECT(:,:)   = PREF_CONSTANT
  
  U_DOWN_RECT(1,:,:) = UX_CONSTANT
  U_DOWN_RECT(2,:,:) = UY_CONSTANT
  U_DOWN_RECT(3,:,:) = UZ_CONSTANT
  P_DOWN_RECT(:,:)   = PREF_CONSTANT
  
  U_FRONT_RECT(1,:,:) = UX_CONSTANT
  U_FRONT_RECT(2,:,:) = UY_CONSTANT
  U_FRONT_RECT(3,:,:) = UZ_CONSTANT
  P_FRONT_RECT(:,:)   = PREF_CONSTANT
  
  U_BACK_RECT(1,:,:) = UX_CONSTANT
  U_BACK_RECT(2,:,:) = UY_CONSTANT
  U_BACK_RECT(3,:,:) = UZ_CONSTANT
  P_BACK_RECT(:,:)   = PREF_CONSTANT
  
ENDIF
    


      
    OPEN(1,FILE='ReadData/PU_RECT_SUR.plt',STATUS='OLD')
    READ(1,*)
    DO K=1,NL
    DO J=1,NL
        READ(1,*) XP,YP,ZP,P_RECT(J,K),U_RECT(1,J,K),U_RECT(2,J,K),U_RECT(3,J,K)
    ENDDO
    ENDDO
    CLOSE(1)  
      
    OPEN(1,FILE='ReadData/PU_CY_SUR.plt',STATUS='OLD')
        READ(1,*)
        DO K=1,NK
        DO J=1,NJ
            READ(1,*) XP,YP,ZP,P_CY(J,K),U_CY(1,J,K),U_CY(2,J,K),U_CY(3,J,K)
        ENDDO
        ENDDO  
    CLOSE(1)   
      
    IF (MJ.LE.PNJ) THEN
        DO KP=1,NK
        K=KP
        DO JP=1,CNJ
        J=JP+CNJ*(MJ-1)
        U_INLET_CY(:,JP,KP)=U_CY(:,J,K)
        ENDDO
        ENDDO
    ELSE
        DO KP=1,CNL
        K=KP+CNL*(MJ2-1)
        DO JP=1,CNL
        J=JP+CNL*(MJ1-1)
        U_INLET_RECT(:,JP,KP)=U_RECT(:,NL-K+1,J)
        ENDDO
        ENDDO
    ENDIF  
      
    If ( iMPI_MyID.EQ.0 ) THEN
      CALL SYSTEM_CLOCK(COUNT3,COUNT_RATE,COUNT_MAX)
      WRITE(*,*) 'CALCULATE_INI_UP_FIELD FINISHED',(COUNT3-COUNT2)/DBLE(COUNT_RATE)
      ENDIF
    CALL MPI_Barrier( MPI_COMM_WORLD, iMPI_ErrorInfo )
    END SUBROUTINE
    
    
    
    
    
