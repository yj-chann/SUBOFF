PROGRAM MAIN

    !----------------------------------------
    
    INCLUDE 'head.fi'
    
    !----------------------------------------
    !  execute statement
    !----------------------------------------
    !并行初始化并划分笛卡尔拓扑：
    CALL MPI_Init( iMPI_ErrorInfo )
    CALL MPI_Comm_Rank(MPI_COMM_WORLD, iMPI_MyID , iMPI_ErrorInfo)
    CALL MPI_Comm_Size(MPI_COMM_WORLD, iMPI_NumProcs , iMPI_ErrorInfo)
    DIMS(0)= NPR
    DIMS(1)= NPC
    PERIODIC(0) = .FALSE.
    PERIODIC(1) = .FALSE.
    CALL MPI_CART_CREATE(MPI_COMM_WORLD, 2, DIMS, PERIODIC,.FALSE.,MPI_COMM_CART,iMPI_ErrorInfo)
    CALL MPI_CART_COORDS(MPI_COMM_CART, iMPI_MyID, 2, COORD,iMPI_ErrorInfo)
    
    
    CALL MPI_CART_SUB(MPI_COMM_CART, (/.TRUE.,.FALSE./), MPI_COMM_COL, iMPI_ErrorInfo)
    CALL MPI_CART_SUB(MPI_COMM_CART, (/.FALSE.,.TRUE./), MPI_COMM_ROW, iMPI_ErrorInfo)
    !Check the number of process
    !!!NP=NPROW*NPCOL
    IF ( NPP .EQ. iMPI_NumProcs ) THEN
    
    !   the number of process is right
    
        CALL MPI_Barrier( MPI_COMM_WORLD, iMPI_ErrorInfo )
    
    !   Main Program
    
        CALL MAINPROGRAM
    
    ENDIF
    
    !Print*, 'FINISHED'
    
    CALL MPI_Barrier( MPI_COMM_WORLD, iMPI_ErrorInfo )
    CALL MPI_Finalize( iMPI_ErrorInfo )
    
    !----------------------------------------
    
    END PROGRAM
    
    !===== END PROGRAM ======================
    
    !-------------------------------------------------------------------------
    !  MAIN PROGRAM
    !-------------------------------------------------------------------------
    
    SUBROUTINE MAINPROGRAM
    
    !----------------------------------------
    
    INCLUDE 'head.fi'
    
    INCLUDE 'parameter.h'
    
    REAL TIME
    REAL*8 X1(N2,N3),Y1(N2,N3),Z1(N2,N3)
    REAL*8 X,Y,Z,PAI
    REAL*8 A1,A2,A3,A4,A5,A6,A7,UX,UY,UZ,R
    REAL*8 UREF(3,N2,N3),RATIO(N2,N3),UINF,PHI_REF,UREF_ST(3,N2,N3),UREF_ED(3,N2,N3),PHI1,PHI_ED,PHI_ST,RATIO1
    REAL*8 INLET_X(N2,N3),INLET_Y(N2,N3),INLET_Z(N2,N3)
    REAL*8 TEMP1,TEMP2,TEMP3,TEMP4,U_RECT(3,N2,N3)
    INTEGER I,J,K,N0,IP,STI,EDI,NT0,NN1,NN2,NN3,TI,II
    INTEGER INTER
    CHARACTER*400 CHAR1,DIR
    CHARACTER*20 CHAR3
    CHARACTER*1000000 CHAR4
    CHARACTER*1 CHAR2
    
    
    
    PAI=ACOS(-1.0)
    ! DIR='../../../../07_mesh_for_turb/postProcessing/probes1/0/'
    DIR='../../../../06_tg/vel/'


    UINF=1.649194
    PHI_REF=5.623942839039030E-002
        
    OPEN(iMPI_MyID+1000,FILE='ReadData/PU_RECT_SUR.plt',STATUS='OLD',ACTION='READ')
    READ(iMPI_MyID+1000,*)
    DO K=1,N3
    DO J=1,N2
        READ(iMPI_MyID+1000,*) INLET_X(J,K),INLET_Y(J,K),INLET_Z(J,K),TEMP4,U_RECT(1,J,K),U_RECT(2,J,K),U_RECT(3,J,K)
    ENDDO
    ENDDO
    CLOSE(iMPI_MyID+1000)
    
    DO K=1,N3
    DO J=1,N2
    R=SQRT((INLET_Y(J,K))**2+(INLET_Z(J,K))**2)
    TEMP1=DELTA_RATIO*0.5*RMAX
    TEMP2=R-( RMAX-TEMP1 )
    IF (TEMP2.LE.-TEMP1) THEN
    RATIO(J,K)=1.0
    ELSEIF (TEMP2.GE.TEMP1) THEN
    RATIO(J,K)=0.0
    ELSE
    RATIO(J,K)=1.0-0.5*(1.0+TEMP2/TEMP1+SIN(PAI*TEMP2/TEMP1)/PAI)
    ENDIF
    ENDDO
    ENDDO
    
    
    
    IF (iMPI_MyID.EQ.0) THEN	
    OPEN(1,FILE='RATIO.plt',STATUS='UNKNOWN',ACTION='WRITE') 
    WRITE(1,104) N2,N3
    DO K=1,N3
    DO J=1,N2  
        WRITE(1,105) INLET_X(J,K),INLET_Y(J,K),INLET_Z(J,K),RATIO(J,K)
    ENDDO
    ENDDO
    CLOSE(1)	
    
    ! 用于Tecplot绘制的版本
    OPEN(1,FILE='Tecplot_InputFiles/ITIweightFunc.plt',STATUS='UNKNOWN',ACTION='WRITE') 
    WRITE(1, '(A)') 'TITLE = "ITI weight function"'
    WRITE(1, '(A)') 'VARIABLES = "X", "Y", "Z" ,"weight"'
    WRITE(1, '(A, I0, A, I0, A)') 'ZONE T="ITI weight function", I=', N2, ', J=', N3, ', F=POINT'
    DO K=1,N3
    DO J=1,N2  
        WRITE(1,105) INLET_X(J,K),INLET_Y(J,K),INLET_Z(J,K),RATIO(J,K)
    ENDDO
    ENDDO
    CLOSE(1)
    ENDIF
    104 FORMAT('Zone  i=',I0,' j=',I0,' f=point') 
    105 FORMAT(4(E22.15,2X)) 
    
    
    
    ! IF (iMPI_MyID.EQ.0) THEN
    !     OPEN(11,FILE=TRIM(ADJUSTL(DIR))//'u000001.plt',STATUS='OLD',ACTION='READ') 
    !     OPEN(1,FILE='../INLET/0/U',STATUS='UNKNOWN',ACTION='WRITE')
    !     WRITE(1,'(I0)') NALL
    !     WRITE(1,*) '('
    !     DO K=1,N3
    !     DO J=1,N2
    !         READ(11,*) TEMP1, TEMP2 ,TEMP3
    !         TEMP1=( TEMP1-UINF)*RATIO(J,K) +U_RECT(1,J,K)
    !         TEMP2=( TEMP2 )*RATIO(J,K) +U_RECT(2,J,K)
    !         TEMP3=( TEMP3 )*RATIO(J,K) +U_RECT(3,J,K)
    !     WRITE(1,230) '( ',TEMP1,' ',TEMP2,' ',TEMP3,' )'
    !     ENDDO
    !     ENDDO
    !     WRITE(1,*) ')'
    !     CLOSE(1)
    !     CLOSE(11)
    ! ENDIF

    ! 0 时刻无脉动
    IF (iMPI_MyID.EQ.0) THEN
        OPEN(1,FILE='../INLET/0/U',STATUS='UNKNOWN',ACTION='WRITE')
        WRITE(1,'(I0)') NALL
        WRITE(1,*) '('
        DO K=1,N3
        DO J=1,N2
        WRITE(1,230) '( ',U_RECT(1,J,K),' ',U_RECT(2,J,K),' ', U_RECT(3,J,K),' )'
        ENDDO
        ENDDO
        WRITE(1,*) ')'
        CLOSE(1)
    ENDIF
    
    
    
    IF (iMPI_MyID .EQ. 0) THEN
    
    OPEN(11,FILE=TRIM(ADJUSTL(DIR))//'u000001.plt',STATUS='OLD',ACTION='READ')
    OPEN(14,FILE='ReadData/surfaceFieldValue.dat',STATUS='OLD',ACTION='READ') 
    DO I=1,5
    READ(14,*)
    ENDDO
    DO K=1,N3
        DO J=1,N2
            READ(11,*)  UREF_ST(1,J,K),UREF_ST(2,J,K),UREF_ST(3,J,K)
        ENDDO
    ENDDO
    CLOSE(11)
    READ(14,*) TIME,PHI_ST
    
    WRITE(CHAR3,'(I6.6)') NT-MID+1
    OPEN(11,FILE=TRIM(ADJUSTL(DIR))//'u'//TRIM(ADJUSTL(CHAR3))//'.plt',STATUS='OLD',ACTION='READ')
    DO I=1,NT-MID-1
        READ(14,*)
    ENDDO
    DO K=1,N3
        DO J=1,N2
            READ(11,*)  UREF_ED(1,J,K),UREF_ED(2,J,K),UREF_ED(3,J,K)
        ENDDO
    ENDDO
    READ(14,*) TIME,PHI_ED
    CLOSE(11)
    CLOSE(14)
    
    DO K=1,N3
        DO J=1,N2
            TEMP1=( UREF_ED(1,J,K)-UINF )*RATIO(J,K)+U_RECT(1,J,K)
            UREF_ED(1,J,K)=TEMP1
            TEMP2=( UREF_ED(2,J,K))*RATIO(J,K)+U_RECT(2,J,K)
            UREF_ED(2,J,K)=TEMP2
            TEMP3=( UREF_ED(3,J,K))*RATIO(J,K)+U_RECT(3,J,K)
            UREF_ED(3,J,K)=TEMP3
            TEMP1=( UREF_ST(1,J,K)-UINF)*RATIO(J,K)+U_RECT(1,J,K)
            UREF_ST(1,J,K)=TEMP1
            TEMP2=( UREF_ST(2,J,K))*RATIO(J,K)+U_RECT(2,J,K)
            UREF_ST(2,J,K)=TEMP2
            TEMP3=( UREF_ST(3,J,K))*RATIO(J,K)+U_RECT(3,J,K)
            UREF_ST(3,J,K)=TEMP3
    ENDDO
    ENDDO
    ENDIF
    
    ! ---------------------------------------------------------
    !  Broadcast the data to ALL other ranks
    ! ---------------------------------------------------------
    ! Make sure all processes wait here before broadcasting
    CALL MPI_BARRIER(MPI_COMM_WORLD, iMPI_ErrorInfo)
    CALL MPI_BCAST(PHI_ED, 1, mpi_real8, 0, MPI_COMM_WORLD, iMPI_ErrorInfo)
    CALL MPI_BCAST(PHI_ST, 1, mpi_real8, 0, MPI_COMM_WORLD, iMPI_ErrorInfo)
    CALL MPI_BCAST(UREF_ED, 3*N2*N3, mpi_real8, 0, MPI_COMM_WORLD, iMPI_ErrorInfo)
    CALL MPI_BCAST(UREF_ST, 3*N2*N3, mpi_real8, 0, MPI_COMM_WORLD, iMPI_ErrorInfo)
    
    
    
    
    NT0=NT/NPP
    STI=iMPI_MyID*NT0+ST ! STRAT FROM 0
    EDI=( iMPI_MyID+1 )*NT0+ST-1
    NN1=NINT(1.0/DT)
    
    OPEN(iMPI_MyID+14,FILE='ReadData/surfaceFieldValue.dat',STATUS='OLD',ACTION='READ') 
    DO I=1,5+STI
        READ(iMPI_MyID+14,*)
    ENDDO
    
    WRITE(CHAR3, '(I0)') iMPI_MyID
    OPEN(iMPI_MyID+20, FILE='Tecplot_InputFiles/UINLET/UINLET_Animate_'//TRIM(ADJUSTL(CHAR3))//'.plt', STATUS='UNKNOWN',ACTION='WRITE')
    WRITE(iMPI_MyID+20, '(A)') 'TITLE = "UInlet Animation"'
    WRITE(iMPI_MyID+20, '(A)') 'VARIABLES = "X", "Y", "Z", "u", "v", "w"'
    
    
    DO I=STI,EDI
        NN2=MOD(I-ST+1,NN1)
        IF (NN2.EQ.0) THEN
        NN3=(I-ST+1)/NN1
        WRITE(CHAR1,'(I0)') NN3
        ELSE 
        TIME=DT*(I-ST+1)
        WRITE(CHAR1,'(F10.7)') TIME
    505     IP=INDEX(CHAR1,'0',BACK=.TRUE.) 
        CHAR1(IP:IP)=''
        CHAR2=CHAR1(IP-1:IP-1)
        READ(CHAR2,'(I1)') N0
        IF (N0.EQ.0) GOTO 505
        ENDIF

        WRITE(CHAR3,'(I6.6)') I+1
        OPEN(iMPI_MyID+11,FILE=TRIM(ADJUSTL(DIR))//'u'//TRIM(ADJUSTL(CHAR3))//'.plt',STATUS='OLD',ACTION='READ')
        DO K=1,N3
            DO J=1,N2
                READ(iMPI_MyID+11,*) UREF(1,J,K),UREF(2,J,K),UREF(3,J,K)
            ENDDO
        ENDDO
        CLOSE(iMPI_MyID+11)

        READ(iMPI_MyID+14,*) TIME,PHI1


        OPEN(iMPI_MyID+1000,FILE='../INLET/'//TRIM(ADJUSTL(CHAR1))//'/U',STATUS='UNKNOWN',ACTION='WRITE')
        WRITE(iMPI_MyID+1000,'(I0)') NALL
        WRITE(iMPI_MyID+1000,*) '('
    
        IF (MOD(I+1, 10) .EQ. 0) THEN
            WRITE(iMPI_MyID+20, '(A, I0, A, I0, A, A)') 'ZONE T="UINLET", I=', N2, ', J=', N3, ', F=POINT, STRANDID=1, SOLUTIONTIME=', TRIM(ADJUSTL(CHAR1))
        ENDIF
    
    
        IF ((I-ST+1).LE.(NT-MID)) THEN !!!
        DO K=1,N3
        DO J=1,N2 
        TEMP1=(( UREF(1,J,K)-UINF )*RATIO(J,K) +U_RECT(1,J,K))*PHI_REF/ABS(PHI1)
        TEMP2=(( UREF(2,J,K) )*RATIO(J,K) +U_RECT(2,J,K))*PHI_REF/ABS(PHI1)
        TEMP3=(( UREF(3,J,K) )*RATIO(J,K) +U_RECT(3,J,K))*PHI_REF/ABS(PHI1)
    
    
        WRITE(iMPI_MyID+1000,230) '( ',TEMP1,' ',TEMP2,' ',TEMP3,' )'
    
            IF (MOD(I+1, 10) .EQ. 0) THEN
                WRITE(iMPI_MyID+20, '(6(E22.15,2X))') INLET_X(J,K), INLET_Y(J,K), INLET_Z(J,K), TEMP1, TEMP2, TEMP3
            ENDIF
    
        ENDDO
        ENDDO
    
        ELSE
            RATIO1=1.0-DBLE(I-ST-NT+MID)/DBLE(MID)
            DO K=1,N3
            DO J=1,N2 
            TEMP1=RATIO1*UREF_ED(1,J,K)*PHI_REF/ABS(PHI_ED)+(1.0-RATIO1)*UREF_ST(1,J,K)*PHI_REF/ABS(PHI_ST)
            TEMP2=RATIO1*UREF_ED(2,J,K)*PHI_REF/ABS(PHI_ED)+(1.0-RATIO1)*UREF_ST(2,J,K)*PHI_REF/ABS(PHI_ST)
            TEMP3=RATIO1*UREF_ED(3,J,K)*PHI_REF/ABS(PHI_ED)+(1.0-RATIO1)*UREF_ST(3,J,K)*PHI_REF/ABS(PHI_ST) 
            WRITE(iMPI_MyID+1000,230) '( ',TEMP1,' ',TEMP2,' ',TEMP3,' )'
                IF (MOD(I+1, 10) .EQ. 0) THEN
                    WRITE(iMPI_MyID+20, '(6(E22.15,2X))') INLET_X(J,K), INLET_Y(J,K), INLET_Z(J,K), TEMP1, TEMP2, TEMP3
                ENDIF
            ENDDO
            ENDDO
        ENDIF
        WRITE(iMPI_MyID+1000,*) ')'
        CLOSE(iMPI_MyID+1000)
    
    ENDDO
    
    CLOSE(iMPI_MyID+14)
    CLOSE(iMPI_MyID+20)
    
    
    ! 重复DT~NT*DT
    
    
    230 FORMAT(A2,E22.15,A1,E22.15,A1,E22.15,A2)  
    
    END SUBROUTINE
    