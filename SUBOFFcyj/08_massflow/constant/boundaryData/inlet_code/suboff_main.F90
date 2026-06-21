!========================================!
!      CHANNEL STAGGERED BI_PARALLEL     !
!========================================!

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
    REAL*8 UREF(3,N2,N3),RATIO(N2,N3)
	REAL*8 INLET_X(N2,N3),INLET_Y(N2,N3),INLET_Z(N2,N3)
	REAL*8 TEMP1,TEMP2,TEMP3,TEMP4,U_RECT(3,N2,N3)
    INTEGER I,J,K,N0,IP,STI,EDI,NT0,NN1,NN2,NN3,TI,II
    INTEGER INTER
    CHARACTER*400 CHAR1,DIR
    CHARACTER*20 CHAR3
    CHARACTER*1000000 CHAR4
    CHARACTER*1 CHAR2
    
    
    PAI=ACOS(-1.0)
    DIR='../../../../07_mesh_for_turb/postProcessing/probes1/0/'


        
    OPEN(iMPI_MyID+1,FILE='PU_RECT_SUR.plt',STATUS='OLD',ACTION='READ')
    READ(iMPI_MyID+1,*)
    DO K=1,N3
    DO J=1,N2
        READ(iMPI_MyID+1,*) INLET_X(J,K),INLET_Y(J,K),INLET_Z(J,K),TEMP4,U_RECT(1,J,K),U_RECT(2,J,K),U_RECT(3,J,K)
    ENDDO
    ENDDO
    CLOSE(iMPI_MyID+1)
    
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
	OPEN(1,FILE='RATIO.plt',STATUS='UNKNOWN') 
    WRITE(1,104) N2,N3
    DO K=1,N3
    DO J=1,N2  
        WRITE(1,105) INLET_X(J,K),INLET_Y(J,K),INLET_Z(J,K),RATIO(J,K)
    ENDDO
    ENDDO
    CLOSE(1)	
    ! 用于Tecplot绘制的版本
    OPEN(1,FILE='Tecplot_InputFiles/ITIweightFunc.plt',STATUS='UNKNOWN') 
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
105 FORMAT(4(E21.14,2X)) 
	

    
    IF (iMPI_MyID.EQ.0) THEN
        OPEN(11,FILE=TRIM(ADJUSTL(DIR))//'Ux',STATUS='OLD') 
		OPEN(12,FILE=TRIM(ADJUSTL(DIR))//'Uy',STATUS='OLD')
		OPEN(13,FILE=TRIM(ADJUSTL(DIR))//'Uz',STATUS='OLD')
        !probe information ends at t=deltaT at line N2*N3, namely, the numder of probes
        !followed by a line of title
        !and then the nmumerical value, so it starts at t=deltaT at line N2*N3+2, namely, the numder of probes plus 2
        !so If ST=0, all the numerical data will be maintained.
        !And ,the first numerical data will be maintained as time=0, while actually not.
        DO I=1,N2*N3+2+ST-1
        READ(11,*)
		READ(12,*)
		READ(13,*)
        ENDDO 
        READ(11,*) TIME,(( UREF(1,J,K),J=1,N2),K=1,N3)
		READ(12,*) TIME,(( UREF(2,J,K),J=1,N2),K=1,N3)
		READ(13,*) TIME,(( UREF(3,J,K),J=1,N2),K=1,N3)
        OPEN(1,FILE='../INLET/0/U',STATUS='UNKNOWN')
        WRITE(1,*) '('
        DO K=1,N3
        DO J=1,N2
		TEMP1=( UREF(1,J,K)-U_RECT(1,J,K) )*RATIO(J,K) +U_RECT(1,J,K)
		TEMP2=( UREF(2,J,K)-U_RECT(2,J,K) )*RATIO(J,K) +U_RECT(2,J,K)
		TEMP3=( UREF(3,J,K)-U_RECT(3,J,K) )*RATIO(J,K) +U_RECT(3,J,K)
        WRITE(1,230) '( ',TEMP1,' ',TEMP2,' ',TEMP3,' )'
        ENDDO
        ENDDO
        WRITE(1,*) ')'
        CLOSE(1)
        CLOSE(11)
		CLOSE(12)
		CLOSE(13)
    ENDIF
    
    
    NT0=NT/NPP
    STI=iMPI_MyID*NT0+ST
    EDI=( iMPI_MyID+1 )*NT0+ST-1
    NN1=NINT(1.0/DT)
    
    OPEN(iMPI_MyID+11,FILE=TRIM(ADJUSTL(DIR))//'Ux',STATUS='OLD') 
	OPEN(iMPI_MyID+12,FILE=TRIM(ADJUSTL(DIR))//'Uy',STATUS='OLD')
	OPEN(iMPI_MyID+13,FILE=TRIM(ADJUSTL(DIR))//'Uz',STATUS='OLD') 
    DO I=1,N2*N3+2+STI-1
        READ(iMPI_MyID+11,*)
        READ(iMPI_MyID+12,*)
        READ(iMPI_MyID+13,*)
    ENDDO
    
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
        READ(iMPI_MyID+11,*) TIME,(( UREF(1,J,K),J=1,N2),K=1,N3)
		READ(iMPI_MyID+12,*) TIME,(( UREF(2,J,K),J=1,N2),K=1,N3)
		READ(iMPI_MyID+13,*) TIME,(( UREF(3,J,K),J=1,N2),K=1,N3)
        OPEN(iMPI_MyID+1,FILE='../INLET/'//TRIM(ADJUSTL(CHAR1))//'/U',STATUS='UNKNOWN')
        WRITE(iMPI_MyID+1,*) '('
        DO K=1,N3
        DO J=1,N2 
		TEMP1=( UREF(1,J,K)-U_RECT(1,J,K) )*RATIO(J,K) +U_RECT(1,J,K)
		TEMP2=( UREF(2,J,K)-U_RECT(2,J,K) )*RATIO(J,K) +U_RECT(2,J,K)
		TEMP3=( UREF(3,J,K)-U_RECT(3,J,K) )*RATIO(J,K) +U_RECT(3,J,K)
        WRITE(iMPI_MyID+1,230) '( ',TEMP1,' ',TEMP2,' ',TEMP3,' )'
        ENDDO
        ENDDO
        WRITE(iMPI_MyID+1,*) ')'
        CLOSE(iMPI_MyID+1)

        !输出中间某个时间步的未加权和加权入口场
        IF(iMPI_MyID+1.EQ.NPP/2.AND.I.EQ.EDI) THEN
            OPEN(iMPI_MyID*10+1,FILE='Tecplot_InputFiles/u_INLET_noWeight.plt',STATUS='UNKNOWN') 
            WRITE(iMPI_MyID*10+1, '(A)') 'TITLE = "u_INLET_noWeight"'
            WRITE(iMPI_MyID*10+1, '(A)') 'VARIABLES = "X", "Y", "Z" ,"u" ,"v" ,"w"'
            WRITE(iMPI_MyID*10+1, '(A, I0, A, I0, A)') 'ZONE T="u_INLET_noWeight", I=', N2, ', J=', N3, ', F=POINT'
            DO K=1,N3
            DO J=1,N2  
                WRITE(iMPI_MyID*10+1,'(6(E21.14,2X))') INLET_X(J,K),INLET_Y(J,K),INLET_Z(J,K),UREF(1,J,K),UREF(2,J,K),UREF(3,J,K)
            ENDDO
            ENDDO
            CLOSE(iMPI_MyID*10+1)

            OPEN(iMPI_MyID*20+1,FILE='Tecplot_InputFiles/u_INLET_Weight.plt',STATUS='UNKNOWN') 
            WRITE(iMPI_MyID*20+1, '(A)') 'TITLE = "u_INLET_Weight"'
            WRITE(iMPI_MyID*20+1, '(A)') 'VARIABLES = "X", "Y", "Z" ,"u" ,"v" ,"w"'
            WRITE(iMPI_MyID*20+1, '(A, I0, A, I0, A)') 'ZONE T="u_INLET_Weight", I=', N2, ', J=', N3, ', F=POINT'
            DO K=1,N3
            DO J=1,N2  
                TEMP1=( UREF(1,J,K)-U_RECT(1,J,K) )*RATIO(J,K) +U_RECT(1,J,K)
                TEMP2=( UREF(2,J,K)-U_RECT(2,J,K) )*RATIO(J,K) +U_RECT(2,J,K)
                TEMP3=( UREF(3,J,K)-U_RECT(3,J,K) )*RATIO(J,K) +U_RECT(3,J,K)
                WRITE(iMPI_MyID*20+1,'(6(E21.14,2X))') INLET_X(J,K),INLET_Y(J,K),INLET_Z(J,K),TEMP1,TEMP2,TEMP3
            ENDDO
            ENDDO
            CLOSE(iMPI_MyID*20+1)
        ENDIF


    ENDDO

    CLOSE(iMPI_MyID+11)
	CLOSE(iMPI_MyID+12)
	CLOSE(iMPI_MyID+13)

    
 
230 FORMAT(A2,E21.14,A1,E21.14,A1,E21.14,A2)  
!-------------------------------------------------------------------------------
!    local variable
!-------------------------------------------------------------------------------
    
!-------------------------------------------------------------------------
!    execute statement
!------------------------------------------------------------------------- 
	

!-------------------------------------------------

    END SUBROUTINE
