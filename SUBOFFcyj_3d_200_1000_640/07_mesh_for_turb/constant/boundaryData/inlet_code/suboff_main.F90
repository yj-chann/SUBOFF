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
    REAL*8 X,Y,Z,PAI
    REAL*8 A1,A2,A3,A4,A5,A6,A7,UX,UY,UZ
    INTEGER I,J,K,N0,IP,STI,EDI,NT0,NN1,NN2,NN3,TI
    CHARACTER*400 CHAR1,DIR
    CHARACTER*20 CHAR3
    CHARACTER*1 CHAR2
    
    
    PAI=ACOS(-1.0)
    DIR='../../../../06_tg/vel/'
    
    !u000001.plt同时用于0文件和DT文件
    IF (iMPI_MyID.EQ.0) THEN
        OPEN(11,FILE=TRIM(ADJUSTL(DIR))//'u000001.plt',STATUS='OLD')   
        OPEN(1,FILE='../INLET/0/U',STATUS='UNKNOWN')
        WRITE(1,'(I0)') NALL
        WRITE(1,*) '('
        DO I=1,NALL
        READ(11,*) UX,UY,UZ
        WRITE(1,230) '( ',UX,' ',UY,' ',UZ,' )'
        ENDDO
        WRITE(1,*) ')'
        CLOSE(1)
        CLOSE(11)
    ENDIF

    NT0=NT/NPP
    STI=iMPI_MyID*NT0+ST
    EDI=( iMPI_MyID+1 )*NT0+ST-1
    NN1=NINT(1.0/DT)
    DO I=STI,EDI
        NN2=MOD(I,NN1)
        IF (NN2.EQ.0) THEN
        NN3=I/NN1
        WRITE(CHAR1,'(I0)') NN3 ! 整数时间，如 "1", "2"
        ELSE 
        TIME=DT*I ! 不管ST取多少，永远从0时刻开始
        WRITE(CHAR1,'(F10.7)') TIME ! 浮点时间，如 "0.4000000"
        !remove all the '0's in the end of the string. eg: 0.4000000 -> 0.4
505     IP=INDEX(CHAR1,'0',BACK=.TRUE.)             
        CHAR1(IP:IP)=''     
	    CHAR2=CHAR1(IP-1:IP-1)  
        READ(CHAR2,'(I1)') N0
        IF (N0.EQ.0) GOTO 505
        ENDIF
        TI=I-ST+1
        WRITE(CHAR3,'(I6.6)') TI

        OPEN(iMPI_MyID*NT0+11,FILE=TRIM(ADJUSTL(DIR))//'u'//TRIM(ADJUSTL(CHAR3))//'.plt',STATUS='OLD')  
        OPEN(iMPI_MyID*NT0+1,FILE='../INLET/'//TRIM(ADJUSTL(CHAR1))//'/U',STATUS='UNKNOWN')
        WRITE(iMPI_MyID*NT0+1,'(I0)') NALL
        WRITE(iMPI_MyID*NT0+1,*) '('
        DO J=1,NALL
        READ(iMPI_MyID*NT0+11,*) UX,UY,UZ
        WRITE(iMPI_MyID*NT0+1,230) '( ',UX,' ',UY,' ',UZ,' )'
        ENDDO
        WRITE(iMPI_MyID*NT0+1,*) ')'
        CLOSE(iMPI_MyID*NT0+1)
        CLOSE(iMPI_MyID*NT0+11)
    ENDDO


230 FORMAT(A2,E21.14,A1,E21.14,A1,E21.14,A2) 
	

!-------------------------------------------------

    END SUBROUTINE
