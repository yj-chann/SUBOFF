!========================================!
!      CHANNEL STAGGERED BI_PARALLEL     !
!========================================!

    PROGRAM MAIN
    
!----------------------------------------
    
    INCLUDE 'head.fi'

    INTEGER MI,MJ,MJ1,MJ2
    INTEGER NIF,NIF1,NIF2
    INTEGER*8 COUNT1,COUNT2,COUNT3,COUNT_RATE,COUNT_MAX
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
    !!!NPP = NPR * NPC
    IF ( NPP .EQ. iMPI_NumProcs ) THEN

    !   the number of process is right

        CALL MPI_Barrier( MPI_COMM_WORLD, iMPI_ErrorInfo )

    !   Main Program
        MI=COORD(0)+1
        MJ=COORD(1)+1
        !WRITE(*,*) 'RENUMBERMESH is starting'
        !CALL RENUMBERMESH
        !WRITE(*,*) 'RENUMBERMESH is ok'
       
        CALL SYSTEM_CLOCK(COUNT1,COUNT_RATE,COUNT_MAX)
        CALL GET_GRID(COUNT1)

        If ( iMPI_MyID.EQ.0 ) THEN
        CALL SYSTEM_CLOCK(COUNT2,COUNT_RATE,COUNT_MAX)
        WRITE(*,*) 'GRID',(COUNT2-COUNT1)/DBLE(COUNT_RATE)
        ENDIF

        CALL MPI_Barrier( MPI_COMM_WORLD, iMPI_ErrorInfo )
        
        CALL MAIN_CODE  

        If ( iMPI_MyID.EQ.0 ) THEN
        CALL SYSTEM_CLOCK(COUNT3,COUNT_RATE,COUNT_MAX)
        WRITE(*,*) 'MAIN',(COUNT3-COUNT2)/DBLE(COUNT_RATE)
        ENDIF

        If ( iMPI_MyID.EQ.0 ) THEN
            Print *, 'FINISHED'
        ENDIF
    
    ENDIF
    
   

    CALL MPI_Barrier( MPI_COMM_WORLD, iMPI_ErrorInfo )
    CALL MPI_Finalize( iMPI_ErrorInfo )

!----------------------------------------

    END PROGRAM

!===== END PROGRAM ======================

    
    !SUBROUTINE RENUMBERMESH
    !INCLUDE 'head.fi'
    !CHARACTER*900 DIR_RENUM
    !INTEGER I
    !INTEGER RENUM0
    !
    !!DIR_RENUM='/es01/paratera/sce4049/zf/20221108/3d_renum_ref/0/'
    !!OPEN(1,FILE=TRIM(ADJUSTL(DIR_RENUM))//'s_renum',STATUS='OLD')
    !!DO I=1,22
    !!    READ(1,*)
    !!ENDDO
    !DO I=0,C_NUM_ALL-1
    !    !READ(1,*) RENUM0
    !    !RENUM(RENUM0-1)=I
    !    !RENUM2(I)=RENUM0-1
    !    RENUM(I)=I
    !    RENUM2(I)=I
    !ENDDO
    !CLOSE(1)
    !
    !END SUBROUTINE
    
    
    SUBROUTINE GET_INDEX(NUM1,I2,J2,K2,IF_CY)
    INCLUDE 'head.fi'
    INTEGER NUM1,I2,J2,K2,REST,IF_CY,TEMP1
    
    IF (NUM1.LT.(NI*NJ*NK)) THEN
    IF_CY=1    
    TEMP1=NUM1
    K2=(TEMP1-MOD(TEMP1,NI*NJ))/NI/NJ+1
    REST=TEMP1-(K2-1)*NI*NJ
    J2=(REST-MOD(REST,NI))/NI+1
    I2=REST-(J2-1)*NI+1 
    ELSE
    IF_CY=0    
    TEMP1=NUM1-NI*NJ*NK
    K2=(TEMP1-MOD(TEMP1,NI*NL))/NI/NL+1
    REST=TEMP1-(K2-1)*NI*NL
    J2=(REST-MOD(REST,NI))/NI+1
    I2=REST-(J2-1)*NI+1    
    ENDIF
    
    END SUBROUTINE
    
    SUBROUTINE GET_INDEX_POINT(NUM1,I2,J2,K2,IF_CY)
    INCLUDE 'head.fi'
    INTEGER NUM1,I2,J2,K2,REST,IF_CY,TEMP1
    
    IF (NUM1.LT.((NI+1)*(NJ+1)*NK)) THEN
    IF_CY=1    
    TEMP1=NUM1
    K2=(TEMP1-MOD(TEMP1,(NI+1)*(NJ+1)))/(NI+1)/(NJ+1)+1
    REST=TEMP1-(K2-1)*(NI+1)*(NJ+1)
    J2=(REST-MOD(REST,NI+1))/(NI+1)+1
    I2=REST-(J2-1)*(NI+1)+1 
    ELSE 
    IF_CY=0    
    TEMP1=NUM1-(NI+1)*(NJ+1)*NK
    K2=(TEMP1-MOD(TEMP1,(NI+1)*(NL-1)))/(NI+1)/(NL-1)+1
    REST=TEMP1-(K2-1)*(NI+1)*(NL-1)
    J2=(REST-MOD(REST,NI+1))/(NI+1)+1
    I2=REST-(J2-1)*(NI+1)+1 
    
    K2=K2+1
    J2=J2+1
    ENDIF
    
    END SUBROUTINE
    
    
    SUBROUTINE GET_CELLID(I,J,K,NUM1,IF_CY)
    INCLUDE 'head.fi'
    INTEGER NUM1,I,J,K,TEMP,IF_CY
    
    IF (IF_CY.EQ.1) THEN
    TEMP=NI*NJ*(K-1)+NI*(J-1)+I-1
    IF (K.EQ.0) TEMP=NI*NJ*(NK-1)+NI*(J-1)+I-1
    IF (K.EQ.(NK+1)) TEMP=NI*(J-1)+I-1
    !NUM1=RENUM(TEMP) 
    NUM1=TEMP
    ELSEIF (IF_CY.EQ.0) THEN
    TEMP=NI*NL*(K-1)+NI*(J-1)+I-1 +NI*NJ*NK
    !NUM1=RENUM(TEMP) 
    NUM1=TEMP
    ENDIF
    
    END SUBROUTINE
    
    
       
    
    
    SUBROUTINE BUBBLE_SORT(A,B)
    INTEGER N
    PARAMETER(N=7)
    INTEGER A(N),B(N)
    INTEGER I,J,TEMP
    
    DO I=1,N
        B(I)=I
    ENDDO
    DO I=N-1,1,-1
    DO J=1,I
    IF (A(J)>A(J+1)) THEN
        TEMP=A(J)
        A(J)=A(J+1)
        A(J+1)=TEMP
        TEMP=B(J)
        B(J)=B(J+1)
        B(J+1)=TEMP
    ENDIF
    ENDDO
    ENDDO
    
    END SUBROUTINE
    
    SUBROUTINE BUBBLE_SORT_POINTS_CY(A,B)
    INCLUDE 'head.fi'
    INTEGER A(P_NUM_CY),B(P_NUM_CY)
    INTEGER I,J,TEMP
    
    DO I=1,P_NUM_CY
        B(I)=I
    ENDDO
    DO I=P_NUM_CY-1,1,-1
    DO J=1,I
    IF (A(J)>A(J+1)) THEN
        TEMP=A(J)
        A(J)=A(J+1)
        A(J+1)=TEMP
        TEMP=B(J)
        B(J)=B(J+1)
        B(J+1)=TEMP
    ENDIF
    ENDDO
    ENDDO
    
    END SUBROUTINE
    
    
    SUBROUTINE BUBBLE_SORT_POINTS_RECT(A,B)
    INCLUDE 'head.fi'
    INTEGER A(P_NUM_RECT),B(P_NUM_RECT)
    INTEGER I,J,TEMP
    
    DO I=1,P_NUM_RECT
        B(I)=I
    ENDDO
    DO I=P_NUM_RECT-1,1,-1
    DO J=1,I
    IF (A(J)>A(J+1)) THEN
        TEMP=A(J)
        A(J)=A(J+1)
        A(J+1)=TEMP
        TEMP=B(J)
        B(J)=B(J+1)
        B(J+1)=TEMP
    ENDIF
    ENDDO
    ENDDO
    
    END SUBROUTINE
    
    
    SUBROUTINE CAL_FIJ_CY(F11,F12,F13,F14,F21,F22,F23,F24)
    INCLUDE 'head.fi'
    INTEGER TEMP
    INTEGER F11,F12,F13,F14,F21,F22,F23,F24
    
    TEMP=F11
    F11=SORT_CY2(TEMP)
    TEMP=F12
    F12=SORT_CY2(TEMP)
    TEMP=F13
    F13=SORT_CY2(TEMP)
    TEMP=F14
    F14=SORT_CY2(TEMP)
    TEMP=F21
    F21=SORT_CY2(TEMP)
    TEMP=F22
    F22=SORT_CY2(TEMP)
    TEMP=F23
    F23=SORT_CY2(TEMP)
    TEMP=F24
    F24=SORT_CY2(TEMP)
    
    END SUBROUTINE
    
    SUBROUTINE CAL_FIJ_RECT(F11,F12,F13,F14,F21,F22,F23,F24)
    INCLUDE 'head.fi'
    INTEGER TEMP
    INTEGER F11,F12,F13,F14,F21,F22,F23,F24
    
    TEMP=F11
    F11=SORT_RECT2(TEMP)
    TEMP=F12
    F12=SORT_RECT2(TEMP)
    TEMP=F13
    F13=SORT_RECT2(TEMP)
    TEMP=F14
    F14=SORT_RECT2(TEMP)
    TEMP=F21
    F21=SORT_RECT2(TEMP)
    TEMP=F22
    F22=SORT_RECT2(TEMP)
    TEMP=F23
    F23=SORT_RECT2(TEMP)
    TEMP=F24
    F24=SORT_RECT2(TEMP)
    
    END SUBROUTINE