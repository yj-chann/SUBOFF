    PROGRAM MKDIR_MAIN
    IMPLICIT NONE
    INCLUDE 'parameter.h'
    
    REAL TIME
    INTEGER I,N0,IP,N1,N2,N3,ST
    CHARACTER*200 CHAR1
    CHARACTER*1 CHAR2
    
    N1=NINT(1.0/DT)
    WRITE(*,'(I0,A)') N1 ,' steps per second'
    OPEN(1,FILE='ofrun.sh',STATUS='UNKNOWN')
    ST=0
    DO I=ST,NT+ST
        N2=MOD(I,N1)
        IF (N2.EQ.0) THEN
        N3=I/N1
        WRITE(CHAR1,'(I0)') N3
        ELSE 
        TIME=DT*I
        WRITE(CHAR1,'(F10.7)') TIME
505     IP=INDEX(CHAR1,'0',BACK=.TRUE.) 
        CHAR1(IP:IP)=''
	    CHAR2=CHAR1(IP-1:IP-1)
        READ(CHAR2,'(I1)') N0
        IF (N0.EQ.0) GOTO 505
        ENDIF
        WRITE(1,*) 'mkdir -p ../INLET/'//TRIM(ADJUSTL(CHAR1))
    ENDDO
    CLOSE(1)
    
    END
