program generate_ofrun
    implicit none

    include "parameter.h"

    integer :: unit_number = 10
    integer :: NN2 , I, NN1,NN3, N0, IALL, IP,NNN2
    REAL*8 :: TIME0,TIME1,TCycle
    character(len=20) :: filename = "ofrun.sh"
    character(len=20) :: CHAR1, CHAR2, CHAR3, CHAR4

    ! Open the file. 'replace' creates it or overwrites it if it exists.
    open(unit=unit_number, file=filename, status='replace', action='write')

    ! Write the bash shebang and the copy command
    write(unit_number, '(A)') '#!/bin/bash'
    write(unit_number, '(A)') 'cp -r ../../../../08_massflow/constant/boundaryData/INLET/0 ../INLET1/0'
    NN1=NINT(1.0/DT)
    NN2=MOD(NT,NN1)
    IF (NN2.EQ.0) THEN
    NN3=NT/NN1
    WRITE(CHAR1,'(I0)') NN3 ! 整数时间，如 "1", "2"
    ELSE 
    TIME0=DT*NNN2
    WRITE(CHAR1,'(F10.7)') TIME0 ! 浮点时间，如 "0.4000000"
    !remove all the '0's in the end of the string. eg: 0.4000000 -> 0.4
504     IP=INDEX(CHAR1,'0',BACK=.TRUE.)             
    CHAR1(IP:IP)=''     
    CHAR2=CHAR1(IP-1:IP-1)  
    READ(CHAR2,'(I1)') N0
    IF (N0.EQ.0) GOTO 504
    ENDIF
    write(unit_number, '(A)') 'cp -r ../../../../08_massflow/constant/boundaryData/INLET/0 ../INLET1/'//trim(ADJUSTL(CHAR1))


!     TCycle=DT*NT
!     IAll=NT*LOOP_N

!     DO I=1, IALL

!     NN1=NINT(1.0/DT)
!     IF (MOD(I,NT).EQ.0) THEN
!         NNN2=NT
!     ELSE
!         NNN2=MOD(I,NT) ! 映射到1~NT
!     ENDIF
!     NN2=MOD(NNN2,NN1)
!     IF (NN2.EQ.0) THEN
!     NN3=NNN2/NN1
!     WRITE(CHAR1,'(I0)') NN3 ! 整数时间，如 "1", "2"
!     ELSE 
!     TIME0=DT*NNN2
!     WRITE(CHAR1,'(F10.7)') TIME0 ! 浮点时间，如 "0.4000000"
!     !remove all the '0's in the end of the string. eg: 0.4000000 -> 0.4
! 505     IP=INDEX(CHAR1,'0',BACK=.TRUE.)             
!     CHAR1(IP:IP)=''     
!     CHAR2=CHAR1(IP-1:IP-1)  
!     READ(CHAR2,'(I1)') N0
!     IF (N0.EQ.0) GOTO 505
!     ENDIF

!     NN2=MOD(I+NT,NN1)
!     IF (NN2.EQ.0) THEN
!     NN3=(I+NT)/NN1
!     WRITE(CHAR3,'(I0)') NN3 ! 整数时间，如 "1", "2"
!     ELSE 
!     TIME1=DT*I+TCycle
!     WRITE(CHAR3,'(F10.7)') TIME1 ! 浮点时间，如 "0.4000000"
!     !remove all the '0's in the end of the string. eg: 0.4000000 -> 0.4
! 506     IP=INDEX(CHAR3,'0',BACK=.TRUE.)             
!     CHAR3(IP:IP)=''     
!     CHAR4=CHAR3(IP-1:IP-1)  
!     READ(CHAR4,'(I1)') N0
!     IF (N0.EQ.0) GOTO 506
!     ENDIF


!     write(unit_number, '(A)') 'cp -r ../../../../08_massflow/constant/boundaryData/INLET/'//trim(ADJUSTL(CHAR1))//' ../INLET1/'//trim(ADJUSTL(CHAR3))

!     ENDDO

!     ! Close the file unit
!     close(unit_number)

    print *, "Successfully generated: ", trim(filename)

end program generate_ofrun