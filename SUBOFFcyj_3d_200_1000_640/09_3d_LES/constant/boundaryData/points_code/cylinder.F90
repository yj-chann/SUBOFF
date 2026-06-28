
    PROGRAM cylinder_main
    IMPLICIT NONE
    INCLUDE 'parameter.h'

   
    REAL*8 INLET_X(NALL),INLET_Y(NALL),INLET_Z(NALL)
    INTEGER I
 
    
      
    

    OPEN(1,FILE='ReadData/INLET1.plt',STATUS='OLD')   
    DO I=1,NALL
        READ(1,*) INLET_X(I),INLET_Y(I),INLET_Z(I)
    ENDDO 
    CLOSE(1)
    

    OPEN(11,FILE='../INLET1/points',STATUS='UNKNOWN')
    WRITE(11,'(I0)') NALL
    WRITE(11,*) '('
    DO I=1,NALL
        WRITE(11,230) '( ',INLET_X(I),' ',INLET_Y(I),' ',INLET_Z(I),' )'
    ENDDO  
    WRITE(11,*) ')'
    CLOSE(11) 

     
230 FORMAT(A2,E22.15,A1,E22.15,A1,E22.15,A2)      
    

    END PROGRAM


