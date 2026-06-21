
    PROGRAM cylinder_main
    IMPLICIT NONE
    INCLUDE 'cylinder.h'

    REAL*8 INLET_X(NALL),INLET_Y(NALL),INLET_Z(NALL)
    INTEGER*8 I,IP
    CHARACTER*200 CHAR1
 
    OPEN(1,FILE='patch_INLET_constant.obj',STATUS='OLD')   
    DO I=1,NALL
        READ(1,'(A)') CHAR1
        IP=INDEX(CHAR1,"v") 
	    CHAR1(IP:IP)=''
        READ(CHAR1,*) INLET_X(I),INLET_Y(I),INLET_Z(I)
    ENDDO 
    CLOSE(1)


    ! 用于Tecplot可视化 z方向优先 z从ZMIN->ZMAX y从YMAX->YMIN
    OPEN(UNIT=20, FILE='Tecplot_InputFiles/Turb_Inlet.plt', STATUS='UNKNOWN', FORM='FORMATTED')
    WRITE(20, '(A)') 'TITLE = "Turb_Inlet"'
    WRITE(20, '(A)') 'VARIABLES = "X", "Y", "Z"'
    WRITE(20, '(A, I0, A, I0, A)') 'ZONE T="Turb_Inlet", I=', N, ', J=', N, ', F=POINT'
    DO I = 1, NALL
        WRITE(20, '(3(ES22.15,2X))') INLET_X(I),INLET_Y(I),INLET_Z(I)
    END DO
    CLOSE(20)
    
    
    
    OPEN(1,FILE='INLET1.plt',STATUS='UNKNOWN')   
    DO I=1,NALL 
        WRITE(1,230) INLET_X(I),INLET_Y(I),INLET_Z(I)
    ENDDO 
    CLOSE(1)
    


230 FORMAT(3(E22.15,2X))
    
    
     

    END PROGRAM


