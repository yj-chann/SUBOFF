MODULE TimeFormatterModule
    IMPLICIT NONE

CONTAINS

    FUNCTION FormatTime(ST, DT, I) RESULT(CHAR1)
        ! Inputs
        REAL(8), INTENT(IN) :: ST
        REAL(8), INTENT(IN) :: DT
        INTEGER, INTENT(IN) :: I
        
        ! Output
        CHARACTER(LEN=20)   :: CHAR1
        
        ! Local Variables
        INTEGER :: N1, N2, N3, IP
        REAL(8) :: TIME

        ! Calculate frequency
        N1 = NINT(1.0_8 / DT)
        N2 = MOD(I, N1)

        IF (N2 == 0) THEN
            ! Whole number step: Format as Integer
            N3 = I / N1
            WRITE(CHAR1, '(I0)') N3
            CHAR1 = ADJUSTL(CHAR1) ! Left-justify the string
            
        ELSE
            ! Fractional step: Format as Float
            TIME = ST + DT * I 
            WRITE(CHAR1, '(F10.7)') TIME
            CHAR1 = ADJUSTL(CHAR1)
            
            ! Safely strip trailing zeros
            DO
                IP = LEN_TRIM(CHAR1)
                IF (IP <= 1) EXIT  ! Stop if string is almost empty
                
                IF (CHAR1(IP:IP) == '0') THEN
                    CHAR1(IP:IP) = ' '  ! Overwrite trailing zero with space
                ELSE
                    EXIT
                END IF
            END DO
            
        END IF
        
    END FUNCTION FormatTime

END MODULE TimeFormatterModule




PROGRAM MAIN
    USE TimeFormatterModule
    IMPLICIT NONE
    include 'parameter.h'

    INTEGER PNJ_ID,PNI_ID,PNK_ID
    INTEGER I,J,K,I_GLOBAL,J_GLOBAL,K_GLOBAL
    CHARACTER(LEN=100) DIR,filename1,filename2,line
    CHARACTER(LEN=20) CHAR1
    REAL UX(NI,NJ+NL/2,NK),UY(NI,NJ+NL/2,NK),UZ(NI,NJ+NL/2,NK),p(NI,NJ+NL/2,NK)
    REAL X(NI,NJ+NL/2,NK),Y(NI,NJ+NL/2,NK),Z(NI,NJ+NL/2,NK)
    REAL X0(NI+1,NJ+NL/2+1),Y0(NI+1,NJ+NL/2+1)
    REAL, PARAMETER :: PI = 4.0 * ATAN(1.0)
    REAL :: THETA(NK)
    REAL :: DTHETA
    REAL :: XC, RC
    INTEGER :: I_TIME,line_cnt,ios

    DIR='../09_3d_startup/hidden.procs/processor'

    UX=0.0; UY=0.0; UZ=0.0;


    DTHETA = PI / DBLE(NK)

    DO K = 1, NK
        THETA(K) = -PI/4.0 + DBLE(K) * DTHETA
    END DO


    OPEN(1,FILE='ReadData/suboff_mesh_2d.plt',STATUS='OLD')   
    READ(1,*)
    DO I=1,NI+1 
    DO J=1,NJ+NL/2+1
        READ(1,*) X0(I,J),Y0(I,J)
    ENDDO
    ENDDO  
    CLOSE(1)


    DO K = 1, NK
        DO J = 1, NJ+NL/2
            DO I = 1, NI

                XC = 0.25 * (X0(I, J)   + X0(I+1, J) + &
                             X0(I, J+1) + X0(I+1, J+1))
                             
                RC = 0.25 * (Y0(I, J)   + Y0(I+1, J) + &
                             Y0(I, J+1) + Y0(I+1, J+1))

                X(I, J, K) = XC
                Y(I, J, K) = RC * COS(THETA(K))
                Z(I, J, K) = RC * SIN(THETA(K))
                
            END DO
        END DO
    END DO

    DO I_TIME = 1,NT


        ! READ CY2+CY3
        DO PNJ_ID=PNJ1+1,PNJ
            DO PNI_ID=1,PNI
                WRITE(CHAR1,'(I0)') (PNJ_ID-1)*PNI+PNI_ID-1
                
                filename1=TRIM(ADJUSTL(DIR))//TRIM(ADJUSTL(CHAR1))//'/'//TRIM(ADJUSTL(FormatTime(ST, DT, I_TIME)))//'/U'
                filename2=TRIM(ADJUSTL(DIR))//TRIM(ADJUSTL(CHAR1))//'/'//TRIM(ADJUSTL(FormatTime(ST, DT, I_TIME)))//'/p'
                ! print* ,filename
                OPEN(10,FILE=filename1,STATUS='OLD',ACTION='READ', IOSTAT=ios)
                OPEN(20,FILE=filename2,STATUS='OLD',ACTION='READ', IOSTAT=ios)
                DO line_cnt=1,23
                    READ(10,*)
                    READ(20,*)
                ENDDO
                DO K=1,NK
                    DO J=1,CNJ
                        DO I=1,CNI
                            I_GLOBAL=(PNI_ID-1)*CNI+I
                            J_GLOBAL=NL/2+(PNJ_ID-1)*CNJ+J
                            K_GLOBAL=K
                            DO line_cnt = 1, LEN_TRIM(line)
                                IF (line(line_cnt:line_cnt) == '(' .OR. line(i:i) == ')') THEN
                                    line(line_cnt:line_cnt) = ' '
                                END IF
                            END DO
                            READ(line, *, IOSTAT=ios) UX(I_GLOBAL,J_GLOBAL,K_GLOBAL),UY(I_GLOBAL,J_GLOBAL,K_GLOBAL),UZ(I_GLOBAL,J_GLOBAL,K_GLOBAL)
                            READ(20,*) p(I_GLOBAL,J_GLOBAL,K_GLOBAL)
                        ENDDO
                    ENDDO
                ENDDO
                CLOSE(10)
                CLOSE(20)


            ENDDO
        ENDDO

            WRITE(CHAR1,'(I6.6)') I_TIME
            filename1='Tecplot_InputFiles/PU/PU'//TRIM(ADJUSTL(CHAR1))//'.plt'
            OPEN(UNIT=40, FILE=filename1, STATUS='REPLACE', ACTION='WRITE', IOSTAT=IOS)
            IF (IOS /= 0) THEN
                PRINT *, "Error: Could not open file ", TRIM(filename1), " for writing."
                STOP
            END IF
            WRITE(40, '(A)') 'TITLE = "3D Flow Field Data"'
            WRITE(40, '(A)') 'VARIABLES = "X", "Y", "Z", "UX", "UY", "UZ", "P"'
            
            ! Note: For F=POINT, Tecplot expects I to vary fastest, then J, then K.
            WRITE(40, '(A, I0, A, I0, A, I0, A)') &
                'ZONE T="3D Flow Field Data", I=', NI, ', J=', NL/2+NJ, ', K=', NK, ', F=POINT'
            DO K = 1, NK
                DO J = 1, NL/2+NJ
                    DO I = 1, NI
                        ! ES15.6 ensures scientific notation with consistent padding
                        WRITE(40, '(7(E22.15, 2X))') &
                            X(I,J,K), Y(I,J,K), Z(I,J,K), &
                            UX(I,J,K), UY(I,J,K), UZ(I,J,K), P(I,J,K)
                    END DO
                END DO
            END DO
            CLOSE(40)

            PRINT *, "Tecplot file successfully written to: ", TRIM(filename1), ' TIME: ', ST+I_TIME*DT
    ENDDO

    
      
    END PROGRAM MAIN




