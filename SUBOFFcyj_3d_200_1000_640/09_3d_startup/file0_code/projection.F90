!--------------------------------------------------------------
!     GRID GENERATION IN CHANNEL FLOW
!     X2 Clustering only
!     JICHOI
!--------------------------------------------------------------

      SUBROUTINE PROJECTION(X1,Y1,ILOC1,JLOC1)
      INCLUDE 'head.fi'
      REAL TIME,T1,T2
      REAL X1,Y1,DOT(REF_N2),L1(REF_N2),DIST(REF_N1)
      REAL TEMP1,TEMP2,TEMP3,TEMP4,A1,A2,B1,B2
      INTEGER I,J,ILOC(1),JLOC(1),K,I1,I2,J1,J2,II,JJ
      INTEGER ILOC1,JLOC1
      
      
      
      
          
      DO K=1,REF_N2
      L1(K)=SQRT((X1-X00_REF00(K))**2+(Y1-Y00_REF00(K))**2)   
      DOT(K)=ABS( (X1-X00_REF00(K))*(X0_REF0(REF_N1+1,K+1)-X0_REF0(REF_N1+1,K))+(Y1-Y00_REF00(K))*(Y0_REF0(REF_N1+1,K+1)-Y0_REF0(REF_N1+1,K))  )/L1(K)/L2_REF(K)
      ENDDO
      JLOC=MINLOC(DOT)
      DO K=1,REF_N1
      DIST(K)=SQRT( (X1-X0_REF0(K,JLOC(1)))**2+(Y1-Y0_REF0(K,JLOC(1)))**2 )
      ENDDO
      ILOC=MINLOC(DIST)
      
      DO II=ILOC(1)-2,ILOC(1)+2
      IF (II.GE.1 .AND. II.LE.REF_N1) THEN
      DO JJ=JLOC(1)-2,JLOC(1)+2
      IF (JJ.GE.1 .AND. JJ.LE.REF_N2) THEN
      I1=II
      I2=II+1
      J1=JJ
      J2=JJ+1
      A1=X1-X0_REF0(I1,J1)
	  A2=Y1-Y0_REF0(I1,J1)
	  B1=X0_REF0(I1,J2)-X0_REF0(I1,J1)
	  B2=Y0_REF0(I1,J2)-Y0_REF0(I1,J1)
  	  TEMP1=A1*B2-A2*B1   
      A1=X1-X0_REF0(I2,J1)
	  A2=Y1-Y0_REF0(I2,J1)
	  B1=X0_REF0(I1,J1)-X0_REF0(I2,J1)
	  B2=Y0_REF0(I1,J1)-Y0_REF0(I2,J1)
  	  TEMP2=A1*B2-A2*B1
      A1=X1-X0_REF0(I2,J2)
	  A2=Y1-Y0_REF0(I2,J2)
	  B1=X0_REF0(I2,J1)-X0_REF0(I2,J2)
	  B2=Y0_REF0(I2,J1)-Y0_REF0(I2,J2)
  	  TEMP3=A1*B2-A2*B1
      A1=X1-X0_REF0(I1,J2)
	  A2=Y1-Y0_REF0(I1,J2)
	  B1=X0_REF0(I2,J2)-X0_REF0(I1,J2)
	  B2=Y0_REF0(I2,J2)-Y0_REF0(I1,J2)
  	  TEMP4=A1*B2-A2*B1  
      IF (TEMP1.GE.0.0 .AND. TEMP2.GE.0.0 .AND. TEMP3.GE.0.0 .AND. TEMP4.GE.0.0) GOTO 101
      ENDIF
      ENDDO
      ENDIF
      ENDDO
101   ILOC1=II
      JLOC1=JJ
      
      
      
212   FORMAT('Zone i= ',I0,'j= ',I0,' f=point')  
213   FORMAT(5(E21.14,2X))  
214   FORMAT(4(I0,2X)) 
      
      END SUBROUTINE
