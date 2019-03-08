      PROGRAM PHOTON2D
      INCLUDE 'COMMON2D.F'
      DOUBLE PRECISION ALARM
      DOUBLE PRECISION XPOINTING(IMAXX,JMAXX), YPOINTING(IMAXX,JMAXX)
      DOUBLE PRECISION RF(JMAXX), TF(JMAXX)
      OPEN(80,FILE='PHOTON2D.IN',STATUS='OLD')
      OPEN(81,FILE='Parameter',STATUS='UNKNOWN')

 501  FORMAT(10000E14.6)
99991 FORMAT(1400(1x,e10.4))
99992 format(i3.3,1x,i3.3,1x,E10.4)
99995 FORMAT(<IMAXX>(1x,e10.4))


 502  FORMAT('************************************************',/
     &'SOLVING THE ',I2,' X',I2,' SUPERCELL PROBLEM',/
     &'EACH CELL HAS ',I3,' X',I3,' GRIDS',/
     &'TOTAL GRIDS:  ',I4,' X',I4,' =',I10,/
     &'PERIODIC DIELECTRIC CONTRAST AND FILLING RATIO',/
     &2F10.5,/
     &'DISORDERS',/
     &5F10.3,/
     &'TOTAL NUMBER OF POINTS TO READ',/
     &I10,/
     &'************************************************')
      
      READ(80,*) IPOLAR
      READ(80,*) IOVERLAP
      READ(80,*) NBANDS
      READ(80,*) NXGRID
      READ(80,*) NX
      READ(80,*) NY
      READ(80,*) A
      READ(80,*) WAVE
      READ(80,*) AMPLITUDE
      READ(80,*) XLENGTH
      READ(80,*) YLENGTH
      READ(80,*) BACKROUND
      READ(80,*) EPS1
      READ(80,*) COND0
      READ(80,*) N0LEVEL0
      READ(80,*) TSTART
      READ(80,*) SC

      LL=10

      PI=4.0*DATAN(1.0D0)
      C0=2.9979D8
      EPS0=8.854D-12
      XMU0=4.0*PI*1.0D-7
      hbar=1.0545887d-34
      electrone=1.6021892d-19
      electronm=9.109534d-31

      
cxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
C	W0=2*PI*3.E15
        W0=C0*2*PI/WAVE
	WI=C0*2*PI/WAVE
	PUMPRATE0=1.0E15
c        PUMPRATE0=0

	
cXXXXXXXXXXXXXXXXXXXXXXXXX

C------TYPICAL PARAMETERS OF ATOMIC SYSTEM DEFINED
CC DELTAW=W0/30 WHEN WAVELENGTH = 5.0E-6 (M)
	T2=1.0E-15

	ttrans31=1.d-3
	ttrans32=1.d-17

	ttrans2=1.0d-15
	ttrans21=1.0d-15
	ttrans10=1.d-17

	ratereal=1.d0/ttrans2
	rateclassic=(electrone**2/electronm)*w0**2/6/Pi/eps0/c0**3
	ki=ratereal/rateclassic*(electrone**2/electronm)
	deltaw= ratereal+1.d0/T2+1.d0/ttrans21
	
        NLINE=NXGRID/2
        NWHITE=4*NLINE
        write(*,*)'nxgrid=',NXGRID
      DO 402 J=1,JMAXX
         HINCR(J)=0.0D0
 402    continue 

Cxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
      ALARM=2.0D0*NLINE
      IF(ALARM.LT.12.0D0.AND.ALARM.GE.6.0D0)THEN
      PRINT*, '*******************************************'
      PRINT*, 'WARNING: FEW GRID POINTS->PROGRAM CONTINUES'
      PRINT*, '*******************************************' 
      ELSEIF(ALARM.LT.6.0D0)THEN
      PRINT*, '*******************************************'
      PRINT*, 'ERROR: VERY FEW GRID POINTS->PROGRAM STOPS' 
      PRINT*, '*******************************************'
      STOP
      ENDIF

	LENGTHX=INT(XLENGTH/A*2*NLINE)
	LENGTHY=INT(YLENGTH/A*2*NLINE)
	LENGTHY=LENGTHX
	RADIUS=XLENGTH/A*2*NLINE
	KDIV=1
C NSTEPS is the total number of time steps
	NSTEPS=(NFFTSTEPS+100)*KDIV
	IKS=1
      
        write (*, *) "OK at 205 line"
      CALL SETUP   

      WRITE(*,*) 'SOLVING THE ',NX,' X',NY,' SUPERCELL PROBLEM'
      WRITE(*,*) 'EACH CELL HAS ',2*NLINE,' X',2*NLINE,' GRIDS'
      WRITE(*,*) 'TOTAL GRIDS:  ',IMAX,' X',JMAX,' =',IMAX*JMAX
      WRITE(*,*) 'DT=',DT,'    DX=',DX
      
	WRITE(1,*)  "SOLVING THE ",NX," X",NY, " SUPERCELL PROBLEM"
      WRITE(1,*) "EACH CELL HAS ",2*NLINE," X",2*NLINE," GRIDS"
      WRITE(1,*) 'TOTAL GRIDS:  ',IMAX,' X',JMAX,' =',IMAX*JMAX
      write(1,*)  'RADIUS=',RADIUS/2/NLINE
	WRITE(1,*) 'DT= ', DT, '    DX= ', DX

 
      NCOUNT=0
      DO 6 JJ = 1, JMAX
         RF(JJ) = 0.0D0
         TF(JJ) = 0.0D0
6     CONTINUE

C NSTEPS is the total number of time steps
      NSTEPS=(NFFTSTEPS+100)
 
      DO 10 IT=1,NSTEPS
         CALL MARCH_FIELDS(IT)
C 	 CALL COLLECT_DATA(IT,NCOUNT)
CXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
     

      IF(MOD(IT,LL).EQ.0)THEN
         write(*,*)
         write(*,*)
         write(*,*)IT
      ENDIF




   
      IF(MOD(IT,LL).EQ.0)THEN

      REWIND(91)

      DO 200 JJ=1,jMAX
      DO 200 II=1,IMAX
      If(IPOLAR.EQ.0)THEN
	   FRR(II,JJ)=F4R(II,JJ)
	   IF(DABS(F1R(II,jJ)).LT.1.0D-100)THEN
	      F1R(II,jJ)=0.0D0
	   ENDIF
	   IF(DABS(FRR(II,jJ)).LT.1.0D-100)THEN
	      FRR(II,jJ)=0.0D0
	   ENDIF
	ENDIF
200	CONTINUE

C---------RECORD FIELD COMPONENTS AND POLARIZATION
      DO 201 JJ=1,jMAX/SC
        WRITE(91,99995)(F1R(SC*II,SC*jJ),II=1,IMAX/SC)     
c        WRITE(94,99995)(FRR(SC*II,SC*jJ)
c     &           /F1R(SC*II,SC*jJ)/EPS0,II=1,IMAX/SC)
                          

201   continue
     
     
        IF(MOD(IT,LL).EQ.0)THEN
C        DO 1993 II=1,IMAX 
         II=IMAX-40
        write(191,*) (F1R(II,JMAX/2))  
         II=400
        write(1911,*) (F1R(II,JMAX/2)) 
C        write(201,*) (N1LEVEL(II,JMAX/2))
C        write(202,*) (N2LEVEL(II,JMAX/2))
        write(204,*) (DELTAN(IMAX/2,JMAX/2))             
C 1993    continue 

 

C          CLOSE(94)
        ENDIF
CCCCC field distribution output CCCC.AND.IT.GE.11000
        IF(MOD(IT,2000).EQ.0)THEN
        DO 444 JJ=1,JMAX
              DO 445 II=1,IMAX
              write(203,*)(F1R(II,JJ))
445           continue
444     continue
        ENDIF
CCCCC

      CLOSE(91)
c      CLOSE(94)

         
c        WRITE(94,99995)(FRR(SC*II,SC*jJ)
c     &           /F1R(SC*II,SC*jJ)/EPS0,II=1,IMAX/SC)

c-----------------------------------------------------------------------
      ENDIF
CXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX


 10   CONTINUE
 

      STOP
      END
*                                                                *
******************************************************************
*                                                                *




*                                                                *
******************************************************************
*                                                                *
	SUBROUTINE SETUP
	INCLUDE 'COMMON2D.F'
	DOUBLE PRECISION EPSIL(IMAXX,JMAXX)
	DOUBLE PRECISION xc,yc,dist1,esum,JZ1,JZ2
	DOUBLE PRECISION ICYL(2,NXMAX),ICYLP(2,NXMAX)
	DOUBLE PRECISION COND(IMAXX, JMAXX)

C-------------------------
C     define constants
C-------------------------
	IF (NX*NY.GT.NXMAX)THEN
		WRITE (*, *) "TOO FEW CELLS AVAILABLE, GIVE A BIGGER NXMAX"
	stop
	ENDIF
	DX=A/(2*NLINE)
	DT=DX/C0/DSQRT(2.0D0)/1.2
	WRITE(*,*)'DX=',DX,NWHITE
C  IF THE STRUCTURE IS SQUARE
C--------------------------------------
	IMAX=NX*2*NLINE+32*2*NLINE+60
	JMAX=NY*2*NLINE+40
C--------------------------------------
	ISTAR=30
c	IENDO=3*NWHITE+2
C	JSTAR=JMAX/2-15*NLINE
C       JENDO=JMAX/2+15*NLINE
        JSTAR=0
	JENDO=JMAX

	IF(IMAX.GT.IMAXX.OR.JMAX.GT.JMAXX)THEN 
	PRINT*,'ERROR: TOO MANY GRID POINTS->PROGRAM STOPS'
	STOP
	ENDIF

 501	format(2000I2)
 505	FORMAT(2000E14.6)

CXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
CCCC	IF(IPOLAR.EQ.0)THEN
C-----------------------------------------------------------------
       IC=0
       DO 20 IY=1,NY
       DO 20 IX=1,NX
          IC=IC+1
          ICYL(1,IC)=NLINE+32*2*NLINE+(IX-1)*2*NLINE
          ICYL(2,IC)=30+(IY-1)*2*NLINE
20     CONTINUE
       NTOTAL=IC
C------------------------------------------------------------------------------
	DO 29 J=1,JMAX
	DO 29 I=1,IMAX
           EPSIL(I,J)=1.0D0
           COND(I,J)=0
           N0LEVEL(I,J)=0.0D0
           N1LEVEL(I,J)=0.0D0
           N2LEVEL(I,J)=0.0D0
           N3LEVEL(I,J)=0.0D0	   
29	CONTINUE

CXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
C        DO 30 I=32*2*NLINE,IMAX-60
C	DO 30 J=20,JMAX-20
        DO 30 I=IMAX/2-100,IMAX/2+100
        DO 30 J=JMAX/2-100,JMAX/2+100
C        DO 30 J=0,JMAX
	  EPSIL(I,J)=BACKROUND
	  N0LEVEL(I,J)=N0LEVEL0
Cxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
C          DO 134 ICK=1,NTotal
C            XC=DFLOAT(I)
C            YC=DFLOAT(J)
C            DIST1=DSQRT((XC-ICYL(1,ICK))**2+(YC-ICYL(2,ICK))**2)   
C            IF(DIST1.LE.RADIUS)THEN
C               EPSIL(I,J)=EPS1
C	       N0LEVEL(I,J)=0.0D0
C            ENDIF
C 134       CONTINUE

30	CONTINUE
  
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C	COEFFICIENTS FOR ABSORPTION BOUNDARY CONDITION

	DO 200 JJ=1,jMAX/SC
	WRITE(88,505) (N0LEVEL(SC*ii,SC*jj),ii=1,iMAX/SC)	
 200	WRITE(99,505) (EPSIL(SC*ii,SC*jj),ii=1,iMAX/SC)
	close(88)
	close(99)  
        
C       STOP
Cxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx


      IF(IPOLAR.EQ.0)THEN
C EZ polarization; EPSIL and EZ are defined on the same grid 
C points. Use CE for updating EZ.
C------------------------------------------------------------------

      DO 50 J=1,JMAX
      DO 50 I=1,IMAX
         NSOURCE(I,J)=0
CCCCCCC***************************************************
         CE1(I,J)=(EPSIL(I,J)*EPS0/DT-COND(I,J)/2)/
     &            (EPSIL(I,J)*EPS0/DT+COND(I,J)/2)
         CE2(I,J)=1/(EPSIL(I,J)*EPS0/DT+COND(I,J)/2)/DX
         CE3(I,J)=-1/(EPSIL(I,J)*EPS0/DT+COND(I,J)/2)/DT

CCCCCCC***************************************************
50	continue 
C------------------------------------------------------------------
      ELSEIF(IPOLAR.EQ.1)THEN
C------------------------------------------------------------------

      ENDIF
C------------------------------------------------------------------
      
      CH=DT/(XMU0*DX)
      
CXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
C-----GAUSSION SOURCE LOCATION DEFINED

	DO 65 J=1,JMAX
	DO 65 I=1,IMAX
	If(I.EQ.(ISTAR+5).AND.J.GT.JSTAR.AND.J.LT.JENDO)THEN
             NSOURCE(I,J)=1
	  ELSE
	     NSOURCE(I,J)=0
	  ENDIF
65	CONTINUE
c-------------------------------------------------------------
C     PARAMETERS IN THE FIELD AND ATOMIC SYSTEM INITIATED
	CPZ1 = 2.0D0*KI*DT**2/(2+DELTAW*DT)
	CPZ2 = (4-2*(W0*DT)**2)/(2+DELTAW*DT)
	CPZ3 = -(2-DELTAW*DT)/(2+DELTAW*DT)
      CN1LEVEL1=1/(2D0*TTRANS31)/(1/DT+1/TTRANS10/2D0)
      CN1LEVEL2=-1/(HBAR*W0*DT)/(1/DT+1/TTRANS10/2D0)
      CN1LEVEL3=1/(2D0*TTRANS21)/(1/DT+1/TTRANS10/2D0)
      CN1LEVEL4=(1/DT-1/TTRANS10/2D0)/(1/DT+1/TTRANS10/2D0)

      CN2LEVEL1=1/(2D0*TTRANS32)/(1/DT+1/TTRANS21/2D0)
      CN2LEVEL2=1/(HBAR*W0*DT)/(1/DT+1/TTRANS21/2D0)
      CN2LEVEL3=(1/DT-1/TTRANS21/2D0)/(1/DT+1/TTRANS21/2D0)

      CN3LEVEL1=1/(1/DT+1/TTRANS32/2D0+1/TTRANS31/2D0)
      CN3LEVEL2=(1/DT-1/TTRANS32/2D0-1/TTRANS31/2D0)
     &/(1/DT+1/TTRANS32/2D0+1/TTRANS31/2D0)
        
      CN0LEVEL1=1
      CN0LEVEL2=-DT
      CN0LEVEL3=DT/TTRANS10/2	
c-------------------------------------------------------------

      SX=C0*DT/DX
      SY=C0*DT/DX
      TX(1,1)=0.5D0*(2-SX)*(1-SX)
      TX(1,2)=SX*(2-SX)
      TX(1,3)=0.5D0*SX*(SX-1)
      TX(1,4)=0.0D0
      TX(1,5)=0.0D0
      TX(1,6)=0.0D0
      TX(1,7)=0.0D0
      TX(2,1)=TX(1,1)**2
      TX(2,2)=2*TX(1,1)*TX(1,2)
      TX(2,3)=2*TX(1,1)*TX(1,3)+TX(1,2)**2
      TX(2,4)=2*TX(1,2)*TX(1,3)
      TX(2,5)=TX(1,3)**2
      TX(2,6)=0.0D0
      TX(2,7)=0.0D0
      TX(3,1)=TX(1,1)*TX(2,1)
      TX(3,2)=TX(1,1)*TX(2,2)+TX(1,2)*TX(2,1)
      TX(3,3)=TX(1,1)*TX(2,3)+TX(1,2)*TX(2,2)+TX(1,3)*TX(2,1)
      TX(3,4)=TX(1,1)*TX(2,4)+TX(1,2)*TX(2,3)+TX(1,3)*TX(2,2)
      TX(3,5)=TX(1,1)*TX(2,5)+TX(1,2)*TX(2,4)+TX(1,3)*TX(2,3)
      TX(3,6)=                 TX(1,2)*TX(2,5)+TX(1,3)*TX(2,4)
      TX(3,7)=                                 TX(1,3)*TX(2,5)

      TY(1,1)=0.5D0*(2-SY)*(1-SY)
      TY(1,2)=SY*(2-SY)
      TY(1,3)=0.5D0*SY*(SY-1)
      TY(1,4)=0.0D0
      TY(1,5)=0.0D0
      TY(1,6)=0.0D0
      TY(1,7)=0.0D0
      TY(2,1)=TY(1,1)**2
      TY(2,2)=2*TY(1,1)*TY(1,2)
      TY(2,3)=2*TY(1,1)*TY(1,3)+TY(1,2)**2
      TY(2,4)=2*TY(1,2)*TY(1,3)
      TY(2,5)=TY(1,3)**2
      TY(2,6)=0.0D0
      TY(2,7)=0.0D0
      TY(3,1)=TY(1,1)*TY(2,1)
      TY(3,2)=TY(1,1)*TY(2,2)+TY(1,2)*TY(2,1)
      TY(3,3)=TY(1,1)*TY(2,3)+TY(1,2)*TY(2,2)+TY(1,3)*TY(2,1)
      TY(3,4)=TY(1,1)*TY(2,4)+TY(1,2)*TY(2,3)+TY(1,3)*TY(2,2)
      TY(3,5)=TY(1,1)*TY(2,5)+TY(1,2)*TY(2,4)+TY(1,3)*TY(2,3)
      TY(3,6)=                TY(1,2)*TY(2,5)+TY(1,3)*TY(2,4)
      TY(3,7)=                                TY(1,3)*TY(2,5)

      DO 70 J=1,7
      DO 70 I=1,3
         TX(I+3,J)=TX(I,J)
 70      TY(I+3,J)=TY(I,J)

      DO 71 J=1,JMAX
      DO 71 I=1,IMAX
CCCCCCC FIELD INITIALIZED
         F1R(I,J)=0.0D0
         F2R(I,J)=0.0D0
         F3R(I,J)=0.0D0
	 F4R(I,J)=0.0D0
	 F4R1(I,J)=0.0D0

 71      continue 

      DO 72 J=1,JMAX
      DO 72 K=1,7
      DO 72 M=1,6
         ABCXR(M,K,J)=0.0D0
 72      continue 
      DO 73 I=1,IMAX
      DO 73 K=1,7
      DO 73 M=1,6
         ABCYR(M,K,I)=0.0D0
 73     continue 

      RETURN
      END
*                                                                *
******************************************************************
*                                                                *


*                                                                *
******************************************************************
*                                                                *
      SUBROUTINE MARCH_FIELDS(IT)
      INCLUDE 'COMMON2D.F'
      DOUBLE PRECISION SUMZ(6),TE,TH,F4RTEMP,N0LEVELTEMP,
     &N1LEVELTEMP,N2LEVELTEMP,N3LEVELTEMP,PUMP
501   FORMAT(1000I2)
502   FORMAT(10000E16.6)
503   FORMAT(20G16.6)

      TE=IT*DT
      TH=(0.50D0+IT)*DT
C*************************************************
      IF(IPOLAR.EQ.0)THEN
C*************************************************
        DO 91 J=1,JMAX
           DO 92 I=1,7
              ABCXR(3,I,J)=ABCXR(2,I,J)
              ABCXR(2,I,J)=ABCXR(1,I,J)
              ABCXR(1,I,J)=F1R(I,J)
              ABCXR(6,I,J)=ABCXR(5,I,J)
              ABCXR(5,I,J)=ABCXR(4,I,J)
              ABCXR(4,I,J)=F1R(IMAX-I+1,J)
92         CONTINUE
91      CONTINUE

        DO 93 I=2,IMAX-1
           DO 94 J=1,7
              ABCYR(3,J,I)=ABCYR(2,J,I)
              ABCYR(2,J,I)=ABCYR(1,J,I)
              ABCYR(1,J,I)=F1R(I,J)
              ABCYR(6,J,I)=ABCYR(5,J,I)
              ABCYR(5,J,I)=ABCYR(4,J,I)
              ABCYR(4,J,I)=F1R(I,JMAX-J+1)
94         CONTINUE
93      CONTINUE
C-------------------------------------------------------
C     UPDATE PZ
        DO 88 J=2,JMAX-1
        DO 88 I=2,IMAX-1
	     F4RTEMP=F4R(I,J)
         DELTAN(I,J)=N1LEVEL(I,J)-N2LEVEL(I,J)
C        DELTAN(I,J)=-N0LEVEL(I,J)
C-----Change Delta N--------------------          
c       DO 303 A=IMAX/2-50,IMAX/2+50
c       DO 303 B=JMAX/2-50,JMAX/2+50
c          DELTAN(A,B)=1.3e25
c303    CONTINUE
C---------------------------
          F4R(I,J)=CPZ1*DELTAN(I,J)*F1R(I,J)
     &     +CPZ2*F4R(I,J)+CPZ3*F4R1(I,J)
          F4R1(I,J)=F4RTEMP

88      CONTINUE
	IF(MOD(IT,LL).EQ.0)THEN
        DO 300 JJ=1,jMAX
	   WRITE(94,502)(F4R(II,jJ)
     &           /F1R(II,jJ)/EPS0,II=1,IMAX)

300      CONTINUE
          CLOSE(94)
        ENDIF

C-------------------------------------------------------
C     UPDATE EZ
        DO 10 J=2,JMAX-1
        DO 10 I=2,IMAX-1
C----------------------------------------------------------
         F1R1(I,J)=F1R(I,J)
         F1R(I,J)=CE1(I,J)*F1R(I,J)+CE2(I,J)
     & *(F3R(I,J)-F3R(I-1,J)-F2R(I,J)+F2R(I,J-1)
     & -NSOURCE(I,J)*HINCR(J)) 
     & +CE3(I,J)*(F4R(I,J)-F4R1(I,J))
C-------------------------------------------------------
10     CONTINUE
C-----
        DO 110 J=1,JMAX
           DO 111 K=1,6
111        SUMZ(K)=0.0D0         
           DO 112 M=1,7
           DO 112 K=1,6
112        SUMZ(K)=SUMZ(K)+TX(K,M)*ABCXR(K,M,J)
           F1R(1,J)=2*SUMZ(1)-SUMZ(2)
           F1R(IMAX,J)=2*SUMZ(4)-SUMZ(5)

110     CONTINUE

C        DO 120 I=2,IMAX-1
C           DO 121 K=1,6
C121        SUMZ(K)=0.0D0         
C           DO 122 M=1,7
C           DO 122 K=1,6
C122        SUMZ(K)=SUMZ(K)+TY(K,M)*ABCYR(K,M,I)
C           F1R(I,1)=2*SUMZ(1)-SUMZ(2)
C           F1R(I,JMAX)=2*SUMZ(4)-SUMZ(5)

C------------change to periodic boundary condition
         DO 120 I=2,IMAX-1
            F1R(I,1)=CE1(I,1)*F1R(I,1)+CE2(I,1)
     & *(F3R(I,1)-F3R(I-1,1)-F2R(I,1)+F2R(I,JMAX)
     & -NSOURCE(I,1)*HINCR(1)) 
     & +CE3(I,1)*(F4R(I,1)-F4R1(I,1))
            F1R(I,JMAX)=CE1(I,JMAX)*F1R(I,JMAX)+CE2(I,JMAX)
     & *(F3R(I,JMAX)-F3R(I-1,JMAX)-F2R(I,JMAX)+F2R(I,JMAX-1)
     & -NSOURCE(I,JMAX)*HINCR(JMAX)) 
     & +CE3(I,JMAX)*(F4R(I,JMAX)-F4R1(I,JMAX))
C            F1R(I,1)=2*SUMZ(1)-SUMZ(2)
C            F1R(I,JMAX)=2*SUMZ(4)-SUMZ(5)

120     CONTINUE
C-------------------------------------------------------
C Update HX
        CALL SOURCE(TE,TH,IT)
        DO 20 J=1,JMAX-1
        DO 20 I=1,IMAX
           F2R(I,J)=F2R(I,J)-CH*(F1R(I,J+1)-F1R(I,J))
20      continue 
C-------------------------------------------------------
	DO 30 J=1,JMAX
        DO 30 I=1,IMAX-1
           F3R(I,J)=F3R(I,J)+CH*(F1R(I+1,J)-F1R(I,J))
c     &     -NSOURCE(I+1,J)*EINCR(J))
 30     CONTINUE

C------------change to periodic boundary condition
        DO 121 I=2,IMAX-1
            F2R(I,JMAX)=F2R(I,JMAX)-CH*(F1R(I,1)-F1R(I,JMAX))
            F3R(I,JMAX)=F3R(I,JMAX)+CH*(F1R(I+1,JMAX)-F1R(I,JMAX))

121     CONTINUE
 
C-------------------------------------------------------
C Update N3,N2,N1,N0
C       PUMP = 1.0E15
C	PUMP = 1.33E9
        PUMP = PUMPRATE0
C    & *(DEXP(-(DFLOAT(IT-TSTART))**2/(2*110.0**2)))**2
CXXXXXXXXXXXXXXXstop pumpingXXXXXXX
C        IF(IT.GT.20000)THEN
C           PUMP=0
C        ENDIF

CXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX


        DO 89 I=2,IMAX-1
        DO 89 J=2,JMAX-1

           N3LEVELTEMP=N3LEVEL(I,J)
           N2LEVELTEMP=N2LEVEL(I,J)
           N1LEVELTEMP=N1LEVEL(I,J)
           N0LEVELTEMP=N0LEVEL(I,J)
C           N3LEVEL(I,J)=CN3LEVEL1*PUMP*N0LEVEL(I,J)
C     &     +CN3LEVEL2*N3LEVEL(I,J)
C           N2LEVEL(I,J)=CN2LEVEL1*(N3LEVEL(I,J)+N3LEVELTEMP)
C     &     +CN2LEVEL2*F1R(I,J)*(F4R(I,J)-F4R1(I,J))
C     &     +CN2LEVEL3*N2LEVEL(I,J)
C           
C           N1LEVEL(I,J)=CN1LEVEL1*(N3LEVEL(I,J)+N3LEVELTEMP)
C     &     +CN1LEVEL2*F1R(I,J)*(F4R(I,J)-F4R1(I,J)) 
C     &     +CN1LEVEL3*(N2LEVEL(I,J)+N2LEVELTEMP)
C     &     +CN1LEVEL4*N1LEVEL(I,J)
C
C           N0LEVEL(I,J)=CN0LEVEL1*N0LEVEL(I,J)
C     &     +CN0LEVEL2*PUMP*N0LEVEL(I,J)
C     &     +CN0LEVEL3*(N1LEVEL(I,J)+N1LEVELTEMP)
           N3LEVEL(I,J)=CN3LEVEL1*PUMP*N0LEVEL(I,J)
     &     +CN3LEVEL2*N3LEVEL(I,J)
           N2LEVEL(I,J)=CN2LEVEL1*(N3LEVEL(I,J)+N3LEVELTEMP)
     &     +CN2LEVEL2*(F1R(I,J)+F1R1(I,J))*(F4R(I,J)-F4R1(I,J))
     &     +CN2LEVEL3*N2LEVEL(I,J)
           
           N1LEVEL(I,J)=CN1LEVEL1*(N3LEVEL(I,J)+N3LEVELTEMP)
     &     +CN1LEVEL2*(F1R(I,J)+F1R1(I,J))*(F4R(I,J)-F4R1(I,J)) 
     &     +CN1LEVEL3*(N2LEVEL(I,J)+N2LEVELTEMP)
     &     +CN1LEVEL4*N1LEVEL(I,J)

           N0LEVEL(I,J)=CN0LEVEL1*N0LEVEL(I,J)
     &     +CN0LEVEL2*PUMP*N0LEVEL(I,J)
     &     +CN0LEVEL3*(N1LEVEL(I,J)+N1LEVELTEMP)
	   
 89      CONTINUE
C*************************************************
C*************************************************
c      ELSEIF(IPOLAR.EQ.1)THEN
C*************************************************

C**********************************************************
      ENDIF
C**********************************************************

      RETURN

      END
*                                                                      *
************************************************************************      
*                                                                      *
      SUBROUTINE COLLECT_DATA(N,NCOUNT)
      INCLUDE 'COMMON2D.F'
      DOUBLE PRECISION COLLECTORR(JMAX),
     &COLLECTOTT(JMAX),COLLECTOPP(IMAX),
     &COLLECTON2(IMAX)
      OPEN(84,FILE='RR.OUT',STATUS='UNKNOWN')
      OPEN(85,FILE='TT.OUT',STATUS='UNKNOWN')
      OPEN(86,FILE='PP.OUT',STATUS='UNKNOWN')
      OPEN(87,FILE='N1.OUT',STATUS='UNKNOWN')

 501  FORMAT(3000E16.6)
      DO 10 JJ=1,JMAX
         IF(ABS(F1R(55,JJ)).LE.1.0D-100)THEN
           COLLECTORR(JJ)=0.0D0		   
         ELSE
           COLLECTORR(JJ)=F1R(55,JJ)
         ENDIF
 
         IF(ABS(F1R(IMAX-30,JJ)).LE.1.0D-100)THEN
           COLLECTOTT(JJ)=0.0D0	   
         ELSE
           COLLECTOTT(JJ)=F1R(IMAX-30,JJ)
         ENDIF
         
	 IF(ABS(F4R(IMAX/2,JJ)).LE.1.0D-100)THEN
           COLLECTOPP(JJ)=0.0D0	   
         ELSE
           COLLECTOPP(JJ)=F4R(IMAX/2,JJ)
     &                 /F1R(IMAX/2,JJ)/EPS0
         ENDIF
	 COLLECTON2(II)=N2LEVEL(II,JMAX/2)-N1LEVEL(II,JMAX/2)
	 
10    CONTINUE
      DO II = 1,IMAX
	 IF(ABS(F4R(II,JMAX/2)).LE.1.0D-100)THEN
           COLLECTOPP(II)=0.0D0	   
         ELSE
           COLLECTOPP(II)=F4R(II,JMAX/2)
     &                 /F1R(II,JMAX/2)/EPS0
         ENDIF
      ENDDO         

      WRITE(84,501)(COLLECTORR(JJ),JJ=1,JMAX)
      WRITE(85,501)(COLLECTOTT(JJ),JJ=1,JMAX)      
      WRITE(86,501)(COLLECTOPP(II),II=1,IMAX)      
      WRITE(87,501)(COLLECTON2(II),II=1,IMAX)
      RETURN
      END
*                                                                      *
************************************************************************      
*                                                                      *
      SUBROUTINE SOURCE(TE,TH,IT)
      INCLUDE 'COMMON2D.F'
      DOUBLE PRECISION TE,TH,AMP1,AMP2,PH1R,PH2R,TMP1,TMP2

         AMP1=AMPLITUDE
         AMP2=AMPLITUDE

         PH1R=DCOS(WI*TE)*AMP1
C    &       *DEXP(-(DFLOAT(IT-400))**2/(2*330.0**2))
         PH2R=DCOS(WI*TH+WI/C0*DX/2)*AMP1
C     &       *DEXP(-(DFLOAT(IT-400))**2/(2*330.0**2))

        AMP2=-1.0D0/XMU0/C0
	TMP1=(JENDO+JSTAR)/2
	TMP2=(JENDO-JSTAR)

	DO 10 J=1,JMAX
         AMP1=1.0d0
             EINCR(J)=PH1R*AMP1
             HINCR(J)=PH2R*AMP1*AMP2
C             EINCR(J)=PH1R*AMP1*DEXP(-4*PI*(J-TMP1)**2/TMP2**2)
C            HINCR(J)=PH2R*AMP1*AMP2*DEXP(-4*PI*(J-TMP1)**2/TMP2**2)
10      CONTINUE
	
	RETURN
      END

*                                                                      *
************************************************************************      
*                                                                      *
