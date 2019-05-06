 
      PROGRAM UNIF
 
* source directed along x direction
* with exponential bias, parameter a
* see MCNP manual p3-64

      REAL         CORN1(3), CORN2(3)
      REAL         SPOS(3),  EPOS(3)
      REAL         D(3)
      LOGICAL      ILOPEN
      CHARACTER*80 FILNAM
      CHARACTER*80 LINEIN
 
*-----------------------------------------------------------------------
*L              1.         Initialise
 
 
*L                  1.1      Problem dependent numbers
      CORN1(1)=-4700.
      CORN1(2)=-600.
      CORN1(3)=-850.
      CORN2(1)=1064.
      CORN2(2)=600.
      CORN2(3)=860.

      SPOS(1)=CORN1(1)+10.
      SPOS(2)=0.5*(CORN1(2)+CORN2(2))
      SPOS(3)=0.5*(CORN1(3)+CORN2(3))
* Offset slightly from centre
      SPOS(2)=0.49*CORN1(2)+0.51*CORN2(2)
      SPOS(3)=0.49*CORN1(3)+0.51*CORN2(3)
 
*L                  1.2      Defaults
      NNUM=100
      NAM=1
      NAE=0
      SINF=1.E+22
      EPS=1.E-4
      EPSM=1.-EPS
      PI=4.*ATAN(1.)
 
*L                  1.3      Admin
      I0=ICHAR('0')
* File units
* the input data is in argument list
* and qry output is on channel NPUNCH
      NREAD=1
      NPUNCH=2
 
*-----------------------------------------------------------------------
*L              2.         argument list
 
      INARGS=IARGC()
      I=0
      DO 200 J=1,INARGS
      CALL GETARG(J,LINEIN)
      IF (LINEIN(1:1).NE.'-') THEN
        I=I+1
        READ(LINEIN,*) IVAL
*       IVAL=ICHAR(LINEIN(1:1))-I0
        IF (I.EQ.1) THEN
          NNUM=IVAL
        ELSE IF (I.EQ.2) THEN
          NAM=IVAL
        ELSE IF (I.EQ.3) THEN
          NAE=IVAL
        END IF
      END IF
  200 CONTINUE
*L                  2.1      open new file
      INQUIRE(UNIT=NPUNCH,OPENED=ILOPEN)
      IF (ILOPEN) CLOSE(UNIT=NPUNCH)
      WRITE(FILNAM,9000) NNUM
 9000 FORMAT('unif',I7.7,'.qry')
      OPEN(UNIT=NPUNCH,FILE=FILNAM(1:16))
      NNUMD=2*NNUM
      WRITE(NPUNCH,*) 'POINTS ',NNUMD
 
* set constants
      A=REAL(NAM)*10.**REAL(-NAE)
      WRITE(*,*) 'A=',A
      IF (ABS(A).LT.1) THEN
        ANG=ASIN(A)*180/PI
        WRITE(*,*) 'ANG(DEG)=',ANG
      ELSE
        WRITE(*,*) 'ABS(A) GE 1 MAKES NO SENSE'
        STOP
      END IF
 
*-----------------------------------------------------------------------
*L              3.         main loop
 
      DO 303 J=1,NNUM
  300 CONTINUE
      XIM=RANDLX(0.)
      ZMU=A*XIM
      ZSMU=SQRT(1.-ZMU**2)
      XIP=RANDLX(0.)
      ARG=2.*PI*XIP
      D(1)=ZSMU
      D(2)=ZMU*SIN(ARG)
      D(3)=ZMU*COS(ARG)
      ZSMIN=SINF
      ZS1=-1.
      ZS2=-1.
 
      DO 301 JI=1,3
      IF (ABS(D(JI)).GT.EPS) ZS1=(CORN1(JI)-SPOS(JI))/D(JI)
      IF (ZS1.GT.0.) THEN
        ZSMIN=MIN(ZS1,ZSMIN)
      END IF
      IF (ABS(D(JI)).GT.EPS) ZS2=(CORN2(JI)-SPOS(JI))/D(JI)
      IF (ZS2.GT.0.) THEN
        ZSMIN=MIN(ZS2,ZSMIN)
      END IF
  301 CONTINUE
 
      IF (ZSMIN.LT.SINF) THEN
        DO 302 JI=1,3
        EPOS(JI)=SPOS(JI)+EPSM*ZSMIN*D(JI)
  302   CONTINUE
      ELSE
        GO TO 300
      END IF
 
      WRITE(NPUNCH,9001) (SPOS(JI),JI=1,3)
      WRITE(NPUNCH,9001) (EPOS(JI),JI=1,3)
 
  303 CONTINUE
 
      STOP
 9001 FORMAT(3(1X,1PG14.7))
      END
