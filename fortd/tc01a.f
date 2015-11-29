      SUBROUTINE TC01A (XN,FN,GN,N,X,Y,KN)
C
C************************************************************
C*  PURPOSE.                                                *
C*    GIVEN A CUBIC SPLINE S(X) DEFINED IN TERMS OF ITS     *
C*    KNOTS X(I),VALUES S(X(I)) AND FIRST DERIVATIVE AT THE *
C*    KNOTS S'(X(I)).THIS ROUTINE FINDS THE VALUE OF X      *
C*    GIVEN THE VALUE OF S(X).AN ESTIMATE OF X MUST BE      *
C*    SUPPLIED AND THE CALCULATED VALUE IS THE VALUE        *
C*    NEAREST THIS ESTIMATE.                                *
C*                                                          *
C*  ARGUMENT LIST.                                          *
C*    XN   IS A REAL ARRAY OF LENGTH N WHICH CONTAINS THE   *
C*         VALUES OF THE KNOTS X(1)<X(2)<...<X(N)           *
C*    FN   IS A REAL ARRAY OF LENGTH N WHICH CONTAINS THE   *
C*         VALUES OF THE SPLINE AT THE KNOTS.SO THAT S(X(I))*
C*         IS FOUND IN FN(I),I=1,2,..,N                     *
C*    GN   IS A REAL ARRAY OF LENGTH N WHICH CONTAINS THE   *
C*         VALUES OF THE FIRST DERIVATIVE AT THE KNOTS,SO   *
C*         THAT S'(X(I)) IS FOUND IN GN(I),I=1,2,...,N      *
C*    N    N IS AN INTEGER SPECIFYING THE NUMBER OF KNOTS   *
C*    X    IS A REAL VARIABLE SET BY USER TO ESTIMATE OF X  *
C*    Y    IS A REAL VARIABLE SET BY USER TO VALUE OF S(X)  *
C*    KN   IS AN INTEGER SET BY THE ROUTINE TO THE INTERVAL *
C*         IN WHICH THE ROOT IS FOUND. (0<KN<N) THAT IS THE *
C*    INTERVAL (X(KN),X(KN+1)).                             *
C*         KN=0 IF X(1)<=X<=X(N) , THERE WAS NO SOLUTION    *
C*         KN=-1 VALUE OF N WAS NON POSITIVE                *
C************************************************************
C     NON-STANDARD:
C     ASSUME THAT ALL INPUTS HAVE SENSIBLE VALUES.
C     AND DIAGNOSTICS REMOVED
C************************************************************
C
      DIMENSION XN(N),FN(N),GN(N)
C
C     EPS IS THE MACHINE PRECISION*10
C
      EPS=EPSILON(1.D0)*10.0
      P5=REAL(1)/REAL(2)
      IR = 0
      IF(N-1)680,100,110
C
C     THERE IS ONLY ONE KNOT POINT
  100 IF(FN(1) .NE. Y)GO TO 670
      X = XN(1)
      GO TO 700
C
  110 KN = 0
      IF(X .LE. XN(1))GO TO 140
      IF(X .GE. XN(N))GO TO 180
C
C     LOOK FOR THE INTERVAL WHICH CONTAINS THE USERS ESTIMATE X
C     THIS WILL BE THE FIRST INTERVAL TO LOOK FOR A ROOT
      IB = 1
      IE = N
  120 IF(IE .LE. IB+1)GO TO 210
      K = (IB+IE)/2
      IF(XN(K) .LT. X)GO TO 130
      IE = K
      GO TO 120
  130 IB = K
      GO TO 120
C
C      THE ESTIMATE OF X SUPPLIED BY THE USER WAS <X(1),SCAN
C      THE INTERVALS FROM (X(1),X(2)),TO X(N-1),X(N) FOR A ROOT
  140 IB = 1
  150 IE = -1
  160 IF(IB-N)210,670,670
  170 IB = IB+1
      GO TO 160
C
C     THE ESTIMATE OF X SUPPLIED BY THE USER WAS >X(N),SCAN THE
C     INTERVALS FROM (X(N-1),X(N)) TO (X(1),X(2)) FOR A ROOT
  180 IB = N
  190 IE = 0
  200 IB = IB-1
      IF(IB-1)670,210,210
C
C     IS THERE AROOT IN THE INTERVAL IN (X(I),X(I+1))
  210 I = IB
  220 A0 = FN(I)-Y
      SB = FN(I+1)-Y
      H = XN(I+1)-XN(I)
      A1 = H*GN(I)
      QA = H*GN(I+1)
      AL1 = SB-A0-QA
      AL0 = SB-A0-A1
      A31H = 3*AL1/H
      A30H = 3*AL0/H
      B2 = -AL0-AL1
      A2 = AL0-B2
      JR = 1
      IF(A0*SB .LT. 0)GO TO 290
      JR = 0
      R1 = P5*A1+A0
      R2 = -P5*QA+SB
      IF(A0)230,250,240
  230 IF(MAX(R1,R2))620,260,260
  240 IF(MIN(R1,R2))260,260,620
C
C      THERE IS A ROOT AT A KNOT XN(I)
  250 TC = 0
      GO TO 390
C
C     MORE TESTS ARE REQUIRED TO DETECT A POSSIBLE ROOT
  260 IF(B2 .NE. 0)GO TO 270
      IF(A2 .EQ. 0)GO TO 620
      R = -P5*A1/A2
      GO TO 280
C
  270 AD = 1/(3*B2)
      U = A2*AD
      AE = A1*AD
      V = U*U-AE
      IF(V .LT. 0)GO TO 620
      R = - SIGN(1.,U)*( ABS(U)+ SQRT(V))
      IF(A0*A2 .GT. 0)R = AE/R
  280 IF( ABS(R-P5) .GT. P5)GO TO 620
      SR = ((B2*R+A2)*R+A1)*R +A0
      IF(A0*SR .GT. 0)GO TO 620
C
C     FIND A ROOT IN THE INTERVAL
C
C     SET UP BRACKET FOR THE ROOT
  290 TC = -1
      IF(X .LT. XN(I))GO TO 320
      IF(X .GT. XN(I+1))GO TO 300
      TC = (X-XN(I))/H
      IF(JR .EQ. 1)GO TO 320
      IF(TC-R)320,310,310
  300 IF(JR .EQ. 1)GO TO 320
  310 TB = R
      FB = SR
      GO TO 330
  320 TB = 0
      FB = A0
      IF(JR .EQ. 1)GO TO 330
      TE = R
      FE = SR
      GO TO 340
  330 TE = 1
      FE = SB
  340 IF(TC .LT. 0)TC = P5*(TE+TB)
      DELTA = EPS*( ABS(A0) + ABS(SB) + 1.481*( ABS(AL1)+ ABS(AL0)))
C
C     USE NEWTON RAPHSON TO FIND A ROOT
  350 TT = 1-TC
      TA = TT*TC
      PT = TC*(TA*AL1 + SB)+TT*(-TA*AL0+A0)
      IF( ABS(PT) .LE. DELTA)GO TO 390
C
C     REVISE THE BRACKET ON THE ROOT
      IF(PT*FB .LT. 0)GO TO 360
      FB = PT
      TB = TC
      GO TO 370
  360 FE = PT
      TE = TC
  370 PTD = TC*(A31H*TT+GN(I+1))+TT*(A30H*TC+GN(I))
      IF(PTD .EQ. 0)GO TO 380
      TC =TC-PT/PTD
C
C     CHECK THAT NEXT VALUE OF TC IS INSIDE THE INTERVAL (TB,TE)
      IF((TC-TB)*(TC-TE) .LT. 0)GO TO 350
C
C     CALCULATE THE NEXT VALUE OF T BY BISECTION
  380 TC = (TB+TE)*P5
      IF((TC-TB)*(TC-TE))350,390,390
C
C     ONE ROOT FOUND,LOOK AT THE ROOTS OF THE DEFLATED POLYNOMIAL
C     B0 + B1*T + B2*T*T
  390 B1 = B2*TC+AL0-B2
      B0 = B1*TC+A1
      TC = H*TC+XN(I)
      IF(B2 .NE. 0)GO TO 400
      IF(B1 .EQ. 0)GO TO 530
      R1 = -B0/B1
      IF( ABS(R1-P5) .GT. P5)GO TO 530
      R1 = H*R1+XN(I)
      GO TO 440
  400 AD = 1/B2
      AF = B0*AD
      U = P5*B1*AD
      V = U*U-AF
      IF(V .LT. 0)GO TO 530
C
C     STORE VALUES OF ROOTS OF QUADRATIC IN R1 .LE. R2
      R2 =  ABS(U)+ SQRT(V)
      R1 = 0
      IF(R2 .EQ. 0)GO TO 420
      IF(U .GE. 0)GO TO 410
      R1 = AF/R2
      GO TO 420
  410 R1 = -R2
      R2 = AF/R1
  420 R1 = H*R1+XN(I)
      R2 = H*R2+XN(I)
      IF(JR .NE. 1)GO TO 430
C
C     THERE ARE EITHER THREE ROOTS OR ONE ROOT IN INTERVAL
      AE = P5- ABS(U+P5)
      IF(AE .LT. 0)GO TO 530
      IF(V-AE*AE)450,450,530
C
C     THERE ARE THREE ROOTS IN THE INTERVAL
  430 IF(U*U .LE. V)R1 = R2
  440 R2 = R1
  450 S2 = TC
C
C    ORDER ROOTS IN THE INTERVAL
      IF(TC .GT. R1)GO TO 460
      S2 = R1
      R1 = TC
      GO TO 470
  460 IF(TC .LE. R2)GO TO 470
      S2 = R2
      R2 = TC
C
C     FIND ROOT NEAREST X
  470 IF(X .LE. S2)GO TO 490
      IF(X .LE. R2)GO TO 510
  480 Q1 = R2
      JJ = 1-IR
      IN1 = I
      GO TO 540
  490 IF(X .LE. R1)GO TO 520
      Q1 = R1
      Q2 = S2
  500 IN1 = I
      IN2 = I
      GO TO 550
  510 Q1 = S2
      Q2 = R2
      GO TO 500
  520 Q2 = R1
      IN2 = I
      JJ = IR-1
      GO TO 540
C
C     THERE IS ONLY ONE ROOT
  530 R1 = TC
      R2 = TC
      IF(X-TC)520,520,480
  540 IF(JJ .NE. 0)GO TO 580
  550 IF((X-Q1) .LE. (Q2-X))GO TO 570
  560 X = Q2
      KN = IN2
      GO TO 700
  570 X = Q1
      KN = IN1
      GO TO 700
  580 IF(IE)560,570,590
  590 IR = 1
      IF(JJ)600,550,610
  600 IF(Q2-X-X+XN(IB))560,560,190
  610 IF(XN(IE)-X-X+Q1)640,570,570
C
C
C     NO ROOT IN INTERVAL X(I)<X(I+1)
  620 IF(IE)170,200,630
  630 IF((X-XN(IB)) .GT. (XN(IE)-X))GO TO 650
      IF(IB .GE. 2)GO TO 200
  640 IB = IE
      GO TO 150
  650 IF(IE .GE. N)GO TO 190
      I = IE
      IE = IE+1
      GO TO 220
  660 IF(JJ)560,550,570
C
C     NO ERROR MESSAGES AND RETURNS
  670 IF(IR .EQ. 1)GO TO 660
      GO TO 700
  680 KN = -1
  700 RETURN
      END
