      SUBROUTINE TG03A(K,MPK,A,T,MP2K,WK,IW,XVALUE,ID,S)
C      **********************************************************
C                            PURPOSE
C      **********************************************************
C
C      LET S(X) BE A SPLINE OF DEGREE K-1 WITH M INTERIOR KNOTS
C      IN A FINITE INTERVAL.  THE PURPOSE OF THIS SUBROUTINE
C      IS TO EVALUATE THE DERIVATIVES S(J-1)(X),J=1,...,DERIV
C      WHERE DERIV=MAX(1,MIN(ID,K)), AT A GIVEN VALUE OF X
C      GIVEN THE M+K COEFFICIENTS IN THE NORMALISED B-SPLINE
C      REPRESENTATION OF S(X).
C
C      **********************************************************
      REAL A(MPK),T(MP2K),WK(IW),S(ID)
      INTEGER DERIV
C      ********************************************************
C      ASSUME THAT K AND MPK HAVE SENSIBLE VALUES
C      ********************************************************
C     IF(K.LT.1) GOTO 150
C     IF(MPK.LT.K) GOTO 160
C      END OF CHECK.
C      *********************************************************
      DERIV=MAX(MIN(ID,K),1)
C      *********************************************************
C      INITIALISING S TO ZERO.
C      **********************************************************
           DO 10 I=1,DERIV
   10      S(I)= 0
C      END OF INITIALISATION.
C      **********************************************************
      MPKP1=MPK+1
      IF(XVALUE.LT.T(K).OR.XVALUE.GT.T(MPKP1)) RETURN
C      *********************************************************
C      COMPUTE THE VALUE OF JINT
C      **********************************************************
      IL=K
      IR=MPKP1
   20 IF(IR-IL.LE.1) GOTO 40
      MIDDLE=(IL+IR)/2
      IF(XVALUE.LT.T(MIDDLE)) GOTO 30
      IL=MIDDLE
      GOTO 20
   30 IR=MIDDLE
      GOTO 20
   40 JINT=IL
C      JINT IS SUCH THAT  T(JINT).LE.XVALUE.LT.T(JINT+1).
C      *********************************************************
      E1=XVALUE-T(JINT)
      E2=T(JINT+1)-XVALUE
      R1=T(JINT+1)-T(JINT)
      KP1=K+1
C      *********************************************************
C      CALCULATE S(J-1)(X),J=1,...,DERIV AT X=XVALUE
C      *********************************************************
           DO 140 IJ=1,DERIV
           KMIJP1=KP1-IJ
           RK=REAL(KMIJP1)
C      COMPUTE THE B-SPLINES N(K+1-IJ)(X) AT X=XVALUE.
           WK(KP1)= 1
           IF(KMIJP1.EQ.1)GOTO 110
           WK(KP1)= 1/R1
           K1=KMIJP1-1
           IF(K1.LT.2) GOTO 60
                DO 50 J=2,K1
                JINTPJ=JINT+J
                KPJ=K+J
   50           WK(KPJ)=E1*WK(KPJ-1)/(T(JINTPJ)-T(JINT))
   60      KPJ=K+KMIJP1
           WK(KPJ)=E1*WK(KPJ-1)
                DO 100 J=1,K1
                JINTMJ=JINT-J
                E3=XVALUE-T(JINTMJ)
                WK(KP1)=E2*WK(KP1)
                IF(J.EQ.K1)GOTO 70
                WK(KP1)=WK(KP1)/(T(JINT+1)-T(JINTMJ))
   70           NKJ=KMIJP1-J
                IF(NKJ.LT.2)GOTO100
                N1=NKJ-1
                IF(N1.LT.2)GOTO 90
                     DO 80 L=2,N1
                     JINTPL=JINT+L
                     KPL=K+L
   80                WK(KPL)=(E3*WK(KPL-1)+(T(JINTPL)-XVALUE)*WK(KPL))/(
     $               T(JINTPL)-T(JINTMJ))
   90           KPNKJ=K+NKJ
                JN=JINT+NKJ
                WK(KPNKJ)=E3*WK(KPNKJ-1)+(T(JN)-XVALUE)*WK(KPNKJ)
  100           CONTINUE
C      END OF COMPUTING THE B-SPLINES.
  110      J1=JINT
           JINTMK=JINT-K
                DO 130 I=1,KMIJP1
                KPI=K+I
                IF(IJ.GT.1)GOTO 120
                IR=JINTMK+I
                WK(I)=A(IR)
                GOTO 130
  120           J1=J1+1
                J2=J1-KMIJP1
                WK(I)=RK*(WK(I+1)-WK(I))/(T(J1)-T(J2))
  130           S(IJ)=S(IJ)+WK(I)*WK(KPI)
  140      CONTINUE
C      END OF CALCULATING THE DERIVATIVES.
C      ********************************************************
      RETURN
      END
