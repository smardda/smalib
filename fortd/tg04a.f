C
         SUBROUTINE TG04A(N,X,K,NORM,XVALUE,JINT,V,VINT)
C
C     ********************************************************
C     PURPOSE
C     ********************************************************
C     GIVEN THE DATA POINTS X(1)<=X(2)<=,...,<=X(N-1)<=X(N)
C     AND A VALUE OF X, X=XVALUE, WHERE X(1)<=XVALUE<=X(N),
C     THIS SUBROUTINE COMPUTES THE VALUES OF THE B-SPLINES OF
C     DEGREE K-1, (MKI(X)) 1<=I<=N-K, AT X=XVALUE, AND THE
C     VALUES OF THE CORRESPONDING INTEGRALS OVER THE RANGE
C     X(I)<=X<=XVALUE.
C     N.B. THE FUNCTION MKI(X) IS A SPLINE OF DEGREE K-1 WITH
C     KNOTS AT THE POINTS X(I),X(I+1),...,X(I+K), AND IS
C     NON-ZERO ONLY OVER THE RANGE X(I)<X<X(I+K).
C     ********************************************************
       REAL        X(*),     V(*),     VINT(*)
C     ********************************************************
C     NON-STANDARD:
C     ASSUME THAT ALL INPUTS HAVE SENSIBLE VALUES.
C     AND DIAGNOSTICS REMOVED
C     ********************************************************
    3    XLEFT=X(1)
         XRITE=X(N)
C     ********************************************************
C     INITIALISING V AND VINT TO ZERO.
C     ********************************************************
         DO 4 I=1,K
         VINT(I)=0
    4    V(I)=0
C     END OF INITIALISATION.
C     ********************************************************
C
         IF(XVALUE.GE.XLEFT) GOTO 5
         JINT=1
         RETURN
C
C     **********************************************************
C     CHECK, AND IF NECESSARY COMPUTE, THE VALUE OF JINT.
C     **********************************************************
    5    IF(XVALUE.LT.XRITE) GOTO 8
C     CHECK TO SEE IF XVALUE(=XRITE) IS A MULTIPLE KNOT
         I1=N
         I2=N-1
         DO 6 I=1,I2
         IRMI=N-I
         IF(X(IRMI).LT.XVALUE)GOTO 7
         I1=IRMI-1
    6    CONTINUE
    7    JINT=I1
         IF(JINT.GT.0) GOTO 13
         GOTO 19
C
    8    IF(JINT.LT.1.OR.JINT.GE.N) GOTO 9
         IF(XVALUE.LT.X(JINT)) GOTO 9
         IF(XVALUE.LT.X(JINT+1)) GOTO 13
    9    IL=1
         IR=N
   10    IF(IR-IL.LE.1) GOTO 12
         MIDDLE=(IL+IR)/2
         IF(XVALUE.LT.X(MIDDLE)) GOTO 11
         IL=MIDDLE
         GOTO 10
   11    IR=MIDDLE
         GOTO 10
   12    JINT=IL
C     **********************************************************
C
C     **********************************************************
C     REMOVE
C     CHECK THAT THE DATA POINTS USED IN THE CALCULATION OF THE
C     B-SPLINES ARE IN ASCENDING ORDER
C     **********************************************************
   13    JINTMK=JINT-K
C     **********************************************************
C
C     **********************************************************
C     REMOVE
C     CHECK THAT NO MORE THAN K OF THESE DATA POINTS COALESCE
C     **********************************************************
C     ********************************************************
   19    JJ1=N-JINT
         IF(JJ1.GT.0) GOTO 20
         IF(NORM.EQ.1) VINT(1)=1/REAL(K)
         IF(NORM.EQ.2) VINT(1)=(X(JINT)-X(JINTMK))/REAL(K)
         RETURN
C
C     ********************************************************
C     COMPUTE M1JINT(XVALUE) AND (XVALUE-X(JINT))*M1JINT(XVALUE)
C     AND TEST FOR K=1
C     ********************************************************
   20    E1=XVALUE-X(JINT)
         E2=X(JINT+1)-XVALUE
         V(1)=1/(X(JINT+1)-X(JINT))
         VINT(1)=E1*V(1)
         IF(K.GT.1)GOTO 21
         IF(NORM.EQ.1) RETURN
         V(1)=1
         VINT(1)=E1
         RETURN
C     ********************************************************
C
C     ********************************************************
C     COMPUTE AND STORE IN V(J) THE VALUE OF MJJINT(XVALUE)
C     AND STORE IN VINT(J) THE CORRESPONDING VALUE OF
C     (XVALUE-X(JINT))*MJJINT(XVALUE) FOR J=2,...,NK
C     ********************************************************
   21    NK=MIN(K,JJ1)
         IF(NK.EQ.1)GOTO 24
         DO 23 J=2,NK
         JINTPJ=JINT+J
         IF(J.EQ.K.AND.NORM.EQ.2)GOTO 22
         V(J)=E1*V(J-1)/(X(JINTPJ)-X(JINT))
         VINT(J)=E1*V(J)
         GOTO 23
   22    V(J)=E1*V(J-1)
         VINT(J)=E1*V(J)/(X(JINTPJ)-X(JINT))
   23    CONTINUE
C     ********************************************************
C
   24    IF(JINT.EQ.1) GOTO 30
C
C     ********************************************************
C     COMPUTE AND STORE IN V(1) THE VALUE OF
C     MJ+1JINT-J(XVALUE) AND STORE IN VINT(1) THE APPROPRIATE
C     SUM FOR THE INTEGRAL
C     ********************************************************
         MJ=MIN(JINT,K)-1
         DO 29 J=1,MJ
         JINTMJ=JINT-J
         E3=XVALUE-X(JINTMJ)
         IF(J+1.EQ.K.AND.NORM.EQ.2)GOTO 25
         V(1)=E2*V(1)/(X(JINT+1)-X(JINTMJ))
         VINT(1)=VINT(1)+E3*V(1)
         GOTO 26
   25    V(1)=E2*V(1)
         VINT(1)=VINT(1)+E3*V(1)/(X(JINT+1)-X(JINTMJ))
C     ********************************************************
C
   26    NKJ=MIN(K-J,JJ1)
         IF(NKJ.LE.1)GOTO 29
C
C     ********************************************************
C     COMPUTE AND STORE IN V(L) THE VALUE MJ+LJINT-J(XVALUE)
C     AND STORE IN VINT(L) THE ACCUMULATED SUM FOR THE
C     INTEGRALS, FOR L=2,...,NKJ=MIN(K-J,JJ1)
C     *********************************************************
         DO 28 L=2,NKJ
         JINTPL=JINT+L
         IF(J+L.EQ.K.AND.NORM.EQ.2)GOTO 27
         V(L)=(E3*V(L-1)+(X(JINTPL)-XVALUE)*V(L))/(X(JINTPL)-X(
     +   JINTMJ))
         VINT(L)=VINT(L)+E3*V(L)
         GOTO 28
   27    V(L)=E3*V(L-1)+(X(JINTPL)-XVALUE)*V(L)
         VINT(L)=VINT(L)+E3*V(L)/(X(JINTPL)-X(JINTMJ))
   28    CONTINUE
   29    CONTINUE
C     *********************************************************
C
C     ********************************************************
C     FORMING THE INTEGRALS AND RESETTING THE APPROPRIATE
C     LOCATIONS OF V AND VINT TO ZERO
C     ********************************************************
   30    DO 32 J=1,K
         JINTPJ=JINT+J
         JTPJMK=JINTPJ-K
         IF(JINTPJ.GT.N)RETURN
         IF(JTPJMK.LT.1)GOTO 31
         IF(NORM.EQ.1)VINT(J)=VINT(J)/REAL(K)
         IF(NORM.EQ.2)VINT(J)=(X(JINTPJ)-X(JTPJMK))*VINT(J)/REAL(K)
         GOTO 32
   31    V(J)=0
         VINT(J)=0
   32    CONTINUE
         RETURN
         END
