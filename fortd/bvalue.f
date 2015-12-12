      SUBROUTINE BVALUE (T,BCOEF,N,K,X,JDERIV,KUNIF,MFLAG,PVALUE)
C  FROM  * A PRACTICAL GUIDE TO SPLINES *  BY C. DE BOOR    
C
CALLS  INTERV AND
CALLS  INTERU IF NODES ARE UNIFORMLY SPACED AS ASSUMED BY BEQ
C
CALCULATES VALUE AT  X  OF  JDERIV-TH DERIVATIVE OF SPLINE FROM B-REPR.
C  THE SPLINE IS TAKEN TO BE CONTINUOUS FROM THE RIGHT, EXCEPT AT THE
C  RIGHTMOST KNOT, WHERE IT IS TAKEN TO BE CONTINUOUS FROM THE LEFT.
C
C******  I N P U T ******
C  T, BCOEF, N, K......FORMS THE B-REPRESENTATION OF THE SPLINE  F  TO
C        BE EVALUATED. SPECIFICALLY,
C  T.....KNOT SEQUENCE, OF LENGTH  N+K, ASSUMED NONDECREASING.
C  BCOEF.....B-COEFFICIENT SEQUENCE, OF LENGTH  N .
C  N.....LENGTH OF  BCOEF  AND DIMENSION OF SPLINE(K,T),
C        A S S U M E D  POSITIVE .
C  K.....ORDER OF THE SPLINE .
C  KUNIF.....UNITY FOR QUASI-UNIFORM KNOTS
C
C  W A R N I N G . . .   THE RESTRICTION  K .LE. KMAX (=20)  IS IMPOSED
C        ARBITRARILY BY THE DIMENSION STATEMENT FOR  AJ, DL, DR  BELOW,
C        BUT IS  N O W H E R E  C H E C K E D  FOR.
C
C  X.....THE POINT AT WHICH TO EVALUATE .
C  JDERIV.....INTEGER GIVING THE ORDER OF THE DERIVATIVE TO BE EVALUATED
C        A S S U M E D  TO BE ZERO OR POSITIVE.
C
C******  O U T P U T  ******
C  PVALUE.....THE VALUE OF THE (JDERIV)-TH DERIVATIVE OF  F  AT  X .
C  MFLAG.....NON-ZERO IMPLIES ERROR
C
C******  M E T H O D  ******
C     THE NONTRIVIAL KNOT INTERVAL  (T(I),T(I+1))  CONTAINING  X  IS LO-
C  CATED WITH THE AID OF  INTERV . THE  K  B-COEFFS OF  F  RELEVANT FOR
C  THIS INTERVAL ARE THEN OBTAINED FROM  BCOEF (OR TAKEN TO BE ZERO IF
C  NOT EXPLICITLY AVAILABLE) AND ARE THEN DIFFERENCED  JDERIV  TIMES TO
C  OBTAIN THE B-COEFFS OF  (D**JDERIV)F  RELEVANT FOR THAT INTERVAL.
C  PRECISELY, WITH  J = JDERIV, WE HAVE FROM X.(12) OF THE TEXT THAT
C
C     (D**J)F  =  SUM ( BCOEF(.,J)*B(.,K-J,T) )
C
C  WHERE
C                   / BCOEF(.),                     ,  J .EQ. 0
C                   /
C    BCOEF(.,J)  =  / BCOEF(.,J-1) - BCOEF(.-1,J-1)
C                   / ----------------------------- ,  J .GT. 0
C                   /    (T(.+K-J) - T(.))/(K-J)
C
C     THEN, WE USE REPEATEDLY THE FACT THAT
C
C    SUM ( A(.)*B(.,M,T)(X) )  =  SUM ( A(.,X)*B(.,M-1,T)(X) )
C  WITH
C                 (X - T(.))*A(.) + (T(.+M-1) - X)*A(.-1)
C    A(.,X)  =    ---------------------------------------
C                 (X - T(.))      + (T(.+M-1) - X)
C
C  TO WRITE  (D**J)F(X)  EVENTUALLY AS A LINEAR COMBINATION OF B-SPLINES
C  OF ORDER  1 , AND THE COEFFICIENT FOR  B(I,1,T)(X)  MUST THEN BE THE
C  DESIRED NUMBER  (D**J)F(X). (SEE X.(17)-(19) OF TEXT).
C
      INTEGER JDERIV,K,N,   I,ILO,IMK,J,JC,JCMIN,JCMAX,JJ,KMAX,KMJ,KM1
     *                     ,MFLAG,NMI,JDRVP1
      PARAMETER (KMAX = 20)
      REAL BCOEF(N),T(N+K),X,   AJ(KMAX),DL(KMAX),DR(KMAX),FKMJ
      PVALUE = 0.
      IF (JDERIV .GE. K)                GO TO 99
C
C  *** FIND  I   S.T.   1 .LE. I .LT. N+K   AND   T(I) .LT. T(I+1)   AND
C      T(I) .LE. X .LT. T(I+1) . IF NO SUCH I CAN BE FOUND,  X  LIES
C      OUTSIDE THE SUPPORT OF  THE SPLINE  F , HENCE  PVALUE = 0.
C      (THE ASYMMETRY IN THIS CHOICE OF  I  MAKES  F  RIGHTCONTINUOUS, EXCEPT
C      AT  T(N+K) WHERE IT IS LEFTCONTINUOUS.)
      IF (KUNIF.EQ.0) CALL INTERV ( T, N+K, X, I, MFLAG )
      IF (KUNIF.EQ.1) THEN
         ZRDT=1/(T(K+2)-T(K+1))
         IOFF=1+(K-1)/2
         CALL INTERU ( ZRDT, N, X-T(K), I, 0, IOFF, K )
         MFLAG=0
      END IF
      IF (MFLAG .NE. 0)                 GO TO 99
C  *** IF K = 1 (AND JDERIV = 0), PVALUE = BCOEF(I).
      KM1 = K - 1
      IF (KM1 .GT. 0)                   GO TO 1
      PVALUE = BCOEF(I)
                                        GO TO 99
C
C  *** STORE THE K B-SPLINE COEFFICIENTS RELEVANT FOR THE KNOT INTERVAL
C     (T(I),T(I+1)) IN AJ(1),...,AJ(K) AND COMPUTE DL(J) = X - T(I+1-J),
C     DR(J) = T(I+J) - X, J=1,...,K-1 . SET ANY OF THE AJ NOT OBTAINABLE
C     FROM INPUT TO ZERO. SET ANY T.S NOT OBTAINABLE EQUAL TO T(1) OR
C     TO T(N+K) APPROPRIATELY.
    1 JCMIN = 1
      IMK = I - K
      IF (IMK .GE. 0)                   GO TO 8
      JCMIN = 1 - IMK
      DO 5 J=1,I
    5    DL(J) = X - T(I+1-J)
      DO 6 J=I,KM1
         AJ(K-J) = 0.
    6    DL(J) = DL(I)
                                        GO TO 10
    8 DO 9 J=1,KM1
    9    DL(J) = X - T(I+1-J)
C
   10 JCMAX = K
      NMI = N - I
      IF (NMI .GE. 0)                   GO TO 18
      JCMAX = K + NMI
      DO 15 J=1,JCMAX
   15    DR(J) = T(I+J) - X
      DO 16 J=JCMAX,KM1
         AJ(J+1) = 0.
   16    DR(J) = DR(JCMAX)
                                        GO TO 20
   18 DO 19 J=1,KM1
   19    DR(J) = T(I+J) - X
C
   20 DO 21 JC=JCMIN,JCMAX
   21    AJ(JC) = BCOEF(IMK + JC)
C
C               *** DIFFERENCE THE COEFFICIENTS  JDERIV  TIMES.
      IF (JDERIV .EQ. 0)                GO TO 30
      DO 23 J=1,JDERIV
         KMJ = K-J
         FKMJ = FLOAT(KMJ)
         ILO = KMJ
         DO 23 JJ=1,KMJ
            AJ(JJ) = ((AJ(JJ+1) - AJ(JJ))/(DL(ILO) + DR(JJ)))*FKMJ
   23       ILO = ILO - 1
C
C  *** COMPUTE VALUE AT  X  IN (T(I),T(I+1)) OF JDERIV-TH DERIVATIVE,
C     GIVEN ITS RELEVANT B-SPLINE COEFFS IN AJ(1),...,AJ(K-JDERIV).
   30 IF (JDERIV .EQ. KM1)              GO TO 39
      JDRVP1 = JDERIV + 1     
      DO 33 J=JDRVP1,KM1
         KMJ = K-J
         ILO = KMJ
         DO 33 JJ=1,KMJ
            AJ(JJ) = (AJ(JJ+1)*DL(ILO) + AJ(JJ)*DR(JJ))/(DL(ILO)+DR(JJ))
   33       ILO = ILO - 1
   39 PVALUE = AJ(1)
C
   99                                   RETURN
      END
