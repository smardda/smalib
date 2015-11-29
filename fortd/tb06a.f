C
         SUBROUTINE TB06A(N,X,F,K,NK,ETA,IL,AN,WK0,ISW,WK,A,KCALL)
C
C     *********************************************************
C     PURPOSE
C     *********************************************************
C     GIVEN FUNCTION VALUES  F(1),...,F(N)  AT THE DATA POINTS
C     X(1)<X(2)<...<X(N-1)<X(N), AND GIVEN N-K KNOTS, ETA(I),
C     I=1,...,N-K IN THE OPEN INTERVAL X(1)<X<X(N),  THIS
C     SUBROUTINE COMPUTES THE COEFFICIENTS, IN THE B-SPLINE
C     REPRESENTATION, OF THE UNIQUE SPLINE FUNCTION S(X)
C     OF DEGREE K-1, (1<=K<=N), WHICH HAS THE KNOTS ETA(I),
C     I=1,...,N-K, AND WHICH SATISFIES THE INTERPOLATION
C     CONDITIONS:
C     S(X(I)) = F(I),  I=1,...,N.
C     N.B. TO ENSURE THE EXISTENCE OF A UNIQUE INTERPOLATING
C     SPLINE, THE KNOTS MUST SATISFY THE INEQUALITIES:
C     X(1)<ETA(1)<=...<=ETA(N-K)<X(N)
C     ETA(J)<ETA(J+K), J=1,...,N-K
C     AND
C     X(J)<ETA(J)<X(J+K), J=1,...,N-K.
C     
C     NON-STANDARD:
C     VERSION ENABLES IL,AN TO BE REUSED BETWEEN CALLS
C     ASSUME THAT ALL INPUTS HAVE SENSIBLE VALUES.
C     AND DIAGNOSTICS REMOVED
C     KCALL MUST BE NON-ZERO INTEGER
C
C     *************************************************************
       REAL
     +   X(N),     F(N),     ETA(NK),  AN(N,K),  WK(ISW),  A(N),
     +   MULT
       REAL        WK0(N,K)
       INTEGER     IL(N)
       SAVE        IFIRST
       DATA        IFIRST/0/
         NORM=2
C
C     **********************************************************
C     INITIALISING THE COEFFICIENTS,OF THE B-SPLINE REPRESENTATION
C     OF S(X),TO ZERO
C     **********************************************************
         DO 1 I=1,N
    1    A(I)=0
C     END OF INITIALISATION.
C     **********************************************************
C
         NMK=N-K
         NPK=N+K
         N1=NPK+1
         N2=N1+K
         N3=N2+K
         NM1=N-1
    2    CONTINUE
C
C     **************************************************************
C     WE INTRODUCE AN ADDITIONAL 2*K KNOTS WK(I),I=1,...,K,N+1,...,N+K
C     SUCH THAT
C     WK(1)=WK(2)=...=WK(K)=X(1)
C     AND
C     WK(N+1)=WK(N+2)=...=WK(N+K)=X(N)
C     AND WE DEFINE INTERMEDIATE VALUES OF WK(I) BY THE EQUATION
C     WK(K+I)=ETA(I),  I=1,...,N-K.
C
C     ***************************************************************
         DO 3 I=1,K
         WK(I)=X(1)
         IPN=I+N
    3    WK(IPN)=X(N)
         IF(NMK.EQ.0) GOTO 5
         DO 4 I=1,NMK
         KPI=K+I
    4    WK(KPI)=ETA(I)
C     *********************************************************
C     THE KNOTS HAVE BEEN SET.
C     *********************************************************
C     **********************************************************
    5    CONTINUE
C
         IF (IFIRST.NE.KCALL) THEN
C     TOGGLE SWITCH, IF KCALL CHANGES, RE-EVALUATE COEFFICIENTS
            IFIRST=KCALL
C     SET UP THE COEFFICIENT MATRIX  AN
C     **********************************************************
C
            IF(N.EQ.2) GOTO 14
            AN(1,1)=1
            DO 99 J=2,K
            AN(1,J)=0
   99       CONTINUE
            IL(1)=K
            DO 7 I =2,NM1
            CALL TG04A(NPK,WK(1),K,NORM,X(I),IL(I),WK(N1),WK(N2))
            DO 6 J =1,K
            NPKPJ=NPK+J
    6       AN(I,J)=WK(NPKPJ)
    7       CONTINUE
C.DBG            WRITE(*,*) IFIRST, KCALL, N, K
C.DBG            WRITE(*,*) IL
C.DBG            WRITE(*,*) AN

C     **********************************************************
C     THE COEFFICIENT MATRIX  AN  HAS BEEN SET UP.
C     **********************************************************
         END IF
C
         DO 8 I =1,NM1
         DO 8 J =1,K
    8    WK0(I,J)=AN(I,J)
C
C     *********************************************************
C     TAKE A COPY OF THE FUNCTION VALUES  F(1),...,F(N)
C     *********************************************************
         DO 9 I=1,N
         N3PI=N3+I
    9    WK(N3PI)=F(I)
C     *********************************************************
C     END OF COPY.
C     *********************************************************
C
C     **********************************************************
C     BEGIN GAUSSIAN ELIMINATION
C     **********************************************************
         NM2=NM1-1
         DO 13 JCOL=1,NM2
         JCOLP1=JCOL+1
         IF(IL(JCOLP1).GE.JCOL+K) GOTO 12
         ILJC=IL(JCOL)
         JINC=K-ILJC
         J=JCOL+JINC
         PIVOT=WK0(JCOL,J)
         N3PJC=N3+JCOL
         DO 11 I=JCOLP1,NM1
         IF(IL(I).GE.JCOL+K)GOTO 12
         N3PI=N3+I
         IINC=K-IL(I)
         J1=JCOL+IINC
         MULT=WK0(I,J1)/PIVOT
         IF(JCOL.EQ.1) GOTO 11
         DO 10 J=JCOLP1,ILJC
         J1=J+JINC
         J2=J+IINC
   10    WK0(I,J2)=WK0(I,J2)-MULT*WK0(JCOL,J1)
   11    WK(N3PI)=WK(N3PI)-MULT*WK(N3PJC)
   12    CONTINUE
   13    CONTINUE
C     *********************************************************
C     END OF GAUSSIAN ELIMINATION.
C     *********************************************************
C
C     *********************************************************
C     BEGIN BACK SUBSTITUTION
C     *********************************************************
C
   14    A(1)=F(1)
         A(N)=F(N)
         IF(N.EQ.2) GOTO 18
         DO 17 IT=2,NM1
         J=N-IT+1
         N3PJ=N3+J
         ILJ=IL(J)
         JINC=K-ILJ
         J1=J+JINC
         PIVOT=WK0(J,J1)
         SUM=0
         JP1=J+1
         DO 15 JPI=JP1,ILJ
         J1=JPI+JINC
   15    SUM=SUM+WK0(J,J1)*A(JPI)
   16    A(J)=(WK(N3PJ)-SUM)/PIVOT
   17    CONTINUE
C     *********************************************************
C     END OF BACK SUBSTITUTION.
C.WA
C     *********************************************************
   18    CONTINUE
         RETURN
C     *********************************************************
C
         END
