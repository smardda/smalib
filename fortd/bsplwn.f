C
         SUBROUTINE BSPLWN(KLEFT,PALPHA,KXT,PCOEFF)
C
C Evaluate B-splines using explicit polynomial
C Rember that for 1st and 5th options in KLEFT order
C INTERW must allow alpha up to 2
C KXT = N1+1 on input
         REAL PCOEFF(*)
C
C 1) 43: part (ws(30),3);
C
      IF (KLEFT.LE.4) THEN
         PCOEFF(1)=(-(PALPHA-2)**3)/8
         PCOEFF(2)=((19*PALPHA**2-90*PALPHA+108)*PALPHA)/72
         PCOEFF(3)=(-(13*PALPHA-36)*PALPHA**2)/72
         PCOEFF(4)=PALPHA**3/24
C
C
C 2) 44:  part(ws(31),3);
C
      ELSE IF (KLEFT.EQ.5) THEN
         PCOEFF(1)=(-(PALPHA-1)**3)/9
         PCOEFF(2)=(23*PALPHA**3-42*PALPHA**2-12*PALPHA+40)/72
         PCOEFF(3)=(-(9*PALPHA**3-6*PALPHA**2-12*PALPHA-8))/24
         PCOEFF(4)=PALPHA**3/6
C
C 3) 45:  part(ws(32),3);
C
      ELSE IF (KLEFT.EQ.6) THEN
         PCOEFF(1)=(-(PALPHA-1)**3)/8
         PCOEFF(2)=(11*PALPHA**3-21*PALPHA**2-3*PALPHA+17)/24
         PCOEFF(3)=(-(3*PALPHA**3-3*PALPHA**2-3*PALPHA-1))/6
         PCOEFF(4)=PALPHA**3/6
C
C 5) 47:  part(ws(34),3);
C
      ELSE IF (KLEFT.EQ.KXT) THEN
         PCOEFF(1)=(-(PALPHA-2)**3)/24
         PCOEFF(2)=((13*PALPHA+10)*(PALPHA-2)**2)/72
         PCOEFF(3)=(-(19*PALPHA**2+14*PALPHA+4)*(PALPHA-2))/72
         PCOEFF(4)=PALPHA**3/8
C
C 6) 48:  part(ws(35),3);
C
      ELSE IF (KLEFT.EQ.KXT-1) THEN
         PCOEFF(1)=(-(PALPHA-1)**3)/6
         PCOEFF(2)=(9*PALPHA**3-21*PALPHA**2+3*PALPHA+17)/24
         PCOEFF(3)=(-(23*PALPHA**3-27*PALPHA**2-27*PALPHA-9))/72
         PCOEFF(4)=PALPHA**3/9

C
C 7) 49:  part(ws(36),3);
C
      ELSE IF (KLEFT.EQ.KXT-2) THEN
         PCOEFF(1)=(-(PALPHA-1)**3)/6
         PCOEFF(2)=(3*PALPHA**3-6*PALPHA**2+4)/6
         PCOEFF(3)=(-(11*PALPHA**3-12*PALPHA**2-12*PALPHA-4))/24
         PCOEFF(4)=PALPHA**3/8
C
C 4) 46:  part(ws(33),3);
C
      ELSE 
C        IF (KLEFT.EQ.7:KXT-2) THEN
         ZALFM1=1-PALPHA
         PCOEFF(1)=ZALFM1**3/6
         PCOEFF(2)=(3*PALPHA**3-6*PALPHA**2+4)/6
         PCOEFF(3)=(3*ZALFM1**3-6*ZALFM1**2+4)/6
         PCOEFF(4)=PALPHA**3/6
C
      END IF
C
      RETURN
      END
