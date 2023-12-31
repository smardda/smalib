C
         SUBROUTINE INTERW(PRDX,KXT,PX,KPAD,KOFF,K,KLEFT,PALPHA)
C
C.WA REPLACES INTERV FOR UNIFORM NODAL DISTRIBUTION
C.WA WITH KOFF SPLINE OFFSET IN RANGE AND KPAD PADDING
C.WA DOUBLE SEPARATION OF REPEATED END NODES FROM REST
C.WA ALSO RETURN FRACTIONAL PART
         ZX=PX*PRDX
         IP=INT(ZX)
         PALPHA=ZX-MIN(KXT-3,IP)
         IF (IP.EQ.1) PALPHA=PALPHA+1
         PALPHA=MIN(2.,PALPHA)
         IP=IP+KPAD
         KLEFT=MIN( MAX(IP-(KOFF-1),0), KXT-K ) + K
         RETURN
         END
