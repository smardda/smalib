C
         PROGRAM STIVTK
C
C     STITCH VTK FILES TOGETHER
C     READ TEXTFILE CONTAINING NAMES ON NLEDGE
C     FINAL OUTPUT TO NEW FILE ON NPUNCH
       CHARACTER   *80 FILENM,         NEWFIL,   FILOUT
       CHARACTER   *80 RECLN,          RECORD,   RECLIN,   RLINE
       CHARACTER   *80 RECHED(5)
       INTEGER
     +   ICOUNT,   IENBL,    IMAXU,    IUNIT,    IPOSN,    IRECL,
     +   ININC ,   ISTAT,    ISTERR
      REAL, DIMENSION(:,:), ALLOCATABLE :: RPTS, ZPTS 
      INTEGER, DIMENSION(:), ALLOCATABLE :: NODL, INODL
      INTEGER, DIMENSION(:), ALLOCATABLE :: TYPES, ITYPES 
      INTEGER, DIMENSION(:), ALLOCATABLE :: NDATA, IDATA 
      REAL, DIMENSION(:), ALLOCATABLE :: RDATA, ZDATA 
       LOGICAL     EXST
         NLEDGE=8
         NPUNCH=7
         IUNIT=NLEDGE+1
         IFILE=0
C
C-----------------------------------------------------------------------
CL              1.         Open file
C
C     Read in number of arguments
         INARGS=IARGC()
         IF (INARGS.GT.0) THEN
            CALL GETARG(1,FILENM)
            CALL LENWD(FILENM,IEND)
            ISTERR=IEND
            OPEN(UNIT=NLEDGE,FILE=FILENM(1:IEND),
     +     STATUS='OLD',IOSTAT=ISTAT)
            IF(ISTAT.NE.0)THEN
               PRINT*,'File ',FILENM(1:IEND),' does not exist'
               GO TO 403
            END IF
            IF (INARGS.EQ.2) THEN
               CALL GETARG(2,FILOUT)
            ELSE
               FILOUT=FILENM(1:IEND)//'.sti.vtk'
            END IF
            INQUIRE(FILE=FILOUT,EXIST=EXST)
            IF(EXST)THEN
               PRINT*,'***Warning output file already exists'
            END IF
            OPEN(UNIT=NPUNCH,FILE=FILOUT,IOSTAT=ISTAT)
            IF(ISTAT.NE.0)THEN
               PRINT*,'Opening file ',FILOUT,' gives error'
               GO TO 403
            END IF
         ELSE
            PRINT*,'***No input file given'
            STOP
         END IF
  100    CONTINUE
C     Read in next line
         READ(NLEDGE,'(A80)',END=400,ERR=501)RECORD
C     Transfer line to buffer
         NEWFIL=RECORD
C
C-----------------------------------------------------------------------
CL              2.         Included files
C
  200    CALL XBLCH(NEWFIL,1,80,IENDP)
         IF(IENDP .EQ. 0)THEN
            PRINT*,'***Warning-No filename found in ',RECORD
            GOTO 403
         ENDIF
C     Inquire if file exists
  201    INQUIRE(FILE=NEWFIL,EXIST=EXST)
         IF(EXST)THEN
            OPEN(IUNIT,FILE=NEWFIL,STATUS='OLD',ERR=500)
            IFILE=IFILE+1
         ELSE
            PRINT*,'File ',NEWFIL(1:IENDP),' does not exist'
            GO TO 100
         ENDIF
C
CL                  2.1      read in data into scratch arrays
  210    CONTINUE
         DO 211 J=1,4
         IF (IFILE.EQ.1) THEN
            READ(IUNIT,'(A80)',END=212,ERR=502)RECHED(J)
         ELSE
            READ(IUNIT,'(A80)',END=212,ERR=502)RECLIN
         END IF
  211    CONTINUE
         CALL VFILE_COPY(IUNIT,NEWFIL,ZPTS,INODL,ITYPES,IDATA,ZDATA)
C
         GOTO 210
  212    CONTINUE
C
C-----------------------------------------------------------------------
CL              3.         Add scratch arrays to main arrays
C
         CALL VFILE_ADDARR(ZPTS,INODL,ITYPES,IDATA,ZDATA,
     +   RPTS,NODL,NTYPES,NDATA,RDATA)
C     get new file
         CLOSE(IUNIT)
         GO TO 100
C
C-----------------------------------------------------------------------
CL              4.         Output arrays in vtk format
C
  400    CONTINUE
         DO 401 J=1,4
         WRITE(NPUNCH,'(A80)',ERR=503)RECHED(J)
  401    CONTINUE
         CALL VFILE_WRTBODY(NPUNCH,FILOUT,RPTS,NODL,NTYPES,NDATA,RDATA)
C
C     Close units
  402    CLOSE(NPUNCH)
  403    CLOSE(NLEDGE)
         STOP
C
C-----------------------------------------------------------------------
CL              5.         Errors
C
  500    PRINT*,'***Error opening file ',NEWFIL
         GO TO 100
  501    PRINT*,'***Error reading in file ',FILENM(1:ISTERR)
         PRINT*,'Processing stopped'
         GO TO 402
  502    PRINT*,'***Error reading in file ',NEWFIL
         PRINT*,'Processing continuing on other files'
         CLOSE(IUNIT)
  503    PRINT*,'***Error writing file ',NEWFIL
         PRINT*,'Run abandoned'
         GO TO 402
         END
C
         SUBROUTINE XBLCH(KCSTR,KSTA,KEND,KEO)
C
C     Compress string, removing blanks
       CHARACTER   *(*) KCSTR
C
C-----------------------------------------------------------------------
CL              1.         Check arguments
C
         ISLEN=LEN(KCSTR)
         IEND=MAX(KSTA,KEND)
         IF (KSTA.LT.1.OR.IEND.GT.ISLEN) THEN
            WRITE(*,*) KSTA,KEND,ISLEN,'*** RUN ABANDONED IN XBLCH'
            STOP
         END IF
C
C-----------------------------------------------------------------------
CL              2.         remove blanks
C
         I=0
         DO 200 IENBL=KSTA,KEND
         IF (KCSTR(IENBL:IENBL).EQ.' ') GO TO 200
         I=I+1
         KCSTR(I:I)=KCSTR(IENBL:IENBL)
  200    CONTINUE
         KEO=I
C-----------------------------------------------------------------------
         RETURN
         END
C
         SUBROUTINE LENWD(KCH,KEND)
C
C     CALCULATE LENGTH OF WORD STARTING AT CHARACTER 1
       CHARACTER   *(*) KCH
         ISLEN=LEN(KCH)
         DO 1 IENBL=1,ISLEN
         IF (KCH(IENBL:IENBL).EQ.' ') GO TO 2
    1    CONTINUE
    2    CONTINUE
         KEND=IENBL-1
         RETURN
         END
