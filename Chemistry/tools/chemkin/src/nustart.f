C     CVS $Revision: 1.1.1.1 $ created $Date: 2006/05/26 19:09:33 $

C///////////////////////////////////////////////////////////////////////
C
C     N U S T A R T
C
C     VERSION 1.00 OF DECEMBER 1997
C
C     CONVERT THE OPPDIF BINARY SOLUTION FILE INTO A RESTART FILE FOR
C     A DIFFERENT CHEMICAL MECHANISM.
C
C     BASED ON A PROGRAM WRITTEN BY FRAN RUPLEY FOR THE PREMIX AND SPIN
C     PROGRAMS.  MODIFIED FOR USE WITH OPPDIF BY ANDY LUTZ.  THIS
C     PROGRAM HAS CHANGE BLOCKS FOR USE WITH PREMIX AND SPIN, BUT ONLY
C     THE OPPDIF VERSION IS KNOWN TO FUNCTION WITH THE CURRENT VERSIONS
C     OF THE PROGRAMS.
C
C     WRITTEN BY: DR. ANDREW E. LUTZ
C                 FRAN M. RUPLEY
C                 SANDIA NATIONAL LABORATORIES
C                 MAIL STOP 9051
C                 LIVERMORE, CALIFORNIA  94551-0969  USA
C
C                 (510) 294-2761
C                 (510) 294-3657
C
C                 aelutz@california.sandia.gov
C                 fran@california.sandia.gov
C
C///////////////////////////////////////////////////////////////////////
C
C     DOCUMENTATION:
C
C     A. E. Lutz, R. J. Kee, J. F. Grcar and F. M. Rupley, "OPPDIF: A
C     Fortran Program for Computing Opposed-Flow Diffusion Flames,"
C     Sandia National Laboratories Report SAND96-8243, Livermore,
C     California, May 1997.
C
C///////////////////////////////////////////////////////////////////////
C
C     CHANGES FROM THE PREVIOUS VERSION:
C
C     1) THIS COMMENT BLOCK ADDED BY JOE GRCAR.  THE VERSION 1.00 WILL
C        NOT BE HONORED BY THE CVS REPOSITORY, OF COURSE.
C
C///////////////////////////////////////////////////////////////////////

      PROGRAM NUSTRT
      PARAMETER (LENICK=5000, LENRCK=5000, LENCCK=500,
     1           LENISK=5000, LENRSK=5000, LENCSK=500,
     2           LENIMC=500, LENRMC=50000, LINKCK=25,
     3           LINKSK=26,   LINKMC=35, KMAX=100, JMAX=500,
     4           NMAX = KMAX*JMAX,
     4           LIN=5, LOUT=6, LSAVE=14, LROUT=15,
C
C     NUSTRT VERSION NUMBER INITIATED WITH THIS VERSION
C     V 1.0 Consistent with Oppdif. Lutz, Grcar (11/9/97)
C
C*****solution > spin
C     5           NH=1, NF=2, NG=3, NT=4, NL=5, NYS=5, MAXP=10)
C*****END solution > spin
C*****solution > oppdif
     5           NT=1, NG=2, NF=3, NH=4, NYS=4, MAXP=1)
C*****END solution > oppdif
C*****solution > premix
C     5           NT=1, NYS=2, MAXP=1)
C*****END solution > premix
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      DIMENSION ICKWRK(LENICK), RCKWRK(LENRCK),
     1          ISKWRK(LENISK), RSKWRK(LENRSK),
     2          IMCWRK(LENIMC), RMCWRK(LENRMC),
     3          X(JMAX), S(NMAX), SNU(NMAX), SFRAC(KMAX),
     4          KFIRST(MAXP), KLAST(MAXP), KKPHAS(MAXP),
     5          NGIVEN(KMAX), XGIVEN(KMAX,JMAX),
     6          SGIVEN(KMAX,JMAX), XMF(KMAX,JMAX),
     7          YMF(KMAX,JMAX), VALUE(5)
C
      CHARACTER*16 KOLD(KMAX), KNEW(KMAX), CCKWRK(LENCCK),
     1             CSKWRK(LENCSK), ICHR, ISOLUT, ICKLNK,
     2             IMCLNK, ISKLNK, VERS, PREC
      CHARACTER*80 LINE
      LOGICAL IERR, KERR
      DATA ISOLUT/'SOLUTION'/,
     1     ICKLNK/'CKLINK'/, IMCLNK/'MCLINK'/,
     2     ISKLNK/'SKLINK'/, KERR/.FALSE./
      DATA X/JMAX*0.0/, S/NMAX*0.0/, SNU/NMAX*0.0/,
     1     KOLD/KMAX*' '/, KNEW/KMAX*' '/, CCKWRK/LENCCK*' '/
      DATA KSOLD,KGOLD,KKSURF,KKGAS,KKTOT,KKBULK/6*0/
C
      OPEN (LSAVE,  FORM='UNFORMATTED', STATUS='UNKNOWN',
     1      FILE='oldsave.bin')
      OPEN (LROUT,  FORM='UNFORMATTED', STATUS='UNKNOWN',
     1      FILE='newsave.bin')
      OPEN (LINKCK, FORM='FORMATTED', STATUS='UNKNOWN',
     1      FILE='chem.asc')
      OPEN (LINKMC, FORM='FORMATTED', STATUS='UNKNOWN',
     1      FILE='tran.asc')
C*****solution > spin
C      OPEN (LINKSK, FORM='UNFORMATTED', STATUS='UNKNOWN',
C     1      FILE='surf.asc')
C*****END solution > spin
C
C     READ SOLUTION FROM TAPE 14
C
      REWIND (LSAVE)
      REWIND (LROUT)
C
  320 CONTINUE
      READ (LSAVE, END=500, ERR=500) ICHR
C
      IF (ICHR .EQ. ICKLNK) THEN
C
C        OLD CHEMKIN LINKING FILE
C
         CALL CKPNT (LSAVE, LOUT, NPOINT, VERS, PREC, LENI, LENR,
     1               LENC, KERR)
         IF (LENI.LE.LENICK .AND. LENR.LE.LENRCK .AND. LENC.LE.
     1       LENCCK) THEN
            READ (LSAVE, END=500, ERR=500) (ICKWRK(L), L = 1, LENI)
            READ (LSAVE, END=500, ERR=500) (RCKWRK(L), L = 1, LENR)
            READ (LSAVE, END=500, ERR=500) (CCKWRK(L), L = 1, LENC)
C
            CALL CKINDX (ICKWRK, RCKWRK, MM, KGOLD, II, NFIT)
            KKOLD = KGOLD
            CALL CKSYMS (CCKWRK, LOUT, KOLD, KERR)
         ELSE
            WRITE (LOUT, *) ' ERROR...CKWRK, ARRAY LENGTHS TOO SMALL'
            STOP
         ENDIF
C
C        NEW CHEMKIN LINKING FILE
C
         CALL CKINIT (LENICK, LENRCK, LENCCK, LINKCK, LOUT,
     1                ICKWRK, RCKWRK, CCKWRK, IFLAG)
         IF (IFLAG .GT. 0) STOP
         CALL CKINDX (ICKWRK, RCKWRK, MM, KKGAS, II, NFIT)
         CALL CKSYMS (CCKWRK, LOUT, KNEW, KERR)
         KKTOT = KKGAS
         WRITE (LROUT) ICKLNK
         CALL CKSAVE (LOUT, LROUT, ICKWRK, RCKWRK, CCKWRK)
C
      ELSEIF (ICHR .EQ. IMCLNK) THEN
C
C        OLD TRANSPORT LINKING FILE
C
         DO 60 N = 1, 3
            READ (LSAVE, END=500, ERR=500)
   60    CONTINUE
C
C        NEW TRANSPORT LINKING FILE
C
         CALL MCINIT (LINKMC, LOUT, LENIMC, LENRMC, IMCWRK, RMCWRK,
     1                IFLAG)
         IF (IFLAG .GT. 0) STOP
         WRITE (LROUT) IMCLNK
         CALL MCSAVE (LOUT, LROUT, IMCWRK, RMCWRK)
C
C*****solution > spin
C      ELSEIF (ICHR .EQ. ISKLNK) THEN
CC
CC        OLD SURFACE LINKING FILE
CC
C         CALL SKPNT (LSAVE, LOUT, VERS, PREC, LENI, LENR,
C     1               LENC, KERR)
C         IF (LENI.LE.LENISK .AND. LENR.LE.LENRSK .AND. LENC.LE.
C     1       LENCSK) THEN
C            READ (LSAVE, END=500, ERR=500) (ISKWRK(L), L = 1, LENI)
C            READ (LSAVE, END=500, ERR=500) (RSKWRK(L), L = 1, LENR)
C            READ (LSAVE, END=500, ERR=500) (CSKWRK(L), L = 1, LENC)
CC
C            CALL SKINDX (ISKWRK, NNS, KGOLD, KSOLD, KBOLD, KKOLD,
C     1                   NNPHAS,
C     1                   NNSURF, NFSURF, NLSURF, NNBULK, NFBULK,
C     2                   NLBULK, IISUR)
C            CALL SKSYMS (ISKWRK, CSKWRK, LOUT, KOLD, KERR)
C         ELSE
C            WRITE (LOUT, *) ' ERROR...SURFACE ARRAY LENGTHS TOO SMALL'
C            STOP
C         ENDIF
CC
CC        NEW SURFACE LINKING FILE
CC
C         CALL SKINIT (LENISK, LENRSK, LENCSK, LINKSK, LOUT, ISKWRK,
C     1                RSKWRK, CSKWRK, IFLAG)
C         IF (IFLAG .GT. 0) STOP
C         CALL SKINDX (ISKWRK, NELEM, KKGAS, KKSURF, KKBULK, KKTOT,
C     1                NNPHAS, NNSURF, NFSURF, NLSURF, NNBULK, NFBULK,
C     2                NLBULK, IISUR)
C         CALL SKSYMS (ISKWRK, CSKWRK, LOUT, KNEW, KERR)
C         CALL SKPKK (ISKWRK, KKPHAS, KFIRST, KLAST)
C         WRITE (LROUT) ISKLNK
C         CALL SKSAVE (LOUT, LROUT, ISKWRK, RSKWRK, CSKWRK)
C*****END solution > spin
C
      ELSEIF (ICHR .EQ. ISOLUT) THEN
C*****solution > spin
C         READ (LSAVE, END=500, ERR=500)
C     1         NEQ, NATJ, KKS, KKB, NSS, JJ, P, PW, CN
C*****END solution > spin
C*****solution > oppdif
         READ (LSAVE, END=500, ERR=500)
     1         NATJ, JJ, P, AFUEL, AOXID, VFUEL, VOXID
         NEQ = NATJ * JJ
C*****END solution > oppdif
C*****solution > premix
C         READ (LSAVE, END=500, ERR=500) NATJ, JJ, P, FLRT
C*****END solution > premix
         READ (LSAVE, END=500, ERR=500) (X(J), J=1,JJ)
         READ (LSAVE, END=500, ERR=500) (S(N), N = 1, NEQ)
         GO TO 300
      ENDIF
      GO TO 320
C
  300 CONTINUE
C
C     COPY OLD SOLUTION TO NEW SOLUTION
C
      NMOLD = KGOLD + 2
      NMNEW = KKGAS + 2
C
      DO 330 J = 1, JJ
C*****solution > premix
C	 SNU(KKSURF+KKBULK +(NYS+KKGAS)*(J-1) +NMNEW) =
C     1     S(KSOLD + (NYS+KGOLD)*(J-1) + NMOLD)
C*****END solution > premix
	 DO 330 K = 1, NYS
            SNU(KKSURF+KKBULK + (NYS+KKGAS)*(J-1) + K) =
     1      S  (KSOLD  + (NYS+KGOLD)*(J-1) + K)
  330 CONTINUE
C
      DO 335 K = 1, KKTOT
         ICHR = ' '
         ICHR = KNEW(K)
C*****solution > spin
C         CALL SKCOMP (ICHR, KOLD, KKOLD, KNUM, NTIMES)
C*****END solution > spin
C*****solution > oppdif
         CALL CKCOMP (ICHR, KOLD, KKOLD, KNUM)
C*****END solution > oppdif
C*****solution > premix
C         CALL CKCOMP (ICHR, KOLD, KKOLD, KNUM)
C*****END solution > premix
         IF (KNUM .GT. 0) THEN
            IF (KNUM .LE. KGOLD) THEN
               DO 325 J = 1, JJ
                  SNU(KKSURF +KKBULK+ (NYS+KKGAS)*(J-1) + NYS + K)
     1              = S(KSOLD + (NYS+KGOLD)*(J-1) + NYS + KNUM)
  325          CONTINUE
            ELSE
               SNU(K-KKGAS-KKBULK) = S(KNUM-KGOLD-KBOLD)
            ENDIF
         ENDIF
  335 CONTINUE
      DO 336 K = 1, KKBULK
         SNU(KKSURF+K) = 1.0
  336 CONTINUE
C
C     READ NEW INPUT FOR SPECIES
C
      DO 550 K = 1, KMAX
         LINE   = ' '
         READ (LIN, 505, END=510) LINE
  505    FORMAT (A)
         IF (LINE(1:3) .EQ. 'END') GO TO 510
C
         CALL CKSNUM (LINE, -2, LOUT, KNEW, KKTOT, KSP,
     1                NVAL, VALUE, IERR)
         IF (IERR .OR. KSP.LE.0) THEN
            KERR = .TRUE.
            WRITE (LOUT, *) ' ERROR IN SPECIES INPUT...',LINE
         ELSE
            IF (KSP .LE. KKGAS) THEN
               IF (NVAL .EQ. 2) THEN
                  NGIVEN(KSP) = NGIVEN(KSP) + 1
                  XGIVEN(KSP,NGIVEN(KSP)) = VALUE(1)
                  IF (NGIVEN(KSP).GT.1 .AND.
     1                XGIVEN(KSP,NGIVEN(KSP)).LT.
     2                XGIVEN(KSP,NGIVEN(KSP)-1)) THEN
                      KERR = .TRUE.
                      WRITE (LOUT, *)
     1                ' ERROR...X,Y PAIRS MUST BE ORDERED...'
                  ENDIF
                  SGIVEN(KSP,NGIVEN(KSP)) = VALUE(2)
               ELSE
                  KERR = .TRUE.
                  WRITE (LOUT, *)
     1            ' ERROR...MUST GIVE X LOCATION FOR SPECIES...',
     2            LINE
               ENDIF
            ELSE
               SNU(KSP-KKGAS) = VALUE(1)
            ENDIF
         ENDIF
  550 CONTINUE
C
  510 CONTINUE
C
C        INTERPOLATE GAS-PHASE SOLUTION
C
      DO 520 J = 1, JJ
C
C        NORMALIZE MASS FRACTIONS, CONVERT TO MOLE FRACTIONS
C
         SUM = 0.0
         DO 515 K = 1, KKGAS
            YMF(K,J) = SNU(KKSURF + KKBULK+
     1                 (NYS+KKGAS)*(J-1) + NYS + K)
            SUM = SUM + YMF(K,J)
  515    CONTINUE
         IF (SUM .NE. 0.0) THEN
            DO 518 K = 1, KKGAS
               YMF(K,J) = YMF(K,J) / SUM
  518       CONTINUE
         ENDIF
C
         CALL CKYTX (YMF(1,J), ICKWRK, RCKWRK, XMF(1,J))
  520 CONTINUE
C
      DO 535 K = 1, KKGAS
         DO 535 N = 1, NGIVEN(K)
            DO 535 J = 1, JJ
               CALL TEMP (NGIVEN(K), X(J), XGIVEN(K,1), SGIVEN(K,1),
     1                      XMF(K,J))
  535 CONTINUE
C
      DO 600 J = 1, JJ
C
C        RE-NORMALIZE MOLE FRACTIONS
C
         SUMJ = 0.0
         DO 555 K = 1, KKGAS
            SUMJ = SUMJ + XMF(K,J)
  555    CONTINUE
         IF (SUMJ .GT. 0.0) THEN
            DO 560 K = 1, KKGAS
               XMF(K,J) = XMF(K,J) / SUMJ
  560       CONTINUE
         ENDIF
C
C        CONVERT TO MASS FRACTIONS
C
         CALL CKXTY (XMF(1,J), ICKWRK, RCKWRK, YMF(1,J))
         DO 570 K = 1, KKGAS
            SNU(KKSURF + KKBULK+
     1                (NYS+KKGAS)*(J-1) + NYS + K) = YMF(K,J)
  570    CONTINUE
  600 CONTINUE
C
C     NORMALIZE SITE FRACTIONS
C
      IF (KKSURF .GT. 0) THEN
         DO 650 N = 2, NNSURF+1
            SUMS = 0.0
            DO 625 K = KFIRST(N), KLAST(N)
               SUMS = SUMS + SNU(K)
  625       CONTINUE
            IF (SUMS .GT. 0.0) THEN
               DO 630 K = KFIRST(N), KLAST(N)
                  SNU(K) = SNU(K) / SUMS
  630          CONTINUE
            ENDIF
  650   CONTINUE
      ENDIF
C
C     WRITE NEW SOLUTION
C
      WRITE (LROUT) ISOLUT
C*****solution > spin
CC      NSSOLV = 0
C      NSSOLV = NSS
C      NATJ = KKGAS + NYS
C      NEQ = KKSURF + 2*KKBULK + NSSOLV + NATJ*JJ
C      WRITE (LROUT) NEQ, NATJ, KKSURF, KKBULK, NSSOLV, JJ, P,
C     1              PW, CN
C*****END solution > spin
C*****solution > oppdif
      NATJ = KKGAS + NYS
      NEQ = NATJ*JJ
      WRITE (LROUT)
     1         NATJ, JJ, P, AFUEL, AOXID, VFUEL, VOXID
C*****END solution > oppdif
C*****solution > premix
C      NATJ = KKGAS + NYS
C      NEQ = NATJ*JJ
C      WRITE (LROUT) NATJ, JJ, P, FLRT
C*****END solution > premix
      WRITE (LROUT) (X(J), J = 1, JJ)
      WRITE (LROUT) (SNU(N), N = 1, NEQ)
      REWIND (LSAVE)
      REWIND (LROUT)
      STOP
C
  500 CONTINUE
  700 FORMAT (1X,A)
      WRITE (LOUT, '(A)')' Error reading solution file, ICHR='//ICHR
      END
      SUBROUTINE TEMP(NPTS, X, XX, TT, T)
C
C*****precision > double
        IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C
C*****precision > single
C        IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      DIMENSION XX(*), TT(*)
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C     THIS SUBROUTINE USES BISECTION TO LINEARLY INTERPOLATE
C     AN ARRAY OF XX,TT PAIRS.  GIVEN AN XX,TT PAIR THIS ROUTINE
C     RETURNS THE INTERPOLATED VALUE OF THE T AT THE POINT X.
C
C INPUT-
C   NPTS   - NUMBER OF XX,TT PAIRS.
C   X      - LOCATION AT WHICH INTERPOLATED T IS DESIRED.
C   XX     - ARRAY OF X POINTS AT WHICH TT ARE GIVEN.
C   TT     - ARRAY OF T VALUES AT THE XX LOCATIONS.
C
C OUTPUT-
C   T     - INTERPOLATED T AT POINT X
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C             check for x outside (1,npts)
C
      IF (X .LE. XX(2)) THEN
         N = 2
         S = (TT(N) - TT(N-1)) / (XX(N) - XX(N-1))
      ELSEIF (X .GE. XX(NPTS-1)) THEN
         N = NPTS-1
         S = (TT(N+1) - TT(N)) / (XX(N+1) - XX(N))
      ELSE
         NLO = 1
         NHI = NPTS
         S   = 0.0
C
C        bisect interval
C
50       CONTINUE
         N = (NLO+NHI)/2
         IF (X .LT. XX(N)) THEN
            IF (X .LT. XX(N-1)) THEN
               NHI = N
               GO TO 50
            ELSEIF (X .EQ. XX(N-1)) THEN
               N = N-1
            ELSE
               S = (TT(N) - TT(N-1)) / (XX(N) - XX(N-1))
            ENDIF
         ELSEIF (X .GT. XX(N)) THEN
            IF (X .GT. XX(N+1)) THEN
               NLO = N
               GO TO 50
            ELSEIF (X .EQ. XX(N+1)) THEN
               N = N + 1
            ELSE
               S = (TT(N+1) - TT(N)) / (XX(N+1) - XX(N))
            ENDIF
         ENDIF
      ENDIF
C
  100 CONTINUE
      T      = TT(N) + S * (X - XX(N))
      RETURN
      END
