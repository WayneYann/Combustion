C     CVS $Revision: 1.1.1.1 $ created $Date: 2006/05/26 19:09:33 $

C///////////////////////////////////////////////////////////////////////
C
C     P O S T O P P
C
C     VERSION 8.2 OF 04/17/98
C
C     READ THE OPPDIF BINARY SOLUTION FILE AND WRITE TEXT FILES FOR
C     PLOTTING.
C
C     BASED ON A PROGRAM WRITTEN BY FRAN RUPLEY.  MODIFIED FOR USE WITH
C     OPPDIF BY ANDY LUTZ.
C
C     WRITTEN BY: DR. ANDREW E. LUTZ
C                 SANDIA NATIONAL LABORATORIES
C                 MAIL STOP 9051
C                 LIVERMORE, CALIFORNIA  94551-0969  USA
C
C                 (510) 294-2761
C                 (510) 294-3657
C
C                 aelutz@california.sandia.gov
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
C     CHANGES FOR v. 8.2
C     Action #0170: (E. Meeks)
C     1) Replaced KERR with IERR (logical) in all occurences in 
C        the mail routine.
C     2) Removed unused variables H, WDOT, X from main routine.
C     3) Fixed open statements for IFORM=2 option in RDSAVE
C
C     CHANGES FROM THE PREVIOUS VERSION:
C     1) THIS COMMENT BLOCK ADDED BY JOE GRCAR.
C
C///////////////////////////////////////////////////////////////////////

C     VERSION 1.1: PLOTS SOLUTION FROM OPPDIF
C     VERSION 1.2: PLOTS FUEL MIXTURE FRACTION
C     VERSION 1.3: INTEGRATES AND PLOTS EMISSION INDEX  4/17/92
C     VERSION 1.4: CALC FLAME THICKNESS 4/22/92
C     VERSION 1.5: INTEGRATE TAKENO'S EI DEFINITIION  4/22/92
C     VERSION 1.6: EFFECTIVE EI DEFINITIION,
C                  PRINT DATA FILES TO UNIT 11. 4/24/92
C     V 2.0: TEXT POST-PROCESSING FOR OPPDIF     8/94
C     V 2.1: READ SAVE OR RECOVER FILES          9/94
C     V 2.2: IMPROVED FILE NAMING, NEW MIXTURE FRACTION DEFINITION 12/95
C     V 2.3: FIXED MINOR BUGS FOR UNDECLARED LOGICALS  11/97
 
C*****single precision
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END single precision
C*****double precision
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END double precision

      PARAMETER (KMAX=60, IMAX=260, JMAX=250, LENRCK=9500, LENICK=9500,
     1           LENCCK=500, LENRMC=60 000, LENIMC=250,
     2           LIN=5, LOUT=6, LSAVE=55, LTEXT=11,
     4           NT=1, NG=2, NF=3, NH=4, NY1=5)

      DIMENSION
     1  ICKWRK(LENICK), IMCWRK(LENIMC),
     2  T(JMAX), XCM(JMAX), S(KMAX,JMAX),
     3  SN(KMAX,JMAX), CONT(JMAX,IMAX), XMF(JMAX,KMAX),
     4  XT(KMAX), Q(IMAX), CIK(IMAX), SUMRT(JMAX),
     5  RCKWRK(LENRCK), RMCWRK(LENRMC),
     6  RHO(JMAX), V(JMAX), EDOTS(JMAX), EDOTN(JMAX),
     7  WT(KMAX), FMIXT(JMAX), EI(JMAX,5), WTM(JMAX), U(JMAX),
     8  KSPEC(KMAX), KPROD(KMAX), ILIST(IMAX)

      CHARACTER*16 KSYM(KMAX), CCKWRK(LENCCK), ELEM(4)
      CHARACTER*48 FILNAM
      CHARACTER*32 IREAC(IMAX)
      CHARACTER*8  PREFIX

      LOGICAL IERR, MOL, EINOX

      COMMON /MIX/ AWT(4), NCF(4,200), NEC, NEH, NEO
      DATA MDIM/4/

      DATA MOL/.TRUE./, EINOX/.FALSE./, NSOL/1/

      WRITE (LOUT, '(//5X,A/5X,A//)' ) 
     1 'POSTOPP: POST-PROCESSING FOR OPPDIF' ,
     2 '         VERSION 8.2, (04/98)'
      WRITE (LOUT, '(//5X,A/5X,A/5X,A)' ) 
     1 'Read file formatted as either:',
     2 '(1) save',
     3 '(2) recover'
      READ (LIN,*) IFORM
      WRITE (LOUT, '(/,5X,A)') 'Enter file name'
      READ (LIN, '(A)') FILNAM
      WRITE (LOUT,'(5X,A//)') 'READING FILE...'
      IF (IFORM .EQ. 1) THEN
        OPEN (LSAVE, STATUS='OLD', FORM='UNFORMATTED', FILE=FILNAM)
      ELSEIF (IFORM .EQ. 2) THEN
        OPEN (LSAVE, STATUS='OLD', FORM='UNFORMATTED', FILE=FILNAM)
      ELSE
        WRITE (LOUT, *) 'Stop - Form not found.'
        STOP
      ENDIF
      WRITE (LOUT,'( /5X,A)') 'Enter prefix for text filenames:'
      READ (LIN,'(A)') PREFIX
      CALL SQUEEZ (ILEN, PREFIX)

C     READ BINARY FILE
      ISEN = 0
      CALL RDSAVE (LSAVE, LOUT, LENICK, ICKWRK, LENRCK, RCKWRK,
     1             LENCCK, CCKWRK, LENIMC, IMCWRK, LENRMC, RMCWRK,
     2             NSOL, ISEN, KMAX, JMAX, IMAX, II, KSYM, IREAC,
     3             KK, JJ, XCM, S, SN, CONT, P, IFORM)

      CALL CKRP (ICKWRK, RCKWRK, RU, RUC, PATM)
      CALL CKWT (ICKWRK, RCKWRK, WT)
      CALL CKAWT (ICKWRK, RCKWRK, AWT)
      CALL CKNCF (MDIM, ICKWRK, RCKWRK, NCF)
      CALL CKSYMS (CCKWRK, LOUT, KSYM, IERR)
      IF (IERR) THEN
        WRITE (LOUT,'(//5X,A)') 'ERROR GETTING SPECIES NAMES!'
        STOP
      ENDIF
      CALL CKSYME (CCKWRK, LOUT, ELEM, IERR)
      IF (IERR) THEN
        WRITE (LOUT,'(//5X,A)') 'ERROR GETTING ELEMENT NAMES!'
        STOP
      ENDIF

C     READ KEYWORD INPUT
      CALL RDKEY (LIN, LOUT, NSOL, KSYM, KK, 
     1            KSPEC, NSPEC, KPROD, NPROD, MOL, EINOX, 
     2            KFUEL)

      CALL CKCOMP ('H', ELEM, 4, NEH)
      CALL CKCOMP ('C', ELEM, 4, NEC)
      CALL CKCOMP ('O', ELEM, 4, NEO)

      WRITE (6,'(/5X,A,1PE13.3)') 'Pressure eigenvalue =', S(NH,1) 

      DO 100 J = 1, JJ

         T(J) = S(NT,J)

C        COMPUTE FUEL MIXTURE FRACTION
         CALL MIXTRF ( S(NY1,J), KK, WT, FMIXT(J) )
C         WRITE (6,'(/,A,I3,1X,1PE15.7,/)') 'FMIXT(J) = ',J, FMIXT(J)

C        COMPUTE DENSITY
         CALL CKRHOY (P, S(NT,J), S(NY1,J), ICKWRK, RCKWRK, RHO(J))
         WTM(J) = RHO(J) *RU *S(NT,J) / P

C        RADIAL STRAIN RATE: V = v/r = G / rho
         V(J) = - S(NG,J)/RHO(J)

C        AXIAL VELOCITY
         U(J) = 2.*S(NF,J)/RHO(J)

C        CONVERT TO MOL FRACTIONS
         CALL CKYTX (S(NY1,J), ICKWRK, RCKWRK, XT)
         DO 50 K = 1, KK
            XMF(J,K) = XT(K)
50       CONTINUE

100   CONTINUE

      DO 120 J = 1, JJ-1
         DX = XCM(J+1) - XCM(J)

C        NORMAL STRAIN RATE: du/dx
         EDOTN(J) = ( U(J+1) - U(J) )/DX

C        SHEAR STRAIN RATE: d(G/rho)/dx
         EDOTS(J) = ( V(J+1) - V(J) )/DX
120   CONTINUE
      EDOTN(JJ) = 0.
      EDOTS(JJ) = 0.

C        MEASURE FLAME THICKNESS

      CALL MINMAX (T, JJ, TMN, TMX, IPK)
      T10 = TMN + (TMX-TMN)*0.1
      X1 = 0.
      X2 = 0.
      DO 200 J = 1, JJ
        IF (X1.EQ.0.) THEN
          IF (T(J).GE.T10) X1 = XCM(J)
        ELSE
          IF (T(J).LE.T10) THEN
            X2 = XCM(J)
            GOTO 220
          ENDIF
        ENDIF
200   CONTINUE
      WRITE (LOUT,'(A)') ' WARNING, FLAME THICKNESS NOT FOUND?'
220   CONTINUE
      FLAMET = X2 - X1
      WRITE (LOUT,'(/,A,1PE14.7)') ' FLAME THICKNESS (CM)=', FLAMET
      WRITE (LOUT,'(A,1PE14.7)')   ' FOR TEMPERATURE (K) =', T10
      WRITE (LOUT,'(A,1PE14.7,/)') ' MAXIMUM TEMP    (K) =', TMX

C     MEASURE AVERAGE NORMAL STRAIN

      X1S = X1 * 0.75
      X2S = X2 * 1.25
      DO 230 J = 2, JJ
        IF ( (X1S .GE. XCM(J-1)) .AND. (X1S .LE. XCM(J)) ) THEN
          U1S = U(J) + (U(J)-U(J-1))*(X1S-XCM(J-1))/(XCM(J)-XCM(J-1))
          GOTO 231
        ENDIF
230   CONTINUE
      WRITE (LOUT,'(//2X,A)') 'X1S NOT FOUND!'    
231   CONTINUE
      DO 240 J = 2, JJ
        IF ( (X2S .GE. XCM(J-1)) .AND. (X2S .LE. XCM(J)) ) THEN
          U2S = U(J) + (U(J)-U(J-1))*(X2S-XCM(J-1))/(XCM(J)-XCM(J-1))
          GO TO 241
        ENDIF
240   CONTINUE
      WRITE (LOUT,'(//2X,A)') 'X2S NOT FOUND!'
241   CONTINUE

      EDOTAV = (U2S - U1S)/(X2S - X1S) 
      WRITE (LOUT,'(//2X,A,1PE15.6/)')
     1  'AVERAGE NORMAL STRAIN (1/S) =', EDOTAV
      

C      WRITE TEMP AND STUFF

      FILNAM = PREFIX(1:ILEN)//'.tmp'
      OPEN (UNIT=LTEXT,FORM='FORMATTED',STATUS='UNKNOWN',FILE=FILNAM)
      WRITE (LTEXT, '(1X,9(1X,A))')  'X(cm)', 'T(K)', 
     1  'MIXTURE_FRACTION',
     1  'RHO(g/cc)',  'Mean_WT(g/mol)', 'u(cm/s)', 
     2  'v/r(1/s)', 
     3  'Normal-strain-rate(1/s)',  'Shear-strain-rate(1/cm-s)'
     
      DO 300 J = 1, JJ
        WRITE (11,'(1X,F9.7, 8(1X,1PE12.5))')
     1  XCM(J), S(NT,J), FMIXT(J), RHO(J), WTM(J), U(J), V(J),
     2  EDOTN(J), EDOTS(J)
300   CONTINUE
      CLOSE (LTEXT)

C     WRITE SPECIES

      FILNAM = PREFIX(1:ILEN)//'.spc'
      OPEN (UNIT=LTEXT,FORM='FORMATTED',STATUS='UNKNOWN',FILE=FILNAM)
      K1 = 1
      NPERL = 9
1000  CONTINUE
        K2 = MIN( K1 + NPERL - 1, NSPEC )
        IF (MOL) THEN
          WRITE (11,'(1X,A)') 'MOL FRACTIONS:'
          WRITE (11,'(10(1X,A))')
     1    'X(cm)', (KSYM(KSPEC(K)), K = K1, K2)
          DO 1100 J = 1, JJ          
            WRITE (11,'(1X,F9.7, 9(1X,1PE12.5))')
     1      XCM(J), (XMF(J,KSPEC(K)), K = K1, K2)
1100      CONTINUE
        ELSE
          WRITE (11,'(1X,A)') 'MASS FRACTIONS:'
          WRITE (11,'(10(1X,A))')
     1    'X(cm)', (KSYM(KSPEC(K)), K = K1, K2)
          DO 1200 J = 1, JJ          
            WRITE (11,'(1X,F9.7, 9(1X,1PE12.5))')
     1      XCM(J), (S(NH+KSPEC(K),J), K = K1, K2)
1200      CONTINUE
        ENDIF
        IF (K2 .LT. NSPEC) THEN
          K1 = K1 + NPERL
          GO TO 1000
        ENDIF

      CLOSE (LTEXT)



      IF (EINOX) THEN

C       EMISSION INDEX

        CALL EMISS (V, EI, FMIXT, JMAX, KK, KSYM, LOUT, JJ, RHO,
     1               U, WT, WTM, XCM, XMF, KFUEL)

        FILNAM = PREFIX(1:ILEN)//'.ei'
        OPEN (UNIT=LTEXT,FORM='FORMATTED',STATUS='UNKNOWN',FILE=FILNAM)
        WRITE (LTEXT, '(4(3X,A))')
     1  'X(cm)', 'NOx(g/cm^2-s)', 'NO2(g/cm^2-s)', 'Fuel(g/cm^2-s)' 
        DO 1300 J = 1, JJ
          WRITE (11,'(1X,F9.7, 6(1X,1PE12.5))')
     1    XCM(J), ( EI(J,I), I = 1, 3 )
1300    CONTINUE
        CLOSE (LTEXT)

      ENDIF
     
      IF (NPROD .GE. 0) THEN

        CALL PDRATE (KPROD, NPROD, XCM, ICKWRK, RCKWRK, P, JJ,
     1                   KK, II, IREAC, ILIST, XT, KMAX, JMAX, KSYM,
     2                   S, NY1, NT, Q, CIK, CONT, SUMRT, LTEXT, WT,
     3                   PREFIX, ILEN )

      ENDIF

      STOP
      END


C----------------------------------------------------------------------C
      SUBROUTINE MINMAX (V, JJ, VMIN, VMAX, IPK)
 
C*****single precision
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END single precision
C*****double precision
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END double precision

      DIMENSION V(*)
      VMAX = V(1)
      VMIN = V(1)
      ILO = 1
      IPK = 1
      DO 100 J=2,JJ
         IF (V(J) .LT. VMIN) THEN
           ILO = J
           VMIN = V(J)
         ENDIF
         IF (V(J) .GT. VMAX) THEN
           IPK = J
           VMAX = V(J)
         ENDIF
  100 CONTINUE
      IF (ABS(VMIN) .GT. ABS(VMAX)) IPK = ILO
      RETURN
      END

C----------------------------------------------------------------------C

      SUBROUTINE RDSAVE (LS, LOUT, LENICK, ICKWRK, LENRCK, RCKWRK,
     1                   LENCCK, CCKWRK, LENIMC, IMCWRK, LENRMC,
     2                   RMCWRK, NSOL, ISEN, KMAX, JMAX, IMAX,
     3                   II, KSYM, IREAC, KK, JJ, X, S, SN, CONT, P,
     4                   IFORM)
C

C*****single precision
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END single precision
C*****double precision
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END double precision

      DIMENSION 
     1 ICKWRK(*), IMCWRK(*), X(*), S(KMAX,*), SN(KMAX,*),
     2 CONT(JMAX,*), RCKWRK(*), RMCWRK(*)

      CHARACTER ICHR*16, VERS*16, PREC*16, CDUM*16, CCKWRK(*)*(*),
     1          KSYM(*)*(*), IREAC(*)*(*)
      LOGICAL KERR, IERR

C     IF READING RECOVER FILE, INITIALIZE CHEMKIN FROM LINKING FILE

      IF (IFORM .EQ. 2) THEN
         LINKCK = 25
C*****linking files > binary
C        OPEN (LINKCK, FORM='UNFORMATTED', FILE='./chem.bin')
C*****end linking files > binary
C*****linking files > ascii
      OPEN(LINKCK,STATUS='OLD',FORM='FORMATTED',FILE='./chem.asc')
C*****end linking files > ascii
        CALL CKLEN (LINKCK, LOUT, LENI, LENR, LENC, IFLAG)
        IF (IFLAG .NE. 0) THEN
           WRITE (LOUT,*) ' STOP. BAD CKLEN!'
           STOP
        ENDIF
        IF (LENI.LE.LENICK .AND. LENR.LE.LENRCK .AND. LENC.LE.LENCCK)
     1     THEN
           CALL CKINIT (LENICK, LENRCK, LENCCK, LINKCK, LOUT, 
     1                  ICKWRK, RCKWRK, CCKWRK, IFLAG)
           IF (IFLAG .NE. 0) THEN
             WRITE (LOUT,*) ' STOP. BAD CKINIT!'
             STOP
           ENDIF
           CALL CKINDX (ICKWRK, RCKWRK, MM, KK, II, NFIT)
C
           IF (KK .LE. KMAX) THEN
              CALL CKSYMS (CCKWRK, LOUT, KSYM, IERR)
           ELSE
              WRITE (LOUT, *) ' Error...species dimension too small,',
     1                        ' KMAX should be at least ', KK
              KERR = .TRUE.
           ENDIF
C
           IF (II .LE. IMAX) THEN
              DO 5 I = 1, II
                 CALL CKSYMR (I, LOUT, ICKWRK, RCKWRK, CCKWRK, LT,
     1                        IREAC(I), IERR)
5             CONTINUE
           ELSE
              WRITE (LOUT, *) ' Error...reaction dimension too small,',
     1                        ' IMAX should be at least ',II
              KERR = .TRUE.
           ENDIF
        ELSE
           WRITE (LOUT,'(//5X,A)') 'STOP. REQUIRE CHEMKIN WORK SPACE:'
           WRITE (LOUT,'(5X,A,2X,I10)') 'INTEGER  ', LENI
           WRITE (LOUT,'(5X,A,2X,I10)') 'REAL     ', LENR
           WRITE (LOUT,'(5X,A,2X,I10)') 'CHARACTER', LENC
           STOP
        ENDIF
      ENDIF

      REWIND (LS)
      NS = 0
      KERR = .FALSE.
   10 CONTINUE
      ICHR = ' '
      READ (LS, ERR=100, END=100) ICHR

      IF (ICHR .EQ. 'CKLINK') THEN
C
C        Initialize Chemkin Library common block
C
         CALL CKPNT (LS, LOUT, NPOINT, VERS, PREC, LENI, LENR, LENC,
     1               IERR)
         KERR = KERR.OR.IERR
C
C*****double precision
         IF (INDEX(PREC,'DOUB') .LE. 0) THEN
C*****END double precision
C*****single precision
C         IF (INDEX(PREC,'SING') .LE. 0) THEN
C*****END single precision
C
            WRITE (LOUT, 600) ICHR
            STOP
        ENDIF
C
C       Initialize Chemkin work arrays
C
        IF (LENI.LE.LENICK .AND. LENR.LE.LENRCK .AND. LENC.LE.LENCCK)
     1     THEN
           READ (LS, ERR=200, END=200) (ICKWRK(L), L = 1, LENI)
           READ (LS, ERR=200, END=200) (RCKWRK(L), L = 1, LENR)
           READ (LS, ERR=200, END=200) (CCKWRK(L), L = 1, LENC)
           CALL CKINDX (ICKWRK, RCKWRK, MM, KK, II, NFIT)
C
           IF (KK .LE. KMAX) THEN
              CALL CKSYMS (CCKWRK, LOUT, KSYM, IERR)
              WRITE (LOUT, '(/,4(2X,A))') (KSYM(K),K=1,KK)
           ELSE
              WRITE (LOUT, *) ' Error...species dimension too small,',
     1                        ' KMAX should be at least ', KK
              KERR = .TRUE.
           ENDIF
C
           IF (II .LE. IMAX) THEN
              DO 20 I = 1, II
                 CALL CKSYMR (I, LOUT, ICKWRK, RCKWRK, CCKWRK, LT,
     1                        IREAC(I), IERR)
   20         CONTINUE
           ELSE
              WRITE (LOUT, *) ' Error...reaction dimension too small,',
     1                        ' IMAX should be at least ',II
              KERR = .TRUE.
           ENDIF
        ELSE
           READ (LS, ERR=200, END=200) (IDUM, L = 1, LENI)
           READ (LS, ERR=200, END=200) (RDUM, L = 1, LENR)
           READ (LS, ERR=200, END=200) (CDUM, L = 1, LENC)
           WRITE (LOUT, 700)
     1        ICHR, LENI, LENICK, LENR, LENRCK, LENC, LENCCK
           KERR = .TRUE.
        ENDIF
C
      ELSEIF (ICHR .EQ. 'MCLINK') THEN
C
C        Initialize Transport Library common block
C
         CALL MCPNT (LS, LOUT, NPOINT, VERS, PREC, LENI, LENR, IERR)
         KERR = KERR.OR.IERR
C
C*****double precision
         IF (INDEX(PREC,'DOUB') .LE. 0) THEN
C*****END double precision
C*****single precision
C         IF (INDEX(PREC,'SING') .LE. 0) THEN
C*****END single precision
C
            WRITE (LOUT, 600) ICHR
            STOP
         ENDIF
C
C        Initialize Transport work arrays
C
        IF (LENI.LE.LENIMC .AND. LENR.LE.LENRMC) THEN
           READ (LS, ERR=200, END=200) (IMCWRK(L), L = 1, LENI)
           READ (LS, ERR=200, END=200) (RMCWRK(L), L = 1, LENR)
        ELSE
           READ (LS, ERR=200, END=200) (IDUM, L = 1, LENI)
           READ (LS, ERR=200, END=200) (RDUM, L = 1, LENR)
           WRITE (LOUT, 701) ICHR, LENI, LENIMC, LENR, LENRMC
           KERR = .TRUE.
        ENDIF
C
      ELSEIF (ICHR .EQ. 'SOLUTION') THEN
         IF (KERR) STOP
C
C*****OLD OPPDIF
C         READ (LS, ERR=300, END=300) NATJ, JJ, P
C         READ (LS, ERR=300, END=300) VFUEL, VOXID
C*****END OLD OPPDIF
C*****NEW OPPDIF
         READ (LS, ERR=300, END=300)
     1   NATJ, JJ, P, AFUEL, AOXID, VFUEL, VOXID
C*****END NEW OPPDIF
C
         IF (JJ .GT. JMAX) THEN
            WRITE (LOUT, *) ' Error...node dimension too small,',
     1                      ' JMAX should be at least ', JJ
            STOP
         ENDIF
C
         IF (NATJ .GT. KMAX) THEN
            WRITE (LOUT, *) ' Error...solution dimension too small,',
     1                      ' KMAX should be at least ',NATJ
            STOP
         ENDIF
C
         READ (LS, ERR=300, END=300) (X(J), J=1,JJ)
         READ (LS, ERR=300, END=300) ((S(N,J),N=1,NATJ),J=1,JJ)
         NS = NS + 1
         WRITE (LOUT, *) ' Read solution ',NS
         IF (NS.EQ.NSOL .AND. ISEN.LE.0) RETURN
C
      ELSEIF (ICHR .EQ. 'SENSITIVITY') THEN
         DO 40 I = 1, II
            READ (LS, ERR=400, END=400) IS,
     1                                  ((SN(N,J),N=1,NATJ),J=1,JJ)
            IF (ISEN .GT. 0) THEN
               DO 35 J = 1, JJ
                  CONT(J,I) = SN(ISEN,J)
   35          CONTINUE
            ENDIF
   40    CONTINUE
         WRITE (LOUT, *) ' Read sensitities for solution ',NS
         IF (NS .EQ. NSOL) RETURN
      ENDIF
      GO TO 10
C
  100 CONTINUE
      WRITE  (LOUT, 110) ICHR
  110 FORMAT (' Error reading ', A ,' record')
      STOP
C
  200 CONTINUE
      WRITE  (LOUT, 210) ICHR
  210 FORMAT (' Error reading ', A, ' data...')
      STOP
C
  300 CONTINUE
      WRITE (LOUT, *) ' Error reading solution...'
      WRITE (LOUT, *) ' NATJ=',NATJ,', JJ=',JJ,', P=',P,
     1                ', VFUEL=',VFUEL,', VOXID=',VOXID
      STOP
C
  400 CONTINUE
      WRITE (LOUT, *) ' Error reading gas-phase sensitivities...'
      WRITE (LOUT, *) ' ISEN=',ISEN,', NSOL=',NSOL,', NS=',NS,
     1               ', II=',II,', I=',I
      STOP
C
  600 FORMAT (' Error...',A,
     1        ' precision incompatible with post-processor')
  700 FORMAT (/,' Error...not enough ',A,' work space provided:',
     1        /,'                PROVIDED        REQUIRED',
     2        /,' INTEGER  ', 2I15,
     3        /,' REAL     ', 2I15,
     4        /,' CHARACTER', 2I15,/)
 701  FORMAT (/,' Error...not enough ',A,' work space provided:',
     1        /,'                PROVIDED        REQUIRED',
     2        /,' INTEGER  ', 2I15,
     3        /,' REAL     ', 2I15,/)
C
      END
C----------------------------------------------------------------------C

      SUBROUTINE PDRATE (KPROD, NPROD, XCM, ICKWRK, RCKWRK, P, JJ,
     1                   KK, II, IREAC, ILIST, XT, KMAX, JMAX, KSYM,
     2                   S, NY1, NT, Q, CIK, CONT, SUMRT, LTEXT, WT,
     3                   PREFIX, ILEN )
C
C*****single precision
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END single precision
C*****double precision
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END double precision

      DIMENSION 
     1 ICKWRK(*), RCKWRK(*), XCM(*), XT(*), S(KMAX,*), KPROD(*),
     2 Q(*), CIK(*), CONT(JMAX,*), SUMRT(*), ILIST(*), WT(*)

      CHARACTER IREAC(*)*(*), KSYM(*)*(*)
      CHARACTER*60 FILNAM
      CHARACTER*8  PREFIX

      DO 1000 N = 1, NPROD

        K = KPROD(N)

C       REACTION CONTRIBUTIONS TO SPECIES K

        DO 100 J = 1, JJ
          CALL CKYTX (S(NY1,J), ICKWRK, RCKWRK, XT)
          CALL CKQXP  (P, S(NT,J), XT, ICKWRK, RCKWRK, Q)
          CALL CKCONT (K, Q, ICKWRK, RCKWRK, CIK)
          CALL CKRHOX (P, S(NT,J), XT, ICKWRK, RCKWRK, RHO)

C         CONT(J,I) = CONTRIBUTION OF REACTION I TO SPECIES K AT PT J
C         SUMRT(J)  = TOTAL PRODUCTION RATE OF SPECIES K AT PT J
C         UNITS: 1/SEC
          SUMRT(J) = 0
          DO 100 I = 1, II
            CONT(J,I) = CIK(I) / RHO * WT(K)
            SUMRT(J) = SUMRT(J) + CONT(J,I)
100     CONTINUE

C       SELECT SIGNIFICANT REACTIONS

        CALL MINMAX (SUMRT, JJ, SMIN, SMAX, IPK)
        FRACT = 0.15
        SMIN = FRACT * ABS(SMIN)
        SMAX = FRACT * ABS(SMAX)
        NLIST = 0
        DO 200 I = 1, II
          CALL MINMAX (CONT(1,I), JJ, CMIN, CMAX, IPK)
          IF ( (ABS(CMIN) .GT. SMIN) .OR. (ABS(CMAX) .GT. SMAX) ) THEN
            NLIST = NLIST + 1
            ILIST(NLIST) = I
            CALL SQUEEZ ( LENGTH, IREAC(I) )
          ENDIF          
200     CONTINUE

C       WRITE SIGNIFICANT REACTIONS

        FILNAM = PREFIX(1:ILEN)//'.'//KSYM(K)
        CALL SQUEEZ (LENGTH, FILNAM)
        OPEN (LTEXT,  STATUS='UNKNOWN', FORM='FORMATTED', FILE=FILNAM)

        I1 = 1
        NPERL = 7
        I2 = MIN( I1 + NPERL - 1, NLIST )
        WRITE (LTEXT,'(10(1X,A))')
     1  'X(cm)', 'Rate(1/s)', ( IREAC(ILIST(I)), I = I1, I2 )
        DO 300 J = 1, JJ          
          WRITE (LTEXT,'(1X,F9.7, 9(1X,1PE12.5))')
     1    XCM(J), SUMRT(J), ( CONT(J,ILIST(I)), I = I1, I2 )
300     CONTINUE

        IF (I2 .LT. NLIST) THEN
          NPERL = 9
          I1 = I1 + NPERL
        ELSE
          RETURN
        ENDIF

500     CONTINUE

        I2 = MIN( I1 + NPERL - 1, NLIST )
        WRITE (LTEXT,'(1X,A)') 'MOL FRACTIONS:'
        WRITE (LTEXT,'(10(1X,A))')
     1  'X(cm)', ( IREAC(ILIST(I)), I = I1, I2 )
        DO 600 J = 1, JJ          
          WRITE (LTEXT,'(1X,F9.7, 9(1X,1PE12.5))')
     1    XCM(J), ( CONT(J,ILIST(I)), I = I1, I2 )
600     CONTINUE

        IF (I2 .LT. NLIST) THEN
          NPERL = 9
          I1 = I1 + NPERL
          GO TO 500
        ENDIF

        CLOSE (LTEXT)  

1000  CONTINUE

      RETURN
      END

C----------------------------------------------------------------------C

      SUBROUTINE RDKEY (LIN, LOUT, NSOL, KSYM, KK,
     1                  KSPEC, NSPEC, KPROD, NPROD, MOL, EINOX,
     2                  KFUEL)

C*****single precision
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END single precision
C*****double precision
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END double precision

      EXTERNAL CKCHUP

      DIMENSION KSPEC(1), KPROD(1), NRAY(10)

      CHARACTER KEY*4, LINE*76, KSYM(*)*(*), FUEL*16, CKCHUP*4

      LOGICAL MOL, EINOX, IERR, KERR

      NPROD = 0
      NSPEC = 0
      NDIM = KK
      FUEL = ' '
      KERR = .FALSE.
      IERR = .FALSE.

   10 CONTINUE

      KEY = ' '
      LINE = ' '
      WRITE (LOUT, '(/2X,A,/5X,A)') 'ENTER KEYWORD:',
     1  'MASS, MOLE, EINO, NSOL, SPEC, PROD, HELP, END'
      READ (LIN, '(A4,A76)', END=50) KEY, LINE
      KEY = CKCHUP(KEY, 4)

      IF (KEY .EQ. 'END') THEN
         IF (KERR) THEN
            WRITE (LOUT, *) ' Error in keyword input '
            STOP
         ENDIF
         GO TO 50

      ELSEIF (KEY .EQ. 'MOLE') THEN
         MOL = .TRUE.

      ELSEIF (KEY .EQ. 'MASS') THEN
         MOL = .FALSE.

      ELSEIF (KEY .EQ. 'NSOL') THEN
         CALL CKXNUM (LINE, 1, LOUT, NVAL, VALUE, IERR)
         KERR = KERR.OR.IERR
         NSOL = INT(VALUE)

      ELSEIF (KEY .EQ. 'SPEC') THEN
         CALL CKCRAY (LINE, KK, KSYM, LOUT, NDIM, KSPEC, NFD, IERR)
         IF (IERR) THEN
            WRITE (LOUT,'(A)')
     1      ' ERROR READING KEYWRD '//KEY
            KERR = .TRUE.
            IF (NFD .GT. KK) THEN
             WRITE (LOUT,'(A,I4,A)') 'FOUND TOO MANY SPECIES ?'
             NSPEC = KK
            ENDIF
         ELSE
           NSPEC = NFD
         ENDIF
         KERR = KERR.OR.IERR


      ELSEIF (KEY .EQ. 'PROD') THEN
         CALL CKCRAY (LINE, KK, KSYM, LOUT, NDIM, KPROD, NFD, IERR)
         IF (IERR) THEN
            WRITE (LOUT,'(A)')
     1      ' ERROR READING KEYWRD '//KEY
            KERR = .TRUE.
            IF (NFD .GT. KK) THEN
             WRITE (LOUT,'(A,I4,A)') 'FOUND TOO MANY SPECIES ?'
             NPROD = KK
            ENDIF
         ELSE
           NPROD = NFD
         ENDIF
         KERR = KERR.OR.IERR

      ELSEIF (KEY .EQ. 'EINO') THEN
         EINOX = .TRUE.
C        IF WANT INDEX, IDENTIFY FUEL
         READ (LIN, '(A4,A76)', END=50) KEY, LINE
         KEY = CKCHUP(KEY, 4)

         IF (KEY .EQ. 'FUEL') THEN
           NDIM = 10
           CALL CKCRAY (LINE, KK, KSYM, LOUT, NDIM, NRAY, NFD, IERR)
           IF (IERR .OR. NFD .NE. 1) THEN
             WRITE (LOUT,*) 'STOP, Could not find FUEL.'
             STOP
           ELSE
             FUEL = KSYM(KFUEL)
           ENDIF
         ELSE
             WRITE (LOUT,*) 'Must identify FUEL for Emission Index.'
             EINOX = .FALSE.
         ENDIF

      ELSEIF (KEY .EQ. 'HELP') THEN
         WRITE (LOUT, '(/2X,A,11(/2X,A))') 
     1 'Keywords available:',
     2 '-------------------',
     3 'MASS: want species mass fractions',
     4 'MOLE: want species mole fractions',
     5 'EINO: want emission indicies',
     6 '      EINO must be followed by keyword FUEL (next line)',
     7 'FUEL: fuel species name',
     8 'NSOL: solution index for multiple solutions',
     9 'SPEC character list : species to be printed',
     * 'PROD character list : production rates for species listed',
     1 'END: last keyword in input'

      ENDIF

      GO TO 10

   50 CONTINUE

C     WHICH FUEL TYPE FOR MIXTURE FRACTION?
      IF (EINOX) THEN
        CALL CKCOMP (FUEL, KSYM, KK, KFUEL)          
        WRITE (LOUT,*) 'FUEL IS ', FUEL
        WRITE (LOUT,*) 'SPECIES # ', KFUEL
      ENDIF

      RETURN
      END


C---------------------------------------------------------------------
C
      SUBROUTINE EMISS (V, EI, FMIX, JDIM, KK, KSYM, LOUT, JJ,
     1                  RHO, U, WT, WTM, X, XMF, KFUEL)
C
C*****single precision
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END single precision
C*****double precision
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END double precision
C
      DIMENSION XMF(JDIM,*), FMIX(*), NRAY(10), WT(*), WTM(*),
     1          EI(JDIM,*), X(*), RHO(*), V(*), U(*)
      CHARACTER KSYM(*)*(*), SPECIES*40
      LOGICAL KERR
      DATA SMALL/1.E-30/, ICALL/0/
C
      IF (ICALL .EQ. 0) THEN
        ICALL = 1
C
C      SEARCH MECHANISM FOR SPECIES IN RATIO
C
        SPECIES = 'NO NO2'
        NDIM = 10
        CALL CKCRAY (SPECIES, KK, KSYM, LOUT, NDIM, NRAY, NF, KERR)
        IF (KERR) THEN
          WRITE (LOUT,'(/,1X,A,/)') 'Computing Emission Index.'
          WRITE (LOUT,*) 'STOP, Could not find species list.'
          WRITE (LOUT,*) SPECIES
          WRITE (LOUT,*) 'NO. FOUND = ', NF
          WRITE (LOUT,*) 'INDICIES:', (NRAY(I),I=1,NF)
          STOP
        ELSEIF (NF .NE. 2) THEN
          WRITE (LOUT,'(/,1X,A,/)') 'Computing Emission Index.'
          WRITE (LOUT,*) 'STOP, wrong no of species in list.'
          WRITE (LOUT,*) SPECIES
          WRITE (LOUT,*) 'NO. FOUND = ', NF
          WRITE (LOUT,*) 'INDICIES:', (NRAY(I),I=1,NF)
          STOP
        ELSE
          KNO = NRAY(1)
          KNO2 = NRAY(2)
          WNO2 = WT(KNO2)
          WRITE (LOUT,'(/,1X,A)') 'Species in EINOx'
          DO 50 I=1,2
            WRITE (LOUT,'(1X,A,1X,I3)')
     1      KSYM(NRAY(I)), NRAY(I)
50        CONTINUE
        ENDIF

        SPECIES = 'CO'
        NDIM = 10
        CALL CKCRAY (SPECIES, KK, KSYM, LOUT, NDIM, NRAY, NF, KERR)
        IF (KERR .OR. NF .NE. 1) THEN
          WRITE (LOUT,*) 'WARNING: Could not find CO.'
        ELSE
          KCO = NRAY(1)
          WCO = WT(KCO)
          WRITE (LOUT,'(/,1X,A)') 'Species'
          WRITE (LOUT,'(1X,A,1X,I3)') KSYM(NRAY(1)), NRAY(1)
        ENDIF

        WFUEL = WT(KFUEL)

      ENDIF

C     INTEGRATE NO+NO2 OVER GRID USING MOLE FRACTIONS
C     NOTE: EI DEFINITION USES W_NO2 FOR NO

      EMFUEL = RHO(1) *XMF(1,KFUEL)*WFUEL/WTM(1) * U(1)

      XNOMX = 0.
      XNO2MX = 0.
      EINO2A = 0.
      EICOA  = 0.
      EINOXA = 0.
      EIFA   = 0.

      DO 200 J = 1, JJ-1
        DX   = X(J+1) - X(J)
        XCO  = XMF(J,KCO)
        XNO  = XMF(J,KNO)
        XNO2 = XMF(J,KNO2)
        XNOX = XNO + XNO2
        XNOMX  = MAX(XNOMX,  XNO)
        XNO2MX = MAX(XNO2MX, XNO2)
        XCOMX  = MAX(XCOMX,  XCO)

C       NOTE: V = v/r

C       NOx mass production rate / area
        EI(J,1) = 2.*RHO(J) *V(J) *(XNO+XNO2) *WNO2/WTM(J)
        EINOXA  = EINOXA + EI(J,1) *DX

C       NO2 mass production rate / area
        EI(J,2) = 2.*RHO(J) *V(J) *XNO2 *WNO2/WTM(J)
        EINO2A = EINO2A + EI(J,2) *DX

C       Fuel mass consumption rate / area
        YFUEL = XMF(J,KFUEL)*WFUEL/WTM(J)
        EI(J,3) = 2.*RHO(J) *V(J) *YFUEL
        EIFA   = EIFA  + EI(J,3) *DX

C       CO mass production rate / area
        EICOA  = EICOA + 2.*RHO(J) *V(J) *XCO *WCO/WTM(J) *DX
200   CONTINUE
C
C     Takeno's eqn 3 definition
      EINOT = EINOXA /EMFUEL
C
C     Divide integrals by fuel consumption rate
      EINOX = EINOXA / ( EMFUEL - EIFA )
      EINO2 = EINO2A / ( EMFUEL - EIFA )
      EICO  = EICOA  / ( EMFUEL - EIFA )
      PNO2  = EINO2 / (EINOX + SMALL) *100.
C
C     Fraction of fuel that exits without burning.
      EIF   = EIFA   /EMFUEL
C
      WRITE (LOUT,'(//2X,A,1PE12.5)') 'NOx Emission Index  =', EINOX
      WRITE (LOUT,'(2X,A,1PE12.5)')   'Percent of NO2      =', PNO2
      WRITE (LOUT,'(2X,A,1PE12.5)')   'CO Emission Index   =', EICO
      WRITE (LOUT,'(2X,A,1PE12.5)')   'Fraction unburnt F  =', EIF
      WRITE (LOUT,'(2X,A,1PE12.5)')   'Takeno NOx EI       =', EINOT
      WRITE (LOUT,'(/2X,A,1PE12.5)')  'NOx EI per area     =', EINOXA
C
      WRITE (LOUT,'(/2X,A,1PE12.5)')  'Max. NO Fracion     =', XNOMX
      WRITE (LOUT,'(2X,A,1PE12.5)')   'Max. NO2 Fracion    =', XNO2MX
      WRITE (LOUT,'(2X,A,1PE12.5/)')  'Max. CO Fracion     =', XCOMX
C
      RETURN
      END

C----------------------------------------------------------------------

      SUBROUTINE SQUEEZ (LENGTH, STRING)

C///////////////////////////////////////////////////////////////////////
C
C     T O O L S
C
C     SQUEEZ
C
C     SQUEEZE LEADING BLANKS AND MULTIPLE BLANKS FROM A CHARACTER
C     STRING.  RETURN THE LENGTH OF THE STRING (OR 1 IF ALL BLANK).
C
C     THIS VERSION MODIFIED BY A. LUTZ TO REPLACE BLANKS WITH "."
C
C///////////////////////////////////////////////////////////////////////

      IMPLICIT COMPLEX (A - P, R - Z), INTEGER (Q)
      CHARACTER CHAR*1, STRING*(*)
      INTEGER J, LENGTH
      LOGICAL BLANK

C///  SQUEEZE THE STRING.

      LENGTH = 0
      BLANK = .TRUE.
      DO 0100 J = 1, LEN (STRING)
         CHAR = STRING (J : J)
         IF (.NOT. BLANK .OR. CHAR .NE. ' ') THEN
            BLANK = CHAR .EQ. ' '
            LENGTH = LENGTH + 1
            STRING (LENGTH : LENGTH) = CHAR
         END IF
0100  CONTINUE

C///  ADJUST THE LENGTH AND PAD THE STRING.

      IF (0 .LT. LENGTH) THEN
         IF (STRING (LENGTH : LENGTH) .EQ. ' ') LENGTH = LENGTH - 1
         IF (LENGTH .LT. LEN (STRING)) STRING (LENGTH + 1 : ) = ' '
      ELSE
         LENGTH = 1
      END IF

C     REPLACE BLANKS WITH PERIODS (ANDY

      DO 0200 J = 1, LENGTH
         CHAR = STRING (J : J)
         IF (CHAR .EQ. ' ') THEN
            STRING (J : J) = "."
         END IF
0200  CONTINUE

C///  EXIT.

      RETURN
      END


C--------------------------------------------------------------------

      SUBROUTINE MIXTRF ( Y, KK, WT, FMIX )

C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C*****precision > double
        IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double

      DIMENSION Y(*), WT(*)
      COMMON /MIX/ AWT(4), NCF(4,200), NEC, NEH, NEO

C     This subroutine calculates the fuel mixture fraction for the
C     2-reactors with fuel defined by:
C     NEC, NEH, NEO = element pointers for C, H, O
C     Y(K)     = species mass fractions
C     FMIX     = mixture fraction of fuel
C     AWT(M)   = element weights
C     NCF(M,K) = coefficient of element-m in species-k
C     WT(K)    = species molecular mass

      FMIX = 0.
      DO 10 K = 1, KK
          TERM = NCF(NEH,K)*AWT(NEH) / WT(K)
          FMIX = FMIX + Y(K)*TERM
10    CONTINUE
      IF (NEC .GT. 0) THEN
        DO 20 K = 1, KK
          TERM = NCF(NEC,K)*AWT(NEC) /WT(K)
          FMIX = FMIX + Y(K)*TERM
20      CONTINUE
      ENDIF
      END

