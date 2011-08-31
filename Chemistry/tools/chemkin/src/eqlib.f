C  CVS $Revision: 1.1.1.1 $  created $Date: 2006/05/26 19:09:33 $
      SUBROUTINE EQABS
C
C  START PROLOGUE
C
C     CHEMLIB IMPLEMENTATION OF STANJAN-III EQUILIBRIUM PROGRAM
C
C     WRITTEN BY:
C         ROBERT J. KEE, ANDREW E. LUTZ, AND FRAN M. RUPLEY
C         COMPUTATIONAL MECHANICS DIVISION
C         SANDIA NATIONAL LABORATORIES
C         LIVERMORE, CA  94550
C         (510) 294-2761        aelutz@ca.sandia.gov
C
C     CONTRIBUTION TO REAL-GAS VERSION:
C         MARC D. RUMMINGER
C         DEPARTMENT OF MECHANICAL ENGINEERING
C         UC BERKELEY
C
C  V.3.12 97/4/17 A. Lutz
C  1. fix bug on constant-V cases by resetting densities to zero 
C     (DO 113).
C  2. put SAVE statement in EQUIL
C  V.3.11 97/12/17 A. Lutz
C  1. move initialization loop in EQUIL to allow continuation
C  V.3.10 97/07/17
C  1. add KMON = 2 if DIAG keyword used
C  V.3.8 97/03/01
C  1. remove "Surface" change blocks and use logical LSURF instead,
C     and make new main "driver" program to set up arrays,
C     to satisfy Reaction Design requirement of providing object
C     files instead of source code.
C  V.3.7, 96/08/07 A. Lutz
C  1. added arrays to call lists for CJ option
C  V.3.6 (4/21/96 F. Rupley)
C  1.  initial CHEMKIN-III version
C  2.  added SUBROUTINE EQABS to hold documentation
C  3.  appended the STANJAN subroutines to the eqlib.f subroutines
C  4.  change blocks to OPEN ascii linkfiles
C  5.  modify error handling for keyword input
C  6.  add START/END PROLOGUE language to subroutine comments
C  7.  SUBROUTINE EQSTRT logical QUIT initialized .FALSE.
C  8.  need SAVE in SUBROUTINEs SJTP and SJEQLB
C
C  END PROLOGUE
C
      END
C
      SUBROUTINE EQINTP (LINCK, LINSK, LIN, LOUT, LSAVE, LSURF,
     1                   LIWORK, I, LRWORK, R, LCWORK, C)
C
C            INTERACTIVE DRIVER FOR STANJAN-III EQUILIBRIUM PROGRAM
C
C     WRITTEN BY:
C         ROBERT J. KEE, ANDREW E. LUTZ, and FRAN M. RUPLEY
C         COMPUTATIONAL MECHANICS DIVISION
C         SANDIA NATIONAL LABORATORIES
C         LIVERMORE, CA  94550
C         (415) 294-3272
C
C*****precision > double
        IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C        IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      PARAMETER (LENIRG=1, LENRRG=1, LENCRG=1)
      DIMENSION I(LIWORK), R(LRWORK), IRGWRK(LENIRG), RGWRK(LENRRG)
      CHARACTER*16 C(LCWORK), CRGWRK(LENCRG)
      LOGICAL LCNTUE, EQST, LPRNT, REVERS, FAZEQ, KERR
      DATA LPRNT, LCNTUE, EQST, REVERS, FAZEQ /.TRUE.,4*.FALSE./
C
      WRITE  (LOUT, 15)
   15 FORMAT (/' EQUIL:  Chemkin interface for Stanjan-III',
     1        /'         CHEMKIN-III Version 3.12, 97/12/23',
C*****precision > double
     2        /'         DOUBLE PRECISION')
C*****END precision > double
C*****precision > single
C     2        /'         SINGLE PRECISION')
C*****END precision > single
C
      KERR = .FALSE.
      CALL CKLEN (LINCK, LOUT, LENICK, LENRCK, LENCCK, IFLAG)
      IF (IFLAG.NE.0) KERR = .TRUE.
C
      IICK = 1
      IRCK = 1
      ICCK = 1
      LENISK = 0
      LENRSK = 0
      LENCSK = 0
      IF (LSURF .GT. 0) THEN
         CALL SKLEN (LINSK, LOUT, LENISK, LENRSK, LENCSK, IFLAG)
         IF (IFLAG.NE.0) KERR = .TRUE.
      ENDIF
C
      IISK = IICK + LENICK
      IRSK = IRCK + LENRCK
      ICSK = ICCK + LENCCK
      IITOT = LENICK + LENISK - 1
      IRTOT = LENRCK + LENRSK - 1
      ICTOT = LENCCK + LENCSK - 1
      IF (LIWORK .LT. IITOT) THEN
         WRITE (LOUT, *) 'EQINTP ERROR: IWORK must be at least ', IITOT
         KERR = .TRUE.
      ENDIF
      IF (LRWORK .LT. IRTOT) THEN
         WRITE (LOUT, *) 'EQINTP ERROR: RWORK must be at least ', IRTOT
         KERR = .TRUE.
      ENDIF
      IF (LCWORK .LT. ICTOT) THEN
         WRITE (LOUT, *) 'EQINTP ERROR: CWORK must be at least ', ICTOT
         KERR = .TRUE.
      ENDIF
      IF (KERR) RETURN
C
      CALL CKINIT (LENICK, LENRCK, LENCCK, LINCK, LOUT, I, R, C, IFLAG)
      KERR = KERR .OR. (IFLAG.NE.0)
      IF (LSURF .GT. 0)
     1   CALL SKINIT (LENISK, LENRSK, LENCSK, LINSK, LOUT,
     2                I(IISK), R(IRSK), C(ICSK), IFLAG)
      KERR = KERR .OR. (IFLAG.NE.0)
      IF (KERR) RETURN
C
      CALL CKINDX (I(IICK), R(IRCK), MM, KKGAS, II, NFIT)
      IF (LSURF .GT. 0) THEN
         CALL SKINDX (I(IISK), MM, KKGAS, KKSURF, KKBULK, KKTOT, 
     1                NPHASE, NNSURF, NFSURF, NLSURF, NNBULK, NFBULK, 
     2                NLBULK, IISUR)
         IKCON = IISK + LENISK
         IXCON = IRSK + LENRSK
         IATOM = ICSK + LENCSK
      ELSE
         KKSURF = 0
         KKBULK = 0
         KKTOT = KKGAS
         NPHASE = 1
         NNSURF = 0
         NNBULK = 0
         IISUR  = 0
         IKCON = IICK + LENICK
         IXCON = IRCK + LENRCK
         IATOM = ICCK + LENCCK
      ENDIF
      CALL EQLEN (MM, KKGAS, KKBULK, KKTOT, NPHASE, NNBULK, 
     1            LENIEQ, LENREQ, LENCEQ)
C
      IKKPHA= IKCON  + KKTOT
      IKFIRS= IKKPHA + NPHASE
      IKLAST= IKFIRS + NPHASE
      IIEQ  = IKLAST + NPHASE
      IITOT = IIEQ   + LENIEQ - 1
      IF (LIWORK .LT. IITOT) THEN
         WRITE (LOUT, *) 'ERROR: IWORK must be at least ', IITOT
         KERR = .TRUE.
      ENDIF
C
      IREAC = IXCON + KKTOT
      IACT  = IREAC + KKTOT
      IPDEN = IACT  + KKTOT
      IVINP = IPDEN + NPHASE
      IREQ  = IVINP + 10
      IRTOT = IREQ  + LENREQ - 1
      IF (LRWORK .LT. IRTOT) THEN
         WRITE (LOUT, *) 'ERROR: RWORK must be at least ', IRTOT
         KERR = .TRUE.
      ENDIF
C
      IKSYM = IATOM + MM
      IPSYM = IKSYM + KKTOT
      ICEQ  = IPSYM + NPHASE
      ICTOT = ICEQ  + LENCEQ - 1
      IF (LCWORK .LT. ICTOT) THEN
         WRITE (LOUT, *) 'ERROR: CWORK must be at least ', ICTOT
         KERR = .TRUE.
      ENDIF
      IF (KERR) RETURN
C
  100 CONTINUE
C
C     INITIALIZE INPUT FOR PROBLEM
      CALL EQSTRT (LIN, LOUT, LENICK, LENRCK, LENCCK, LINCK, LENISK,
     1             LENRSK, LENCSK, LINSK, MM, KKTOT, NPHASE, I(IICK),
     2             R(IRCK), C(ICCK), I(IISK), R(IRSK), C(ICSK),
     3             I(IKKPHA), I(IKFIRS), I(IKLAST), NPHASE, NNBULK,
     4             NFBULK, NLBULK, R(IACT), MM, KKGAS, KKSURF, KKBULK,
     5             KKTOT, C(IATOM), C(IKSYM), C(IPSYM), NOP, KMON, KFRZ,
     7             R(IREAC), TEST, PEST, LCNTUE, NCON, I(IKCON), 
     8             R(IXCON), R(IVINP), REVERS, FAZEQ, LENIRG, LENRRG, 
     9             LENCRG, IRGWRK, RGWRK, CRGWRK, LINKRG, LSURF, KERR)
      IF (KERR) RETURN
C
C     INTERFACE TO STANJAN
      CALL EQUIL (LOUT, LPRNT, LSAVE, EQST, LCNTUE, I(IICK), R(IRCK),
     1            I(IISK), R(IRSK), I(IKKPHA), I(IKFIRS), I(IKLAST),
     2            NPHASE, NNBULK, NFBULK, NLBULK, R(IACT), R(IPDEN),
     3            LENIEQ, I(IIEQ), LENREQ, R(IREQ), LENCEQ, C(ICEQ), 
     4            MM, KKGAS, KKSURF, KKBULK, KKTOT, C(IATOM), C(IKSYM),
     5            NOP, KMON, KFRZ, R(IREAC), TEST, PEST, NCON, I(IKCON), 
     6            R(IXCON), R(IVINP), REVERS, FAZEQ, IRGWRK, RGWRK,
     7            LSURF, KERR)
C
      IF (LCNTUE) GO TO 100
      RETURN
      END
C-------------------------------------------------------------------
C
      SUBROUTINE EQLEN (MM, KKGAS, KKBULK, KKTOT,
     1                 NNPHAS, NNBULK, NITOT, NRTOT, NCTOT)
C
C  START PROLOGUE
C
C     SET POINTERS AND COMPUTE NECESSARY WORK SPACE
C     INPUT:
C     MM    - NUMBER OF ELEMENTS
C     KKGAS - NUMBER GAS SPECIES
C     KKBULK- NUMBER BULK SPECIES
C     KKTOT - TOTAL SPECIES
C     NNPHAS- MAXIMUM NUMBER OF PHASES (INCLUDING SURFACE PHASES)
C     NNBULK- MAXIMUM NUMBER OF BULK PHASES (CONDENSED)
C     OUTPUT:
C     NITOT, NRTOT, NCTOT - INTEGER, REAL, CHARACTER WORKSPACE REQUIRED
C
C  END PROLOGUE
C
C*****precision > double
        IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C        IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      COMMON /IPAR/ NIWORK, NICMP, NIKNT, NRWORK, NRADD, NSML,
     2              NHML, NWT, NDEN, NXCON, NKCON, NAMAX,
     3              NPHASE, NX1,   NX2,   NY1,    NY2,   NT1,   NT2,
     4              NP1,    NP2,   NV1,   NV2,    NWM1,  NWM2,
     5              NS1,    NS2,   NU1,   NU2,    NH1,   NH2,
     6              NC1,    NC2,   NCDET, NTEST,  NPEST, NSMOL1,
     7              NSMOL2, NYG1,  NXP1,  NXP2,   NWMG1, NWMG2,
     8              NVG1,   NVG2, KKSURF, LSURF
C
C    If calling SJROP:
C
      NAMAX = KKTOT
      NPHASE = 1 + NNBULK
      NPMAX = NPHASE
      NSMAX = KKGAS + KKBULK
      NWORK = MAX(2*NAMAX, NAMAX+NPMAX)
C
C    If using a species data file or calling SJSET:
      NIWORK = 22 + 14*NAMAX + 4*NPMAX + 8*NSMAX + 2*NAMAX*NSMAX
      NRWORK = 24 + 16*NAMAX + 12*NAMAX*NAMAX + 3*NPMAX*NAMAX + 6*NPMAX
     1            + 18*NSMAX + NWORK*NWORK + NWORK
C
      NCTOT = NSMAX
C
C     pointers for additional integer arrays
C
      NICMP = NIWORK + 1
      NIKNT = NICMP + NAMAX * KKTOT
      NKCON = NIKNT + KKTOT
      NITOT = NKCON + KKTOT - 1
C
C     pointers for additional real arrays
C
      NSML  = NRWORK + 1
      NHML  = NSML  + KKTOT
      NWT   = NHML  + KKTOT
      NDEN  = NWT   + KKTOT
      NXCON = NDEN  + KKTOT
      NX1   = NXCON + KKTOT
      NX2   = NX1   + KKTOT
      NY1   = NX2   + KKTOT
      NY2   = NY1   + KKTOT
      NYG1  = NY2   + KKTOT
      NXP1  = NYG1  + KKTOT
      NXP2  = NXP1  + KKTOT
      NSMOL1= NXP2  + KKTOT
      NSMOL2= NSMOL1+ NNPHAS
      NT1   = NSMOL2+ NNPHAS
      NT2   = NT1   + 1
      NP1   = NT2   + 1
      NP2   = NP1   + 1
      NS1   = NP2   + 1
      NS2   = NS1   + 1
      NU1   = NS2   + 1
      NU2   = NU1   + 1
      NV1   = NU2   + 1
      NV2   = NV1   + 1
      NVG1  = NV2   + 1
      NVG2  = NVG1  + 1
      NWM1  = NVG2  + 1
      NWM2  = NWM1  + 1
      NWMG1 = NWM2  + 1
      NWMG2 = NWMG1 + 1
      NH1   = NWMG2 + 1
      NH2   = NH1   + 1
      NC1   = NH2   + 1
      NC2   = NC1   + 1
      NCDET = NC2   + 1
      NTEST = NCDET + 1
      NPEST = NTEST + 1
C
      NRTOT = NPEST
C     NRADD is added space for Stanjan's SW
      NRADD = NRTOT - NRWORK
C
C     end of SUBROUTINE EQLEN
      RETURN
      END
C
C-------------------------------------------------------------------
C
      SUBROUTINE EQPRNT (KKGAS, KKBULK, KKTOT, NNBULK,
     1                   NFBULK, NLBULK, KFIRST, KLAST, NOP, KFRZ,
     2                   LOUT, MM, LSAVE, LCNTUE, KSYM, REQWRK)
C
C  START PROLOGUE
C  END PROLOGUE
C
C*****precision > double
        IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C        IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      DIMENSION REQWRK(*), KFIRST(*), KLAST(*)
      CHARACTER KSYM(*)*(*)
      LOGICAL LCNTUE
      COMMON /IPAR/ NIWORK, NICMP, NIKNT, NRWORK, NRADD, NSML,
     2              NHML, NWT, NDEN, NXCON, NKCON, NAMAX,
     3              NPHASE, NX1,   NX2,   NY1,    NY2,   NT1,   NT2,
     4              NP1,    NP2,   NV1,   NV2,    NWM1,  NWM2,
     5              NS1,    NS2,   NU1,   NU2,    NH1,   NH2,
     6              NC1,    NC2,   NCDET, NTEST,  NPEST, NSMOL1,
     7              NSMOL2, NYG1,  NXP1,  NXP2,   NWMG1, NWMG2,
     8              NVG1,   NVG2, KKSURF, LSURF

      DATA XMIN/1.E-12/
C
      IF (NNBULK .GT. 0) THEN
         KFBULK = KFIRST(NFBULK)
         KLBULK = KLAST (NLBULK)
      ENDIF
C
      IF (KFRZ .NE. 0) THEN
        WRITE (LOUT, '(/1X,A/)')
     1 ' CHEMICAL COMPOSITION IS FROZEN.'
        WRITE (LOUT, 7201)
      ELSE
        WRITE (LOUT, 7200)
      ENDIF
C
      WRITE (LOUT,7203)
     1   REQWRK(NP1), REQWRK(NP2), REQWRK(NT1), REQWRK(NT2),
     2   REQWRK(NV1), REQWRK(NV2), REQWRK(NH1), REQWRK(NH2),
     3   REQWRK(NU1), REQWRK(NU2), REQWRK(NS1), REQWRK(NS2),
     4   REQWRK(NWM1), REQWRK(NWM2)
         WRITE (LOUT, 7210)
         DO 50 K = 1, KKGAS
           IF( MAX(REQWRK(NX1+K-1),REQWRK(NX2+K-1)) .GT. XMIN )
     1     WRITE (LOUT, 7215) KSYM(K), REQWRK(NX1+K-1), REQWRK(NX2+K-1)
50       CONTINUE
        IF (LSURF .GT. 0) THEN
        IF (KKBULK .GT. 0) THEN
         WRITE (LOUT, 7215) (KSYM(K), REQWRK(NX1+K-1), REQWRK(NX2+K-1),
     2                          K = KFIRST(NFBULK), KLAST(NLBULK))
        ENDIF
        WRITE (LOUT, '(/,2X,A)') ' GAS PHASE'
        WRITE (LOUT, 7220) REQWRK(NSMOL1), REQWRK(NSMOL2)
        WRITE (LOUT, 7230) REQWRK(NWMG1), REQWRK(NWMG2),
     1                    REQWRK(NVG1),  REQWRK(NVG2)
        WRITE (LOUT, 7210)
        WRITE (LOUT, 7215) ( KSYM(K), REQWRK(NXP1+K-1),
     1                      REQWRK(NXP2+K-1), K = 1, KKGAS )
        IF (NNBULK .GT. 0) THEN
          DO 100 N = NFBULK, NLBULK
            WRITE (LOUT, '(/,2X,A,I3)') ' BULK PHASE #',N-1
            WRITE (LOUT, 7220) REQWRK(NSMOL1+N-1), REQWRK(NSMOL2+N-1)
            WRITE (LOUT, 7210)
            WRITE (LOUT, 7215) (KSYM(K), REQWRK(NXP1+K-1),
     1                   REQWRK(NXP2+K-1), K=KFIRST(N), KLAST(N))
  100     CONTINUE
        ENDIF
        ENDIF
C
      IF (NOP .EQ. 10)
     1   WRITE (LOUT, 7100) REQWRK(NC1), REQWRK(NC2), REQWRK(NCDET),
     2                      REQWRK(NCDET) / REQWRK(NC1)
C
      IF (LSAVE .GT. 0) THEN
C
C     WRITE BINARY RECORD OF STARTING AND EQUILIBRIUM VALUES
C
      IF (LSURF .LE. 0) THEN
      WRITE (LSAVE) REQWRK(NT1), REQWRK(NP1), REQWRK(NU1), REQWRK(NH1),
     1              REQWRK(NS1), REQWRK(NV1), REQWRK(NC1),
     2              (REQWRK(NX1 + K - 1), K = 1, KKGAS),
     4              REQWRK(NT2), REQWRK(NP2), REQWRK(NU2), REQWRK(NH2),
     5              REQWRK(NS2), REQWRK(NV2), REQWRK(NC2),
     6              (REQWRK(NX2 + K - 1), K = 1, KKGAS)
      ELSE
      WRITE (LSAVE) REQWRK(NT1), REQWRK(NP1), REQWRK(NU1), REQWRK(NH1),
     1              REQWRK(NS1), REQWRK(NV1), REQWRK(NC1),
     2              (REQWRK(NX1 + K - 1), K = 1, KKGAS),
     3              (REQWRK(NX1 + K - 1), K = KFBULK, KLBULK),
     4              REQWRK(NT2), REQWRK(NP2), REQWRK(NU2), REQWRK(NH2),
     5              REQWRK(NS2), REQWRK(NV2), REQWRK(NC2),
     6              (REQWRK(NX2 + K - 1), K = 1, KKGAS),
     7              (REQWRK(NX2 + K - 1), K = KFBULK, KLBULK)
      ENDIF
      ENDIF
C
C
      IF (LCNTUE) THEN
         WRITE (LOUT, *)
         WRITE (LOUT, *) '*************************',
     1                   'CONTINUING TO NEW PROBLEM',
     2                   '*************************'
         WRITE (LOUT, *)
      ENDIF
C
C      FORMATS
C
 7200 FORMAT (/'  MIXTURE:     ',T20,
     *         'INITIAL STATE:',T40,'EQUILIBRIUM STATE:')
 7201 FORMAT (/'  MIXTURE:     ',T20,
     *         'INITIAL STATE:',T40,'FROZEN STATE:')
 7203 FORMAT (/' P (atm)       ',T20,1PE14.4,T40,1PE14.4,
     2        /' T (K)         ',T20,1PE14.4,T40,1PE14.4,
     3        /' V (cm3/gm)    ',T20,1PE14.4,T40,1PE14.4,
     4        /' H (erg/gm)    ',T20,1PE14.4,T40,1PE14.4,
     5        /' U (erg/gm)    ',T20,1PE14.4,T40,1PE14.4,
     6        /' S (erg/gm-K)  ',T20,1PE14.4,T40,1PE14.4,
     7        /' W (gm/mol)    ',T20,1PE14.4,T40,1PE14.4)
 7210 FORMAT(  ' Mol Fractions')
 7215 FORMAT(  ' ', A16,            T20,1PE14.4,T40,1PE14.4)
 7220 FORMAT(  ' Mols          ',T20,1PE14.4,T40,1PE14.4)
 7230 FORMAT(  ' W (gm/mol)    ',T20,1PE14.4,T40,1PE14.4,
     1        /' V (cm3/gm)    ',T20,1PE14.4,T40,1PE14.4)
 7100 FORMAT (/' C-J DETONATION PROPERTIES',
     1        /' C (cm/s)      ',T20,1PE14.4,T40,1PE14.4,
     2        /' CDET (cm/s)   ',T40, 1PE14.4,
     3        /' MACH          ',T40, 1PE14.4 )
C
 7216 FORMAT(  ' ', A16,            T20,1PE14.4)
 7221 FORMAT(  ' Mols          ',T20,1PE14.4)
 7231 FORMAT(  ' W (gm/mol)    ',T20,1PE14.4,
     1        /' V (cm3/gm)    ',T20,1PE14.4)
C
C     end of SUBROUTINE EQPRNT
      RETURN
      END
C
C------------------------------------------------------------------
C
      SUBROUTINE EQRUN (KKGAS, NNBULK, NFBULK, NLBULK, KKPHAS, KFIRST,
     1                  KLAST, NOP, KMON, LOUT, MM, NCON, ATOM, CEQWRK,
     3                  IEQWRK, REQWRK, ICKWRK, RCKWRK, ISKWRK,
     4                  RSKWRK, NSMAX, KERR, HUGE, KFRZ, FAZEQ,
     5                  IRGWRK, RRGWRK)
C
C  START PROLOGUE
C
C    THIS ROUTINE DRIVES THE CHEMKIN INTERFACE TO THE STANJAN-3
C    EQUILIBRIUM CODE.  IT CAN ONLY BE CALLED AFTER SUBROUTINE
C    EQINIT IS CALLED TO SET UP THE INTERNAL WORKING AND STORAGE SPACE.
C    ALSO, THE INCLUDED ROUTINES SJTP,SJTPRP **MUST** BE USED INSTEAD OF
C    THE ORIGINAL ROUTINES PROVIDED IN THE STANJAN LIBRARY.
C
C  INPUT
C    NOP   - INTERGER DESIGNATING THE PROBLEM TYPE:
C            1.  SPECIFIED T AND P
C            2.  SPECIFIED T AND V
C            3.  SPECIFIED T AND S
C            4.  SPECIFIED P AND V
C            5.  SPECIFIED P AND H
C            6.  SPECIFIED P AND S
C            7.  SPECIFIED V AND U
C            8.  SPECIFIED V AND H
C            9.  SPECIFIED V AND S
C            10. CHAPMAN-JOUGUET DETONATION (H, S, V, T CONTAIN THE
C                UNBURNED STATE, AND TE IS THE BURNED ESTIMATE.
C    KMON  - INTEGER COTROLLING PRINTED OUTPUT FROM STANJAN.  FOR
C            NO RUN-TIME INFORMATION USE KMON=0.
C    LOUT  - OUTPUT UNIT ONTO WHICH STANJAN DIAGNOSTICS ARE WRITTEN.
C    NPHASE- NUMBER OF PHASES IN THE SYSTEM.
C    KFRZ  - 0 EQUILIBRIUM COMPOSITION, 1 FROZEN COMPOSITION
C    REAC  - VECTOR OF REACTANT MOLE FRACTIONS, FROM WHICH THE ATOM
C            POPULATIONS ARE COMPUTED.  DIMENSION REAC(*) AT LEAST KK,
C            WHERE KK IS THE NUMBER OF SPECIES.
C    DCS   - VECTOR OF THE CONDENSED PHASE SPECIES DENSITIES.
C            DCS(*) MUST BE DIMENSIONED AT LEAST KK IN THE CALLING
C            PROGRAM, BUT ALL ENTRIES CAN BE 0.
C    ATOM  - ARRAY OF CHARACTER*16 NAMES OF THE ELEMENTS IN THE PROBLEM.
C            DECLARE AND DIMENSION CHAR*16 ATOM(MM), WHERE MM IS THE
C             NUMBER OF ELEMENTS.
C    CEQWRK- ARRAY OF CHARACTER*16 NAMES OF THE SPECIES IN THE PROBLEM.
C            DECLARE AND DIMENSION CHAR*16 CEQWRK(KK).
C    IEQWRK- INTEGER EQUILIBRIUM WORK SPACE.
C            IEQWRK(*) MUST BE DIMENSIONED AT LEAST LENIEQ.
C    REQWRK- REAL EQUILIBRIUM WORK SPACE.
C            REQWRK(*) MUST BE DIMENSIONED AT LEAST LENREQ.
C    LENIEQ- ACTUAL DIMENSION DECLARED FOR IEQWRK(*)
C            LENIEQ MUST BE AT LEAST
C            22 + 14*MM + 4*NPHASE + 8*KK + 2*MM*KK
C            WHERE MM IS THE NUMBER OF ELEMENTS AND KK IS THE NUMBER OF
C            SPECIES.
C    LENREQ- ACTUAL DIMENSION DECLARED FOR REQWRK(*)
C            LENREQ MUST BE AT LEAST
C            24 + 16*MM + 12*MM*MM + 3*NPHASE*MM + 6*NPHASE
C                 + 18*KK + NWORK*NWORK + NWORK,  WHERE NWORK IS
C            MAX(2*MM, MM+NPHASE).
C    ICKWRK- INTEGER CHEMKIN WORK SPACE
C    RCKWRK - REAL CHEMKIN WORK SPACE.
C    LENRCK- ACTUAL DECLARED DIMENSION OF THE REAL CHEMKIN WORK SPACE.
C
C  INPUT AND OUTPUT
C    (ONLY THE PROPERTIES REQUIRED FOR THE GIVEN OPTION NEED TO BE GIVEN
C    P     - PRESSURE (DYNES/CM**2)
C    T     - TEMPERATURE (K)
C    H     - MIXTURE ENTHALPY (J/KG)
C    S     - MXTURE ENTROPY (J/KG-K)
C    U     - MIXTURE INTERNAL ENERGY (J/KG)
C    V     - MIXTURE SPECIFIC VOLUME (M**3/KG)
C    TE    - ESTIMATED BURNED TEMPERTURE FOR THE C-J DETONATION OPTION.
C    TP    - TP(N) IS THE TEMPERATURE OF THE N-TH PHASE.
C
C  OUTPUT
C    C     - SOUND SPEED (M/S) IF THE OTION CALLS FOR COMPUTING IT
C    CDET  - C-J DETONATION WAVESPEED (M/S)
C    WM    - EQUILIBRIUM MEAN MOLECULAR WEIGHT (G/MOLE)
C    PMOL  - ARRAY OF RELATIVE MOLES OF THE N-TH PHASE.
C    SMOL  - ARRAY OF RELATIVE MOLES OF THE K-TH SPECIES.
C    XP    - PHASE MOLE FRACTION OF THE K-TH SPECIES.
C    XM    - MIXTURE MOLE FRACTION OF THE K-TH SPECIES.
C            DIMENSION XM(*) AT LEAST KK.
C    Y     - MIXTURE MASS FRACTION.  DIMENSION Y(*) AT LEAST KK.
C
C  END PROLOGUE
C
C*****precision > double
        IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C        IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      DIMENSION IEQWRK(*), REQWRK(*), ICKWRK(*), RCKWRK(*),
     1          ISKWRK(*), RSKWRK(*), KKPHAS(*), KFIRST(*),
     2          KLAST(*), IRGWRK(*), RRGWRK(*)
      CHARACTER ATOM(*)*16, CEQWRK(*)*16
      LOGICAL KERR, FAZEQ
C
C 80 pointers required by SJEQLB
       COMMON /SJEPTR/
     ;   IoKERR,IoKMON,IoKTRE,IoKUMO,IoNA,IoNB,IoNP,IoNS,IoIB,IoIBO,
     ;   IoJB,IoJBAL,IoJBA,IoJBB,IoJBO,IoJBX,IoJS2,IoKB,IoKB2,IoKBA,
     ;   IoKBB,IoKBO,IoKPC,IoKPCX,IoLB2,IoMPA,IoMPJ,IoN,LoN,IoNSP,
     ;   IoFRND,IoHUGE,IoR1,IoR2,IoR3,IoA,LoA,IoB,LoB,IoBBAL,
     ;   IoCM,LoCM,IoD,LoD,IoDC,LoDC,IoDPML,IoDLAM,IoDLAY,IoE,
     ;   LoE,IoEEQN,IoELAM,IoELMA,IoELMB,IoF,IoG,IoHA,IoHC,IoPA,
     ;   IoPC,IoPMOL,IoQ,LoQ,IoQC,LoQC,IoRC,LoRC,IoRL,IoRP,
     ;   IoSMOA,IoSMOB,IoSMOO,IoSMOL,IoSMUL,IoW,IoX,IoXO,IoY,IoZ
C
C  20 additional pointers required by  SJTP
       COMMON /SJTPTR/
     ;   IoKFRZ,IoCVCJ,IoPATM,IoRGAS,IoP,IoT,IoH,IoS,IoU,IoV,
     ;   IoWM,IoTP,IoDCS,IoDHF0,IoWMS,IoHMH0,IoS0,IoWMP,IoXM,IoYM
C
C  10 additional pointers required by SJSET
       COMMON /SJSPTR/
     ;   IoKNFS,IoKSME,IoKUFL,IoNAF,IoNSF,IoJFS,IoNF,LoNF,IoITHL,IoITHM
C
C  20 additional pointers required by SJROP
       COMMON /SJRPTR/
     ; IoKDET,IoKRCT,IoKSND,IoKTRI,IoKTRT,IoI2,IoI3,IoNCON,IoTE,IoC,
     ; IoCDET,IoH1,IoP1,IoT1,IoV1,IoR5,IoR6,IoSMFA,IoSMFB,IoSMOZ
C
C--------------------------------------------------------------------
C      Species are loaded by phases:
C	   phase 1:  species  1, 2, ... , NSP(1)
C	   phase 2:  species  NSP(1)+1, ... , NSP(2)
C	   ...
C	   phase NP: species  NSP(NP-1)+1, ... , NSP(NP)
C
C	   CHEM(J): CHARACTER*16 name of species J, for J=1, ... ,NS
C
C      Atoms appearing in the species are numbered sequentially:
C	   ATOM(I): CHARACTER*16 name of the atom, for I=1, ..., NA
C
C      N(I,J) is the number of Ith atoms in a Jth molecule, which may
C      be negative for electrons in positive ions.
C
C      Species may be cross-referenced to a data file for property
C      look-up: JFS(J) is the file reference of species
C      (JFS le NSMAX).
C------------------------------------------------------------------
C    If calling SJROP:
C
C    Load the following into the work arrays:
C
C	   Mixture specification parameters:
C	       NA	   number of atom types in the system
C	       NP	   number of phases in the system
C	       NS	   number of species in the system
C	       NSP(M)	   number of species in Mth phase
C	       N(I,J)	   number of Ith atoms in molecule of Jth species
C	       PA(I)	   (relative) mols of Ith atoms
C
C	   Pointers to species data in the file:
C              JFS(J)	   index of Jth species in the species data file
C
C	   Species properties required for solution:
C	       DCS(J)	   density of Jth species, KG/M**3 (0 for gas)
C	       WMS(J)	   molal mass of Jth species, KG/KG-MOL
C
C	   Properties specified (load two specified properties):
C	       P	   pressure, Pa    (trial value if not specified)
C	       T	   temperature, K  (trial value if not specified)
C	       H	   mixture enthalpy, J/kg
C	       S	   mixture entropy, J/kg-K
C	       U	   mixture internal energy, J/kg
C	       V	   mixture specific volume, m**3/kg
C
C	   Other property specifications for special options:
C	       TE	   estimate of burned gas T (C-J det. run only)
C	       TP(K)	   T of Kth reactant phase (reactants run only)
C
C	   Control parameters:
C	       KFRZ 0 equilibrium composition
C                  -1 frozen composition with same T all phases
C                  -2 frozen composition, different phase temperatures
C	       KMON	   runtime monitor control (0 none)
C	       KUMO	   output unit for runtime monitor
C
C	   Select the option:
C	       NOP	   option selection:
C			 1 specified T and P
C			 2 specified T and V
C			 3 specified T and S
C			 4 specified P and V
C			 5 specified P and H
C			 6 specified P and S
C			 7 specified V and U
C			 8 specified V and H
C			 9 specified V and S
C			10 Chapman-Jouguet detonation
C		       (H,S,V,T contain unburned state,
C                        TE is burned estimate)
C
C       CALL SJROP (NAMAX,NSMAX,NIWORK,NRWORK,NSWORK,
C     ;		   ATOM,CEQWRK,IWORK,RWORK,SWORK,NOP)
C
C      Output in the work arrays:
C
C	   Calculated temperature-dependent properties of the species:
C	       G(J)	   g(T,P)/RT for Jth species
C	       HMH0(J)	   H(T)-H(298.15) for Jth species, kcal/mol
C	       S0(J)	   S0(T) for Jth species, cal/mol-K
C
C	   Properties of the entire mixture:
C	       C	   sound speed, m/s (if calculated)
C	       H	   mixture enthalpy, J/kg
C	       S	   mixture entropy, J/kg-K
C	       U	   mixture internal energy, J/kg
C	       V	   mixture specific volume, m**3/kg
C	       WM	   molal mass of mixture, g/mol
C
C	   Distributed properties of the mixture:
C	       PMOL(M)	   (relative) mols of the Mth phase
C	       SMOL(J)	   (relative) mols of Jth species
C	       WMP(M)	   molal mass of Mth phase, g/mol
C	       X(J)	   phase mol fraction of Jth species
C	       XM(J)	   mixture mol fraction Jth species
C	       YM(J)	   mixture mass fraction Jth species
C
C	   Detonation properties (if C-J option).
C	       CDET	   C-J detonation wavespeed, m/s (if calculated)
C	       H1	   Chap. Jouguet unburned enthalpy, J/kg
C	       P1	   Chap. Jouguet unburned pressure, Pa
C	       T1	   Chap. Jouguet unburned temperature, K
C	       V1	   Chap. Jouguet volume, m**3/kg
C	       T2	   Chap. Jouguet detonation temperature, K ????
C
C	   Flags:
C	       KERR	   error flag
C
C	   Other data of possible interest:
C	       NB	   number of independent atoms
C	       IB(K)	   atom index of Kth indpendent atom
C	       ELAM(K)	   element potential of Kth independent atom
C-------------------------------------------------------------------
C    Check KERR = 0 for successful calculation.
C         LOAD THE STANJAN COMMON BLOCKS
C
      COMMON /IPAR/ NIWORK, NICMP, NIKNT, NRWORK, NRADD, NSML,
     2              NHML, NWT, NDEN, NXCON, NKCON, NAMAX,
     3              NPHASE, NX1,   NX2,   NY1,    NY2,   NT1,   NT2,
     4              NP1,    NP2,   NV1,   NV2,    NWM1,  NWM2,
     5              NS1,    NS2,   NU1,   NU2,    NH1,   NH2,
     6              NC1,    NC2,   NCDET, NTEST,  NPEST, NSMOL1,
     7              NSMOL2, NYG1,  NXP1,  NXP2,   NWMG1, NWMG2,
     8              NVG1,   NVG2, KKSURF, LSURF
C
      REQWRK(2) = HUGE
      KERR = .FALSE.
      CALL SJSPTS( NAMAX, NPHASE, NSMAX, NIWORK, NRWORK, NRADD,
     1             REQWRK, LOUT)
C
C         SET THE RUN-TIME MONITOR
C
      IEQWRK(IoKMON) = KMON
      IF (KMON .NE. 0) THEN
         IEQWRK(IoKUMO) = LOUT
      ELSE
         IEQWRK(IoKUMO) = 40
      ENDIF
C
      DO 260 M = 1, MM
         REQWRK(IoPA+M) = 0.0
         DO 255 K = 1, KKGAS
C           SET THE RELATIVE ATOM POPULATIONS FROM REAC(K)
            RCMP = FLOAT (IEQWRK(NICMP + (K-1)*MM + M - 1))
            X    = REQWRK(NX1 + K - 1)
            REQWRK(IoPA+M) = REQWRK(IoPA+M) + X * RCMP
            IEQWRK(IoN+M+LoN*K) = IEQWRK(NICMP+(K-1)*MM+M-1)
  255    CONTINUE
         IF (NNBULK .GT. 0) THEN
            KNUM = KKGAS
            DO 256 K = KFIRST(NFBULK), KLAST(NLBULK)
               KNUM = KNUM + 1
               RCMP = FLOAT (IEQWRK(NICMP + (K-1)*MM + M - 1))
               X    = REQWRK(NX1 + K - 1)
               REQWRK(IoPA+M) = REQWRK(IoPA+M) + X * RCMP
               IEQWRK(IoN+M+LoN*KNUM) = IEQWRK(NICMP+(K-1)*MM+M-1)
  256       CONTINUE
         ENDIF
  260 CONTINUE
C
C       MODIFY NA FOR CONSTRAINTS
C       IEQWRK(IoNA) = NAMAX
       IF (FAZEQ) THEN
         IEQWRK(IoNA) = MM + 2
       ELSE
         IEQWRK(IoNA) = MM + NCON
       ENDIF
       IF (IEQWRK(IoNA) .GT. NSMAX) THEN
          WRITE (LOUT,'(/5X,A,/,3(5X,A,I6)/)')
     1    'WARNING: TOO MANY CONSTRAINTS.',
     2    '  CONSTRAINTS = ', NCON,
     3    '  ELEMENTS    = ', MM,
     4    '  SPECIES     = ', NSMAX
          KERR = .TRUE.
          RETURN
       ENDIF
       IEQWRK(IoNP) = NPHASE
       IEQWRK(IoNS) = NSMAX
       IEQWRK(IoNSP + 1) = KKGAS
       IF (NNBULK .GT. 0) THEN
          NP = 1
          DO 275 N = NFBULK, NLBULK
             NP = NP + 1
             IEQWRK(IoNSP + NP) = KKPHAS(N)
275       CONTINUE
       ENDIF
C
C      ITERATION LOOP FOR VARIABLE MEAN MOLECULAR WEIGHT
C
      ITER = 0
      TMOL = 1.
310   CONTINUE
C
      IF (FAZEQ) THEN
C
C        PHASE EQUIL: CONSTRAIN THE TOTAL MOLE FRACTION OF
C        SPECIES THAT DO NOT HAVE A CONDENSED PHASE
C
         M1 = MM + 1
         M2 = MM + 2
cm         M3 = MM + 3
         IF (ITER .EQ. 0) THEN
           XCTOT = 0.
           DO 315 N = 1, NCON
             XCON = REQWRK(NXCON + N - 1)
             XCTOT = XCTOT + XCON
315        CONTINUE
         ENDIF
         REQWRK(IoPA+M1) = XCTOT * TMOL
         REQWRK(IoPA+M2) = 0.
C
C        SET "ATOM POPULATION" MATRIX FOR CONSTRAINT
C
         DO 320 K = 1, NSMAX
           IEQWRK(IoN+M1+LoN*K) = 0
           IEQWRK(IoN+M2+LoN*K) = 0
320      CONTINUE
         DO 322 N = 1, NSMAX
           XCON = REQWRK(NXCON + N - 1)
           KCON = IEQWRK(NKCON + N - 1)
           IF (XCON .GT. 1.E-20) THEN
             IEQWRK(IoN+M1+LoN*KCON) = 1
           ELSE
             IEQWRK(IoN+M2+LoN*KCON) = 1
           ENDIF
322      CONTINUE
C
      ELSE
       DO 350 N = 1, NCON
         M = MM + N
         XCON = REQWRK(NXCON + N - 1)
         KCON = IEQWRK(NKCON + N - 1)
C
C        SET CONSTRAINTS ON INDIVIDUAL SPECIES MOLE FRACTIONS
C
         REQWRK(IoPA+M) = XCON * TMOL
C
C        SET "ATOM POPULATION" MATRIX FOR CONSTRAINTS
C
         DO 325 K = 1, NSMAX
            IEQWRK(IoN+M+LoN*K) = 0
  325    CONTINUE
         IEQWRK(IoN+M+LoN*KCON) = 1
  350  CONTINUE
      ENDIF
c
c   check by printing population matrix
c
C*****checkprint
C      if (iter .eq. 0) then
C        write (6,'(//a/)') 'Population matrix:'
C        do 360 m=1,ieqwrk(iona)
C          write (6,'(11(1x,i3),/)') (ieqwrk(ion+m+lon*k),k=1,nsmax)
C360     continue
C      endif
C*****END checkprint
C
C     FROZEN COMPOSITION OPTION
C
      IEQWRK(IoKFRZ) = - KFRZ
C
C     IoJFS = species data pointers
C     IoDCS = condensed phase densities
C     IoWMS = molar masses
C
      DO 500 K = 1, KKGAS
         IEQWRK(IoJFS + K) = K
         REQWRK(IoDCS + K) = REQWRK(NDEN + K - 1) *1.E+3
         REQWRK(IoWMS + K) = REQWRK(NWT  + K - 1)
         REQWRK(IoSMOL + K) = REQWRK(NX1 + K - 1)
500   CONTINUE
C
      IF (NNBULK .GT. 0) THEN
         KNUM = KKGAS
         DO 505 K = KFIRST(NFBULK), KLAST(NLBULK)
            KNUM = KNUM + 1
            IEQWRK(IoJFS + KNUM) = KNUM
            REQWRK(IoDCS + KNUM) = REQWRK(NDEN + K - 1) *1.E+3
            REQWRK(IoWMS + KNUM) = REQWRK(NWT  + K - 1)
            REQWRK(IoSMOL + KNUM) = REQWRK(NX1 + K - 1)
  505   CONTINUE
      ENDIF
C
C         SET THE SELECTED PROPERTIES
C
      REQWRK(IoP) = REQWRK(NP1)
      IF ((NOP.EQ.2).OR.(NOP.EQ.3).OR.(NOP.EQ.7).OR.(NOP.EQ.9))
     1    REQWRK(IoP) = REQWRK(NPEST)
C
      REQWRK(IoT) = REQWRK(NT1)
      IF ((NOP.GE.4) .AND. (NOP.LE.9)) REQWRK(IoT) = REQWRK(NTEST)
C
      REQWRK(IoH) = REQWRK(NH1)
      REQWRK(IoS) = REQWRK(NS1)
      REQWRK(IoU) = REQWRK(NU1)
      REQWRK(IoV) = REQWRK(NV1)
      REQWRK(IoTE) = REQWRK(NTEST)
C
C      CALL STANJAN:
C      NRADD is NSW
C      REQWRK(NRWORK+1) is SW(1), used as extra space in SJTPRP
C
       CALL SJROP( NAMAX, NSMAX, NIWORK, NRWORK, NRADD,
     1             ATOM, CEQWRK, IEQWRK, REQWRK, REQWRK(NRWORK+1), NOP,
     2             ICKWRK, RCKWRK, ISKWRK, RSKWRK,
     3             IRGWRK, RRGWRK )
C
       IF (IEQWRK(IoKERR) .NE. 0) THEN
           WRITE (LOUT,'(/A)') ' ERROR RETURNING FROM STANJAN !'
           WRITE (LOUT,'(A,I3)') ' KERR =', IEQWRK(IoKERR)
           KERR = .TRUE.
       ENDIF
C
C       CHECK CHANGE IN MEAN MOLECULAR WEIGHT
C
       IF (NCON .GT. 0) THEN
         TMOLL = TMOL
         TMOL = 0.
         DO 700 K = 1, KKGAS
           TMOL = TMOL + REQWRK(IoSMOL+K)
700      CONTINUE
         IF (NNBULK .GT. 0) THEN
           KNUM = KKGAS
           DO 710 K = KFIRST(NFBULK), KLAST(NLBULK)
             KNUM = KNUM + 1
             TMOL = TMOL + REQWRK(IoSMOL+KNUM)
710        CONTINUE
         ENDIF
         DEL = ABS(TMOL - TMOLL)/TMOL
         ITER = ITER + 1
         IF (ITER .GT. 20) THEN
           WRITE (6,*) 'ITERATION NOT CONVERGED!'
           WRITE (6,*) 'TMOL =', TMOL
           KERR = .TRUE.
           RETURN
         ENDIF
         IF (DEL .GT. 1.E-4) GO TO 310
      ENDIF
C
C           RETRIEVE THE STANJAN RESULTS
C
      DO 800 K = 1, KKGAS
         REQWRK(NX2 + K - 1) = REQWRK(IoXM+K)
         REQWRK(NY2 + K - 1) = REQWRK(IoYM+K)
800   CONTINUE
      IF (NNBULK .GT. 0) THEN
         KNUM = KKGAS
         DO 805 K = KFIRST(NFBULK), KLAST(NLBULK)
            KNUM = KNUM + 1
            REQWRK(NX2 + K - 1) = REQWRK(IoXM+KNUM)
            REQWRK(NY2 + K - 1) = REQWRK(IoYM+KNUM)
  805    CONTINUE
      ENDIF
C
      REQWRK(NP2) = REQWRK(IoP)
      REQWRK(NT2) = REQWRK(IoT)
      REQWRK(NH2) = REQWRK(IoH)
      REQWRK(NU2) = REQWRK(IoU)
      REQWRK(NS2) = REQWRK(IoS)
      REQWRK(NV2) = REQWRK(IoV)
      REQWRK(NC2) = REQWRK(IoC)
      REQWRK(NWM2)= REQWRK(IoWM)
      REQWRK(NCDET) = REQWRK(IoCDET)
C
C     end of SUBROUTINE EQRUN
      RETURN
      END
C
C-------------------------------------------------------------------
C
      SUBROUTINE EQSOL (KKGAS, REQWRK, X, Y, T, P, H, V, S, WM, C, CDET,
     1                  KFIRST, KLAST, NNBULK, NFBULK, NLBULK  )
C
C  START PROLOGUE
C
C  RETURNS THE EQUILIBRIUM SOLUTION IN CGS UNITS
C
C  END PROLOGUE
C
C*****precision > double
        IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C        IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      DIMENSION REQWRK(*), X(*), Y(*), KFIRST(*), KLAST(*)
      COMMON /IPAR/ NIWORK, NICMP, NIKNT, NRWORK, NRADD, NSML,
     2              NHML, NWT, NDEN, NXCON, NKCON, NAMAX,
     3              NPHASE, NX1,   NX2,   NY1,    NY2,   NT1,   NT2,
     4              NP1,    NP2,   NV1,   NV2,    NWM1,  NWM2,
     5              NS1,    NS2,   NU1,   NU2,    NH1,   NH2,
     6              NC1,    NC2,   NCDET, NTEST,  NPEST, NSMOL1,
     7              NSMOL2, NYG1,  NXP1,  NXP2,   NWMG1, NWMG2,
     8              NVG1,   NVG2, KKSURF, LSURF
C
      DO 50 K = 1, KKGAS
         X(K) = REQWRK(NX2 + K - 1)
         Y(K) = REQWRK(NY2 + K - 1)
   50 CONTINUE
      IF (NNBULK .GT. 0) THEN
        DO 60 K = KFIRST(NFBULK), KLAST(NLBULK)
          X(K) = REQWRK(NX2 + K - 1)
          Y(K) = REQWRK(NY2 + K - 1)
   60   CONTINUE
      ENDIF
      P = REQWRK(NP2)
      T = REQWRK(NT2)
      V = REQWRK(NV2)
      WM = REQWRK(NWM2)
      S = REQWRK(NS2)
      H = REQWRK(NH2)
C      U = REQWRK(NU2)
      C = REQWRK(NC2)
      CDET = REQWRK(NCDET)
C
C     end of SUBROUTINE EQSOL
      RETURN
      END
C
C-------------------------------------------------------------------
C
      SUBROUTINE EQSTRT (LIN, LOUT, LENICK, LENRCK, LENCCK, LINKCK,
     *                    LENISK, LENRSK, LENCSK, LINKSK, MDIM,
     1                    KDIM, MXPHAS, ICKWRK, RCKWRK, CCKWRK,
     2                    ISKWRK, RSKWRK, CSKWRK, KKPHAS, KFIRST,
     3                    KLAST, NNPHAS, NNBULK, NFBULK, NLBULK, ACT,
     4                    MM, KKGAS, KKSURF, KKBULK, KKTOT, ATOM, KSYM,
     5                    PSYM, NOP, KMON, KFRZ, REAC, TEST, PEST,
     6                    LCNTUE, NCON, KCON, XCON, VINPUT, REVERS,
     7                    FAZEQ,
     8                    LENIRG, LENRRG, LENCRG, IRGWRK, RRGWRK,CRGWRK,
     9                    LINKRG, LSURF, KERR)
C
C  START PROLOGUE
C
C  Input:  LIN  - unit number for input data
C          LOUT - unit number for output messages
C          LENICK - length of CHEMKIN integer array
C          LENRCK - length of CHEMKIN real array
C          LENCCK - length of CHEMKIN character array
C          LINKCK - unit number for CHEMKIN linking file
C          LENISK - length of SURFACE integer array
C          LENRSK - length of SURFACE real array
C          LENCSK - length of SURFACE character array
C          LINKSK - unit number for SURFACE linking file
C          MDIM   - dimension for atomic elements
C          KDIM   - dimension for molecular species
C  Output: ICKWRK - CHEMKIN integer array
C          RCKWRK - CHEMKIN real array
C          CCKWRK - CHEMKIN character array
C          ISKWRK - SURFACE integer array
C          RSKWRK - SURFACE real array
C          CSKWRK - SURFACE character array
C          MM     - total number of atomic elements
C          KKTOT  - total number of molecular species
C          ATOM   - character array of element names
C          KSYM   - character array of species names
C          NOP    - integer equilibrium option number
C          KMON   - integer monitor flag
C          REAC   - real array of input mole fractions
C                   of the mixture
C          TEMP   - starting temperature for the problem
C          TEST   - estimated equilibrium temperature
C          PRES   - starting pressure for the problem
C          PEST   - estimated equilibrium pressure
C          LCNTUE - logical flag for continuation (initialize
C                   in calling program)
C          NCON   - integer number of species to be held constant
C          KCON   - integer array of index numbers of the species
C                   to be held constant
C          XCON   - real mole fractions input for the species to
C                   be held constant
C          FAZEQ  - logical to constrain species that do not have
C                   a bulk phase 'brother'
C          EOS    - character that defines which equation of state
C                   used by RG-CHEMKIN.  The default is the
C                   Ideal gas EOS ('IGA')
C          IRGWRK - RG-CHEMKIN integer array
C          RRGWRK - RG-CHEMKIN real array
C          CRGWRK - RG-CHEMKIN character array
C
C  END PROLOGUE
C
C*****precision > double
        IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C        IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      DIMENSION REAC(*), XCON(*), KCON(*), RCKWRK(*), ICKWRK(*),
     1          RSKWRK(*), ISKWRK(*), KKPHAS(*), KFIRST(*),
     2          KLAST(*), ACT(*), VINPUT(*), VLAST(6),
     3          IRGWRK(*), RRGWRK(*)
      LOGICAL LCNTUE, IERR, REVERS, BADLIN, QUIT, LIVE, LPREV, FAZEQ
      LOGICAL KERR
      CHARACTER KEYWRD*4, LINE*76, EOS*15, CKCHUP*4
      EXTERNAL CKCHUP
      CHARACTER*3 EOSDAT(7)
      CHARACTER*(*) CCKWRK(*), CSKWRK(*), ATOM(*), KSYM(*), PSYM(*),
     1              CRGWRK(*)
      DATA NVAR/6/, NTIN/1/, NPIN/2/, NVIN/3/, NHIN/4/, NUIN/5/,
     1     NSIN/6/, EOS/'IGA'/
      DATA EOSDAT /'IGA', 'VDW', 'RED', 'SOV', 'PEN', 'BKW', 'NBA'/
      DATA IEOS /0/
C
      QUIT = .FALSE.
      KERR = .FALSE.
      IF (.NOT. LCNTUE) THEN
C
C            INITIALIZE VARIABLES
C
         DO 10 K = 1, KDIM
            REAC(K) = 0.
            ACT(K)  = 0.
            XCON(K) = 0.
10       CONTINUE
         DO 12 K = 1, NVAR
            VINPUT(K) = 0.
12       CONTINUE
         TEST = 0.0
         PEST = 0.0
         NCON = 0
         NOP  = 0
         KMON = 0
         KFRZ = 0
         NNBULK = 0
         REVERS = .FALSE.
         LIVE   = .FALSE.
         LPREV  = .FALSE.
         FAZEQ  = .FALSE.
C
         CALL CKINIT (LENICK, LENRCK, LENCCK, LINKCK, LOUT,
     1                ICKWRK, RCKWRK, CCKWRK, IFLAG)
         IF (IFLAG .NE. 0) THEN
              WRITE (LOUT, *) ' STOP. ERROR FROM CKINIT.'
              KERR = .TRUE.
              RETURN
         ENDIF
C
C              THE MANUAL GIVES THE FOLLOWING CALL:
C         CALL RGCOMBLK (LOUT, IRGWRK, RRGWRK, IFLAG)
C
          IF (LSURF .LE. 0) THEN
             CALL CKINDX (ICKWRK, RCKWRK, MM, KKGAS, II, NFIT)
C
C           CHECK ELEMENT SPACE
            IF (MDIM .GE. MM) THEN
               CALL CKSYME (CCKWRK, LOUT, ATOM, IERR)
            ELSE
               WRITE (LOUT, *)
     1         ' ERROR...ELEMENT DIMENSION MUST BE AT LEAST ',MM
               QUIT = .TRUE.
            ENDIF
C           CHECK SPECIES SPACE
            KKTOT = KKGAS
            IF (KDIM .GE. KKTOT) THEN
               CALL CKSYMS (CCKWRK, LOUT, KSYM, IERR)
            ELSE
               WRITE (LOUT, *)
     1         ' ERROR...SPECIES DIMENSION MUST BE AT LEAST ',KKTOT
               QUIT = .TRUE.
            ENDIF
C
            NNPHAS = 1
            NNBULK = 0
            NNSURF = 0
         ELSE
            CALL SKINIT (LENISK, LENRSK, LENCSK, LINKSK, LOUT,
     1                    ISKWRK, RSKWRK, CSKWRK, IFLAG)
            IF (IFLAG .NE. 0) THEN
               WRITE (LOUT, *) ' STOP. ERROR FROM SKINIT!'
               KERR = .TRUE.
               RETURN
            ENDIF
            CALL SKINDX (ISKWRK, MM, KKGAS, KKSURF, KKBULK, KKTOT,
     1                   NNPHAS, NNSURF, NFSURF, NLSURF, NNBULK, NFBULK,
     2                   NLBULK, IISUR)
C
C           CHECK ELEMENT SPACE
            IF (MDIM .GE. MM) THEN
               CALL SKSYME (ISKWRK, CSKWRK, LOUT, ATOM, IERR)
            ELSE
               WRITE (LOUT, *)
     1         ' ERROR...ELEMENT DIMENSION MUST BE AT LEAST ',MM
               QUIT = .TRUE.
            ENDIF
C
C           CHECK SPECIES SPACE
            IF (KDIM .GE. KKTOT) THEN
               CALL SKSYMS (ISKWRK, CSKWRK, LOUT, KSYM, IERR)
            ELSE
               WRITE (LOUT, *)
     1         ' ERROR...SPECIES DIMENSION MUST BE AT LEAST ',KKTOT
               QUIT = .TRUE.
            ENDIF
C
C           CHECK PHASE SPACE
            IF (NNPHAS .GT. MXPHAS) THEN
               WRITE (LOUT, *) ' Error...MXPHAS must be at least ',
     1         NNPHAS
               QUIT = .TRUE.
            ELSE
               CALL SKSYMP (ISKWRK, CSKWRK, LOUT, PSYM, IERR)
               CALL SKPKK (ISKWRK, KKPHAS, KFIRST, KLAST)
               IF (KDIM .LT. KKTOT) THEN
                  WRITE (LOUT, *) ' Error...KDIM must be at least ',
     1            KKTOT
                  QUIT = .TRUE.
C               ELSE
C                  IF (NNBULK .GT. 0) THEN
C                     DO 111 N = NFBULK, NLBULK
C                        DO 111 K = KFIRST(N), KLAST(N)
C                           ACT(K) = 1.0 / KKPHAS(N)
C  111                CONTINUE
C                  ENDIF
                ENDIF
             ENDIF
          ENDIF
C
          IF (QUIT) THEN
             KERR = .TRUE.
             RETURN
          ENDIF
C
      ELSE
C       NOT FIRST TIME THROUGH
        LPREV = .TRUE.
      ENDIF
C
      QUIT   = .FALSE.
      BADLIN = .FALSE.
      LCNTUE = .FALSE.
      DO 5 K = 1, NVAR
        VLAST(K)  = VINPUT(K)
        VINPUT(K) = 0.
5     CONTINUE
C
      WRITE (LOUT,'(/A/)') '           KEYWORD INPUT '
C
90    CONTINUE
      KEYWRD = ' '
      LINE = ' '
      IERR = .FALSE.
      READ  (LIN,  '(A)')     LINE
      WRITE (LOUT, '(1X, A)') LINE
      CALL CKDTAB (LINE)
      KEYWRD = CKCHUP(LINE(1:4), 4)
      LINE(1:4) = ' '
C
C               IS THIS A KEYWORD COMMENT?
C
      IF (KEYWRD(1:1) .EQ. '.' .OR. KEYWRD(1:1) .EQ. '/' .OR.
     1    KEYWRD(1:1) .EQ. '!') GO TO 90
      IND = MAX (INDEX(LINE,'!'), INDEX(LINE,'/'))
      IF (IND .GT. 0) LINE(IND:) = ' '
C
      IF (KEYWRD .EQ. 'DIAG') THEN
         KMON = 2
C
C
      ELSEIF (KEYWRD .EQ. 'EOS') THEN
C*****CHEMKIN
         WRITE (LOUT,*) ' WARNING:  TO USE A NON-IDEAL EOS, YOU MUST ',
     $                  ' USE THE REAL GAS VERSION. '
C*****END CHEMKIN
C*****RG CHEMKIN
C         CALL CKCOMP(LINE, EOSDAT, 7, IEOS)
C         IF (IEOS .EQ. 0) THEN
C           WRITE (LOUT,*)
C     1      ' ERROR.  EQUATION OF STATE NOT FOUND. '
C           IERR = .TRUE.
CC
CC          IF THE USER REQUESTS THE 'BKW' EOS, REJECT THE REQUEST.
CC
C        ELSEIF (IEOS .EQ. 6) THEN
C           WRITE (LOUT,*)
C     1      ' ERROR.  BKW EOS NOT ALLOWED. '
C           IERR = .TRUE.
C        ELSE
C           EOS = EOSDAT(IEOS)
CC
CC              IF THE USER REQUESTS A NEW EOS, RE-INITIALIZE
CC                 REAL GAS CHEMKIN
CC
C           CALL RGINIT (LENIRG, LENRRG, LENCRG, LINKRG, LOUT,
C     1                  IRGWRK, RRGWRK, CRGWRK, EOS)
C           CALL RGCOMBLK (LOUT)
C
C        ENDIF
CC
C*****END RG CHEMKIN
C
      ELSEIF (KEYWRD .EQ. 'FROZ') THEN
C        Flag specifying: equilibrium, =0; frozen composition, =1.
         KFRZ = 1
C
      ELSEIF (KEYWRD .EQ. 'FREE') THEN
C        Return KFRZ selection to equilibrium
         KFRZ = 0
         FAZEQ = .FALSE.
         NCON = 0
C
      ELSEIF (KEYWRD .EQ. 'FAZE') THEN
C        Flag to constrain gas-only species
         FAZEQ = .TRUE.
C
      ELSEIF (KEYWRD .EQ. 'TP  ' .OR. KEYWRD .EQ. 'PT  ') THEN
         NOP = 1
      ELSEIF (KEYWRD .EQ. 'TV  ' .OR. KEYWRD .EQ. 'VT  ') THEN
         NOP = 2
      ELSEIF (KEYWRD .EQ. 'TS  ' .OR. KEYWRD .EQ. 'ST  ') THEN
         NOP = 3
      ELSEIF (KEYWRD .EQ. 'PV  ' .OR. KEYWRD .EQ. 'VP  ') THEN
         NOP = 4
      ELSEIF (KEYWRD .EQ. 'HP  ' .OR. KEYWRD .EQ. 'PH  ') THEN
         NOP = 5
      ELSEIF (KEYWRD .EQ. 'PS  ' .OR. KEYWRD .EQ. 'SP  ') THEN
         NOP = 6
      ELSEIF (KEYWRD .EQ. 'VU  ' .OR. KEYWRD .EQ. 'UV  ') THEN
         NOP = 7
      ELSEIF (KEYWRD .EQ. 'VH  ' .OR. KEYWRD .EQ. 'HV  ') THEN
         NOP = 8
      ELSEIF (KEYWRD .EQ. 'VS  ' .OR. KEYWRD .EQ. 'SV  ') THEN
         NOP = 9
      ELSEIF (KEYWRD .EQ. 'CJ  ') THEN
         NOP = 10

      ELSEIF (KEYWRD .EQ. 'CONX') THEN
C
C        Constant mole fraction (X)
        CALL CKSNUM (LINE, 1, LOUT, KSYM, KKGAS, KSPEC, NVAL,
     1                VALUE, IERR)
        NCON = NCON + 1
        IF (IERR .OR. KSPEC.LE.0) THEN
           WRITE (LOUT,'(A)')
     1      ' ERROR READING DATA FOR KEYWORD '//KEYWRD
           IERR = .TRUE.
        ELSEIF (NCON .GT. KKTOT-MM) THEN
           WRITE (LOUT,'(A,I3,A)')
     1      ' ERROR, CAN HAVE NO MORE THAN ', KKTOT-MM, ' SPECIES'
           IERR = .TRUE.
        ELSE
           XCON(NCON) = VALUE
           KCON(NCON) = KSPEC
        ENDIF
C
      ELSEIF (KEYWRD .EQ. 'TEMP') THEN
C        Starting temperature (K)
         CALL CKXNUM (LINE, 1, LOUT, NVAL, VINPUT(NTIN), IERR)
C
      ELSEIF (KEYWRD .EQ. 'PRES') THEN
C        Starting pressure (atm)
         CALL CKXNUM (LINE, 1, LOUT, NVAL, VINPUT(NPIN), IERR)
C
      ELSEIF (KEYWRD .EQ. 'VOL ') THEN
C        Starting volume (cm**3/gm)
         CALL CKXNUM (LINE, 1, LOUT, NVAL, VINPUT(NVIN), IERR)
C
      ELSEIF (KEYWRD .EQ. 'ENTH') THEN
C        Starting enthalpy (erg/gm)
         CALL CKXNUM (LINE, 1, LOUT, NVAL, VINPUT(NHIN), IERR)
C
      ELSEIF (KEYWRD .EQ. 'ENTR') THEN
C        Starting enthalpy (erg/gm)
         CALL CKXNUM (LINE, 1, LOUT, NVAL, VINPUT(NSIN), IERR)
C
      ELSEIF (KEYWRD .EQ. 'ENGY') THEN
C        Starting energy (erg/gm)
         CALL CKXNUM (LINE, 1, LOUT, NVAL, VINPUT(NUIN), IERR)
C
      ELSEIF (KEYWRD .EQ. 'TEST') THEN
C        Estimate of equilibrium temperature (K)
         CALL CKXNUM (LINE, 1, LOUT, NVAL, TEST, IERR)
C
      ELSEIF (KEYWRD .EQ. 'PEST') THEN
C        Estimate of equilibrium pressure (atm)
         CALL CKXNUM (LINE, 1, LOUT, NVAL, PEST, IERR)
C
      ELSEIF (KEYWRD .EQ. 'REAC') THEN
C        Reactant species
         IF (LSURF .LE. 0) THEN
            CALL CKSNUM (LINE, 1, LOUT, KSYM, KKGAS, KSPEC, NVAL,
     1                   VALUE, IERR)
         ELSE
            CALL SKSNUM (LINE, 1, LOUT, KSYM, KKTOT, PSYM, NNPHAS,
     1                   KKPHAS, KSPEC, NT, NVAL, VALUE, IERR)
         ENDIF
         IF (IERR .OR. KSPEC.LE.0) THEN
            WRITE (LOUT,'(A)')
     1      ' ERROR READING DATA FOR KEYWORD '//KEYWRD
            IERR = .TRUE.
         ELSE
            REAC(KSPEC) = VALUE
         ENDIF
C
      ELSEIF (KEYWRD .EQ. 'CNTN') THEN
C        More jobs to follow
         LCNTUE = .TRUE.
C
      ELSEIF (KEYWRD .EQ. 'LIVE') THEN
C        User is live and can respond to bad lines
         LIVE = .TRUE.
C
      ELSEIF (KEYWRD(:3) .EQ. 'END') THEN
         GO TO 1000
C
      ELSE
C        Unknown keyword
         WRITE (LOUT, *) ' Unrecognized input...',KEYWRD
         IERR = .TRUE.
C
      ENDIF
      BADLIN = BADLIN.OR.IERR
      GO TO 90
C
 1000 CONTINUE
C
C     Check for completeness
C
      SUM = 0.
      DO 1100 K = 1, NVAR
        SUM = SUM + VINPUT(K)
1100  CONTINUE
      IF ( VINPUT(NTIN)*VINPUT(NPIN) .GT. 0.) THEN
C       INPUT T,P AS USUAL
        REVERS = .FALSE.
        DO 1105 K = 1, NVAR
          VLAST(K) = 0.
1105    CONTINUE
      ELSEIF (SUM .NE. 0.) THEN
C       INPUT OTHER VARIABLES, REVERSE PROBLEM
        REVERS = .TRUE.
        DO 1107 K = 1, NVAR
          VLAST(K) = 0.
1107    CONTINUE
      ELSEIF (SUM .EQ. 0) THEN
C       NO NEW INPUTS, CHECK FOR CONTINUATION
        IF (LPREV) THEN
          DO 1110 K = 1, NVAR
            VINPUT(K) = VLAST(K)
1110      CONTINUE
        ELSE
          WRITE (LOUT, '(//,2X,A//)')
     1    ' Error...must provide two properties: T, P, V, H, U, or S'
          BADLIN = .TRUE.
        ENDIF
      ENDIF
C
C*****RG CHEMKIN
C      IF (IEOS .EQ. 0) THEN
C
C         EOS = 'PEN'
C         IEOS = 5
C         WRITE (LOUT,*) 'USING DEFAULT EOS:  IDEAL GAS', eos
CC
CC           INITIALIZE REAL GAS CHEMKIN WITH DEFAULT EOS
CC
C         CALL RGINIT (LENIRG, LENRRG, LENCRG, LINKRG, LOUT,
C     1                IRGWRK, RRGWRK, CRGWRK, EOS)
C         CALL RGCOMBLK (LOUT)
C
C      ENDIF
CC
C*****END RG CHEMKIN
C
      IF ((NOP.EQ.2 .OR. NOP.EQ.3 .OR. NOP.EQ.7 .OR. NOP.EQ.8 .OR.
     1    NOP.EQ.9) .AND. PEST.LE.0.0) THEN
        IF (VINPUT(NPIN).GT.0.) THEN
          PEST = VINPUT(NPIN)
        ELSE
          PEST = 1.
        ENDIF
      ENDIF
C
      IF (NOP.GT.3 .AND. TEST.LE.0.0) THEN
        IF (VINPUT(NTIN).GT.0.) THEN
          TEST = VINPUT(NTIN)
        ELSE
          TEST = 300.
        ENDIF
      ENDIF
C
      IF (BADLIN.AND.(.NOT.LIVE) ) THEN
         WRITE (LOUT, '(//5X,A//)') ' STOP due to errors in input.'
         KERR = .TRUE.
         RETURN
      ENDIF
C
      WRITE (LOUT, *)
      IF (NOP .EQ. 1) THEN
         WRITE (LOUT, *) ' Constant temperature and pressure problem. '
      ELSEIF (NOP .EQ. 2) THEN
         WRITE (LOUT, *) ' Constant temperature and volume problem. '
      ELSEIF (NOP .EQ. 3) THEN
         WRITE (LOUT, *) ' Constant temperature and entropy problem. '
      ELSEIF (NOP .EQ. 4) THEN
         WRITE (LOUT, *) ' Constant pressure and volume problem. '
      ELSEIF (NOP .EQ. 5) THEN
         WRITE (LOUT, *) ' Constant pressure and enthalpy problem. '
      ELSEIF (NOP .EQ. 6) THEN
         WRITE (LOUT, *) ' Constant pressure and entropy problem. '
      ELSEIF (NOP .EQ. 7) THEN
         WRITE (LOUT, *) ' Constant volume and energy problem. '
      ELSEIF (NOP .EQ. 8) THEN
         WRITE (LOUT, *) ' Constant volume and enthalpy problem. '
      ELSEIF (NOP .EQ. 9) THEN
         WRITE (LOUT, *) ' Constant volume and entropy problem. '
      ELSEIF (NOP .EQ. 10) THEN
         WRITE (LOUT, *) ' Chapman=Jouguet detonation. '
      ENDIF
      IF (FAZEQ) WRITE (LOUT,'(/,2X,A/)') 'Phase equilibrium only'
C
C     end of SUBROUTINE EQSTRT
      RETURN
      END
C
      SUBROUTINE EQUIL (LOUT, LPRNT, LSAVE, EQST, LCNTUE, ICKWRK,
     1           RCKWRK, ISKWRK, RSKWRK, KKPHAS, KFIRST, KLAST,
     2           NNPHAS, NNBULK, NFBULK, NLBULK, ACT, PDEN,
     3           LENIEQ, IEQWRK, LENREQ, REQWRK, LENCEQ, CEQWRK,
     3           MM, KKGAS, KKSUR, KKBULK, KKTOT, ATOM, KSYM, NOP,
     4           KMON, KFRZ, REAC, TEST, PEST, NCON, KCON, XCON,
     5           VINPUT, REVERS, FAZEQ, IRGWRK, RRGWRK, LSUR, 
     6           KERR)
C
C  START PROLOGUE
C
C  INPUT:  LOUT - UNIT NUMBER FOR OUTPUT MESSAGES
C          LPRNT  - LOGICAL, PRINT TEXT IF TRUE
C          LSAVE  - UNIT NUMBER FOR BINARY SOLUTION FILE
C          EQST   - LOGICAL, INITIALIZE STANJAN IF TRUE
C          LCNTUE - LOGICAL, TRUE FOR CONTINUATION
C          ICKWRK - CHEMKIN INTEGER ARRAY
C          RCKWRK - CHEMKIN REAL ARRAY
C          ISKWRK - SURFACE CHEMKIN INTEGER ARRAY
C          RSKWRK - SURFACE CHEMKIN REAL ARRAY
C          KKPHAS - PHASE SPECIES INDICIES, DIMENSION NPHASE
C          KFIRST - 1ST SPECIES INDEX IN PHASE, DIMENSION NPHASE
C          KLAST  - LAST SPECIES INDEX IN PHASE, DIMENSION NPHASE
C          NNPHAS - TOTAL NUMBER OF PHASES
C          NNBULK - TOTAL NUMBER OF BULK SPECIES
C          NFBULK - FIRST BULK SPECIES INDEX
C          NLBULK - LAST BULK SPECIES INDEX
C          ACT    - SPECIES ACTIVITIES, DIMENSION KKTOT
C          PDEN   - BULK SPECIES DENSITY, DIMENSION KKTOT
C          IEQWRK - INTEGER EQUIL WORK SPACE, DIMENSION LENIEQ
C          REQWRK - REAL EQUIL WORK SPACE, DIMENSION LENREQ
C          CEQWRK - CHARACTER EQUIL WORK SPACE, DIMENSION LENCEQ
C          MM     - TOTAL NUMBER OF ATOMIC ELEMENTS
C          KKGAS  - NUMBER OF GAS SPECIES
C          KKSUR  - NUMBER OF SURFACE SPECIES
C          KKBULK - NUMBER OF BULK SPECIES
C          KKTOT  - TOTAL NUMBER OF SPECIES
C          ATOM   - CHARACTER ARRAY OF ELEMENT NAMES, DIMENSION KKTOT
C          KSYM   - CHARACTER ARRAY OF SPECIES NAMES, DIMENSION KKTOT
C          NOP    - INTEGER EQUILIBRIUM OPTION NUMBER
C          KMON   - INTEGER MONITOR FLAG
C          KFRZ   - INTEGER FLAG:
C                   = 0 FOR CHEMICAL EQUILIBRIUM (DEFAULT),
C                   = 1 FOR FROZEN COMPOSITION
C          REAC   - ARRAY, INPUT MOLE FRACTIONS, DIMENSION KKTOT
C          TEST   - ESTIMATED EQUILIBRIUM TEMPERATURE
C          PEST   - ESTIMATED EQUILIBRIUM PRESSURE
C          NCON   - INTEGER NUMBER OF SPECIES TO BE HELD CONSTANT
C          KCON   - ARRAY, INDICIES OF SPECIES TO BE HELD CONSTANT
C          XCON   - ARRAY, MOLE FRACTIONS FOR CONSTANT SPECIES
C          VINPUT - STARTING VARIABLES FOR THE PROBLEM, DIMENSION 6
C          VINPUT(1) - STARTING TEMPERATURE (K)
C          VINPUT(2) - STARTING PRESSURE (ATM)
C          VINPUT(3) - STARTING VOLUME (CM**3/G)
C          VINPUT(4) - STARTING ENTHALPY (ERG/G)
C          VINPUT(5) - STARTING ENERGY (ERG/G)
C          VINPUT(6) - STARTING ENTROPY (ERG/G)
C          REVERS - LOGICAL, TEMP IS UNKNOWN IF TRUE
C          FAZEQ  - LOGICAL, CONSTRAIN GAS SPECIES IF TRUE
C  END PROLOGUE
C
C*****precision > double
        IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C        IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      SAVE

      DIMENSION ICKWRK(*), RCKWRK(*), ISKWRK(*), RSKWRK(*),
     1          IEQWRK(LENIEQ), REQWRK(LENREQ),
     2          REAC(*), KCON(*), XCON(*), KKPHAS(*),
     3          KFIRST(*), KLAST(*), ACT(*), PDEN(*), VINPUT(*),
     4          IRGWRK(*), RRGWRK(*)
      CHARACTER*(*) KSYM(KKTOT), ATOM(MM), CEQWRK(LENCEQ)
      LOGICAL LPRNT, EQST, LCNTUE, KERR, REVERS, FAZEQ, NOFAZ, FAZCON
C
      COMMON /IPAR/ NIWORK, NICMP, NIKNT, NRWORK, NRADD, NSML,
     2              NHML, NWT, NDEN, NXCON, NKCON, NAMAX,
     3              NPHASE, NX1,   NX2,   NY1,    NY2,   NT1,   NT2,
     4              NP1,    NP2,   NV1,   NV2,    NWM1,  NWM2,
     5              NS1,    NS2,   NU1,   NU2,    NH1,   NH2,
     6              NC1,    NC2,   NCDET, NTEST,  NPEST, NSMOL1,
     7              NSMOL2, NYG1,  NXP1,  NXP2,   NWMG1, NWMG2,
     8              NVG1,   NVG2, KKSURF, LSURF
C
      DATA HUGE/1.0E+33/, SMALL/1.E-20/
      DATA NVAR/6/, NTIN/1/, NPIN/2/, NVIN/3/, NHIN/4/, NUIN/5/, NSIN/6/
C
      LSURF = LSUR
      KKSURF = KKSUR
      KERR   = .FALSE.
      FAZCON = .FALSE.
      NOFAZ  = .FALSE.
      NCONEQ = 0

      IF (NOP.LT.1 .OR. NOP.GT.10) THEN
          WRITE (LOUT,'(//5X,A//)') 'MUST HAVE NOP BETWEEN 1 AND 10.'
          RETURN
      ENDIF

      IF (.NOT. EQST) THEN

C        MOVE THIS LOOP INSIDE 1ST CALL
         DO 5 I = 1, LENREQ
           REQWRK(I) = 0.
  5     CONTINUE

C
C        INITIALIZE POINTERS
C
         CALL EQLEN (MM, KKGAS, KKBULK, KKTOT,
     1               NNPHAS, NNBULK, NITOT, NRTOT, NCTOT)
C
         WRITE  (LOUT, 15)
   15    FORMAT (/' EQUIL:  Chemkin interface for Stanjan-III',
     1           /'         CHEMKIN-III Version 3.12, 97/12/23',
C*****precision > double
     2           /'         DOUBLE PRECISION')
C*****END precision > double
C*****precision > single
C     2           /'         SINGLE PRECISION')
C*****END precision > single
C
         lencdum = lenceq
         WRITE (LOUT, 7010) LENIEQ, NITOT, LENREQ, NRTOT, LENCEQ, NCTOT
 7010    FORMAT (/'         WORKING SPACE REQUIREMENTS',
     1           /,'             PROVIDED        REQUIRED ',
     2           /,' INTEGER ' , 2I15,
     3           /,' REAL    ' , 2I15,
     4           /,' CHAR    ' , 2I15)
C
         IF (NITOT.GT.LENIEQ .OR. NRTOT.GT.LENREQ .OR. NCTOT.GT.LENCEQ)
     1      THEN
            WRITE (LOUT, *)
     1      '  FATAL ERROR, NOT ENOUGH WORK SPACE PROVIDED'
            KERR = .TRUE.
            RETURN
         ENDIF
C
         WRITE (LOUT, *) ' STANJAN:  Version 3.95, September 1993 '
         WRITE (LOUT, *) '           W. C. Reynolds, Stanford Univ. '
C
C        INITIALIZE COMMON BLOCK and CHEMKIN VARIABLES
C
         CALL CKRP (ICKWRK, RCKWRK, RU, RUC, PATM)
         DO 50 K = 1, KKGAS
            CEQWRK(K) = KSYM(K)
   50    CONTINUE
         IF (KKBULK .GT. 0) THEN
            KNUM = KKGAS
            DO 100 K = KFIRST(NFBULK), KLAST(NLBULK)
               KNUM = KNUM + 1
               CEQWRK(KNUM) = KSYM(K)
  100       CONTINUE
         ENDIF
C
         IF (LSURF .LE. 0) THEN
            CALL CKNCF (MM, ICKWRK, RCKWRK, IEQWRK(NICMP))
            CALL CKWT (ICKWRK, RCKWRK, REQWRK(NWT))
         ELSE
            CALL SKNCF (MM, ISKWRK, IEQWRK(NICMP))
            CALL SKWT (ISKWRK, RSKWRK, REQWRK(NWT))
         ENDIF
         EQST = .TRUE.
      ENDIF
C
C     CHECK KFRZ FLAG
C
      IF ((KFRZ .LT. 0) .OR. (KFRZ .GT.1)) THEN
          WRITE (LOUT,'(/3(/5X,A))')
     1     'MUST SPECIFY:',
     2     'KFRZ=0,   EQUILIBRIUM COMPOSITION',
     3     'KFRZ=1,   FROZEN COMPOSITION'
          RETURN
      ELSE
          KFRZEQ = KFRZ
      ENDIF
C
C     STACK REACTANT MOLS & MASS INTO WORK ARRAY
C
      REQWRK(NSMOL1) = 0.
      GASMAS = 0.
      DO 250 K = 1, KKGAS
         REQWRK(NSMOL1)  = REQWRK(NSMOL1)  + REAC(K)
         REQWRK(NX1 + K - 1) = REAC(K)
         REQWRK(NY1 + K - 1) = REAC(K) * REQWRK(NWT + K - 1)
         GASMAS = GASMAS + REQWRK(NY1 + K - 1)
  250 CONTINUE
      TOTMOL = REQWRK(NSMOL1)
      TOTMAS = GASMAS
      IF (KKBULK .GT. 0) THEN
         DO 300 N = NFBULK, NLBULK
            REQWRK(NSMOL1+N-1) = 0.
            DO 275 K = KFIRST(N), KLAST(N)
               REQWRK(NX1 + K - 1) = REAC(K)
               REQWRK(NY1 + K - 1) = REAC(K) * REQWRK(NWT + K - 1)
               REQWRK(NSMOL1+N-1) = REQWRK(NSMOL1+N-1) + REAC(K)
               TOTMAS = TOTMAS + REQWRK(NY1 + K - 1)
  275       CONTINUE
            TOTMOL = TOTMOL + REQWRK(NSMOL1+N-1)
  300    CONTINUE
      ENDIF
C
C     NORMALIZE GAS PHASE
C
      IF (REQWRK(NSMOL1) .GT. 0.) THEN
        DO 350 K = 1, KKGAS
          REQWRK(NXP1 + K - 1) = REQWRK(NX1 + K - 1) / REQWRK(NSMOL1)
  350   CONTINUE
      ELSE
        DO 355 K = 1, KKGAS
          REQWRK(NXP1 + K - 1) = 0.
  355   CONTINUE
      ENDIF
C
C     NORMALIZE BULK PHASE(S)
C
      IF (KKBULK .GT. 0) THEN
         DO 370 N = NFBULK, NLBULK
            IF (REQWRK(NSMOL1+N-1) .GT. 0.) THEN
              DO 360 K = KFIRST(N), KLAST(N)
               REQWRK(NXP1 + K - 1) = REQWRK(NX1 + K - 1) /
     1                                   REQWRK(NSMOL1+N-1)
  360         CONTINUE
            ELSE
              DO 365 K = KFIRST(N), KLAST(N)
                 REQWRK(NXP1 + K - 1) = 0.
  365         CONTINUE
            ENDIF
  370    CONTINUE
      ENDIF
C
C     NORMALIZE OVER ALL MIXTURE
C
      DO 450 K = 1, KKGAS
         REQWRK(NX1 + K - 1) = REQWRK(NX1 + K - 1) / TOTMOL
         REQWRK(NY1 + K - 1) = REQWRK(NY1 + K - 1) / TOTMAS
  450 CONTINUE
      IF (KKBULK .GT. 0) THEN
         DO 460 K = KFIRST(NFBULK), KLAST(NLBULK)
           REQWRK(NX1 + K - 1) = REQWRK(NX1 + K - 1) / TOTMOL
           REQWRK(NY1 + K - 1) = REQWRK(NY1 + K - 1) / TOTMAS
  460    CONTINUE
      ENDIF
C
      REQWRK(NWM1) = TOTMAS / TOTMOL
      REQWRK(NPEST) = PEST * PATM
      REQWRK(NTEST) = TEST
C
C     ---- SKIP THERMODYNAMIC PROPERTY EVALUATION IF REVERSE ---
C
      IF (.NOT. REVERS) THEN
        REQWRK(NT1) = VINPUT(NTIN)
        IF (VINPUT(NPIN) .GT. 0.) THEN
          REQWRK(NP1) = VINPUT(NPIN) * PATM
        ELSE
          WRITE (LOUT,'(/3X,A/)') 'WARNING: PRESSURE MUST BE > 0!'
          RETURN
        ENDIF
C
C       THERMO PROPERTIES.  NOTE THAT WE CALL SK BEFORE RG.  THIS
C       PREVENTS THE REAL GAS PROPERTIES FROM BEING OVERWRITTEN BY
C       SURFACE CHEMKIN.
C
       IF (LSURF .LE. 0) THEN
           CALL CKHML (REQWRK(NT1), ICKWRK, RCKWRK, REQWRK(NHML))
           CALL CKSML (REQWRK(NT1), ICKWRK, RCKWRK, REQWRK(NSML))
       ELSE
          CALL SKSML (REQWRK(NT1), ISKWRK, RSKWRK, REQWRK(NSML))
          CALL SKHML (REQWRK(NT1), ISKWRK, RSKWRK, REQWRK(NHML))
C
C         REQWRK(NDEN) ARE BULK DENSITIES IN GM/CM**3
          CALL SKSDEN (ISKWRK, RSKWRK, PDEN)
          CALL CKCOPY (KKGAS, REQWRK(NXP1), ACT)
          IF (NNBULK .GT. 0) THEN
             DO 112 N = NFBULK, NLBULK
                PACT = 1.0 / KKPHAS(N)
                DO 112 K = KFIRST(N), KLAST(N)
                   ACT(K) = PACT
  112        CONTINUE
          ENDIF
          CALL SKDEN  (REQWRK(NP1), REQWRK(NT1), ACT, PDEN,
     1                   ISKWRK, RSKWRK, REQWRK(NDEN))
C         RESET DENSITY TO ZERO FOR GAS PHASE SPECIES FOR STANJAN.
C         OLD VERSION USED ACT(K)=0 FOR GAS, BUT SURF-CHEM-III -> NaN
          IF (NNBULK .GT. 0) THEN
            DO 113 K = 1, KKGAS
               REQWRK(NDEN+K-1) = 0.
113         CONTINUE
          ENDIF
        ENDIF

       REQWRK(NC1) = 0.
        REQWRK(NVG1) = 0.
        IF (REQWRK(NSMOL1) .GT. 0.0) THEN
C*****RG CHEMKIN
C          CALL RGRHOX (REQWRK(NP1), REQWRK(NT1), REQWRK(NXP1), ICKWRK,
C     1                RCKWRK, IRGWRK, RRGWRK, RHO)
C*****END RG CHEMKIN
C*****CHEMKIN
          CALL CKRHOX (REQWRK(NP1), REQWRK(NT1), REQWRK(NXP1), ICKWRK,
     1                RCKWRK, RHO)
C*****END CHEMKIN

          IF (RHO .NE. 0.0) REQWRK(NVG1) = 1.0 / RHO
C         VOLUME (CM**3)
          VGAS1 = REQWRK(NVG1)* GASMAS
          CALL CKMMWX (REQWRK(NXP1), ICKWRK, RCKWRK, REQWRK(NWMG1))
C
          IF (NOP .EQ. 10) THEN
            CALL CKXTY  (REQWRK(NXP1), ICKWRK, RCKWRK, REQWRK(NYG1))
            CALL CKCPBS (REQWRK(NT1), REQWRK(NYG1), ICKWRK, RCKWRK, CP)
            CALL CKCVBS (REQWRK(NT1), REQWRK(NYG1), ICKWRK, RCKWRK, CV)
            GAMMA = CP/CV
            REQWRK(NC1) = SQRT (GAMMA *RU *REQWRK(NT1) /REQWRK(NWMG1))
          ENDIF
C
C       CALL RG-CHEMKIN TO GET PARTIAL MOLAR PROPERTIES
C
C*****RG CHEMKIN
C          CALL RGHPML (REQWRK(NP1), REQWRK(NT1), REQWRK(NXP1),
C     $                 ICKWRK, RCKWRK, IRGWRK, RRGWRK, REQWRK(NHML))
C          CALL RGSPML (REQWRK(NP1), REQWRK(NT1), REQWRK(NXP1),
C     $                 ICKWRK, RCKWRK, IRGWRK, RRGWRK, REQWRK(NSML))
C*****END RG CHEMKIN

        ENDIF
C
C       SUM PARTIAL MOLAR PROPERTIES OVER GAS & BULK PHASES
C
        REQWRK(NH1) = 0.0
        REQWRK(NU1) = 0.0
        REQWRK(NS1) = 0.0
C
        DO 750 K = 1, KKGAS
           REQWRK(NH1) = REQWRK(NH1) + REQWRK(NX1+K-1)*REQWRK(NHML+K-1)
C
C*****RG CHEMKIN
C          REQWRK(NS1) = REQWRK(NS1)
C     1                 + REQWRK(NX1+K-1) * REQWRK(NSML+K-1)
C*****END RG CHEMKIN
C*****CHEMKIN
          IF (REQWRK(NXP1+K-1) .GT. SMALL) THEN
C           ENTROPY OF MIXING & PRESSURE CORRECTION FOR IDEAL GAS
            REQWRK(NS1) = REQWRK(NS1)
     1                 + REQWRK(NX1+K-1) *( REQWRK(NSML+K-1)
     2                  - RU *LOG(REQWRK(NXP1+K-1)*REQWRK(NP1)/PATM) )
          ENDIF
C*****END CHEMKIN
  750   CONTINUE

C
        VBULK1 = 0.0
        TOTBLK = 0.0
        IF (KKBULK .GT. 0) THEN
          DO 800 K = KFIRST(NFBULK), KLAST(NLBULK)
            TOTBLK = TOTBLK + REAC(K)
            VBULK1 = VBULK1 + REAC(K) *REQWRK(NWT+K-1) /REQWRK(NDEN+K-1)
            REQWRK(NH1) = REQWRK(NH1) + REQWRK(NX1+K-1)*REQWRK(NHML+K-1)
            IF (REQWRK(NXP1+K-1) .GT. SMALL) THEN
C             ENTROPY OF MIXING TERMS FOR BULK PHASES
              REQWRK(NS1) = REQWRK(NS1)
     1                 + REQWRK(NX1+K-1) *( REQWRK(NSML+K-1)
     2                  - RU *LOG(REQWRK(NXP1+K-1)) )
            ENDIF
  800     CONTINUE
C         ENTHALPY (PER MOL); PRESSURE CORRECTION FOR BULK PHASES
          IF (TOTBLK .GT. 0.) REQWRK(NH1) = REQWRK(NH1)
     1              + (REQWRK(NP1) - PATM) *VBULK1 /TOTBLK
        ENDIF
C       TOTAL VOLUME = GAS + BULK PHASES
        REQWRK(NV1) = VGAS1 + VBULK1
C       NORMALIZE TOTAL VOLUME BY MASS
        REQWRK(NV1) = REQWRK(NV1) / TOTMAS
C       ENERGY PROPERTIES PER MASS
        REQWRK(NH1) = REQWRK(NH1) / REQWRK(NWM1)
        REQWRK(NS1) = REQWRK(NS1) / REQWRK(NWM1)
C       USE H-PV TO COMPUTE ENERGY.
        REQWRK(NU1) = REQWRK(NH1) - REQWRK(NP1) * REQWRK(NV1)
C
      ELSE
C
C       REVERSE PROBLEM: SET VARIABLES FROM VINPUT
C
        DO 825 N = NT1, NVG2
          REQWRK(N) = 0.
  825   CONTINUE
        IF (NOP .EQ. 2) THEN
          REQWRK(NT1) = VINPUT(NTIN)
          REQWRK(NV1) = VINPUT(NVIN)
          REQWRK(NP1) = REQWRK(NPEST)
        ELSEIF (NOP .EQ. 3) THEN
          REQWRK(NT1) = VINPUT(NTIN)
          REQWRK(NS1) = VINPUT(NSIN)
          REQWRK(NP1) = REQWRK(NPEST)
        ELSEIF (NOP .EQ. 4) THEN
          REQWRK(NP1) = VINPUT(NPIN) * PATM
          REQWRK(NV1) = VINPUT(NVIN)
          REQWRK(NT1) = REQWRK(NTEST)
        ELSEIF (NOP .EQ. 5) THEN
          REQWRK(NP1) = VINPUT(NPIN) * PATM
          REQWRK(NH1) = VINPUT(NHIN)
          REQWRK(NT1) = REQWRK(NTEST)
        ELSEIF (NOP .EQ. 6) THEN
          REQWRK(NP1) = VINPUT(NPIN) * PATM
          REQWRK(NS1) = VINPUT(NSIN)
          REQWRK(NT1) = REQWRK(NTEST)
        ELSEIF (NOP .EQ. 7) THEN
          REQWRK(NV1) = VINPUT(NVIN)
          REQWRK(NU1) = VINPUT(NUIN)
          REQWRK(NT1) = REQWRK(NTEST)
          REQWRK(NP1) = REQWRK(NPEST)
        ELSEIF (NOP .EQ. 8) THEN
          REQWRK(NV1) = VINPUT(NVIN)
          REQWRK(NH1) = VINPUT(NHIN)
          REQWRK(NT1) = REQWRK(NTEST)
          REQWRK(NP1) = REQWRK(NPEST)
        ELSEIF (NOP .EQ. 9) THEN
          REQWRK(NV1) = VINPUT(NVIN)
          REQWRK(NS1) = VINPUT(NSIN)
          REQWRK(NT1) = REQWRK(NTEST)
          REQWRK(NP1) = REQWRK(NPEST)
        ELSEIF ((NOP .EQ. 1) .OR. (NOP .EQ. 10) ) THEN
          WRITE (LOUT,*) ' CANNOT HAVE REVERS=.TRUE. FOR NOP=', NOP
          RETURN
        ENDIF
C
        IF (LSURF .GE. 0) THEN
C          REQWRK(NDEN) ARE BULK DENSITIES IN GM/CM**3
           CALL SKSDEN (ISKWRK, RSKWRK, PDEN)
           CALL SKDEN  (REQWRK(NP1), REQWRK(NT1), ACT, PDEN,
     1                ISKWRK, RSKWRK, REQWRK(NDEN))
        ENDIF
C
      ENDIF
C
C      CONSTRAIN GAS SPECIES IF WANT PHASE EQUIL ONLY
C
      IF (NCON .GT. 0) THEN
        NCONEQ = NCON
      ENDIF
      IF (FAZEQ) THEN
        FAZCON = .TRUE.
        IF (KKBULK .EQ. 0) THEN
C*****checkprint
C          WRITE (LOUT, 9000)
C9000  FORMAT (//5X,' NO NEED FOR PHASE CONSTRAINT:')
C          WRITE (LOUT, '(5X,A//)')
C     1    ' NO BULK SPECIES IN MECHANISM.'
C*****END checkprint
          NOFAZ = .TRUE.
        ELSE
          NCONEQ = 0
          NCONZ  = 0
          SUMXC = 0.
          DO 845 K = 1, KKGAS
C         eliminate phase-brothers using molecular weight
            DO 840 J = KFIRST(NFBULK), KLAST(NLBULK)
              IF ( ABS( REQWRK(NWT+J-1) - REQWRK(NWT+K-1) )
     1           /REQWRK(NWT+K-1) .LE. 1.E-6 ) GO TO 845
840         CONTINUE
C           here if no match
            NCONEQ = NCONEQ +1
            KCON(NCONEQ) = K
            XCON(NCONEQ) = REQWRK(NX1+K-1)
            SUMXC = SUMXC + XCON(NCONEQ)
845       CONTINUE
          IF (NCONEQ .EQ. 0) THEN
C*****checkprint
C            WRITE (LOUT, 9000)
C            WRITE (LOUT,'(5X,A//)')
C     1    ' ALL SPECIES MAY EXIST IN BULK PHASE'
C*****END checkprint
            FAZCON = .FALSE.
          ELSEIF ((NCONEQ.GT.0) .AND. (ABS(SUMXC-1.).LE.1.E-4)) THEN
C*****checkprint
C            WRITE (LOUT, 9000)
C            WRITE (LOUT,'(5X,A//)')
C     1    ' NO BULK SPECIES PRESENT IN MIXTURE.'
C*****END checkprint
            NOFAZ = .TRUE.
          ENDIF
        ENDIF
      ENDIF
C
C     CHECK CONSTRAINTS
C
C*****checkprint
C      IF (NCONEQ .GT. 0) THEN
C        WRITE (LOUT, '(//2X,A)') ' CONSTRAINED SPECIES MOLE FRACTIONS:'
C        DO 846 K = 1, NCONEQ
C          WRITE (LOUT, '(2X,A,2X,1PE12.5)') KSYM(KCON(K)), XCON(K)
C846     CONTINUE
C        WRITE (LOUT,'(//)')
C      ENDIF
C*****END checkprint
C
C
C        CONVERT TO MKS UNITS:
C           P: pressure in Pa (N/M**2)  1 Pa = 10 dynes/cm2
C           T: Kelvins
C           S: J/Kg-K
C           H: J/Kg
C           U: J/Kg
C           V: M**3/Kg
C
      REQWRK(NP1) = REQWRK(NP1) * 1.0E-1
      REQWRK(NPEST) = REQWRK(NPEST) * 1.0E-1
      REQWRK(NH1) = REQWRK(NH1) * 1.0E-4
      REQWRK(NS1) = REQWRK(NS1) * 1.0E-4
      REQWRK(NU1) = REQWRK(NU1) * 1.0E-4
      REQWRK(NV1) = REQWRK(NV1) * 1.0E-3
      REQWRK(NC1) = REQWRK(NC1) * 1.0E-2
      REQWRK(NCDET) = REQWRK(NCDET) * 1.0E-2
C
      DO 850 N = 1, NCONEQ
         REQWRK(NXCON + N - 1) = XCON(N)
         IEQWRK(NKCON + N - 1) = KCON(N)
  850 CONTINUE
C
C     CATCH PHASE-ONLY WHEN NOT NECESSARY
C
      IF (FAZEQ .AND. NOFAZ) THEN
        IF (REVERS) THEN
          KFRZEQ = 1
          NCONEQ = 0
          FAZCON = .FALSE.
        ELSE
C         CONVERT UNITS & COPY STATE 1 TO 2
          REQWRK(NT2) = REQWRK(NT1)
          REQWRK(NP1) = REQWRK(NP1) * 1.0E+1 / PATM
          REQWRK(NP2) = REQWRK(NP1)
          REQWRK(NV2) = REQWRK(NV1)
          REQWRK(NH1) = REQWRK(NH1) * 1.0E+4
          REQWRK(NH2) = REQWRK(NH1)
          REQWRK(NS1) = REQWRK(NS1) * 1.0E+4
          REQWRK(NS2) = REQWRK(NS1)
          REQWRK(NU1) = REQWRK(NU1) * 1.0E+4
          REQWRK(NU2) = REQWRK(NU1)
          REQWRK(NV1) = REQWRK(NV1) * 1.0E+3
          REQWRK(NV2) = REQWRK(NV1)
          REQWRK(NC1) = REQWRK(NC1) * 1.0E+2
          REQWRK(NC2) = REQWRK(NC1)
          REQWRK(NCDET) = REQWRK(NCDET) * 1.0E+2
          REQWRK(NWM2) = REQWRK(NWM1)
          REQWRK(NWMG2) = REQWRK(NWMG1)
          REQWRK(NSMOL2) = REQWRK(NSMOL1)
          DO 900 K = 1, KKGAS
            REQWRK(NX2+K-1) = REQWRK(NX1+K-1)
            REQWRK(NXP2+K-1) = REQWRK(NXP1+K-1)
900       CONTINUE
          IF (KKBULK .GT. 0) THEN
             DO 910 K = KFIRST(NFBULK), KLAST(NLBULK)
               REQWRK(NX2+K-1) = REQWRK(NX1+K-1)
               REQWRK(NXP2+K-1) = REQWRK(NXP1+K-1)
910          CONTINUE
            DO 1000 N = NFBULK, NLBULK
              REQWRK(NSMOL2+N-1) = REQWRK(NSMOL1+N-1)
              DO 940 K = KFIRST(N), KLAST(N)
                REQWRK(NXP2+K-1) = REQWRK(NXP1+K-1)
940           CONTINUE
              IF (REQWRK(NSMOL2+N-1) .GT. 0.) THEN
              DO 960 K = KFIRST(N), KLAST(N)
                REQWRK(NXP2+K-1) = REQWRK(NXP1+K-1)
960           CONTINUE
              ENDIF
1000        CONTINUE
          ENDIF
C         SKIP EQUIL CALC
          GO TO 3000
        ENDIF
      ENDIF
C
C          COMPUTE EQUILIBRIUM COMPOSITION
C
      NSMAX = KKGAS + KKBULK
      CALL EQRUN (KKGAS, NNBULK, NFBULK, NLBULK, KKPHAS, KFIRST, KLAST,
     1            NOP, KMON, LOUT, MM, NCONEQ, ATOM, CEQWRK, IEQWRK,
     2            REQWRK, ICKWRK, RCKWRK, ISKWRK, RSKWRK, NSMAX, KERR,
     3            HUGE, KFRZEQ, FAZCON,
     4            IRGWRK, RRGWRK)
C
C     CONVERT BACK TO CGS
C
      REQWRK(NP1) = REQWRK(NP1) * 1.0E+1 / PATM
      REQWRK(NP2) = REQWRK(NP2) * 1.0E+1 / PATM
      REQWRK(NH1) = REQWRK(NH1) * 1.0E+4
      REQWRK(NH2) = REQWRK(NH2) * 1.0E+4
      REQWRK(NS1) = REQWRK(NS1) * 1.0E+4
      REQWRK(NS2) = REQWRK(NS2) * 1.0E+4
      REQWRK(NU1) = REQWRK(NU1) * 1.0E+4
      REQWRK(NU2) = REQWRK(NU2) * 1.0E+4
      REQWRK(NV1) = REQWRK(NV1) * 1.0E+3
      REQWRK(NV2) = REQWRK(NV2) * 1.0E+3
      REQWRK(NC1) = REQWRK(NC1) * 1.0E+2
      REQWRK(NC2) = REQWRK(NC2) * 1.0E+2
      REQWRK(NCDET) = REQWRK(NCDET) * 1.0E+2
C
C     COMPUTE ABSOLUTE MOLES USING MASS, MEAN WEIGHT
C     AND MOLE FRACTIONS IN NX2,
C
      TOTMOL = TOTMAS / REQWRK(NWM2)
      DO 2000 K = 1, KKGAS
        REQWRK(NX2+K-1) = REQWRK(NX2+K-1) *TOTMOL
2000  CONTINUE
      IF (NNBULK .GT. 0) THEN
         DO 2010 K = KFIRST(NFBULK), KLAST(NLBULK)
           REQWRK(NX2+K-1) = REQWRK(NX2+K-1) *TOTMOL
2010     CONTINUE
      ENDIF
C
C     NOW NORMALIZE WITHIN EACH PHASE
C
      REQWRK(NSMOL2) = 0.
      DO 2020 K = 1, KKGAS
         REQWRK(NXP2+K-1) = REQWRK(NX2+K-1)
         REQWRK(NSMOL2) = REQWRK(NSMOL2) + REQWRK(NXP2+K-1)
2020  CONTINUE
      IF (REQWRK(NSMOL2) .GT. 0.) THEN
        DO 2030 K = 1, KKGAS
          REQWRK(NXP2+K-1) = REQWRK(NXP2+K-1) / REQWRK(NSMOL2)
2030    CONTINUE
      ENDIF
      CALL CKMMWX (REQWRK(NXP2), ICKWRK, RCKWRK, REQWRK(NWMG2))
      P2 = REQWRK(NP2)*PATM
C*****CHEMKIN
      CALL CKRHOX (P2, REQWRK(NT2), REQWRK(NXP2), ICKWRK, RCKWRK, RHO)
C*****END CHEMKIN
C*****RG CHEMKIN
C      CALL RGRHOX (P2, REQWRK(NT2), REQWRK(NXP2), ICKWRK, RCKWRK,
C     $             IRGWRK, RRGWRK, RHO)
C*****END RG CHEMKIN
      IF (RHO .GT. 0.) REQWRK(NVG2) = 1./RHO
      IF (NNBULK .GT. 0) THEN
        DO 2100 N = NFBULK, NLBULK
          REQWRK(NSMOL2+N-1) = 0.
          DO 2040 K = kFIRST(N), KLAST(N)
            REQWRK(NXP2+K-1) = REQWRK(NX2+K-1)
            REQWRK(NSMOL2+N-1) = REQWRK(NSMOL2+N-1) + REQWRK(NXP2+K-1)
2040      CONTINUE
          IF (REQWRK(NSMOL2+N-1) .GT. 0.) THEN
            DO 2060 K = KFIRST(N), KLAST(N)
              REQWRK(NXP2+K-1) = REQWRK(NXP2+K-1)/ REQWRK(NSMOL2+N-1)
2060        CONTINUE
          ENDIF
2100    CONTINUE
      ENDIF
C
C     RENORMALIZE OVER ALL MIXTURE
C
      DO 2110 K = 1, KKGAS
        REQWRK(NX2+K-1) = REQWRK(NX2+K-1) /TOTMOL
2110  CONTINUE
      IF (NNBULK .GT. 0) THEN
         DO 2120 K = KFIRST(NFBULK), KLAST(NLBULK)
           REQWRK(NX2+K-1) = REQWRK(NX2+K-1) /TOTMOL
2120     CONTINUE
      ENDIF
C
3000  CONTINUE
      IF (LPRNT) CALL EQPRNT (KKGAS, KKBULK, KKTOT,
     1           NNBULK, NFBULK, NLBULK, KFIRST, KLAST, NOP, KFRZ,
     2           LOUT, MM, LSAVE, LCNTUE, KSYM, REQWRK)
C
C     end of SUBROUTINE EQUIL
      IF (KERR)  WRITE (LOUT, *) ' Failed due to STANJAN error...'
      RETURN
      END
