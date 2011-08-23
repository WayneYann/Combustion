      PROGRAM CKDRIV
      implicit none
      CHARACTER*72 files, CKINP, THINP, TRINP, CKLOG, CKLIN, MCLIN
C
      integer LIN, LINCK, LOUT, LTHRM, LINMC, LTRAN
      PARAMETER (LIN = 34, LINCK = 35, LOUT = 33, LTHRM = 17, LINMC = 37, LTRAN = 36)
C
      integer IDIM, KDIM, MDIM, MAXORD, MAXSP, MAXTB, MAXTP
      PARAMETER
     +  (IDIM = 1000,
     +   KDIM = 200,
     +   MDIM = 20,
     +   MAXORD = 10,
     +   MAXSP = 12,
     +   MAXTB = 10,
     +   MAXTP = 3)
C
      integer LCWORK, LIWORK, LLWORK, LRWORK
      PARAMETER
     +  (LCWORK = KDIM + MDIM,
     +   LIWORK = IDIM * (27 + MAXORD + 2 * MAXSP + MAXTB)
     +      + KDIM * (4 + MDIM),
     +   LLWORK = KDIM,
     +   LRWORK = IDIM * (33 + MAXORD + MAXSP + MAXTB)
     +      + 2 * KDIM * (4 * MAXTP - 3) + MDIM)
C
      integer LCMCWK, LIMCWK, LLMCWK, LRMCWK
      PARAMETER
     +  (LCMCWK = LCWORK + KDIM,
     +   LIMCWK = LIWORK + 3 * KDIM,
     +   LLMCWK = LLWORK,
     +   LRMCWK = LRWORK + 300 + KDIM * (48 + 4 * KDIM + 8 * MAXTP))


      CHARACTER*16 C(LCMCWK)
      INTEGER I(LIMCWK)
      LOGICAL L(LLMCWK)
      DOUBLE PRECISION R(LRMCWK)
C
      read(5,*) CKINP
      read(5,*) THINP
      read(5,*) TRINP
      read(5,*) CKLIN
      read(5,*) MCLIN
      read(5,*) CKLOG

      OPEN (LIN, STATUS = 'OLD', FORM = 'FORMATTED',
     +   FILE = trim(CKINP))
      OPEN (LTHRM, STATUS = 'OLD', FORM = 'FORMATTED',
     +   FILE = trim(THINP))
      OPEN (LTRAN, STATUS = 'OLD', FORM = 'FORMATTED',
     +   FILE = trim(TRINP))
      OPEN (LINCK, STATUS = 'UNKNOWN', FORM = 'FORMATTED',
     +   FILE = trim(CKLIN))
      OPEN (LOUT, STATUS = 'UNKNOWN', FORM = 'FORMATTED',
     +   FILE = trim(CKLOG))
      OPEN (LINMC, STATUS = 'UNKNOWN', FORM = 'FORMATTED',
     +   FILE = trim(MCLIN))
C
      CALL CKINTP
     +  (MDIM, KDIM, IDIM, MAXTP, MAXSP, MAXTB, LIN, LOUT, LTHRM, LINCK,
     +   MAXORD, LIWORK, I, LRWORK, R, LCWORK, C, L)

      CALL TRANFT
     +  (LINMC, LINCK, LTRAN, MAXTP, LIMCWK, LRMCWK, LCMCWK, I, R,
     +   C)


      CLOSE (LIN)
      CLOSE (LTRAN)
      CLOSE (LTHRM)
      CLOSE (LINCK)
      CLOSE (LINMC)
      CLOSE (LOUT)

      END
