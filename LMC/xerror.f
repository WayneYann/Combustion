      INTEGER FUNCTION LUNSAV (IVALUE, ISET)
      LOGICAL ISET
      INTEGER IVALUE
C-----------------------------------------------------------------------
C LUNSAV saves and recalls the parameter LUNIT which is the logical
C unit number to which error messages are printed.
C
C Saved local variable..
C
C  LUNIT   = Logical unit number for messages.
C            The default is 6 (machine-dependent).
C
C On input..
C
C   IVALUE = The value to be set for the LUNIT parameter,
C            if ISET is .TRUE. .
C
C   ISET   = Logical flag to indicate whether to read or write.
C            If ISET=.TRUE., the LUNIT parameter will be given
C            the value IVALUE.  If ISET=.FALSE., the LUNIT
C            parameter will be unchanged, and IVALUE is a dummy
C            parameter.
C
C On return..
C
C   The (old) value of the LUNIT parameter will be returned
C   in the function value, LUNSAV.
C
C This is a modification of the SLATEC library routine J4SAVE.
C
C Subroutines/functions called by LUNSAV.. None
C-----------------------------------------------------------------------
      INTEGER LUNIT
C-----------------------------------------------------------------------
C The following Fortran-77 declaration is to cause the values of the
C listed (local) variables to be saved between calls to this integrator.
C-----------------------------------------------------------------------
      SAVE LUNIT
      DATA LUNIT /6/
C
      LUNSAV = LUNIT
      IF (ISET) LUNIT = IVALUE
      END

      INTEGER FUNCTION MFLGSV (IVALUE, ISET)
      LOGICAL ISET
      INTEGER IVALUE
C-----------------------------------------------------------------------
C MFLGSV saves and recalls the parameter MESFLG which controls the
C printing of the error messages.
C
C Saved local variable..
C
C   MESFLG = Print control flag..
C            1 means print all messages (the default).
C            0 means no printing.
C
C On input..
C
C   IVALUE = The value to be set for the MESFLG parameter,
C            if ISET is .TRUE. .
C
C   ISET   = Logical flag to indicate whether to read or write.
C            If ISET=.TRUE., the MESFLG parameter will be given
C            the value IVALUE.  If ISET=.FALSE., the MESFLG
C            parameter will be unchanged, and IVALUE is a dummy
C            parameter.
C
C On return..
C
C   The (old) value of the MESFLG parameter will be returned
C   in the function value, MFLGSV.
C
C This is a modification of the SLATEC library routine J4SAVE.
C
C Subroutines/functions called by MFLGSV.. None
C-----------------------------------------------------------------------
      INTEGER MESFLG
C-----------------------------------------------------------------------
C The following Fortran-77 declaration is to cause the values of the
C listed (local) variables to be saved between calls to this integrator.
C-----------------------------------------------------------------------
      SAVE MESFLG
      DATA MESFLG /1/
C
      MFLGSV = MESFLG
      IF (ISET) MESFLG = IVALUE
      END
