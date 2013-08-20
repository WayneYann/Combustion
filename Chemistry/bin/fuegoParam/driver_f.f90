      MODULE procmod
        IMPLICIT NONE
        
        ! Define interface of call-back routine.
        ABSTRACT INTERFACE
           SUBROUTINE callback (n, x, y)
             USE, INTRINSIC :: ISO_C_BINDING
             INTEGER(KIND=C_INT), INTENT(IN), VALUE :: n
             REAL(KIND=C_DOUBLE), INTENT(IN) :: x(n)
             REAL(KIND=C_DOUBLE), INTENT(OUT) :: y
           END SUBROUTINE callback
        END INTERFACE
        
      CONTAINS
        
        ! Define C-bound procedure.
        SUBROUTINE test_observation (cproc, n, x, y) BIND(C)
          USE, INTRINSIC :: ISO_C_BINDING
          TYPE(C_FUNPTR), INTENT(IN), VALUE :: cproc
          INTEGER(KIND=C_INT), INTENT(IN), VALUE :: n
          REAL(KIND=C_DOUBLE), INTENT(IN) :: x(n)
          REAL(KIND=C_DOUBLE), INTENT(OUT) :: y
          
          PROCEDURE(callback), POINTER :: proc

          ! Convert C to Fortran procedure pointer.
          CALL C_F_PROCPOINTER (cproc, proc)
          
          ! Call it.
          CALL proc(n, x, y)
        END SUBROUTINE test_observation
        
      END MODULE procmod
