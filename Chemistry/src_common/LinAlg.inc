!
! F90 interface definitions for (selected) routines in LinAlg.f
!

interface
  subroutine dgefa(a, lda, n, ipvt, info)
    integer, intent(in) :: lda, n
    integer, intent(out) :: info, ipvt(n)
    double precision, intent(inout) :: a(lda, n)
  end subroutine dgefa
  subroutine dgesl(a, lda, n, ipvt, b, job)
    integer, intent(in)             :: lda, n, job, ipvt(n)
    double precision, intent(inout) :: a(lda,n), b(n)
  end subroutine dgesl
end interface
