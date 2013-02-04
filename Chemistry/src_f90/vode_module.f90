module vode_module

  implicit none

  integer, save :: itol, neq
  double precision, save :: rtol, atol

  integer, parameter :: MF_NOSTIFF = 10, MF_NUMERICAL_JAC = 22

  double precision, allocatable, save :: voderwork(:), voderpar(:)
  integer, allocatable, save :: vodeiwork(:), vodeipar(:)
  integer, save :: lvoderwork, lvodeiwork

!$omp threadprivate(voderwork,vodeiwork)

contains

  subroutine vode_init(neq_in)

    implicit none

    integer, intent(in) :: neq_in

    neq = neq_in

    lvoderwork = max(22 + 9*NEQ + 2*NEQ**2, 20+16*NEQ)
    lvodeiwork = 30 + NEQ

    !$omp parallel
    allocate(voderwork(lvoderwork))
    allocate(vodeiwork(lvodeiwork))

    voderwork = 0.d0
    vodeiwork = 0

    allocate(voderpar(2))
    allocate(vodeipar(1))
    !$omp end parallel

  end subroutine vode_init


  subroutine vode_close()
    neq = -1
    !$omp parallel
    deallocate(voderwork)
    deallocate(vodeiwork)
    deallocate(voderpar)
    deallocate(vodeipar)
    !$omp end parallel
  end subroutine vode_close

end module vode_module

