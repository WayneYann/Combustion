module vode_module

  implicit none

  integer, save :: verbose, itol, neq, order, use_ajac
  double precision, save :: rtol, atol

  integer, save :: MF

  double precision, allocatable, save :: voderwork(:), voderpar(:)
  integer, allocatable, save :: vodeiwork(:), vodeipar(:)
  integer, save :: lvoderwork, lvodeiwork

!$omp threadprivate(voderwork,vodeiwork)

contains

  subroutine vode_init(neq_in)

    implicit none

    integer, intent(in) :: neq_in

    neq = neq_in

    if (use_ajac .eq. 0) then

       lvoderwork = max(22 + 9*NEQ + 2*NEQ**2, 20+16*NEQ)

       MF = -21
       
    else

       lvoderwork = max(22 + 9*NEQ +   NEQ**2, 20+16*NEQ)

       MF = 22

    end if

    lvodeiwork = 30 + NEQ

    !$omp parallel
    allocate(voderwork(lvoderwork))
    allocate(vodeiwork(lvodeiwork))

    voderwork = 0.d0
    vodeiwork = 0
    vodeiwork(5) = order

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

