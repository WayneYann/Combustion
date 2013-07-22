module vode_module

  implicit none

  integer, save :: verbose, itol, neq, order
  integer, save :: MF, use_ajac, save_ajac, stiff

  double precision, save :: rtol, atol

  double precision, allocatable, save :: voderwork(:), voderpar(:)
  integer, allocatable, save :: vodeiwork(:), vodeipar(:)
  integer, save :: lvoderwork, lvodeiwork

!$omp threadprivate(voderwork,vodeiwork,voderpar,vodeipar)

contains

  subroutine vode_init(neq_in)

    implicit none

    integer, intent(in) :: neq_in

    neq = neq_in

    if (stiff .eq. 0) then

       MF = 10
       lvoderwork = 20+16*NEQ
       lvodeiwork = 30

    else 

       if (use_ajac .ne. 0) then
          if (save_ajac .ne. 0) then
             MF = 21
             lvoderwork = 22 + 9*NEQ + 2*NEQ**2
          else
             MF = -21
             lvoderwork = 22 + 9*NEQ +   NEQ**2
          end if
       else
          MF = 22
          lvoderwork = 22 + 9*NEQ + 2*NEQ**2
       end if

       lvodeiwork = 30 + NEQ

    end if

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

