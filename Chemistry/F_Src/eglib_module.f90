module eglib_module

  implicit none

  integer, save :: NP
  integer, parameter :: ITLS  = 1
  integer, parameter :: IFLAG = 5
  integer, save :: legwork, legiwork
  double precision, allocatable, save :: egwork(:) 
  integer, allocatable, save :: egiwork(:) 

!$omp threadprivate(egwork,egiwork)

  private

  public egwork, egiwork, legwork, legiwork, eglib_init, eglib_close 

contains

  subroutine eglib_init(nspecies, np_in)
    use omp_module
    integer, intent(in) :: nspecies, np_in

    NP = np_in

    ! for ITLS=1
    if (ITLS .eq. 1) then
       legwork = 23 + 14*nspecies + 32*nspecies**2 + 13*NP &
            + 14*NP*nspecies + NP*nspecies**2
    else if (ITLS .eq. 2) then
       legwork = 23 + 14*nspecies + 32*nspecies**2 + 13*NP &
            + 21*NP*nspecies + NP*(2*nspecies**2+(nspecies*(nspecies+1))/2)
    else
       legwork = 23 + 14*nspecies + 32*nspecies**2 + 13*NP  & 
            + 30*NP*nspecies + 5*NP*nspecies**2
    end if

    legiwork = nspecies
    
    !$omp parallel 
    allocate(egwork(legwork))
    allocate(egiwork(legiwork))
    call egini(NP, 6, IFLAG, ITLS, egwork, legwork, egiwork, legiwork)
    !$omp end parallel

  end subroutine eglib_init


  subroutine eglib_close()
    !$omp parallel 
    deallocate(egwork,egiwork)
    !$omp end parallel
  end subroutine eglib_close

end module eglib_module
