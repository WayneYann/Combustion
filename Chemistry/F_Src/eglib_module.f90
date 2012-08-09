module eglib_module

  implicit none

  integer, parameter :: NP    = 1
  integer, parameter :: ITLS  = 3
  integer, parameter :: IFLAG = 7
  integer, save :: legwork, legiwork
  double precision, allocatable, save :: egwork(:) 
  integer, allocatable, save :: egiwork(:) 

!$omp threadprivate(egwork,egiwork)

  private

  public egwork, egiwork, legwork, legiwork, eglib_init, eglib_close

contains

  subroutine eglib_init(nspecies)
    use omp_module
    integer, intent(in) :: nspecies

    legwork = 23 + 14*nspecies + 32*nspecies**2 + 13*NP  & 
         + 30*NP*nspecies + 5*NP*nspecies**2
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
