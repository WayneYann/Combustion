module eglib_module

  implicit none

  integer, parameter :: NP    = 1
  integer, parameter :: ITLS  = 3
  integer, parameter :: IFLAG = 7
  integer, save :: legwork, legiwork
  double precision, allocatable, save :: egwork(:) 
  integer, allocatable, save :: egiwork(:) 

  private

  public egwork, egiwork, legwork, legiwork, eglib_init, eglib_close

contains

  subroutine eglib_init(nspecies)
    integer, intent(in) :: nspecies

    legwork = 23 + 14*nspecies + 32*nspecies**2 + 13*NP  & 
         + 30*NP*nspecies + 5*NP*nspecies**2
    legiwork = nspecies

    allocate(egwork(legwork))
    allocate(egiwork(legiwork))

    call egini(NP, 6, IFLAG, ITLS, egwork, legwork, egiwork, legiwork)

  end subroutine eglib_init


  subroutine eglib_close()
    deallocate(egwork,egiwork)
  end subroutine eglib_close

end module eglib_module
