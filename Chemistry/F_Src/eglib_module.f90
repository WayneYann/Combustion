module eglib_module

  implicit none

  double precision, allocatable, save :: egwork(:) 
  integer, allocatable, save :: egiwork(:) 

!$omp threadprivate(egwork,egiwork)

  private

  public egwork, egiwork, eglib_init, eglib_close

contains

  subroutine eglib_init(nspecies, NP, ITLS, IFLAG)
    use omp_module
    integer, intent(in) :: nspecies, NP, ITLS, IFLAG

    integer :: legwork, legiwork

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
       
    allocate(egwork(legwork))
    allocate(egiwork(legiwork))
    call egini(NP, 6, IFLAG, ITLS, egwork, legwork, egiwork, legiwork)
  end subroutine eglib_init

  subroutine eglib_close()
    deallocate(egwork)
    deallocate(egiwork)
  end subroutine eglib_close

end module eglib_module
