module eglib_module

  implicit none

  integer, parameter :: NP    = 1
  integer, parameter :: ITLS  = 3
  integer, parameter :: IFLAG = 7
  integer, save :: legwork, legiwork
  double precision, allocatable, save :: egwork(:,:) 
  integer, allocatable, save :: egiwork(:,:) 

  private

  public egwork, egiwork, legwork, legiwork, eglib_init, eglib_close

contains

  subroutine eglib_init(nspecies)
    use omp_module
    integer, intent(in) :: nspecies

    integer :: TID, nmaxt

    legwork = 23 + 14*nspecies + 32*nspecies**2 + 13*NP  & 
         + 30*NP*nspecies + 5*NP*nspecies**2
    legiwork = nspecies
    
    !$omp parallel
    !$omp master
    nmaxt = omp_get_num_threads()
    !$omp end master
    !$omp end parallel

    allocate(egwork(legwork,   0:nmaxt-1))
    allocate(egiwork(legiwork, 0:nmaxt-1))

    !$omp parallel private(TID)
    TID = omp_get_thread_num()
    call egini(NP, 6, IFLAG, ITLS, egwork(:,TID), legwork, egiwork(:,TID), legiwork)
    !$omp end parallel

  end subroutine eglib_init


  subroutine eglib_close()
    deallocate(egwork,egiwork)
  end subroutine eglib_close

end module eglib_module
