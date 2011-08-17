module castro_burner_module

  use network

contains

  subroutine burner(dens, temp, Xin, ein, dt, time, Xout, eout)

    implicit none
    
    double precision, intent(in) :: dens, temp, Xin(nspec), ein, dt, time
    double precision, intent(out) :: Xout(nspec), eout
    
    integer :: n
    double precision :: enuc, dX
    
    Xout(:) = Xin(:)
    
    enuc = 0.0d0
    do n = 1, nspec
       dX = Xout(n)-Xin(n) 
       enuc = enuc - ebin(n) * dX
    enddo
  
    eout = ein + enuc
  
  end subroutine burner

end module castro_burner_module
