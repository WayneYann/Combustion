module eos_module 

  implicit none

  double precision, save, public :: gamma_const = 5.d0/3.d0

  double precision, save, private :: smalld = 1.d-50
  double precision, save, private :: smallt = 1.d-50
  double precision, save, private :: smallp = 1.d-50
  double precision, save, private :: smallc = 1.d-50

  logical, save, private :: initialized = .false.

contains

  subroutine eos_init(small_temp, small_dens, gamma_in)

    implicit none
 
    double precision, intent(in), optional :: small_temp
    double precision, intent(in), optional :: small_dens
    double precision, intent(in), optional :: gamma_in
    
    if (present(small_temp)) then
       if (small_temp > 0.d0) then
          smallt = small_temp
       end if
    endif
 
    if (present(small_dens)) then
       if (small_dens > 0.d0) then
          smalld = small_dens
       end if
    endif
 
    if (present(gamma_in)) then
       gamma_const = gamma_in
    end if

    initialized = .true.
    
  end subroutine eos_init

  subroutine eos_get_small_temp(small_temp_out)
    
    double precision, intent(out) :: small_temp_out
 
    small_temp_out = smallt
 
  end subroutine eos_get_small_temp
 
  subroutine eos_get_small_dens(small_dens_out)
    
    double precision, intent(out) :: small_dens_out
    
    small_dens_out = smalld
    
  end subroutine eos_get_small_dens


  subroutine eos_get_pcg(p,c,g,rho,e,T,xn)
    double precision, intent(out) :: p,c,g
    double precision, intent(in) :: rho, e, T, xn(0)
    p = (gamma_const - 1.0) * rho * e
    p = max(p, smallp)
    c = sqrt(gamma_const*p/rho)
    g = gamma_const
  end subroutine eos_get_pcg

end module eos_module
