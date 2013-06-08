module eos_module 

  use chemistry_module, only : Ru

  implicit none

  double precision, save, private :: smalld
  double precision, save, private :: smallt
  double precision, save, private :: smallp
  double precision, save, private :: smallc

contains

  subroutine eos_init(small_temp, small_dens)

    implicit none
 
    double precision, intent(in), optional :: small_temp
    double precision, intent(in), optional :: small_dens
    
    if (present(small_temp)) then
       if (small_temp > 0.d0) then
          smallt = small_temp
       else
          smallt = 250.d0
       end if
    else
       smallt = 250.d0
    endif
 
    if (present(small_dens)) then
       if (small_dens > 0.d0) then
          smalld = small_dens
       else
          smalld = 1.d-50
       end if
    else
       smalld = 1.d-50
    endif
 
    smallp = 1.d-50
    smallc = 1.d-50
    
  end subroutine eos_init


end module eos_module
