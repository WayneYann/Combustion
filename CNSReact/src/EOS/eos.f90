module eos_module 

  use chemistry_module, only : nspecies, Ru

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

  subroutine eos_get_small_temp(small_temp_out)

    double precision, intent(out) :: small_temp_out
    
    small_temp_out = smallt
    
  end subroutine eos_get_small_temp
  
  subroutine eos_get_small_dens(small_dens_out)
    
    double precision, intent(out) :: small_dens_out
    
    small_dens_out = smalld
    
  end subroutine eos_get_small_dens


  subroutine eos_get_c(c, rho, T, Y, pt_index)

    implicit none

    double precision, intent(out) :: c
    double precision, intent(in) :: rho, T, Y(nspecies)
    integer, optional, intent(in) :: pt_index(:)

    integer :: iwrk
    double precision :: rwrk, Cv, gamma, p

    call ckcvbs(T, Y, iwrk, rwrk, Cv)

    gamma = (Cv + Ru) / Cv

    call ckpy(rho, T, Y, iwrk, rwrk, p)

    c = sqrt(gamma*p/rho)
    c = max(c, smallc)

  end subroutine eos_get_c


  subroutine eos_get_T(T, e, Y, pt_index)

    implicit none

    double precision, intent(out) :: T
    double precision, intent(in) :: e, Y(nspecies)
    integer, optional, intent(in) :: pt_index(:)

    integer :: iwrk
    double precision :: rwrk

    call feeytt(e, Y, iwrk, rwrk, T)
    T = max(T, smallt)

  end subroutine eos_get_T


  subroutine eos_get_p(p, rho, T, Y, pt_index)

    implicit none

    double precision, intent(out) :: p
    double precision, intent(in) :: rho, T, Y(nspecies)
    integer, optional, intent(in) :: pt_index(:)

    integer :: iwrk
    double precision :: rwrk

    call ckpy(rho, T, Y, iwrk, rwrk, p)
    p = max(p,smallp)

  end subroutine eos_get_p


  subroutine eos_get_TP(T, P, rho, e, Y, pt_index)

    implicit none

    double precision, intent(out) :: T, P
    double precision, intent(in ) :: rho, e, Y(:)
    integer, optional, intent(in   ) :: pt_index(:)
    
    integer :: iwrk
    double precision :: rwrk

    call feeytt(e, Y, iwrk, rwrk, T)
    T = max(T, smallt)

    call ckpy(rho, T, Y, iwrk, rwrk, p)
    p = max(p,smallp)    

  end subroutine eos_get_TP


  subroutine eos_get_eP(e, P, rho, T, Y, pt_index)
    
    implicit none
    
    double precision, intent(  out) :: e, P
    double precision, intent(in   ) :: rho, T, Y(:)
    integer, optional, intent(in   ) :: pt_index(:)
    
    integer :: iwrk
    double precision :: rwrk

    call ckpy(rho, T, Y, iwrk, rwrk, P)
    P = max(P,smallp)

    call ckubms(T,Y,iwrk,rwrk,e)

  end subroutine eos_get_eP


  subroutine eos_given_RTY(e, P, C, G, dpdr, dpde, rho, T, Y, pt_index)
    
    implicit none
    
    double precision, intent(  out) :: e, P, C, G, dpdr, dpde
    double precision, intent(in   ) :: rho, T, Y(:)
    integer, optional, intent(in   ) :: pt_index(:)
    
    integer :: iwrk
    double precision :: rwrk, Cv

    call ckpy(rho, T, Y, iwrk, rwrk, P)
    P = max(P,smallp)

    call ckubms(T,Y,iwrk,rwrk,e)

    call ckcvbs(T, Y, iwrk, rwrk, Cv)
    G = (Cv + Ru) / Cv

    C = sqrt(G*P/rho)
    C = max(C, smallc)

    dpdr = P/rho
    dpde = P/(Cv*T)

  end subroutine eos_given_RTY


  subroutine eos_given_ReY(G, P, C, T, dpdr, dpde, rho, e, Y, pt_index)
 
    implicit none
    
    double precision, intent(  out) :: G, P, C, T, dpdr, dpde
    double precision, intent(in   ) :: rho, e, Y(:)
    integer, optional, intent(in   ) :: pt_index(:)

    integer :: iwrk
    double precision :: rwrk, Cv

    call feeytt(e, Y, iwrk, rwrk, T)
    T = max(T, smallt)
    
    call ckcvbs(T, Y, iwrk, rwrk, Cv)
    G = (Cv + Ru) / Cv

    call ckpy(rho, T, Y, iwrk, rwrk, P)

    C = sqrt(G*P/rho)
    C = max(C, smallc)

    dpdr = P/rho
    dpde = P/(Cv*T)
    
  end subroutine eos_given_ReY

end module eos_module 
