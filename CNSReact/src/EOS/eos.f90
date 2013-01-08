module eos_module 
 
! use cdwrk_module
  use chemsolv_module

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
         smalld = 1.d-3
       end if
    else
       smalld = 1.d-3
    endif
 
    smallp = 1.d-6
    smallc = 1.d-6
 
  end subroutine eos_init

  subroutine eos_get_small_temp(small_temp_out)

    double precision, intent(out) :: small_temp_out

    small_temp_out = smallt

  end subroutine eos_get_small_temp

  subroutine eos_get_small_dens(small_dens_out)

    double precision, intent(out) :: small_dens_out

    small_dens_out = smalld

  end subroutine eos_get_small_dens

  subroutine eos_given_ReY(G, P, C, T, dpdr, dpde, R, e, Y, pt_index)
 
     implicit none

     include "cdwrk.H"
 
     double precision, parameter :: SCALP = 0.1d0, SCALR = 0.001d0
     double precision, intent(  out) :: G, P, C, T, dpdr, dpde
     double precision, intent(in   ) :: R, e, Y(:)
     integer, optional, intent(in   ) :: pt_index(:)

     double precision :: cv,cp

     integer NiterCheck

     NiterCheck = T_from_eY(T,Y,e)
     call CKCVBS(T,Y,iwrk,rwrk,cv)
     call CKCPBS(T,Y,iwrk,rwrk,cp)
     G = cp / cv
     call CKPY(R,T,Y,iwrk,rwrk,P)
     P = max(P,smallp)
     C = sqrt(G*ABS(P)/MAX(R,smalld))
     C = max(C,smallc)

!  next two lines merit checking
     dpdr = P/R
     dpde = 1.d0 / cv
 
  end subroutine eos_given_ReY

  subroutine eos_S_given_ReY(S, R, e, T, Y, pt_index)
 
     implicit none
 
     double precision, intent(  out) :: S
     double precision, intent(in   ) :: R, e, T, Y(:)
     integer, optional, intent(in   ) :: pt_index(:)

     S = 0.d0

  end subroutine eos_S_given_ReY

  subroutine eos_given_RTY(e, P, R, T, Y, pt_index)

     implicit none
     include "cdwrk.H"

     ! in/out variables
     double precision, intent(  out) :: e, P
     double precision, intent(in   ) :: R, T, Y(:)
     integer, optional, intent(in   ) :: pt_index(:)

     double precision :: c_v, h

     call CKPY(R,T,Y,iwrk,rwrk,P)
     P = max(P,smallp)

!  hack . . . Fuego doesn't appear to make a CKUBMS function  
    
     call CKSBMS(P,T,Y,iwrk,rwrk,h)

     e = h - P/R

 
  end subroutine eos_given_RTY

  subroutine eos_get_cv(cv, R, T, Y)

     implicit none
     include "cdwrk.H"

! input/output variables
    double precision, intent(out) :: cv
    double precision, intent(in)  :: R, T, Y(:)

     call CKCVBS(T,Y,iwrk,rwrk,cv)

  end subroutine eos_get_cv

end module eos_module 
