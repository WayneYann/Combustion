module eos_module 
 
  implicit none

  ! Avogradro's Number
  double precision, parameter :: n_A    = 6.02214129d23   ! mol^-1

  ! Boltzmann's constant
  double precision, parameter :: k_B    = 1.3806488d-16   ! erg/K

  double precision, save, private :: smalld
  double precision, save, private :: smallt
  double precision, save, private :: smallp
  double precision, save, private :: smallc
  double precision, save          :: gamma_const 

  ! abar is the mean molecular weight.  Ideally, this should
  ! be computed from X_k, using aion and zion from the network

  ! we use a value of 0.6.  This is what you get for a fully
  ! ionized gas consisting of 0.75 H and 0.25 He-4 by mass
  ! (roughly a solar composition)
  double precision, parameter :: Abar = 0.6d0

contains

  subroutine eos_init(small_temp, small_dens, gamma_in)

    implicit none
 
    double precision, intent(in), optional :: small_temp
    double precision, intent(in), optional :: small_dens
    double precision, intent(in), optional :: gamma_in

    if (present(gamma_in)) then
       gamma_const = gamma_in
    else
       gamma_const = 1.4d0
    end if
 
    if (present(small_temp)) then
      if (small_temp > 0.d0) then
       smallt = small_temp
      else
       smallt = 5.d6
      end if
    else
       smallt = 5.d6
    endif
 
    if (present(small_dens)) then
       if (small_dens > 0.d0) then
         smalld = small_dens
       else
         smalld = 1.d-8
       end if
    else
       smalld = 1.d-8
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

  subroutine eos_given_ReX(G, P, C, T, dpdr, dpde, R, e, X, pt_index)
 
     implicit none
 
     double precision, intent(  out) :: G, P, C, T, dpdr, dpde
     double precision, intent(in   ) :: R, e, X(:)
     integer, optional, intent(in   ) :: pt_index(:)

     double precision :: c_v

     G = gamma_const
     P = (gamma_const - 1.d0) * R * e
     P = max(P,smallp)
     C = sqrt(gamma_const*ABS(P)/MAX(R,smalld))
     C = max(C,smallc)

     ! specific heat
     c_v = k_B*n_A/ (abar * (gamma_const - 1.d0))
     T = e / c_v

     dpdr = (gamma_const-1.d0)*e
     dpde = (gamma_const-1.d0)*R
 
  end subroutine eos_given_ReX

  subroutine eos_S_given_ReX(S, R, e, T, X, pt_index)
 
     implicit none
 
     double precision, intent(  out) :: S
     double precision, intent(in   ) :: R, e, T, X(:)
     integer, optional, intent(in   ) :: pt_index(:)

     S = 0.d0

  end subroutine eos_S_given_ReX

  subroutine eos_given_RTX(e, P, R, T, X, pt_index)

     implicit none

     ! in/out variables
     double precision, intent(  out) :: e, P
     double precision, intent(in   ) :: R, T, X(:)
     integer, optional, intent(in   ) :: pt_index(:)

     double precision :: c_v

     c_v = k_B*n_A/ (abar * (gamma_const - 1.d0))
 
     e = T * c_v
     P = (gamma_const - 1.d0) * R * e
     P = max(P,smallp)
 
  end subroutine eos_given_RTX

  subroutine eos_get_cv(cv, R, T, Y)

! input/output variables
    double precision, intent(out) :: cv
    double precision, intent(in)  :: R, T, Y(:)

    cv = k_B*n_A/ (abar * (gamma_const - 1.d0))

  end subroutine eos_get_cv

end module eos_module 
