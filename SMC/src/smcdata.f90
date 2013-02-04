module smcdata_module

  use multifab_module

  implicit none

  type(multifab),save :: Uprime, Unew
  type(multifab),save :: Q, Fdif
  type(multifab),save :: mu, xi ! viscosity
  type(multifab),save :: lam ! partial thermal conductivity
  type(multifab),save :: Ddiag ! diagonal components of rho * Y_k * D

  private

  public :: Uprime, Unew, Q, Fdif, mu, xi, lam, Ddiag
  public :: build_smcdata, destroy_smcdata

contains

  subroutine build_smcdata(la)

    use variables_module, only : ncons, nprim
    use chemistry_module, only : nspecies
    use derivative_stencil_module, only : stencil_ng
    use probin_module, only : advance_method

    implicit none

    type(layout), intent(in) :: la

    if (advance_method .eq. 1) then ! RK
       call multifab_build(Uprime, la, ncons, 0)
       call multifab_build(Unew,   la, ncons, stencil_ng)
    end if

    call multifab_build(Q, la, nprim, stencil_ng)

    call multifab_build(Fdif, la, ncons, 0)

    call multifab_build(mu , la, 1, stencil_ng)
    call multifab_build(xi , la, 1, stencil_ng)
    call multifab_build(lam, la, 1, stencil_ng)
    call multifab_build(Ddiag, la, nspecies, stencil_ng)

  end subroutine build_smcdata


  subroutine destroy_smcdata()
    use probin_module, only : advance_method
    if (advance_method .eq. 1) then ! RK
       call destroy(Unew)
       call destroy(Uprime)
    end if
    call destroy(Q)
    call destroy(Fdif)
    call destroy(mu)
    call destroy(xi)
    call destroy(lam)
    call destroy(Ddiag)
  end subroutine destroy_smcdata

end module smcdata_module

