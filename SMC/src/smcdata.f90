module smcdata_module

  use multifab_module
  use sdcquad_module

  implicit none

  type(multifab),save :: Uprime, Unew
  type(multifab),save :: Q, Fdif, Upchem
  type(multifab),save :: mu, xi ! viscosity
  type(multifab),save :: lam ! partial thermal conductivity
  type(multifab),save :: Ddiag ! diagonal components of rho * Y_k * D

  private

  public :: Uprime, Unew, Q, Fdif, Upchem, mu, xi, lam, Ddiag
  public :: build_smcdata, destroy_smcdata

contains

  subroutine build_smcdata(la,sdc)
    use variables_module, only : ncons, nprim
    use chemistry_module, only : nspecies
    use derivative_stencil_module, only : stencil_ng
    use probin_module, only : advance_method
    use threadbox_module, only : tb_multifab_setval

    implicit none

    type(layout), intent(in) :: la
    type(sdc_ctx_t),  intent(inout) :: sdc
    integer :: err

    if (advance_method .eq. 1) then ! RK
       call multifab_build(Uprime, la, ncons, 0)
       call multifab_build(Unew,   la, ncons, stencil_ng)
    else
       if (sdc%single_rate) then
          call sdc_imex_allocate(sdc%imex, err)
       end if
       if (sdc%multi_rate) then
          call sdc_mrex_allocate(sdc%mrex, err)
       end if
    end if

    call multifab_build(Q, la, nprim, stencil_ng)
    call tb_multifab_setval(Q, 0.d0, .true.)

    call multifab_build(Fdif, la, ncons, 0)
    call multifab_build(Upchem, la, ncons, 0)
    call tb_multifab_setval(Upchem, 0.d0)

    call multifab_build(mu , la, 1, stencil_ng)
    call multifab_build(xi , la, 1, stencil_ng)
    call multifab_build(lam, la, 1, stencil_ng)
    call multifab_build(Ddiag, la, nspecies, stencil_ng)

  end subroutine build_smcdata


  subroutine destroy_smcdata(sdc)
    use probin_module, only : advance_method
    type(sdc_ctx_t), intent(inout) :: sdc

    if (advance_method .eq. 1) then ! RK
       call destroy(Unew)
       call destroy(Uprime)
    else 
       if (sdc%single_rate) then
          call sdc_imex_deallocate(sdc%imex)
       end if
       if (sdc%multi_rate) then
          call sdc_mrex_deallocate(sdc%mrex)
       end if
    end if

    call destroy(Q)
    call destroy(Fdif)
    call destroy(Upchem)
    call destroy(mu)
    call destroy(xi)
    call destroy(lam)
    call destroy(Ddiag)
  end subroutine destroy_smcdata

end module smcdata_module
