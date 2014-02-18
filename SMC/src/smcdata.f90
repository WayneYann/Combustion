module smcdata_module

  use multifab_module
  use sdcquad_module
  use derivative_stencil_module, only : stencil, narrow, s3d, stencil_ng

  implicit none

  type(multifab),save :: Uprime, Unew
  type(multifab),save :: Q, Fdif, Upchem
  type(multifab),save :: mu, xi ! viscosity
  type(multifab),save :: lam ! partial thermal conductivity
  type(multifab),save :: Ddiag ! diagonal components of rho * Y_k * D
  type(multifab),save :: qx, qy, qz  ! for S3D to store first-derivatives

  private

  public :: Uprime, Unew, Q, Fdif, Upchem, mu, xi, lam, Ddiag, qx, qy, qz
  public :: build_smcdata, destroy_smcdata

contains

  subroutine build_smcdata(la,sdc)
    use variables_module, only : ncons, nprim, ndq
    use chemistry_module, only : nspecies
    use probin_module, only : advance_sdc
    use threadbox_module, only : tb_multifab_setval

    implicit none

    type(layout),  intent(in) :: la
    type(sdc_ctx), intent(inout) :: sdc
    integer :: err

    if (.not. advance_sdc) then
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

    if (stencil .eq. s3d) then
       call multifab_build(qx, la, ndq, stencil_ng)
       call multifab_build(qy, la, ndq, stencil_ng)
       call multifab_build(qz, la, ndq, stencil_ng)
    end if

  end subroutine build_smcdata


  subroutine destroy_smcdata(sdc)
    use probin_module, only : advance_sdc
    type(sdc_ctx), intent(inout) :: sdc

    if (.not. advance_sdc) then
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

    if (stencil .eq. s3d) then
       call destroy(qx)
       call destroy(qy)
       call destroy(qz)
    end if

  end subroutine destroy_smcdata

end module smcdata_module
