module sdcquad_module

  use bl_types
  use sdclib
  use sdclib_multifab

  implicit none

  type :: ctx_t
     double precision :: dx(3)
  end type ctx_t

  type :: sdc_t
     logical    :: single_rate, multi_rate
     integer    :: iters = 4
     real(dp_t) :: tol_residual = -1.0d0

     type(mf_encap_t), pointer  :: mfencap
     type(sdc_encap_t), pointer :: encap

     type(sdc_nset_t), pointer  :: nset_adr, nset_ad, nset_r
     type(sdc_exp_t), pointer   :: exp_adr, exp_ad, exp_r
     type(sdc_mrset_t), pointer :: mrset
     type(sdc_srset_t), pointer :: srset
  end type sdc_t

contains


  !
  ! Build/create single-rate SDC object.
  !
  subroutine sdc_build_single_rate(sdc, qtype, nnodes, ctx, feval, post)
    type(sdc_t),    intent(out) :: sdc
    integer,        intent(in)  :: qtype, nnodes
    type(ctx_t),    intent(in), target :: ctx
    type(c_funptr), intent(in), value :: feval, post

    integer :: err

    sdc%single_rate = .true.
    sdc%multi_rate  = .false.

    allocate(sdc%mfencap, sdc%encap, sdc%nset_adr, sdc%exp_adr, sdc%srset)

    call sdc_multifab_build(sdc%encap, c_loc(sdc%mfencap), err)
    call sdc_nset_build(sdc%nset_adr, nnodes, qtype, "ADR", err)
    call sdc_exp_build(sdc%exp_adr, feval, "ADR", err)
    call sdc_hook_add(sdc%exp_adr%hooks, SDC_HOOK_POST_STEP, post, err)
    call sdc_srset_build(sdc%srset, sdc%nset_adr, c_loc(sdc%exp_adr), sdc%encap, c_loc(ctx), "ADR", err)

  end subroutine sdc_build_single_rate


  !
  ! Build/create multi-rate SDC object.
  !
  subroutine sdc_build_multi_rate(sdc, qtype, nnodes, ctx, f1eval, f2eval, post)
    use probin_module, only: sdc_multirate_type

    type(sdc_t),    intent(out) :: sdc
    integer,        intent(in)  :: qtype, nnodes(2)
    type(ctx_t),    intent(in), target :: ctx
    type(c_funptr), intent(in), value :: f1eval, f2eval, post

    integer :: err

    sdc%single_rate = .false.
    sdc%multi_rate  = .true.

    allocate(sdc%mfencap, sdc%encap, sdc%nset_ad, sdc%nset_r, sdc%exp_ad, sdc%exp_r, sdc%mrset)

    call sdc_multifab_build(sdc%encap, c_loc(sdc%mfencap), err)
    call sdc_nset_build(sdc%nset_ad, nnodes(1), qtype, "AD", err)
    call sdc_nset_build(sdc%nset_r, nnodes(2), qtype, "R", err)
    call sdc_exp_build(sdc%exp_ad, f1eval, "AD", err)
    call sdc_exp_build(sdc%exp_r, f2eval, "R", err)

    call sdc_mrset_build(sdc%mrset, 2, "ADR", err)
    call sdc_mrset_add_nset(sdc%mrset, sdc%nset_ad, c_loc(sdc%exp_ad), sdc%encap, c_loc(ctx), 0, err)

    call sdc_hook_add(sdc%mrset%hooks, SDC_HOOK_PRE_UPDATE, post, err)

    select case(sdc_multirate_type)
    case ("local")
       call sdc_mrset_add_nset(sdc%mrset, sdc%nset_r, c_loc(sdc%exp_r), sdc%encap, c_loc(ctx), SDC_MR_LOCAL, err)
    case ("global")
       call sdc_mrset_add_nset(sdc%mrset, sdc%nset_r, c_loc(sdc%exp_r), sdc%encap, c_loc(ctx), SDC_MR_GLOBAL, err)
    case ("repeated")
       call sdc_mrset_add_nset(sdc%mrset, sdc%nset_r, c_loc(sdc%exp_r), sdc%encap, c_loc(ctx), SDC_MR_REPEATED, err)
    case default
       stop "UNKNOWN MULTIRATE TYPE: should be one of 'local', 'global', or 'repeated'"
    end select
  end subroutine sdc_build_multi_rate


  !
  ! Set multifab layout, feval context etc
  !
  subroutine sdc_set_layout(sdc, la, nc, ng)
    type(sdc_t),  intent(inout) :: sdc
    type(layout), intent(in)    :: la
    integer,      intent(in)    :: nc, ng

    sdc%mfencap%nc = nc
    sdc%mfencap%ng = ng
    sdc%mfencap%la = la
  end subroutine sdc_set_layout


  !
  ! Setup and allocate SDC object.
  !
  subroutine sdc_setup(sdc, la, nc, ng)
    type(sdc_t), intent(inout) :: sdc
    type(layout), intent(in)   :: la
    integer,      intent(in)   :: nc, ng

    integer :: err

    call sdc_set_layout(sdc, la, nc, ng)

    if (sdc%single_rate) then
       call sdc_srset_setup(sdc%srset, err)
    else 
       call sdc_mrset_setup(sdc%mrset, err)
    end if
    
    if (err .ne. 0) then
       stop "FAILED TO SETUP SDC SETS"
    end if
  end subroutine sdc_setup


  !
  ! Destroy SDC object.
  !
  subroutine sdc_destroy(sdc)
    type(sdc_t), intent(inout) :: sdc

    if (sdc%single_rate) then
       call sdc_srset_destroy(sdc%srset)
       call sdc_nset_destroy(sdc%nset_adr)
       call sdc_exp_destroy(sdc%exp_adr)
       deallocate(sdc%exp_adr, sdc%nset_adr, sdc%srset)
    end if

    if (sdc%multi_rate) then
       call sdc_mrset_destroy(sdc%mrset)
       call sdc_nset_destroy(sdc%nset_ad)
       call sdc_nset_destroy(sdc%nset_r)
       call sdc_exp_destroy(sdc%exp_ad)
       call sdc_exp_destroy(sdc%exp_r)
       deallocate(sdc%exp_ad, sdc%exp_r, sdc%nset_ad, sdc%nset_r, sdc%mrset)
    end if

    call sdc_multifab_destroy(sdc%encap)
    deallocate(sdc%mfencap)
  end subroutine sdc_destroy

end module sdcquad_module
