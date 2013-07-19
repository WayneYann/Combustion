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
     type(sdc_mrex_t), pointer  :: mrex
     type(sdc_imex_t), pointer  :: imex
     type(ctx_t),      pointer  :: ctx
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
    sdc%ctx => ctx

    allocate(sdc%mfencap, sdc%encap, sdc%nset_adr, sdc%imex)

    call sdc_multifab_build(sdc%encap, c_loc(sdc%mfencap), err)
    call sdc_nset_build(sdc%nset_adr, nnodes, qtype, 0, err)
    call sdc_imex_build(sdc%imex, sdc%nset_adr, feval, c_null_funptr, c_null_funptr, err)
    call sdc_hooks_add(sdc%imex%hooks, SDC_HOOK_POST_STEP, post, err)

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
    sdc%ctx => ctx

    allocate(sdc%mfencap, sdc%encap, sdc%nset_ad, sdc%nset_r, sdc%mrex)

    call sdc_multifab_build(sdc%encap, c_loc(sdc%mfencap), err)
    call sdc_nset_build(sdc%nset_ad, nnodes(1), qtype, 0, err)
    call sdc_nset_build(sdc%nset_r, nnodes(2), qtype, 0, err)
    call sdc_mrex_build(sdc%mrex, 2, err)
    call sdc_mrex_add_nset(sdc%mrex, sdc%nset_ad, f1eval, sdc%encap, c_loc(ctx), 0, err)

    call sdc_hooks_add(sdc%mrex%hooks, SDC_HOOK_POST_STEP, post, err)

    select case(sdc_multirate_type)
    case ("local")
       call sdc_mrex_add_nset(sdc%mrex, sdc%nset_r, f2eval, sdc%encap, &
            c_loc(ctx), SDC_MR_LOCAL, err)
    case ("global")
       call sdc_mrex_add_nset(sdc%mrex, sdc%nset_r, f2eval, sdc%encap, &
            c_loc(ctx), SDC_MR_GLOBAL, err)
       if (err .ne. 0) then
          call sdc_nset_print(sdc%nset_ad, 2)
          call sdc_nset_print(sdc%nset_r, 2)
          stop "NODES DO NOT NEST PROPERLY"
       end if
    case ("repeated")
       call sdc_mrex_add_nset(sdc%mrex, sdc%nset_r, f2eval, sdc%encap, &
            c_loc(ctx), SDC_MR_REPEATED, err)
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
       call sdc_imex_setup(sdc%imex, sdc%encap, c_loc(sdc%ctx), err)
    else 
       call sdc_mrex_setup(sdc%mrex, err)
       ! call sdc_mrex_print(sdc%mrex, 2)
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
       call sdc_imex_destroy(sdc%imex)
       call sdc_nset_destroy(sdc%nset_adr)
       deallocate(sdc%nset_adr, sdc%imex)
    end if

    if (sdc%multi_rate) then
       call sdc_mrex_destroy(sdc%mrex)
       call sdc_nset_destroy(sdc%nset_ad)
       call sdc_nset_destroy(sdc%nset_r)
       deallocate(sdc%nset_ad, sdc%nset_r, sdc%mrex)
    end if

    call sdc_multifab_destroy(sdc%encap)
    deallocate(sdc%encap, sdc%mfencap)
  end subroutine sdc_destroy

end module sdcquad_module
