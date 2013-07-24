module sdcquad_module

  use bl_types
  use sdclib
  use sdclib_multifab

  implicit none

  type :: sdc_ctx_t
     real(dp_t) :: dx(3)

     logical    :: single_rate, multi_rate
     integer    :: iters = 4
     real(dp_t) :: tol_residual = -1.0d0

     type(mf_encap_t),  pointer :: mfencap
     type(sdc_encap_t), pointer :: encap
     type(sdc_nset_t),  pointer :: nset1, nset2
     type(sdc_mrex_t),  pointer :: mrex
     type(sdc_imex_t),  pointer :: imex
  end type sdc_ctx_t

contains


  !
  ! Build/create single-rate SDC object.
  !
  subroutine sdc_build_single_rate(sdc, qtype, nnodes, feval, post)
    type(sdc_ctx_t), intent(out), target :: sdc
    integer,         intent(in)          :: qtype, nnodes
    type(c_funptr),  intent(in), value   :: feval, post

    integer :: err

    sdc%single_rate = .true.
    sdc%multi_rate  = .false.

    allocate(sdc%mfencap, sdc%encap, sdc%nset1, sdc%imex)

    call sdc_multifab_build(sdc%encap, c_loc(sdc%mfencap), err)
    call sdc_nset_build(sdc%nset1, nnodes, qtype, 0, err)
    call sdc_imex_build(sdc%imex, sdc%nset1, feval, c_null_funptr, c_null_funptr, err)
    call sdc_hooks_add(sdc%imex%hooks, SDC_HOOK_POST_STEP, post, err)

  end subroutine sdc_build_single_rate


  !
  ! Build/create multi-rate SDC object.
  !
  subroutine sdc_build_multi_rate(sdc, qtype, nnodes, f1eval, f2eval, post)
    use probin_module, only: sdc_multirate_type, sdc_multirate_repeat, advance_method

    type(sdc_ctx_t), intent(out), target :: sdc
    integer,         intent(in)          :: qtype, nnodes(2)
    type(c_funptr),  intent(in), value   :: f1eval, f2eval, post
    
    integer :: err

    sdc%single_rate = .false.
    sdc%multi_rate  = .true.

    allocate(sdc%mfencap, sdc%encap, sdc%nset1, sdc%nset2, sdc%mrex)

    call sdc_multifab_build(sdc%encap, c_loc(sdc%mfencap), err)
    call sdc_nset_build(sdc%nset1, nnodes(1), qtype, 0, err)
    call sdc_nset_build(sdc%nset2, nnodes(2), qtype, 0, err)
    call sdc_mrex_build(sdc%mrex, 2, err)
    call sdc_mrex_add_nset(sdc%mrex, sdc%nset1, f1eval, sdc%encap, c_loc(sdc), 1, SDC_MR_ROOT, err)

    call sdc_hooks_add(sdc%mrex%hooks, SDC_HOOK_POST_STEP, post, err)

    select case(sdc_multirate_type)
    case ("local")
       if (sdc_multirate_repeat > 1) then
          call sdc_mrex_add_nset(sdc%mrex, sdc%nset2, f2eval, sdc%encap, &
               c_loc(sdc), sdc_multirate_repeat, SDC_MR_LREP, err)
       else
          call sdc_mrex_add_nset(sdc%mrex, sdc%nset2, f2eval, sdc%encap, &
               c_loc(sdc), 1, SDC_MR_LOCAL, err)
       end if
    case ("global")
       call sdc_mrex_add_nset(sdc%mrex, sdc%nset2, f2eval, sdc%encap, &
            c_loc(sdc), 1, SDC_MR_GLOBAL, err)
       if (err .ne. 0) then
          call sdc_nset_print(sdc%nset1, 2)
          call sdc_nset_print(sdc%nset2, 2)
          stop "NODES DO NOT NEST PROPERLY"
       end if
    case ("repeated")
       call sdc_mrex_add_nset(sdc%mrex, sdc%nset2, f2eval, sdc%encap, &
            c_loc(sdc), 1, SDC_MR_GREP, err)
    case default
       stop "UNKNOWN MULTIRATE TYPE: should be one of 'local', 'global', or 'repeated'"
    end select

    ! this controls the order in which the different multi-rate
    ! components are evaluated.  we always want chemistry to be
    ! evaluated first (so that up-to-date chemistry can be used when
    ! computing the boundary conditions for the adv/diff term)
    if (advance_method == 3) then
       ! chemistry is on the "fast" component
       sdc%mrex%order = 1
    else
       ! chemistry is on the "slow" component
       sdc%mrex%order = -1
    end if

  end subroutine sdc_build_multi_rate


  !
  ! Set multifab layout, feval context etc
  !
  subroutine sdc_set_layout(sdc, la, nc, ng)
    type(sdc_ctx_t),  intent(inout) :: sdc
    type(layout),     intent(in)    :: la
    integer,          intent(in)    :: nc, ng

    sdc%mfencap%nc = nc
    sdc%mfencap%ng = ng
    sdc%mfencap%la = la
  end subroutine sdc_set_layout


  !
  ! Setup and allocate SDC object.
  !
  subroutine sdc_setup(sdc, la, nc, ng)
    type(sdc_ctx_t), intent(inout), target :: sdc
    type(layout),    intent(in)            :: la
    integer,         intent(in)            :: nc, ng

    integer :: err

    call sdc_set_layout(sdc, la, nc, ng)

    if (sdc%single_rate) then
       call sdc_imex_setup(sdc%imex, sdc%encap, c_loc(sdc), err)
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
    type(sdc_ctx_t), intent(inout) :: sdc

    if (sdc%single_rate) then
       call sdc_imex_destroy(sdc%imex)
       call sdc_nset_destroy(sdc%nset1)
       deallocate(sdc%nset1, sdc%imex)
    end if

    if (sdc%multi_rate) then
       call sdc_mrex_destroy(sdc%mrex)
       call sdc_nset_destroy(sdc%nset1)
       call sdc_nset_destroy(sdc%nset2)
       deallocate(sdc%nset1, sdc%nset2, sdc%mrex)
    end if

    call sdc_multifab_destroy(sdc%encap)
    deallocate(sdc%encap, sdc%mfencap)
  end subroutine sdc_destroy

end module sdcquad_module
