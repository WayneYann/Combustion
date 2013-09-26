module sdcquad_module

  use bl_types
  use sdclib
  use sdclib_multifab

  implicit none

  type :: sdc_ctx
     real(dp_t) :: dx(3)

     logical    :: single_rate, multi_rate
     integer    :: iters = 4
     real(dp_t) :: tol_residual = -1.0d0

     type(mf_encap),  pointer :: mfencap
     type(sdc_encap), pointer :: encap
     type(sdc_nodes), pointer :: nodes1, nodes2
     type(sdc_mrex),  pointer :: mrex
     type(sdc_imex),  pointer :: imex
  end type sdc_ctx

contains


  !
  ! Build/create single-rate SDC object.
  !
  subroutine sdc_build_single_rate(sdc, qtype, nnodes, feval, post)
    type(sdc_ctx), intent(out), target :: sdc
    integer,         intent(in)          :: qtype, nnodes
    type(c_funptr),  intent(in), value   :: feval, post

    integer :: err

    sdc%single_rate = .true.
    sdc%multi_rate  = .false.

    allocate(sdc%mfencap, sdc%encap, sdc%nodes1, sdc%imex)

    call sdc_multifab_build(sdc%encap, c_loc(sdc%mfencap), err)
    call sdc_nodes_build(sdc%nodes1, nnodes, qtype, err)
    call sdc_imex_build(sdc%imex, sdc%nodes1, feval, c_null_funptr, c_null_funptr, err)
    call sdc_hooks_add(sdc%imex%hooks, SDC_HOOK_POST_STEP, post, err)

  end subroutine sdc_build_single_rate


  !
  ! Build/create multi-rate SDC object.
  !
  subroutine sdc_build_multi_rate(sdc, qtype, nnodes, f1eval, f2eval, post)
    use probin_module, only: sdc_multirate_type, sdc_multirate_repeat, advance_method

    type(sdc_ctx), intent(out), target :: sdc
    integer,         intent(in)          :: qtype, nnodes(2)
    type(c_funptr),  intent(in), value   :: f1eval, f2eval, post
    
    integer :: err

    sdc%single_rate = .false.
    sdc%multi_rate  = .true.

    allocate(sdc%mfencap, sdc%encap, sdc%nodes1, sdc%nodes2, sdc%mrex)

    call sdc_multifab_build(sdc%encap, c_loc(sdc%mfencap), err)
    call sdc_nodes_build(sdc%nodes1, nnodes(1), qtype, err)
    call sdc_nodes_build(sdc%nodes2, nnodes(2), qtype, err)
    call sdc_mrex_build(sdc%mrex, 2, err)
    call sdc_mrex_add_comp(sdc%mrex, sdc%nodes1, f1eval, sdc%encap, c_loc(sdc), 1, SDC_MR_ROOT, err)

    call sdc_hooks_add(sdc%mrex%hooks, SDC_HOOK_POST_STEP, post, err)

    select case(sdc_multirate_type)
    case ("local")
       if (sdc_multirate_repeat > 1) then
          call sdc_mrex_add_comp(sdc%mrex, sdc%nodes2, f2eval, sdc%encap, &
               c_loc(sdc), sdc_multirate_repeat, SDC_MR_LREP, err)
       else
          call sdc_mrex_add_comp(sdc%mrex, sdc%nodes2, f2eval, sdc%encap, &
               c_loc(sdc), 1, SDC_MR_LOCAL, err)
       end if
    case ("global")
       call sdc_mrex_add_comp(sdc%mrex, sdc%nodes2, f2eval, sdc%encap, &
            c_loc(sdc), 1, SDC_MR_GLOBAL, err)
       if (err .ne. 0) then
          call sdc_nodes_print(sdc%nodes1, 2)
          call sdc_nodes_print(sdc%nodes2, 2)
          stop "NODES DO NOT NEST PROPERLY"
       end if
    case ("repeated")
       call sdc_mrex_add_comp(sdc%mrex, sdc%nodes2, f2eval, sdc%encap, &
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

  function sdc_get_chemterm(ctx, node_advdif) result(r)
    use probin_module, only : advance_method

    type(sdc_ctx), intent(in) :: ctx
    integer,         intent(in) :: node_advdif
    type(multifab),  pointer    :: r

    integer :: trat, node_chem

    type(sdc_nset), pointer :: nset1, nset2
    type(c_ptr),      pointer :: p(:)

    call c_f_pointer(ctx%mrex%nsets, p, [ ctx%mrex%ncomps ])
    call c_f_pointer(p(1), nset1)
    call c_f_pointer(p(2), nset2)

    trat = (nset2%nnodes - 1) / (nset1%nnodes - 1)

    if (advance_method == 3) then
       node_chem = node_advdif * trat
       call c_f_pointer(nset2%F, p, [ nset2%nnodes ])
       call c_f_pointer(p(node_chem+1), r)
    else
       node_chem = node_advdif / trat
       call c_f_pointer(nset1%F, p, [ nset1%nnodes ])
       call c_f_pointer(p(node_chem+1), r)
    end if

  end function sdc_get_chemterm


  !
  ! Set multifab layout, feval context etc
  !
  subroutine sdc_set_layout(sdc, la, nc, ng)
    type(sdc_ctx),  intent(inout) :: sdc
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
    type(sdc_ctx), intent(inout), target :: sdc
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
    type(sdc_ctx), intent(inout) :: sdc

    if (sdc%single_rate) then
       call sdc_imex_destroy(sdc%imex)
       call sdc_nodes_destroy(sdc%nodes1)
       deallocate(sdc%nodes1, sdc%imex)
    end if

    if (sdc%multi_rate) then
       call sdc_mrex_destroy(sdc%mrex)
       call sdc_nodes_destroy(sdc%nodes1)
       call sdc_nodes_destroy(sdc%nodes2)
       deallocate(sdc%nodes1, sdc%nodes2, sdc%mrex)
    end if

    call sdc_multifab_destroy(sdc%encap)
    deallocate(sdc%encap, sdc%mfencap)
  end subroutine sdc_destroy

end module sdcquad_module
