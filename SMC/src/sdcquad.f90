module sdcquad_module

  use bl_types
  use sdclib
  use sdclib_multifab

  implicit none

  type :: ctx_t
     double precision :: dx(3)
  end type ctx_t

  type :: sdc_t

     integer :: &
          iters   = 4                       ! Number of SDC iterations
          
     real(dp_t) :: &
          tol_residual = -1.0d0             ! Residual tolerance (ignored if negative)
     
     type(c_ptr) :: nset1, nset2, exp1, exp2
     type(c_ptr) :: srset, mrset

     type(ctx_t), pointer      :: ctx
     type(mf_encap_t), pointer :: mfencap
     type(c_ptr)               :: encap
     
  end type sdc_t

contains


  !
  ! Build/create single-rate SDC object.
  !
  subroutine sdc_build_single_rate(sdc, qtype, nnodes)
    type(sdc_t),   intent(out) :: sdc
    integer,       intent(in)  :: qtype, nnodes

    sdc%nset1 = sdc_nset_create(nnodes, qtype, "ADR" // c_null_char)
    sdc%srset = sdc_srset_create(sdc%nset1, "ADR" // c_null_char)
    sdc%exp1  = sdc_exp_create("ADR" // c_null_char)

    sdc%mrset = c_null_ptr
    sdc%nset2 = c_null_ptr
    sdc%exp2  = c_null_ptr

  end subroutine sdc_build_single_rate


  !
  ! Build/create multi-rate SDC object.
  !
  subroutine sdc_build_multi_rate(sdc, qtype, nnodes)
    type(sdc_t),   intent(out) :: sdc
    integer,       intent(in)  :: qtype, nnodes(2)

    integer :: err

    sdc%mrset = sdc_mrset_create("ADR" // c_null_char)

    sdc%nset1 = sdc_nset_create(nnodes(1), qtype, "AD" // c_null_char)
    sdc%exp1  = sdc_exp_create("AD" // c_null_char)

    sdc%nset2 = sdc_nset_create(nnodes(2), qtype, "R" // c_null_char)
    sdc%exp1  = sdc_exp_create("R" // c_null_char)

    err = sdc_mrset_add_nset(sdc%mrset, sdc%nset1, 0)
    err = sdc_mrset_add_nset(sdc%mrset, sdc%nset2, 0)

    sdc%srset = c_null_ptr

  end subroutine sdc_build_multi_rate

  
  !
  ! Set multifab layout, feval context etc
  !
  subroutine sdc_set_layout(sdc, la, nc, ng)
    type(sdc_t),  intent(inout) :: sdc
    type(layout), intent(in)    :: la
    integer,      intent(in)    :: nc, ng

    ! create encapsulation
    allocate(sdc%mfencap)
    
    sdc%mfencap%nc = nc
    sdc%mfencap%ng = ng
    sdc%mfencap%la = la

    sdc%encap = sdc_encap_multifab(c_loc(sdc%mfencap))
  end subroutine sdc_set_layout

  subroutine sdc_set_context(sdc, ctx)
    type(sdc_t), intent(inout) :: sdc
    type(ctx_t), intent(in), target :: ctx

    sdc%ctx => ctx
  end subroutine sdc_set_context


  !
  ! Setup and allocate SDC object.
  !
  subroutine sdc_setup(sdc)
    type(sdc_t), intent(inout) :: sdc

    integer :: err

    if (c_associated(sdc%srset)) then
       call sdc_exp_attach(sdc%exp1, sdc_srset_stepper(sdc%srset), sdc%encap, c_loc(sdc%ctx))
       err = sdc_srset_setup(sdc%srset)
    end if

    if (c_associated(sdc%mrset)) then
       call sdc_exp_attach(sdc%exp1, sdc_mrset_stepper(sdc%mrset, 0), sdc%encap, c_loc(sdc%ctx))
       call sdc_exp_attach(sdc%exp2, sdc_mrset_stepper(sdc%mrset, 1), sdc%encap, c_loc(sdc%ctx))
       err = sdc_mrset_setup(sdc%srset)
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

    call sdc_srset_destroy(sdc%srset)
    call sdc_mrset_destroy(sdc%mrset)
    call sdc_nset_destroy(sdc%nset1)
    call sdc_nset_destroy(sdc%nset2)
  end subroutine sdc_destroy

end module sdcquad_module
