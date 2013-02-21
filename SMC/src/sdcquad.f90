module sdcquad_module

  use bl_types
  use sdclib
  use sdclib_multifab

  implicit none

  type :: ctx_t
     double precision :: dx(3)
  end type ctx_t

  type :: sdc_t

     logical :: single_rate, multi_rate

     integer :: &
          iters   = 4                       ! Number of SDC iterations

     real(dp_t) :: &
          tol_residual = -1.0d0             ! Residual tolerance (ignored if negative)

     type(c_ptr) :: nset_adr, exp_adr, srset
     type(c_ptr) :: nset_ad, nset_r, exp_ad, exp_r, mrset

     type(mf_encap_t), pointer :: mfencap
     type(c_ptr)               :: encap

  end type sdc_t

contains


  !
  ! Build/create single-rate SDC object.
  !
  subroutine sdc_build_single_rate(sdc, qtype, nnodes, ctx, feval)
    type(sdc_t),    intent(out) :: sdc
    integer,        intent(in)  :: qtype, nnodes
    type(ctx_t),    intent(in), target :: ctx
    type(c_funptr), intent(in), value :: feval

    sdc%single_rate = .true.
    sdc%multi_rate  = .false.

    allocate(sdc%mfencap)
    sdc%encap = sdc_encap_multifab(c_loc(sdc%mfencap))

    sdc%nset_adr = sdc_nset_create(nnodes, qtype, "ADR" // c_null_char)
    sdc%exp_adr  = sdc_exp_create(feval, "ADR" // c_null_char)
    sdc%srset    = sdc_srset_create(sdc%nset_adr, sdc%exp_adr, sdc%encap, &
         c_loc(ctx), "ADR" // c_null_char)

  end subroutine sdc_build_single_rate


  !
  ! Build/create multi-rate SDC object.
  !
  subroutine sdc_build_multi_rate(sdc, qtype, nnodes, ctx, f1eval, f2eval)
    type(sdc_t),    intent(out) :: sdc
    integer,        intent(in)  :: qtype, nnodes(2)
    type(ctx_t),    intent(in), target :: ctx
    type(c_funptr), intent(in), value :: f1eval, f2eval

    integer :: err

    sdc%single_rate = .false.
    sdc%multi_rate  = .true.

    allocate(sdc%mfencap)
    sdc%encap = sdc_encap_multifab(c_loc(sdc%mfencap))

    sdc%nset_ad = sdc_nset_create(nnodes(1), qtype, "AD" // c_null_char)
    sdc%nset_r  = sdc_nset_create(nnodes(2), qtype, "R" // c_null_char)
    sdc%exp_ad  = sdc_exp_create(f1eval, "AD" // c_null_char)
    sdc%exp_r   = sdc_exp_create(f2eval, "R" // c_null_char)

    sdc%mrset = sdc_mrset_create("ADR" // c_null_char)
    err = sdc_mrset_add_nset(sdc%mrset, sdc%nset_ad, sdc%exp_ad, sdc%encap, c_loc(ctx), 0)
    err = sdc_mrset_add_nset(sdc%mrset, sdc%nset_r,  sdc%exp_r,  sdc%encap, c_loc(ctx), 0)

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
  subroutine sdc_setup(sdc)
    type(sdc_t), intent(inout) :: sdc

    integer :: err
    err = 0

    if (sdc%single_rate) then
       err = sdc_srset_setup(sdc%srset)
    else 
       err = sdc_mrset_setup(sdc%mrset)
       ! call sdc_mrset_print(sdc%mrset, 0)
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
    end if

    if (sdc%multi_rate) then
       call sdc_mrset_destroy(sdc%mrset)
       call sdc_nset_destroy(sdc%nset_ad)
       call sdc_nset_destroy(sdc%nset_r)
       call sdc_exp_destroy(sdc%exp_ad)
       call sdc_exp_destroy(sdc%exp_r)
    end if

    call sdc_encap_multifab_destroy(sdc%encap)
    deallocate(sdc%mfencap)
  end subroutine sdc_destroy

end module sdcquad_module
