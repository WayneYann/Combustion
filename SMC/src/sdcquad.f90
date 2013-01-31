module sdcquad_module

  use bl_types
  use sdclib
  use sdclib_multifab

  implicit none

  type :: ctx_t
     double precision :: dx(3)
  end type ctx_t

  type :: sdcquad

     integer :: &
          qtype   = SDC_GAUSS_LOBATTO, &    ! SDC quadrature type
          nnodes  = 3, &                    ! Number of SDC nodes
          iters   = 4                       ! Number of SDC iterations

     real(dp_t) :: &
          tol_residual = -1.0d0  ! Residual tolerance (ignored if negative)

     type(c_ptr) :: nset1, nset2, exp1, exp2
     type(c_ptr) :: srset, mrset

     type(mf_encap_t), pointer :: mfencap
     type(c_ptr)               :: encap
     
  end type sdcquad

  interface build
     module procedure sdcquad_build
  end interface

  interface create
     module procedure sdcquad_create
  end interface create

  interface destroy
     module procedure sdcquad_destroy
  end interface

contains


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  ! Create (but do not build) SDC object.
  !
  subroutine sdcquad_create(sdc, qtype, nnodes, la, ng, nc, ctx)
    integer,       intent(in ) :: qtype, nnodes
    type(sdcquad), intent(out) :: sdc
    type(layout),  intent(in ) :: la
    integer,       intent(in ) :: ng, nc
    type(ctx_t),   intent(in ), target :: ctx

    integer :: err

    sdc%nnodes = nnodes
    sdc%qtype  = qtype
    sdc%iters  = 2*nnodes - 1

    ! create encapsulation
    allocate(sdc%mfencap)
    
    sdc%mfencap%nc = nc
    sdc%mfencap%ng = ng
    sdc%mfencap%la = la

    sdc%encap = sdc_encap_multifab(c_loc(sdc%mfencap))

    ! setup one node-set for single rate sdc
    sdc%nset1 = sdc_nset_create(sdc%nnodes, sdc%qtype, "ADR" // c_null_char)
    sdc%srset = sdc_srset_create(sdc%nset1, "ADR" // c_null_char)
    sdc%exp1  = sdc_exp_create("ADR" // c_null_char)

    ! attach explicit steppers
    call sdc_exp_attach(sdc%exp1, &
         sdc_srset_stepper(sdc%srset, 0), sdc%encap, c_loc(ctx))

  end subroutine sdcquad_create


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  ! Build SDC object.
  !
  subroutine sdcquad_build(sdc)
    type(sdcquad), intent(inout) :: sdc

    integer :: err

    err = sdc_srset_setup(sdc%srset)

    if (err .ne. 0) then
       stop "FAILED TO BULID SDC SETS"
    end if
  end subroutine sdcquad_build


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  ! Destroy SDC object.
  !
  subroutine sdcquad_destroy(sdc)
    type(sdcquad), intent(inout) :: sdc

    call sdc_srset_destroy(sdc%srset)
    call sdc_nset_destroy(sdc%nset1)

  end subroutine sdcquad_destroy

end module sdcquad_module
