module RNS_boundary_module

  implicit none

contains

  subroutine get_lo_bc_vfac(dir,idx,dx,vfac)

    use prob_params_module, only : physbc_lo,Symmetry,SlipWall,NoSlipWall

    integer, intent(in) :: dir, idx(:)
    double precision, intent(in) :: dx(:)
    double precision, intent(out) :: vfac

    vfac = 1.d0

  end subroutine get_lo_bc_vfac

end module RNS_boundary_module

