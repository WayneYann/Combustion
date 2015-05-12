module RNS_boundary_module

  implicit none

contains

  subroutine get_lo_bc_vfac(dir,idx,dx,vfac)

    integer, intent(in) :: dir, idx(:)
    double precision, intent(in) :: dx(:)
    double precision, intent(out) :: vfac

    vfac = 1.d0

    if (dir .eq. 1 .and. idx(1) .eq. 0) then
       vfac = 0.d0
    end if

  end subroutine get_lo_bc_vfac

end module RNS_boundary_module

