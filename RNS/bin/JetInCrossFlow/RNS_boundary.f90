module RNS_boundary_module

  implicit none

contains

  subroutine get_lo_bc_vfac(dir,idx,dx,vfac)

    use prob_params_module, only : phys_prob_lo
    use probdata_module

    integer, intent(in) :: dir, idx(:)
    double precision, intent(in) :: dx(:)
    double precision, intent(out) :: vfac

    double precision :: x, z

    vfac = 1.d0

    if (dir .eq. 1 .and. idx(1) .eq. 0) then
       x = phys_prob_lo(1) + (idx(1)+0.5d0)*dx(1)
       z = phys_prob_lo(3) + (idx(3)+0.5d0)*dx(3)

       if (x**2+z**2 .gt. r_jet**2) then
          vfac = 0.d0
       end if
    end if

  end subroutine get_lo_bc_vfac

end module RNS_boundary_module

