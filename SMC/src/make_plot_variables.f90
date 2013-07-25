module make_plot_variables_module

  use multifab_module

  use smc_bc_module
  use variables_module
  use make_plotvar_2d_module
  use make_plotvar_3d_module

  implicit none

  integer, save :: icomp_h, icomp_divu, icomp_magvort, icomp_omegadot, &
       icomp_dYdt, icomp_heatRelease, icomp_fuelConsumption, ifuel

  private

  public :: make_plotvar_init, make_plotvar

contains

  subroutine make_plotvar_init(icomp_h_in, icomp_divu_in, icomp_magvort_in, icomp_omegadot_in, &
       icomp_dYdt_in, icomp_heatRelease_in, icomp_fuelConsumption_in, ifuel_in)
    integer, intent(in) :: icomp_h_in, icomp_divu_in, icomp_magvort_in, icomp_omegadot_in, &
         icomp_dYdt_in, icomp_heatRelease_in, icomp_fuelConsumption_in, ifuel_in
    icomp_h               = icomp_h_in
    icomp_divu            = icomp_divu_in 
    icomp_magvort         = icomp_magvort_in
    icomp_omegadot        = icomp_omegadot_in
    icomp_dYdt            = icomp_dYdt_in
    icomp_heatRelease     = icomp_heatRelease_in
    icomp_fuelConsumption = icomp_fuelConsumption_in
    ifuel                 = ifuel_in
  end subroutine make_plotvar_init

  subroutine make_plotvar(plotdata,icomp,Q,dx)
    use threadbox_module
    type(multifab), intent(inout) :: plotdata
    type(multifab), intent(in   ) :: Q
    integer, intent(in) :: icomp
    double precision, intent(in) :: dx(3)

    integer :: i, dm, iblock
    double precision, pointer :: qp(:,:,:,:), pdp(:,:,:,:)
    integer :: qlo(4), qhi(4), pdlo(4), pdhi(4)
    integer ::  lo(Q%dim),  hi(Q%dim)
    integer :: dlo(Q%dim), dhi(Q%dim)

    dm = Q%dim

    !$omp parallel private(i,iblock,lo,hi,qp,pdp,qlo,qhi,pdlo,pdhi,dlo,dhi)
    do i = 1, nfabs(Q)

       if (.not.tb_worktodo(i)) cycle

       qp => dataptr(Q, i)
       pdp => dataptr(plotdata, i)

       qlo = lbound(qp)
       qhi = ubound(qp)
       pdlo = lbound(pdp)
       pdhi = ubound(pdp)

       call get_data_lo_hi(i,dlo,dhi)
       
       do iblock = 1, tb_get_nblocks(i)
          lo = tb_get_block_lo(iblock,i)
          hi = tb_get_block_hi(iblock,i)

          if (dm .eq. 2) then
             if (icomp .eq. icomp_h) then
                call make_h_2d(lo,hi,pdp(:,:,1,icomp),pdlo(1:2),pdhi(1:2), &
                     qp,qlo(1:2),qhi(1:2))
             else if (icomp .eq. icomp_divu) then
                call make_divu_2d(lo,hi,pdp(:,:,1,icomp),pdlo(1:2),pdhi(1:2), &
                     qp,qlo(1:2),qhi(1:2), dx, dlo,dhi)
             else if (icomp .eq. icomp_magvort) then
                call make_magvort_2d(lo,hi,pdp(:,:,1,icomp),pdlo(1:2),pdhi(1:2), &
                     qp,qlo(1:2),qhi(1:2), dx, dlo,dhi)
             else if (icomp .eq. icomp_omegadot) then
                call make_omegadot_2d(lo,hi,pdp(:,:,1,icomp:icomp+nspecies-1),pdlo(1:2),pdhi(1:2), &
                     qp,qlo(1:2),qhi(1:2))
             else if (icomp .eq. icomp_dYdt) then
                call make_dYdt_2d(lo,hi,pdp(:,:,1,icomp:icomp+nspecies-1),pdlo(1:2),pdhi(1:2), &
                     qp,qlo(1:2),qhi(1:2))
             else if (icomp .eq. icomp_heatRelease) then
                call make_heatRelease_2d(lo,hi,pdp(:,:,1,icomp),pdlo(1:2),pdhi(1:2), &
                     qp,qlo(1:2),qhi(1:2))
             else if (icomp .eq. icomp_fuelConsumption) then
                call make_fuelCsmp_2d(lo,hi,pdp(:,:,1,icomp),pdlo(1:2),pdhi(1:2), &
                     qp,qlo(1:2),qhi(1:2), ifuel)
             else
                call bl_error("make_plot_variables_module: unknown icomp")          
             end if
          else
             if (icomp .eq. icomp_h) then
                call make_h_3d(lo,hi,pdp(:,:,:,icomp),pdlo(1:3),pdhi(1:3), &
                     qp,qlo(1:3),qhi(1:3))
             else if (icomp .eq. icomp_divu) then
                call make_divu_3d(lo,hi,pdp(:,:,:,icomp),pdlo(1:3),pdhi(1:3), &
                     qp,qlo(1:3),qhi(1:3), dx, dlo,dhi)
             else if (icomp .eq. icomp_magvort) then
                call make_magvort_3d(lo,hi,pdp(:,:,:,icomp),pdlo(1:3),pdhi(1:3), &
                     qp,qlo(1:3),qhi(1:3), dx, dlo,dhi)
             else if (icomp .eq. icomp_omegadot) then
                call make_omegadot_3d(lo,hi,pdp(:,:,:,icomp:icomp+nspecies-1),pdlo(1:3),pdhi(1:3), &
                     qp,qlo(1:3),qhi(1:3))
             else if (icomp .eq. icomp_dYdt) then
                call make_dYdt_3d(lo,hi,pdp(:,:,:,icomp:icomp+nspecies-1),pdlo(1:3),pdhi(1:3), &
                     qp,qlo(1:3),qhi(1:3))
             else if (icomp .eq. icomp_heatRelease) then
                call make_heatRelease_3d(lo,hi,pdp(:,:,:,icomp),pdlo(1:3),pdhi(1:3), &
                     qp,qlo(1:3),qhi(1:3))
             else if (icomp .eq. icomp_fuelConsumption) then
                call make_fuelCsmp_3d(lo,hi,pdp(:,:,:,icomp),pdlo(1:3),pdhi(1:3), &
                     qp,qlo(1:3),qhi(1:3), ifuel)
             else
                call bl_error("make_plot_variables_module: unknown icomp")          
             end if
          end if
       end do

    end do
    !$omp end parallel

  end subroutine make_plotvar

end module make_plot_variables_module
