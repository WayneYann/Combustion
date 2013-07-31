module make_plot_variables_module

  use multifab_module

  use smc_bc_module
  use variables_module
  use plotvar_index_module
  use make_plotvar_2d_module
  use make_plotvar_3d_module

  implicit none

  private

  public :: make_plotvar

contains

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
             if (icomp .eq. icomp_wbar) then
                call make_wbar_2d(lo,hi,pdp(:,:,1,icomp),pdlo(1:2),pdhi(1:2), &
                     qp,qlo(1:2),qhi(1:2))
             else if (icomp .eq. icomp_h) then
                call make_h_2d(lo,hi,pdp(:,:,1,icomp),pdlo(1:2),pdhi(1:2), &
                     qp,qlo(1:2),qhi(1:2))
             else if (icomp .eq. icomp_rhoh) then
                call make_rhoh_2d(lo,hi,pdp(:,:,1,icomp),pdlo(1:2),pdhi(1:2), &
                     qp,qlo(1:2),qhi(1:2))
             else if (icomp .eq. icomp_divu) then
                call make_divu_2d(lo,hi,pdp(:,:,1,icomp),pdlo(1:2),pdhi(1:2), &
                     qp,qlo(1:2),qhi(1:2), dx, dlo,dhi)
             else if (icomp .eq. icomp_magvort) then
                call make_magvort_2d(lo,hi,pdp(:,:,1,icomp),pdlo(1:2),pdhi(1:2), &
                     qp,qlo(1:2),qhi(1:2), dx, dlo,dhi)
             else if (icomp .eq. icomp_burn) then
                call make_burn_2d(lo,hi,pdp(:,:,1,icomp:icomp+nburn-1),pdlo(1:2),pdhi(1:2), &
                     qp,qlo(1:2),qhi(1:2))
             else
                call bl_error("make_plot_variables_module: unknown icomp")          
             end if
          else
             if (icomp .eq. icomp_wbar) then
                call make_wbar_3d(lo,hi,pdp(:,:,:,icomp),pdlo(1:3),pdhi(1:3), &
                     qp,qlo(1:3),qhi(1:3))
             else if (icomp .eq. icomp_h) then
                call make_h_3d(lo,hi,pdp(:,:,:,icomp),pdlo(1:3),pdhi(1:3), &
                     qp,qlo(1:3),qhi(1:3))
             else if (icomp .eq. icomp_rhoh) then
                call make_rhoh_3d(lo,hi,pdp(:,:,:,icomp),pdlo(1:3),pdhi(1:3), &
                     qp,qlo(1:3),qhi(1:3))
             else if (icomp .eq. icomp_divu) then
                call make_divu_3d(lo,hi,pdp(:,:,:,icomp),pdlo(1:3),pdhi(1:3), &
                     qp,qlo(1:3),qhi(1:3), dx, dlo,dhi)
             else if (icomp .eq. icomp_magvort) then
                call make_magvort_3d(lo,hi,pdp(:,:,:,icomp),pdlo(1:3),pdhi(1:3), &
                     qp,qlo(1:3),qhi(1:3), dx, dlo,dhi)
             else if (icomp .eq. icomp_burn) then
                call make_burn_3d(lo,hi,pdp(:,:,:,icomp:icomp+nburn-1),pdlo(1:3),pdhi(1:3), &
                     qp,qlo(1:3),qhi(1:3))
             else
                call bl_error("make_plot_variables_module: unknown icomp")          
             end if
          end if
       end do

    end do
    !$omp end parallel

  end subroutine make_plotvar

end module make_plot_variables_module
