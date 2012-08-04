module make_plot_variables_module

  use variables
  use multifab_module

  private

  public make_divu, make_omegadot

contains

  subroutine make_divu(plotdata,icomp,Q,dx)
    type(multifab), intent(inout) :: plotdata
    type(multifab), intent(in   ) :: Q
    integer, intent(in) :: icomp
    double precision, intent(in) :: dx(Q%dim)

    integer :: i, dm, ngpd, ngq
    double precision, pointer :: qp(:,:,:,:), pdp(:,:,:,:)
    integer :: lo(Q%dim), hi(Q%dim)

    dm = Q%dim
    ngq = nghost(Q)
    ngpd = nghost(plotdata)

    do i = 1, nboxes(Q)
       if (multifab_remote(Q, i)) cycle

       qp => dataptr(Q, i)
       pdp => dataptr(plotdata, i)
       
       lo = lwb(get_box(Q, i))
       hi = upb(get_box(Q, i))

       if (dm .ne. 3) then
          call bl_error("Only 3D make_plot_variables_module:make_divu is support")
       else
          call make_divu_3d(pdp(:,:,:,icomp), qp, lo, hi, ngpd, ngq, dx)
       end if
    end do

  end subroutine make_divu

  subroutine make_divu_3d(divu, Q, lo, hi, ngd, ngq, dx)

    use derivative_stencil_module, only : ALP, BET, GAM, DEL

    integer, intent(in) :: lo(3), hi(3), ngd, ngq
    double precision, intent(out) :: divu(-ngd+lo(1):hi(1)+ngd,-ngd+lo(2):hi(2)+ngd,-ngd+lo(3):hi(3)+ngd)
    double precision, intent(in ) ::    Q(-ngq+lo(1):hi(1)+ngq,-ngq+lo(2):hi(2)+ngq,-ngq+lo(3):hi(3)+ngq,nprim)
    double precision, intent(in) :: dx(3)

    integer :: i,j,k
    double precision :: ux, vy, wz
    double precision :: dxinv(3)

    do i=1,3
       dxinv(i) = 1.0d0 / dx(i)
    end do

    do k=lo(3),hi(3)
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)

             ux =  (ALP*(q(i+1,j,k,qu)-q(i-1,j,k,qu)) &
                  + BET*(q(i+2,j,k,qu)-q(i-2,j,k,qu)) &
                  + GAM*(q(i+3,j,k,qu)-q(i-3,j,k,qu)) &
                  + DEL*(q(i+4,j,k,qu)-q(i-4,j,k,qu)))*dxinv(1)

             vy =  (ALP*(q(i,j+1,k,qv)-q(i,j-1,k,qv)) &
                  + BET*(q(i,j+2,k,qv)-q(i,j-2,k,qv)) &
                  + GAM*(q(i,j+3,k,qv)-q(i,j-3,k,qv)) &
                  + DEL*(q(i,j+4,k,qv)-q(i,j-4,k,qv)))*dxinv(2)

             wz =  (ALP*(q(i,j,k+1,qw)-q(i,j,k-1,qw)) &
                  + BET*(q(i,j,k+2,qw)-q(i,j,k-2,qw)) &
                  + GAM*(q(i,j,k+3,qw)-q(i,j,k-3,qw)) &
                  + DEL*(q(i,j,k+4,qw)-q(i,j,k-4,qw)))*dxinv(3)

             divu(i,j,k) = ux + vy + wz

          enddo
       enddo
    enddo

  end subroutine make_divu_3d


  subroutine make_omegadot(plotdata,icomp,Q)
    type(multifab), intent(inout) :: plotdata
    type(multifab), intent(in   ) :: Q
    integer, intent(in) :: icomp
    
    integer :: i, dm, ngpd, ngq
    double precision, pointer :: qp(:,:,:,:), pdp(:,:,:,:)
    integer :: lo(Q%dim), hi(Q%dim)

    dm = Q%dim
    ngq = nghost(Q)
    ngpd = nghost(plotdata)

    do i = 1, nboxes(Q)
       if (multifab_remote(Q, i)) cycle

       qp => dataptr(Q, i)
       pdp => dataptr(plotdata, i)
       
       lo = lwb(get_box(Q, i))
       hi = upb(get_box(Q, i))

       if (dm .ne. 3) then
          call bl_error("Only 3D make_plot_variables_module:make_omegadot is support")
       else
          call make_omegadot_3d(pdp(:,:,:,icomp:icomp+nspecies-1), qp, lo, hi, ngpd, ngq)
       end if
    end do

  end subroutine make_omegadot

  subroutine make_omegadot_3d(odot, Q, lo, hi, ngo, ngq)
    use chemistry_module, only : molecular_weight

    integer, intent(in) :: lo(3), hi(3), ngo, ngq
    double precision, intent(out) :: odot(-ngo+lo(1):hi(1)+ngo,-ngo+lo(2):hi(2)+ngo,-ngo+lo(3):hi(3)+ngo,nspecies)
    double precision, intent(in ) ::    Q(-ngq+lo(1):hi(1)+ngq,-ngq+lo(2):hi(2)+ngq,-ngq+lo(3):hi(3)+ngq,nprim)

    integer :: i,j,k, iwrk
    double precision :: rwrk, wdot(nspecies), Xt(nspecies)

    do k=lo(3),hi(3)
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)

             Xt = q(i,j,k,qx1:qx1+nspecies-1)
             call ckwxr(q(i,j,k,qrho), q(i,j,k,qtemp), Xt, iwrk, rwrk, wdot)
             odot(i,j,k,:) = wdot * molecular_weight

          enddo
       enddo
    enddo

  end subroutine make_omegadot_3d

end module make_plot_variables_module
