module make_plot_variables_module

  use multifab_module
  use variables_module

  implicit none

  integer, save :: icomp_h, icomp_divu, icomp_omegadot

  private

  public :: make_plotvar_init, make_plotvar

contains

  subroutine make_plotvar_init(icomp_h_in, icomp_divu_in, icomp_omegadot_in)
    integer, intent(in) :: icomp_h_in, icomp_divu_in, icomp_omegadot_in
    icomp_h        = icomp_h_in
    icomp_divu     = icomp_divu_in 
    icomp_omegadot = icomp_omegadot_in
  end subroutine make_plotvar_init

  subroutine make_plotvar(plotdata,icomp,Q,dx)
    type(multifab), intent(inout) :: plotdata
    type(multifab), intent(in   ) :: Q
    integer, intent(in) :: icomp
    double precision, intent(in) :: dx(Q%dim)

    integer :: i, dm, ngpd, ngq
    double precision, pointer :: qp(:,:,:,:), pdp(:,:,:,:)
    integer ::  lo(Q%dim),  hi(Q%dim)

    dm = Q%dim
    ngq = nghost(Q)
    ngpd = nghost(plotdata)

    do i = 1, nfabs(Q)
       qp => dataptr(Q, i)
       pdp => dataptr(plotdata, i)
       
       lo = lwb(get_box(Q, i))
       hi = upb(get_box(Q, i))

       if (dm .ne. 3) then
          call bl_error("Only 3D make_plotvar is support")
       end if

       if (icomp .eq. icomp_h) then
          call make_h_3d(pdp(:,:,:,icomp), qp, lo, hi, ngpd, ngq)                
       else if (icomp .eq. icomp_divu) then
          call make_divu_3d(pdp(:,:,:,icomp), qp, lo, hi, ngpd, ngq, dx)
       else if (icomp .eq. icomp_omegadot) then
          call make_omegadot_3d(pdp(:,:,:,icomp:nspecies-1), qp, lo, hi, ngpd, ngq)
       else
          call bl_error("make_plot_variables_module: unknown icomp")          
       end if
    end do

  end subroutine make_plotvar


  subroutine make_h_3d(h, Q, lo, hi, ngh, ngq)
    integer, intent(in) :: lo(3), hi(3), ngh, ngq
    double precision, intent(out) :: h(-ngh+lo(1):hi(1)+ngh,-ngh+lo(2):hi(2)+ngh,-ngh+lo(3):hi(3)+ngh)
    double precision, intent(in ) :: Q(-ngq+lo(1):hi(1)+ngq,-ngq+lo(2):hi(2)+ngq,-ngq+lo(3):hi(3)+ngq,nprim)

    integer :: i,j,k, n, qyn, qhn

    !$omp parallel do private(i,j,k,n,qyn,qhn)
    do k=lo(3),hi(3)
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)
             h(i,j,k) = 0.d0
             do n=1, nspecies
                qhn = qh1+n-1
                qyn = qy1+n-1
                h(i,j,k) = h(i,j,k) + q(i,j,k,qyn)*q(i,j,k,qhn)
             end do
          enddo
       enddo
    enddo
    !$omp end parallel do

  end subroutine make_h_3d


  subroutine make_divu_3d(divu, Q, lo, hi, ngd, ngq, dx)

    use derivative_stencil_module, only : stencil_ng, first_deriv_8

    integer, intent(in) :: lo(3), hi(3), ngd, ngq
    double precision, intent(out) :: divu(-ngd+lo(1):hi(1)+ngd,-ngd+lo(2):hi(2)+ngd,-ngd+lo(3):hi(3)+ngd)
    double precision, intent(in ) ::    Q(-ngq+lo(1):hi(1)+ngq,-ngq+lo(2):hi(2)+ngq,-ngq+lo(3):hi(3)+ngq,nprim)
    double precision, intent(in) :: dx(3)

    integer :: i,j,k
    double precision, allocatable :: ux(:,:,:), vy(:,:,:), wz(:,:,:)
    double precision :: dxinv(3)

    do i=1,3
       dxinv(i) = 1.0d0 / dx(i)
    end do

    allocate(ux(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)))
    allocate(vy(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)))
    allocate(wz(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)))

    !$omp parallel private(i,j,k)

    !$omp do
    do k=lo(3),hi(3)
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)
             ux(i,j,k) = dxinv(1)*first_deriv_8(q(i-4:i+4,j,k,qu))
          enddo
       end do
    end do
    !$omp end do nowait

    !$omp do
    do k=lo(3),hi(3)
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)
             vy(i,j,k) = dxinv(2)*first_deriv_8(q(i,j-4:j+4,k,qv))
          enddo
       end do
    end do
    !$omp end do nowait

    !$omp do
    do k=lo(3),hi(3)
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)
             wz(i,j,k) = dxinv(3)*first_deriv_8(q(i,j,k-4:k+4,qw))
          enddo
       end do
    end do
    !$omp end do 

    !$omp workshare
    divu = ux + vy + wz
    !$omp end workshare

    !$omp end parallel

    deallocate(ux,vy,wz)

  end subroutine make_divu_3d


  subroutine make_omegadot_3d(odot, Q, lo, hi, ngo, ngq)
    use chemistry_module, only : molecular_weight

    integer, intent(in) :: lo(3), hi(3), ngo, ngq
    double precision, intent(out) :: odot(-ngo+lo(1):hi(1)+ngo,-ngo+lo(2):hi(2)+ngo,-ngo+lo(3):hi(3)+ngo,nspecies)
    double precision, intent(in ) ::    Q(-ngq+lo(1):hi(1)+ngq,-ngq+lo(2):hi(2)+ngq,-ngq+lo(3):hi(3)+ngq,nprim)

    integer :: i,j,k,iwrk
    double precision :: rwrk, wdot(nspecies), Xt(nspecies)

    !$omp parallel do private(i,j,k,iwrk,rwrk,wdot,Xt)
    do k=lo(3),hi(3)
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)
             Xt(:) = q(i,j,k,qx1:qx1+nspecies-1)
             call ckwxr(q(i,j,k,qrho), q(i,j,k,qtemp), Xt, iwrk, rwrk, wdot)
             odot(i,j,k,:) = wdot * molecular_weight
          enddo
       enddo
    enddo
    !$omp end parallel do

  end subroutine make_omegadot_3d

end module make_plot_variables_module
