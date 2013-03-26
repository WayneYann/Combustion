module make_plot_variables_module

  use multifab_module

  use smc_bc_module
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
    integer :: dlo(Q%dim), dhi(Q%dim)

    dm = Q%dim
    ngq = nghost(Q)
    ngpd = nghost(plotdata)

    do i = 1, nfabs(Q)
       qp => dataptr(Q, i)
       pdp => dataptr(plotdata, i)
       
       lo = lwb(get_box(Q, i))
       hi = upb(get_box(Q, i))

       call get_data_lo_hi(i,dlo,dhi)

       if (dm .ne. 3) then
          call bl_error("Only 3D make_plotvar is support")
       end if

       if (icomp .eq. icomp_h) then
          call make_h_3d(pdp(:,:,:,icomp), qp, lo, hi, ngpd, ngq)                
       else if (icomp .eq. icomp_divu) then
          call make_divu_3d(pdp(:,:,:,icomp), qp, lo, hi, ngpd, ngq, dx, dlo,dhi)
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


  subroutine make_divu_3d(divu, Q, lo, hi, ngd, ngq, dx, dlo, dhi)

    use derivative_stencil_module, only : stencil_ng, first_deriv_8, first_deriv_6, &
                first_deriv_4, first_deriv_l3, first_deriv_r3, first_deriv_rb, first_deriv_lb

    integer, intent(in) :: lo(3), hi(3), ngd, ngq, dlo(3), dhi(3)
    double precision, intent(out) :: divu(-ngd+lo(1):hi(1)+ngd,-ngd+lo(2):hi(2)+ngd,-ngd+lo(3):hi(3)+ngd)
    double precision, intent(in ) ::    Q(-ngq+lo(1):hi(1)+ngq,-ngq+lo(2):hi(2)+ngq,-ngq+lo(3):hi(3)+ngq,nprim)
    double precision, intent(in) :: dx(3)

    integer :: i,j,k
    double precision, allocatable :: ux(:,:,:), vy(:,:,:), wz(:,:,:)
    double precision :: dxinv(3)
    integer :: slo(3), shi(3) 

    ! Only the region bounded by [dlo,dhi] contains good data.
    ! [slo,shi] will be safe for 8th-order stencil
    slo = dlo + stencil_ng
    shi = dhi - stencil_ng

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

          do i=slo(1),shi(1)
             ux(i,j,k) = dxinv(1)*first_deriv_8(q(i-4:i+4,j,k,qu))
          enddo
          
          ! lo-x boundary
          if (dlo(1) .eq. lo(1)) then
             i = lo(1)
             ! use completely right-biased stencil
             ux(i,j,k) = dxinv(1)*first_deriv_rb(q(i:i+3,j,k,qu))

             i = lo(1)+1
             ! use 3rd-order slightly right-biased stencil
             ux(i,j,k) = dxinv(1)*first_deriv_r3(q(i-1:i+2,j,k,qu))

             i = lo(1)+2
             ! use 4th-order stencil
             ux(i,j,k) = dxinv(1)*first_deriv_4(q(i-2:i+2,j,k,qu))

             i = lo(1)+3
             ! use 6th-order stencil
             ux(i,j,k) = dxinv(1)*first_deriv_6(q(i-3:i+3,j,k,qu))
          end if

          ! hi-x boundary
          if (dhi(1) .eq. hi(1)) then
             i = hi(1)-3
             ! use 6th-order stencil
             ux(i,j,k) = dxinv(1)*first_deriv_6(q(i-3:i+3,j,k,qu))

             i = hi(1)-2
             ! use 4th-order stencil
             ux(i,j,k) = dxinv(1)*first_deriv_4(q(i-2:i+2,j,k,qu))

             i = hi(1)-1
             ! use 3rd-order slightly left-biased stencil
             ux(i,j,k) = dxinv(1)*first_deriv_l3(q(i-2:i+1,j,k,qu))

             i = hi(1)
             ! use completely left-biased stencil
             ux(i,j,k) = dxinv(1)*first_deriv_lb(q(i-3:i,j,k,qu))
          end if

       end do
    end do
    !$omp end do nowait

    !$omp do
    do k=lo(3),hi(3)

       do j=slo(2),shi(2)
          do i=lo(1),hi(1)
             vy(i,j,k) = dxinv(2)*first_deriv_8(q(i,j-4:j+4,k,qv))
          enddo
       end do

       ! lo-y boundary
       if (dlo(2) .eq. lo(2)) then
          j = lo(2)
          ! use completely right-biased stencil
          do i=lo(1),hi(1)
             vy(i,j,k) = dxinv(2)*first_deriv_rb(q(i,j:j+3,k,qv))
          enddo

          j = lo(2)+1
          ! use 3rd-order slightly right-biased stencil
          do i=lo(1),hi(1)
             vy(i,j,k) = dxinv(2)*first_deriv_r3(q(i,j-1:j+2,k,qv))
          enddo

          j = lo(2)+2
          ! use 4th-order stencil
          do i=lo(1),hi(1)
             vy(i,j,k) = dxinv(2)*first_deriv_4(q(i,j-2:j+2,k,qv))
          enddo

          j = lo(2)+3
          ! use 6th-order stencil
          do i=lo(1),hi(1)
             vy(i,j,k) = dxinv(2)*first_deriv_6(q(i,j-3:j+3,k,qv))
          enddo
       end if

       ! hi-y boundary
       if (dhi(2) .eq. hi(2)) then
          j = hi(2)-3
          ! use 6th-order stencil
          do i=lo(1),hi(1)
             vy(i,j,k) = dxinv(2)*first_deriv_6(q(i,j-3:j+3,k,qv))
          enddo

          j = hi(2)-2
          ! use 4th-order stencil
          do i=lo(1),hi(1)
             vy(i,j,k) = dxinv(2)*first_deriv_4(q(i,j-2:j+2,k,qv))
          enddo

          j = hi(2)-1
          ! use 3rd-order slightly left-biased stencil
          do i=lo(1),hi(1)
             vy(i,j,k) = dxinv(2)*first_deriv_l3(q(i,j-2:j+1,k,qv))
          enddo

          j = hi(2)
          ! use completely left-biased stencil
          do i=lo(1),hi(1)
             vy(i,j,k) = dxinv(2)*first_deriv_lb(q(i,j-3:j,k,qv))
          enddo
       end if

    end do
    !$omp end do nowait

    !$omp do
    do k=slo(3),shi(3)
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)
             wz(i,j,k) = dxinv(3)*first_deriv_8(q(i,j,k-4:k+4,qw))
          enddo
       end do
    end do
    !$omp end do nowait

    ! lo-z boundary
    if (dlo(3) .eq. lo(3)) then
       k = lo(3)
       ! use completely right-biased stencil
       !$OMP DO
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)
             wz(i,j,k) = dxinv(3)*first_deriv_rb(q(i,j,k:k+3,qw))
          enddo
       enddo
       !$OMP END DO NOWAIT

       k = lo(3)+1
       ! use 3rd-order slightly right-biased stencil
       !$OMP DO
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)
             wz(i,j,k) = dxinv(3)*first_deriv_r3(q(i,j,k-1:k+2,qw))
          enddo
       enddo
       !$OMP END DO NOWAIT

       k = lo(3)+2
       ! use 4th-order stencil
       !$OMP DO
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)
             wz(i,j,k) = dxinv(3)*first_deriv_4(q(i,j,k-2:k+2,qw))
          enddo
       enddo
       !$OMP END DO NOWAIT

       k = lo(3)+3
       ! use 6th-order stencil
       !$OMP DO
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)
             wz(i,j,k) = dxinv(3)*first_deriv_6(q(i,j,k-3:k+3,qw))
          enddo
       enddo
       !$OMP END DO NOWAIT
    end if

    ! hi-z boundary
    if (dhi(3) .eq. hi(3)) then
       k = hi(3)-3
       ! use 6th-order stencil
       !$OMP DO
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)
             wz(i,j,k) = dxinv(3)*first_deriv_6(q(i,j,k-3:k+3,qw))
          enddo
       enddo
       !$OMP END DO NOWAIT

       k = hi(3)-2
       ! use 4th-order stencil
       !$OMP DO
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)
             wz(i,j,k) = dxinv(3)*first_deriv_4(q(i,j,k-2:k+2,qw))
          enddo
       enddo
       !$OMP END DO NOWAIT

       k = hi(3)-1
       ! use 3rd-order slightly left-biased stencil
       !$OMP DO
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)
             wz(i,j,k) = dxinv(3)*first_deriv_l3(q(i,j,k-2:k+1,qw))
          enddo
       enddo
       !$OMP END DO NOWAIT

       k = hi(3)
       ! use completely left-biased stencil
       !$OMP DO
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)
             wz(i,j,k) = dxinv(3)*first_deriv_lb(q(i,j,k-3:k,qw))
          enddo
       enddo
       !$OMP END DO NOWAIT
    end if

    !$omp barrier

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

    integer :: i,j,k,n,np,iwrk
    double precision :: Yt(lo(1):hi(1),nspecies), wdot(lo(1):hi(1),nspecies), rwrk

    np = hi(1) - lo(1) + 1

    !$omp parallel do private(i,j,k,n,iwrk,rwrk,wdot,Yt)
    do k=lo(3),hi(3)
       do j=lo(2),hi(2)

          do n=1, nspecies
             do i=lo(1),hi(1)
                Yt(i,n) = q(i,j,k,qy1+n-1)
             end do
          end do
          
          call vckwyr(np, q(lo(1),j,k,qrho), q(lo(1),j,k,qtemp), Yt, iwrk, rwrk, wdot)

          do n=1, nspecies
             do i=lo(1),hi(1)
                odot(i,j,k,n) = wdot(i,n) * molecular_weight(n)
             end do
          end do

       enddo
    enddo
    !$omp end parallel do

  end subroutine make_omegadot_3d

end module make_plot_variables_module
