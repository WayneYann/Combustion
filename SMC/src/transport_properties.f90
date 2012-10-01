module transport_properties

  use chemistry_module
  use eglib_module
  use multifab_module
  use variables_module

  implicit none

  private

  public get_transport_properties

contains

  subroutine get_transport_properties(Q, mu, xi, lam, Ddiag)

    use smc_bc_module, only : get_data_lo_hi

    type(multifab), intent(in   ) :: Q
    type(multifab), intent(inout) :: mu, xi, lam, Ddiag
 
    integer :: ng, n, dm, lo(Q%dim), hi(Q%dim), dlo(Q%dim), dhi(Q%dim)
    double precision, pointer, dimension(:,:,:,:) :: qp, mup, xip, lamp, dp

    dm = Q%dim
    ng = nghost(Q)

    do n=1,nfabs(Q)
       
       qp => dataptr(Q,n)
       mup => dataptr(mu,n)
       xip => dataptr(xi,n)
       lamp => dataptr(lam,n)
       dp => dataptr(Ddiag,n)

       lo = lwb(get_box(Q,n))
       hi = upb(get_box(Q,n))

       call get_data_lo_hi(n,dlo,dhi)

       if (dm .ne. 3) then
          call bl_error("Only 3D is supported in get_transport_properties")
       else
          call get_trans_prop_3d(lo,hi,ng,qp,mup,xip,lamp,dp,dlo,dhi)
       end if

    end do

  end subroutine get_transport_properties
   
  subroutine get_trans_prop_3d(lo,hi,ng,q,mu,xi,lam,Ddiag,dlo,dhi)
    use omp_module
    integer, intent(in) :: lo(3), hi(3), ng, dlo(3), dhi(3)
    double precision,intent(in )::    q(lo(1)-ng:hi(1)+ng,lo(2)-ng:hi(2)+ng,lo(3)-ng:hi(3)+ng,nprim)
    double precision,intent(out)::   mu(lo(1)-ng:hi(1)+ng,lo(2)-ng:hi(2)+ng,lo(3)-ng:hi(3)+ng)
    double precision,intent(out)::   xi(lo(1)-ng:hi(1)+ng,lo(2)-ng:hi(2)+ng,lo(3)-ng:hi(3)+ng)
    double precision,intent(out)::  lam(lo(1)-ng:hi(1)+ng,lo(2)-ng:hi(2)+ng,lo(3)-ng:hi(3)+ng)
    double precision,intent(out)::Ddiag(lo(1)-ng:hi(1)+ng,lo(2)-ng:hi(2)+ng,lo(3)-ng:hi(3)+ng,nspecies)

    integer :: i, j, k, n, iwrk
    double precision :: rwrk
    integer :: np, ii
    double precision :: Cpck(nspecies)
    double precision, allocatable :: Tt(:), Xt(:,:), Yt(:,:), Cpt(:,:), Wtm(:), D(:,:)
    double precision, allocatable :: ME(:), MK(:), L1(:), L2(:)

    ! eglib parameters
    integer, parameter :: ITLS=1, IFLAG=5

    np = dhi(1) - dlo(1) + 1
    call eglib_init(nspecies, np, ITLS, IFLAG)

    allocate(Tt(np))
    allocate(Xt(nspecies,np))
    allocate(Yt(nspecies,np))
    allocate(Cpt(nspecies,np))
    allocate(D(nspecies,np))

    allocate(Wtm(np))
    allocate(ME(np))
    allocate(MK(np))
    allocate(L1(np))
    allocate(L2(np))

    !$omp parallel do private(i,j,k,n,iwrk,rwrk,ii,Cpck) &
    !$omp private(Tt,Wtm,Xt,Yt,Cpt,D,ME,MK,L1,L2)
    do k=dlo(3),dhi(3)
    do j=dlo(2),dhi(2)

       do i=dlo(1), dhi(1)
          ii = i-dlo(1)+1
          Tt(  ii) = q(i,j,k,qtemp)
          Xt(:,ii) = q(i,j,k,qx1:qx1+nspecies-1)
          Yt(:,ii) = q(i,j,k,qy1:qy1+nspecies-1)
          CALL CKCPMS(Tt(ii), iwrk, rwrk, Cpck)
          Cpt(:,ii) = Cpck
          CALL CKMMWY(Yt(:,ii), iwrk, rwrk, Wtm(ii))
       end do
       
       CALL EGMPAR(np, Tt, Xt, Yt, Cpt, egwork, egiwork)
       
       CALL EGME3(np, Tt, Yt, egwork, ME) 
       mu(dlo(1):dhi(1),j,k) = ME

       CALL EGMK3(np, Tt, Yt, egwork, MK) 
       xi(dlo(1):dhi(1),j,k) = MK
       
       CALL EGMVR1(np, Tt, Yt, egwork, D)
       do n=1,nspecies
          do i=dlo(1), dhi(1)
             ii = i-dlo(1)+1
             Ddiag(i,j,k,n) = D(n,ii) * Wtm(ii) / molecular_weight(n)
          end do
       end do
       
       CALL EGML1(np,  1.d0, Tt, Xt, egwork, L1)
       CALL EGML1(np, -1.d0, Tt, Xt, egwork, L2)
       lam(dlo(1):dhi(1),j,k) = 0.5d0*(L1+L2)

    end do
    end do
    !$omp end parallel do

    deallocate(Tt, Xt, Yt, Cpt, Wtm, D, ME, MK, L1, L2)

    call eglib_close()

  end subroutine get_trans_prop_3d

end module transport_properties
