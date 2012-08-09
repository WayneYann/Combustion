module transport_properties

  use chemistry_module
  use eglib_module
  use multifab_module
  use variables

  implicit none

  private

  public get_transport_properties

contains

  subroutine get_transport_properties(Q, mu, xi, lam, Ddiag)

    use multifab_module
    use variables

    type(multifab), intent(in   ) :: Q
    type(multifab), intent(inout) :: mu, xi, lam, Ddiag
 
    integer :: ng, n, dm, lo(Q%dim), hi(Q%dim)
    double precision, pointer, dimension(:,:,:,:) :: qp, mup, xip, lamp, dp

    dm = Q%dim
    ng = nghost(Q)

    do n=1,nboxes(Q)
       if ( remote(Q,n) ) cycle
       
       qp => dataptr(Q,n)
       mup => dataptr(mu,n)
       xip => dataptr(xi,n)
       lamp => dataptr(lam,n)
       dp => dataptr(Ddiag,n)

       lo = lwb(get_box(Q,n))
       hi = upb(get_box(Q,n))

       if (dm .ne. 3) then
          call bl_error("Only 3D is supported in get_transport_properties")
       else
          call get_trans_prop_3d(lo,hi,ng,qp,mup,xip,lamp,dp)
       end if

    end do

  end subroutine get_transport_properties
   
  subroutine get_trans_prop_3d(lo,hi,ng,q,mu,lam,xi,Ddiag)
    use omp_module
    integer, intent(in) :: lo(3), hi(3), ng
    double precision,intent(in )::    q(lo(1)-ng:hi(1)+ng,lo(2)-ng:hi(2)+ng,lo(3)-ng:hi(3)+ng,nprim)
    double precision,intent(out)::   mu(lo(1)-ng:hi(1)+ng,lo(2)-ng:hi(2)+ng,lo(3)-ng:hi(3)+ng)
    double precision,intent(out)::   xi(lo(1)-ng:hi(1)+ng,lo(2)-ng:hi(2)+ng,lo(3)-ng:hi(3)+ng)
    double precision,intent(out)::  lam(lo(1)-ng:hi(1)+ng,lo(2)-ng:hi(2)+ng,lo(3)-ng:hi(3)+ng)
    double precision,intent(out)::Ddiag(lo(1)-ng:hi(1)+ng,lo(2)-ng:hi(2)+ng,lo(3)-ng:hi(3)+ng,nspecies)

    integer :: i, j, k, n, iwrk, tid
    double precision :: rwrk, Tt, Wtm
    double precision, dimension(nspecies) :: Xt, Yt, Cpt, D

    !$omp parallel do private(tid,i,j,k,n,iwrk,rwrk,Tt,Wtm,Xt,Yt,Cpt,D)
    do k=lo(3)-ng,hi(3)+ng
    do j=lo(2)-ng,hi(2)+ng
    do i=lo(1)-ng,hi(1)+ng

       tid = omp_get_thread_num()

       Tt = q(i,j,k,qtemp)
       Xt = q(i,j,k,qx1:qx1+nspecies-1)
       Yt = q(i,j,k,qy1:qy1+nspecies-1)

       CALL CKCPMS(Tt, iwrk, rwrk, Cpt)
       CALL CKMMWY(Yt, iwrk, rwrk, Wtm)

       CALL EGSPAR(Tt, Xt, Yt, Cpt, egwork(:,tid), egiwork(:,tid))

       CALL EGSE3(Tt, Yt, egwork(:,tid), mu(i,j,k)) 
       CALL EGSK3(Tt, Yt, egwork(:,tid), xi(i,j,k)) 
       CALL EGSVR1(Tt, Yt, egwork(:,tid), D)
       CALL EGSPTC2(Tt, egwork(:,tid), egiwork(:,tid), lam(i,j,k))

       Ddiag(i,j,k,:) = D(:) * Wtm / molecular_weight(:)

    end do
    end do
    end do
    !$omp end parallel do

  end subroutine get_trans_prop_3d

end module transport_properties
