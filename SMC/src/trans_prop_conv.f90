! For convergence study only

module transport_properties

  use chemistry_module
  use eglib_module
  use multifab_module
  use variables_module

  implicit none

  ! eglib parameters
  integer, save :: ITLS=-1, IFLAG=-1

  private

  public get_transport_properties

contains

  subroutine get_transport_properties(Q, mu, xi, lam, Ddiag, ng, ghostcells_only)

    use probin_module, only : use_bulk_viscosity
    use smc_bc_module, only : get_data_lo_hi

    type(multifab), intent(in   ) :: Q
    type(multifab), intent(inout) :: mu, xi, lam, Ddiag
    integer, intent(in), optional :: ng
    logical, intent(in), optional :: ghostcells_only
 
    integer :: ngwork, idim
    logical :: lgco
    integer :: ngq, n, dm, lo(Q%dim), hi(Q%dim), wlo(Q%dim), whi(Q%dim)
    double precision, pointer, dimension(:,:,:,:) :: qp, mup, xip, lamp, dp

    logical, save :: first_call = .true.

    if (first_call) then
       first_call = .false.
       if (use_bulk_viscosity) then
          ITLS  = 1 
          IFLAG = 5
       else
          ITLS  = 1
          IFLAG = 3
       end if
    end if

    dm = Q%dim
    ngq = nghost(Q)

    ngwork = ngq
    if (present(ng)) then
       ngwork = min(ngwork, ng)
    end if

    lgco = .false.
    if (present(ghostcells_only)) then
       lgco = ghostcells_only
    end if

    do n=1,nfabs(Q)
       
       qp => dataptr(Q,n)
       mup => dataptr(mu,n)
       xip => dataptr(xi,n)
       lamp => dataptr(lam,n)
       dp => dataptr(Ddiag,n)

       lo = lwb(get_box(Q,n))
       hi = upb(get_box(Q,n))

       call get_data_lo_hi(n,wlo,whi)

       do idim=1,dm
          wlo(idim) = max(wlo(idim), lo(idim)-ngwork)
          whi(idim) = min(whi(idim), hi(idim)+ngwork)
       end do

       if (dm .ne. 3) then
          call bl_error("Only 3D is supported in get_transport_properties")
       else
          call get_trans_prop_3d(lo,hi,ngq,qp,mup,xip,lamp,dp,wlo,whi,lgco)
       end if

    end do

  end subroutine get_transport_properties

  subroutine get_trans_prop_3d(lo,hi,ng,q,mu,xi,lam,Ddiag,wlo,whi,gco)

    use probin_module, only : use_bulk_viscosity

    logical, intent(in) :: gco  ! ghost cells only
    integer, intent(in) :: lo(3), hi(3), ng, wlo(3), whi(3)
    double precision,intent(in )::    q(lo(1)-ng:hi(1)+ng,lo(2)-ng:hi(2)+ng,lo(3)-ng:hi(3)+ng,nprim)
    double precision,intent(out)::   mu(lo(1)-ng:hi(1)+ng,lo(2)-ng:hi(2)+ng,lo(3)-ng:hi(3)+ng)
    double precision,intent(out)::   xi(lo(1)-ng:hi(1)+ng,lo(2)-ng:hi(2)+ng,lo(3)-ng:hi(3)+ng)
    double precision,intent(out)::  lam(lo(1)-ng:hi(1)+ng,lo(2)-ng:hi(2)+ng,lo(3)-ng:hi(3)+ng)
    double precision,intent(out)::Ddiag(lo(1)-ng:hi(1)+ng,lo(2)-ng:hi(2)+ng,lo(3)-ng:hi(3)+ng,nspecies)

    integer :: i, j, k, n, iwrk
    integer :: np, ii, jj, kk, iisize, jisize, kisize 
    integer :: iindex(whi(1)-wlo(1)-hi(1)+lo(1))
    integer :: jindex(whi(2)-wlo(2)-hi(2)+lo(2))
    integer :: kindex(whi(3)-wlo(3)-hi(3)+lo(3))
    double precision :: rwrk
    double precision, allocatable :: Tt(:), Xt(:,:), Yt(:,:), Cpt(:,:), D(:,:)
    double precision, allocatable :: ME1(:), ME2(:), L1(:), L2(:)

    ! eglib parameters
    integer, parameter :: ITLS=0, IFLAG=2

    if (.not. gco) then

       np = whi(1) - wlo(1) + 1
       call eglib_init(nspecies, np, ITLS, IFLAG)
       
       !$omp parallel private(i,j,k,n,iwrk,rwrk,ii) &
       !$omp private(Tt,Xt,Yt,Cpt,D,ME1,ME2,L1,L2)
       
       allocate(Tt(np))
       allocate(ME1(np))
       allocate(ME2(np))
       allocate(L1(np))
       allocate(L2(np))
       
       allocate(Xt(nspecies,np))
       allocate(Yt(nspecies,np))
       allocate(Cpt(nspecies,np))
       allocate(D(nspecies,np))
       
       !$omp do
       do k=wlo(3),whi(3)
          do j=wlo(2),whi(2)
             
             do i=wlo(1), whi(1)
                ii = i-wlo(1)+1
                Tt(  ii) = q(i,j,k,qtemp)
                CALL CKCPMS(Tt(ii), iwrk, rwrk, Cpt(:,ii))
                Yt(:,ii) = q(i,j,k,qy1:qy1+nspecies-1)
                Xt(:,ii) = q(i,j,k,qx1:qx1+nspecies-1)
             end do
          
             CALL EGMPAR(np, Tt, Xt, Yt, Cpt, egwork, egiwork)
          
!             CALL EGME3(np, Tt, Yt, egwork, ME) 
             CALL EGME1(np,  1.d0, Tt, Xt, egwork, ME1) 
             CALL EGME1(np, -1.d0, Tt, Xt, egwork, ME2) 
             mu(wlo(1):whi(1),j,k) = 0.5d0*(ME1+ME2)
          
!             CALL EGMK3(np, Tt, Yt, egwork, MK) 
             xi(wlo(1):whi(1),j,k) = 0.d0
             
             CALL EGMVR1(np, Tt, Yt, egwork, D)
             do n=1,nspecies
                do i=wlo(1), whi(1)
                   ii = i-wlo(1)+1
                   Ddiag(i,j,k,n) = D(n,ii)
                end do
             end do
             
             CALL EGML1(np,  1.d0, Tt, Xt, egwork, L1)
             CALL EGML1(np, -1.d0, Tt, Xt, egwork, L2)
             lam(wlo(1):whi(1),j,k) = 0.5d0*(L1+L2)
             
          end do
       end do
       !$omp end do
       
       deallocate(Tt, Xt, Yt, Cpt, D, ME1, ME2, L1, L2)
       !$omp end parallel

    else ! ghost cells only 

       call bl_error("transport_properties: Should not be there")

    end if

  end subroutine get_trans_prop_3d

end module transport_properties
