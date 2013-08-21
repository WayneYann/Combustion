! For convergence study only

module transport_properties

  use chemistry_module
  use multifab_module
  use variables_module

  use egz_module

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

    integer :: i, j, k, n, np, qxn, iwrk
    double precision :: rwrk, Cp(nspecies)
    double precision, allocatable :: L1Z(:), L2Z(:), DZ(:,:), XZ(:,:), CPZ(:,:), &
         E1Z(:), E2Z(:)

    ! eglib parameters
    integer, parameter :: ITLS_local=0, IFLAG_local=2

    if (.not. gco) then

       np = whi(1) - wlo(1) + 1

       call egzini(np, ITLS_local, IFLAG_local)
       
       !$omp parallel private(i,j,k,n,qxn,iwrk) &
       !$omp private(rwrk,Cp,L1Z,L2Z,DZ,XZ,CPZ,E1Z,E2Z)

       allocate(L1Z(wlo(1):whi(1)))
       allocate(L2Z(wlo(1):whi(1)))

       allocate(DZ(wlo(1):whi(1),nspecies))
       allocate(XZ(wlo(1):whi(1),nspecies))
       allocate(CPZ(wlo(1):whi(1),nspecies))

       allocate(E1Z(wlo(1):whi(1)))
       allocate(E2Z(wlo(1):whi(1)))
       
       !$omp do
       do k=wlo(3),whi(3)
          do j=wlo(2),whi(2)

             do n=1,nspecies
                qxn = qx1+n-1
                do i=wlo(1),whi(1)
                   XZ(i,n) = q(i,j,k,qxn)
                end do
             end do

             if (iflag > 3) then
                do i=wlo(1),whi(1)
                   call ckcpms(q(i,j,k,qtemp), iwrk, rwrk, Cp)
                   CPZ(i,:) = Cp
                end do
             else
                CPZ = 0.d0
             end if

             call egzpar(q(wlo(1):whi(1),j,k,qtemp), XZ, CPZ)

             call egze1( 1.d0, XZ, E1Z)
             call egze1(-1.d0, XZ, E2Z)
             mu(wlo(1):whi(1),j,k) = 0.5d0*(E1Z+E2Z)

             xi(wlo(1):whi(1),j,k) = 0.d0

             call egzl1( 1.d0, XZ, L1Z)
             call egzl1(-1.d0, XZ, L2Z)
             lam(wlo(1):whi(1),j,k) = 0.5d0*(L1Z+L2Z)

             call EGZVR1(q(wlo(1):whi(1),j,k,qtemp), DZ)
             do n=1,nspecies
                do i=wlo(1),whi(1)
                   Ddiag(i,j,k,n) = DZ(i,n)
                end do
             end do

          end do
       end do
       !$omp end do
       
       deallocate(L1Z, L2Z, DZ, XZ, CPZ, E1Z, E2Z)
       !$omp end parallel

    else ! ghost cells only 

       call bl_error("transport_properties: Should not be there")

    end if

  end subroutine get_trans_prop_3d

end module transport_properties
