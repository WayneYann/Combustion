! For convergence study only

module transport_properties

  use chemistry_module
  use multifab_module
  use variables_module

  use egz_module

  implicit none

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

       if (dm .eq. 1) then
          call get_trans_prop_1d(lo,hi,ngq,qp,mup,xip,lamp,dp,wlo,whi,lgco)
       else if (dm .eq. 2) then
          call get_trans_prop_2d(lo,hi,ngq,qp,mup,xip,lamp,dp,wlo,whi,lgco)
       else
          call get_trans_prop_3d(lo,hi,ngq,qp,mup,xip,lamp,dp,wlo,whi,lgco)
       end if

    end do

  end subroutine get_transport_properties

  subroutine get_trans_prop_1d(lo,hi,ng,q,mu,xi,lam,Ddiag,wlo,whi,gco)
    use probin_module, only : use_bulk_viscosity
    logical, intent(in) :: gco  ! ghost cells only
    integer, intent(in) :: lo(1), hi(1), ng, wlo(1), whi(1)
    double precision,intent(in )::    q(lo(1)-ng:hi(1)+ng,nprim)
    double precision,intent(out)::   mu(lo(1)-ng:hi(1)+ng)
    double precision,intent(out)::   xi(lo(1)-ng:hi(1)+ng)
    double precision,intent(out)::  lam(lo(1)-ng:hi(1)+ng)
    double precision,intent(out)::Ddiag(lo(1)-ng:hi(1)+ng,nspecies)

    integer, parameter :: np = 4

    integer :: i, n, nptot, qxn, iwrk, ib, nb, istart, iend
    double precision :: rwrk, Cp(nspecies)
    double precision :: L1Z(np), L2Z(np), DZ(np,nspecies), XZ(np,nspecies), &
         CPZ(np,nspecies), E1Z(np), E2Z(np)

    if (gco) call bl_error("get_trans_prop_1d: ghost-cells-only not supported")

    nptot = whi(1) - wlo(1) + 1
    nb = nptot / np
    if (nb*np .ne. nptot) call bl_error("get_trans_prop_1d: grid size not supported")

    !$omp parallel private(i,n,qxn,iwrk,ib,istart,iend,rwrk) &
    !$omp private(Cp,E1Z,E2Z,L1Z,L2Z,DZ,XZ,CPZ)

    call egzini(np)
              
    !$omp do
    do ib=0,nb-1
       
       istart = wlo(1) + ib*np
       iend = istart + np - 1

       do n=1,nspecies
          qxn = qx1+n-1
          do i=istart,iend
             XZ(i-start+1,n) = q(i,qxn)
          end do
       end do
          
       if (iflag > 3) then
          do i=istart,iend
             call ckcpms(q(i,qtemp), iwrk, rwrk, Cp)
             CPZ(i-istart+1,:) = Cp
          end do
       else
          CPZ = 0.d0
       end if
       
       call egzpar(q(istart:iend,qtemp), XZ, CPZ)
       
       call egze1( 1.d0, XZ, E1Z)
       call egze1(-1.d0, XZ, E2Z)
       mu(istart:iend) = 0.5d0*(E1Z+E2Z)

       xi(istart:iend) = 0.d0
       
       call egzl1( 1.d0, XZ, L1Z)
       call egzl1(-1.d0, XZ, L2Z)
       lam(istart:iend) = 0.5d0*(L1Z+L2Z)
       
       call EGZVR1(q(istart:iend,qtemp), DZ)
       do n=1,nspecies
          do i=istart,iend
             Ddiag(i,n) = DZ(i-istart+1,n)
          end do
       end do
       
    end do
    !$omp end do
    !$omp end parallel

  end subroutine get_trans_prop_1d

  subroutine get_trans_prop_2d(lo,hi,ng,q,mu,xi,lam,Ddiag,wlo,whi,gco)
    logical, intent(in) :: gco  ! ghost cells only
    integer, intent(in) :: lo(3), hi(3), ng, wlo(3), whi(3)
    double precision,intent(in )::    q(lo(1)-ng:hi(1)+ng,lo(2)-ng:hi(2)+ng,nprim)
    double precision,intent(out)::   mu(lo(1)-ng:hi(1)+ng,lo(2)-ng:hi(2)+ng)
    double precision,intent(out)::   xi(lo(1)-ng:hi(1)+ng,lo(2)-ng:hi(2)+ng)
    double precision,intent(out)::  lam(lo(1)-ng:hi(1)+ng,lo(2)-ng:hi(2)+ng)
    double precision,intent(out)::Ddiag(lo(1)-ng:hi(1)+ng,lo(2)-ng:hi(2)+ng,nspecies)

    integer :: i, j, n, np, qxn, iwrk
    double precision :: rwrk, Cp(nspecies)
    double precision, allocatable :: L1Z(:), L2Z(:), DZ(:,:), XZ(:,:), CPZ(:,:), &
         E1Z(:), E2Z(:)

    if (.not. gco) then

       np = whi(1) - wlo(1) + 1
       
       !$omp parallel private(i,j,n,qxn,iwrk) &
       !$omp private(rwrk,Cp,L1Z,L2Z,DZ,XZ,CPZ,E1Z,E2Z)

       call egzini(np)

       allocate(L1Z(wlo(1):whi(1)))
       allocate(L2Z(wlo(1):whi(1)))

       allocate(DZ(wlo(1):whi(1),nspecies))
       allocate(XZ(wlo(1):whi(1),nspecies))
       allocate(CPZ(wlo(1):whi(1),nspecies))

       allocate(E1Z(wlo(1):whi(1)))
       allocate(E2Z(wlo(1):whi(1)))
       
       !$omp do
       do j=wlo(2),whi(2)

          do n=1,nspecies
             qxn = qx1+n-1
             do i=wlo(1),whi(1)
                XZ(i,n) = q(i,j,qxn)
             end do
          end do

          if (iflag > 3) then
             do i=wlo(1),whi(1)
                call ckcpms(q(i,j,qtemp), iwrk, rwrk, Cp)
                CPZ(i,:) = Cp
             end do
          else
             CPZ = 0.d0
          end if
          
          call egzpar(q(wlo(1):whi(1),j,qtemp), XZ, CPZ)
          
          call egze1( 1.d0, XZ, E1Z)
          call egze1(-1.d0, XZ, E2Z)
          mu(wlo(1):whi(1),j) = 0.5d0*(E1Z+E2Z)
          
          xi(wlo(1):whi(1),j) = 0.d0
          
          call egzl1( 1.d0, XZ, L1Z)
          call egzl1(-1.d0, XZ, L2Z)
          lam(wlo(1):whi(1),j) = 0.5d0*(L1Z+L2Z)
          
          call EGZVR1(q(wlo(1):whi(1),j,qtemp), DZ)
          do n=1,nspecies
             do i=wlo(1),whi(1)
                Ddiag(i,j,n) = DZ(i,n)
             end do
          end do
          
       end do
       !$omp end do
       
       deallocate(L1Z, L2Z, DZ, XZ, CPZ, E1Z, E2Z)
       !$omp end parallel

    else ! ghost cells only 

       call bl_error("transport_properties: Should not be there")

    end if

  end subroutine get_trans_prop_2d


  subroutine get_trans_prop_3d(lo,hi,ng,q,mu,xi,lam,Ddiag,wlo,whi,gco)
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

    if (.not. gco) then

       np = whi(1) - wlo(1) + 1
       
       !$omp parallel private(i,j,k,n,qxn,iwrk) &
       !$omp private(rwrk,Cp,L1Z,L2Z,DZ,XZ,CPZ,E1Z,E2Z)

       call egzini(np)

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
