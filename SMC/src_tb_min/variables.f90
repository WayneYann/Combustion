module variables_module

  use chemistry_module, only : nspecies
  use multifab_module
  use omp_module
  use threadbox_module

  implicit none

  ! Indices
  integer, save :: irho, imx, imy, imz, iene, iry1
  integer, save :: qrho, qu, qv, qw, qpres, qtemp, qe, qy1, qx1, qh1

  ! Number of conserved and primitive variables
  integer, save :: ncons, nprim

  ! Indices for S3D style first-derivatives
  integer, parameter :: idu=1, idv=2, idw=3, idT=4, idp=5, idX1=6 

  ! Arithmetic constants 
  double precision, parameter :: Zero          = 0.d0
  double precision, parameter :: One           = 1.d0
  double precision, parameter :: OneThird      = 1.d0/3.d0
  double precision, parameter :: TwoThirds     = 2.d0/3.d0
  double precision, parameter :: FourThirds    = 4.d0/3.d0
  double precision, parameter :: OneQuarter    = 1.d0/4.d0
  double precision, parameter :: ThreeQuarters = 3.d0/4.d0

contains

  ! 
  ! Initialize various indices
  !
  subroutine init_variables()

    use probin_module, only: dm_in

    irho = 1
    imx = 2
    imy = 3
    if (dm_in .eq. 3) then
       imz = 4
       iene = 5
    else
       imz = -1
       iene = 4
    end if
    iry1 = iene+1

    ncons = iry1-1 + nspecies

    qrho = 1
    qu = 2
    qv = 3
    if (dm_in .eq. 3) then
       qw = 4
       qpres = 5
       qtemp = 6
       qe    = 7
       qy1   = 8
    else
       qw = -1
       qpres = 4
       qtemp = 5
       qe    = 6
       qy1   = 7
    end if
    qx1 = qy1 + nspecies
    qh1 = qx1 + nspecies

    nprim = qh1-1 + nspecies
    
  end subroutine init_variables

  !
  ! Convert conserved variables U to primitive variables Q
  !
  subroutine ctoprim(U, Q, ng, ghostcells_only)
    type(multifab), intent(in   ) :: U
    type(multifab), intent(inout) :: Q
    integer, optional, intent(in) :: ng
    logical, optional, intent(in) :: ghostcells_only

    logical :: lgco
    integer :: ngu, ngq, ngto
    integer :: n, lo(U%dim), hi(U%dim), tid, ulo(4), uhi(4), qlo(4), qhi(4)
    double precision, pointer, dimension(:,:,:,:) :: up, qp

    ngu = nghost(U)
    ngq = nghost(Q)

    if (present(ng)) then
       ngto = ng
    else
       ngto = min(ngu, ngq)
    end if

    lgco = .false.
    if (present(ghostcells_only)) then
       lgco = .true.
    end if

    !$omp parallel private(tid,n,up,qp,lo,hi,qlo,qhi,ulo,uhi)
    tid = omp_get_thread_num()
    do n=1,nfabs(Q)
       up => dataptr(U,n)
       qp => dataptr(Q,n)

       ulo = lbound(up)
       uhi = ubound(up)

       qlo = lbound(qp)
       qhi = ubound(qp)

       if (ngto .eq. 0) then
          lo = tb_get_valid_lo(tid,n)
          hi = tb_get_valid_hi(tid,n)
       else
          lo = tb_get_grown_lo(tid,n)
          hi = tb_get_grown_hi(tid,n)
       end if

       call ctoprim_3d(lo,hi,up,ulo(1:3),uhi(1:3),qp,qlo(1:3),qhi(1:3),lgco)
    end do
    !$omp end parallel

  end subroutine ctoprim

  subroutine ctoprim_3d(lo, hi, u, ulo, uhi, q, qlo, qhi, gco)
    logical, intent(in) :: gco  ! ghost cells only?
    integer, intent(in) :: lo(3), hi(3), ulo(3), uhi(3), qlo(3), qhi(3)
    double precision,intent(in) :: u(ulo(1):uhi(1),ulo(2):uhi(2),ulo(3):uhi(3),ncons)
    double precision            :: q(qlo(1):qhi(1),qlo(2):qhi(2),qlo(3):qhi(3),nprim)
    
    integer :: i, j, k, n, iwrk
    double precision :: rho, rhoinv, rwrk, X(nspecies), Y(nspecies), h(nspecies), ei, Tt, Pt

    do k = lo(3),hi(3)
       do j = lo(2),hi(2)
          do i = lo(1),hi(1)
             
             if (gco) then
                call bl_error("ctoprim_3d: gco needs to be fixed.")
                if ( (i.ge.lo(1) .and. i.le.hi(1)) .and. &
                     (j.ge.lo(2) .and. j.le.hi(2)) .and. &
                     (k.ge.lo(3) .and. k.le.hi(3)) ) then
                   cycle
                end if
             end if

             rho = u(i,j,k,irho)
             rhoinv = 1.d0/rho
             q(i,j,k,qrho) = rho
             q(i,j,k,qu) = u(i,j,k,imx) * rhoinv
             q(i,j,k,qv) = u(i,j,k,imy) * rhoinv
             q(i,j,k,qw) = u(i,j,k,imz) * rhoinv

             do n=1,nspecies
                Y(n) = u(i,j,k,iry1+n-1) * rhoinv
                q(i,j,k,qy1+n-1) = Y(n)
             end do

             call ckytx(Y, iwrk, rwrk, X)

             do n=1,nspecies
                q(i,j,k,qx1+n-1) = X(n)
             end do

             ei = rhoinv*u(i,j,k,iene) - 0.5d0*(q(i,j,k,qu)**2+q(i,j,k,qv)**2+q(i,j,k,qw)**2)
             q(i,j,k,qe) = ei

             call feeytt(ei, Y, iwrk, rwrk, Tt)
             q(i,j,k,qtemp) = Tt

             call CKPY(rho, Tt, Y, iwrk, rwrk, Pt)
             q(i,j,k,qpres) = Pt

             call ckhms(Tt, iwrk, rwrk, h)

             do n=1,nspecies
                q(i,j,k,qh1+n-1) = h(n)
             end do
          enddo
       enddo
    enddo

  end subroutine ctoprim_3d

  !
  ! Compute total density
  !
  subroutine reset_density(U)
    type(multifab), intent(inout) :: U
    
    integer :: n, lo(U%dim), hi(U%dim), ulo(4), uhi(4), tid
    double precision, pointer, dimension(:,:,:,:) :: up
    
    !$omp parallel private(tid,n,up,lo,hi,ulo,uhi)
    tid = omp_get_thread_num()
    do n=1,nfabs(U)
       up => dataptr(U, n)

       ulo = lbound(up)
       uhi = ubound(up)
       
       lo = tb_get_valid_lo(tid,n)
       hi = tb_get_valid_hi(tid,n)

       call reset_rho_3d(lo,hi,up,ulo(1:3),uhi(1:3))
    end do
    !$omp end parallel

  end subroutine reset_density

  subroutine reset_rho_3d(lo, hi, u, ulo, uhi)
    integer, intent(in) :: lo(3), hi(3), ulo(3), uhi(3)
    double precision, intent(inout) :: u(ulo(1):uhi(1),ulo(2):uhi(2),ulo(3):uhi(3),ncons)

    integer :: i, j, k, n
    double precision :: rho

    do k = lo(3),hi(3)
       do j = lo(2),hi(2)
          do i = lo(1),hi(1)
             rho = 0.d0
             do n=1, nspecies
                rho = rho + U(i,j,k,iry1+n-1)
             end do
             U(i,j,k,irho) = rho
          end do
       end do
    end do

  end subroutine reset_rho_3d

end module variables_module
