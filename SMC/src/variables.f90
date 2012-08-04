module variables

  use bl_types
  use chemistry_module, only : nspecies
  use multifab_module

  implicit none

  integer, save :: irho, imx, imy, imz, iene, iry1
  integer, save :: qrho, qu, qv, qw, qpres, qtemp, qe, qy1, qx1, qh1

  integer, save :: ncons, nprim

contains

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


  subroutine ctoprim(U, Q, ng)
    type(multifab), intent(in   ) :: U
    type(multifab), intent(inout) :: Q
    integer, optional, intent(in) :: ng

    integer :: ngu, ngq, ngto
    integer :: n, lo(U%dim), hi(U%dim)
    double precision, pointer, dimension(:,:,:,:) :: up, qp

    ngu = nghost(U)
    ngq = nghost(Q)

    if (present(ng)) then
       ngto = ng
    else
       ngto = min(ngu, ngq)
    end if

    do n=1,nboxes(Q)
       if ( remote(Q,n) ) cycle

       up => dataptr(U,n)
       qp => dataptr(Q,n)

       lo = lwb(get_box(Q,n))
       hi = upb(get_box(Q,n))

       if (U%dim .eq. 2) then
          call bl_error("2D not supported in variables::ctoprim")
       else
          call ctoprim_3d(lo,hi,up,qp,ngu,ngq,ngto)
       end if
    end do

  end subroutine ctoprim

  subroutine ctoprim_3d(lo, hi, u, q, ngu, ngq, ngto)
    integer, intent(in) :: lo(3), hi(3), ngu, ngq, ngto
    double precision, intent(in ) :: u(lo(1)-ngu:hi(1)+ngu,lo(2)-ngu:hi(2)+ngu,lo(3)-ngu:hi(3)+ngu,ncons)
    double precision, intent(out) :: q(lo(1)-ngq:hi(1)+ngq,lo(2)-ngq:hi(2)+ngq,lo(3)-ngq:hi(3)+ngq,nprim)
    
    integer :: i, j, k, n, iwrk
    double precision :: rho, rhoinv, rwrk, X(nspecies), Y(nspecies), h(nspecies), ei, Tt, Pt

    do k = lo(3)-ngto,hi(3)+ngto
       do j = lo(2)-ngto,hi(2)+ngto
          do i = lo(1)-ngto,hi(1)+ngto
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

end module variables
