module variables

  use bl_types
  use chemistry_module, only : nspecies
  use multifab_module

  implicit none

  integer, save :: irho, imx, imy, imz, iene, iry1
  integer, save :: qrho, qu, qv, qw, qpres, qtemp, qy1, qx1, qh1

  ! the total number of plot components
  integer, save :: n_plot_comps = 0
  integer, save :: icomp_rho, icomp_vel, icomp_pres, icomp_temp, icomp_Y, icomp_X, icomp_h

  integer, save :: ncons, nprim

contains

  function get_next_plot_index(num) result (next)

    ! return the next starting index for a plotfile quantity,
    ! and increment the counter of plotfile quantities by num
    integer :: num, next

    next = n_plot_comps + 1
    n_plot_comps = n_plot_comps + num

    return
  end function get_next_plot_index

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
       ncons = 5
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
       qy1 = 7
    else
       qw = -1
       qpres = 4
       qtemp = 5
       qy1 = 6
    end if
    qx1 = qy1 + nspecies
    qh1 = qx1 + nspecies

    nprim = qh1-1 + nspecies
    
  end subroutine init_variables


  subroutine init_plot_variables()

    use probin_module, only : dm_in, plot_X, plot_Y, plot_h

    icomp_rho  = get_next_plot_index(1)
    icomp_vel  = get_next_plot_index(dm_in)
    icomp_pres = get_next_plot_index(1)
    icomp_temp = get_next_plot_index(1)

    if (plot_Y) then
       icomp_Y = get_next_plot_index(nspecies)
    end if

    if (plot_X) then
       icomp_X = get_next_plot_index(nspecies)
    end if

    if (plot_h) then
       icomp_h = get_next_plot_index(nspecies)
    end if

  end subroutine init_plot_variables


  subroutine ctoprim(U, Q, ng)
    type(multifab), intent(in   ) :: U
    type(multifab), intent(inout) :: Q
    integer, optional, intent(in) :: ng

    integer :: ngu, ngto
    integer :: n, lo(U%dim), hi(U%dim)
    double precision, pointer, dimension(:,:,:,:) :: up, qp

    ngu = nghost(U)

    if (present(ng)) then
       ngto = ng
    else
       ngto = ngu
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
          call ctoprim_3d(lo,hi,up,qp,ngu,ngto)
       end if
    end do

  end subroutine ctoprim

  subroutine ctoprim_3d(lo, hi, u, q, ngu, ngto)
    integer, intent(in) :: lo(3), hi(3), ngu, ngto
    double precision, intent(in ) :: u(lo(1)-ngu:hi(1)+ngu,lo(2)-ngu:hi(2)+ngu,lo(3)-ngu:hi(3)+ngu,ncons)
    double precision, intent(out) :: q(lo(1)-ngu:hi(1)+ngu,lo(2)-ngu:hi(2)+ngu,lo(3)-ngu:hi(3)+ngu,nprim)
    
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
                q(i,j,k,qy1+n-1) = u(i,j,k,iry1+n-1) * rhoinv
                Y(n) = q(i,j,k,qy1+n-1)
             end do

             call ckytx(Y, iwrk, rwrk, X)

             do n=1,nspecies
                q(i,j,k,qx1+n-1) = X(n)
             end do

             ei = rhoinv*u(i,j,k,iene) - 0.5d0*(q(i,j,k,qu)**2+q(i,j,k,qv)**2+q(i,j,k,qw)**2)
             call feeytt(ei, Y, iwrk, rwrk, Tt)
             q(i,j,k,qtemp) = Tt



             if (Tt < 0.d0) then
                print *, 'xxxxx ', i, j, k, rho, Tt, ei, q(i,j,k,qu), q(i,j,k,qv), q(i,j,k,qw)
                print *, '      ', X, Y
                call flush()
             end if



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
