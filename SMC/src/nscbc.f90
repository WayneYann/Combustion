module nscbc_module

  use bc_module
  use bl_error_module
  use multifab_module
  use parallel
  
  use chemistry_module
  use derivative_stencil_module, only : stencil_ng, first_deriv_8, first_deriv_6, &
       first_deriv_4, first_deriv_l3, first_deriv_r3, first_deriv_rb, first_deriv_lb
  use physbndry_reg_module
  use probin_module, only : bcx_lo,bcx_hi,bcy_lo,bcy_hi,bcz_lo,bcz_hi, &
       prob_lo_x, prob_lo_y, prob_lo_z, prob_hi_x, prob_hi_y, prob_hi_z
  use smc_bc_module
  use variables_module

  implicit none

  type(physbndry_reg), save :: aux_xlo, aux_xhi, aux_ylo, aux_yhi, aux_zlo, aux_zhi
  integer, parameter :: igamma=1, iWbar=2, icv=3, ics=4, iwdot1=5
  integer, save :: naux

  double precision, parameter :: sigma=0.25d0, Pinfty=1.01325d6
  double precision, save :: Lxdomain, Lydomain, Lzdomain
  double precision, save :: Ma2_xlo, Ma2_xhi, Ma2_ylo, Ma2_yhi, Ma2_zlo, Ma2_zhi

  private

  public :: nscbc, nscbc_init, nscbc_close

contains

  subroutine nscbc_init(la)
    use layout_module
    type(layout), intent(in) :: la

    Lxdomain = prob_hi_x - prob_lo_x
    Lydomain = prob_hi_y - prob_lo_y
    Lzdomain = prob_hi_z - prob_lo_z

    naux = iwdot1 + nspecies - 1

    if (bcx_lo .eq. OUTLET) then
       call physbndry_reg_build(aux_xlo, la, naux, 1, -1, .true.)
    end if

    if (bcy_lo .eq. OUTLET) then
       call physbndry_reg_build(aux_ylo, la, naux, 2, -1, .true.)
    end if

    if (bcz_lo .eq. OUTLET) then
       call physbndry_reg_build(aux_zlo, la, naux, 3, -1, .true.)
    end if

    if (bcx_hi .eq. OUTLET) then
       call physbndry_reg_build(aux_xhi, la, naux, 1, +1, .true.)
    end if

    if (bcy_hi .eq. OUTLET) then
       call physbndry_reg_build(aux_yhi, la, naux, 2, +1, .true.)
    end if

    if (bcz_hi .eq. OUTLET) then
       call physbndry_reg_build(aux_zhi, la, naux, 3, +1, .true.)
    end if
  end subroutine nscbc_init

  subroutine nscbc_close()
    if (aux_xlo%nc .le. 0) then
       call physbndry_reg_destroy(aux_xlo)
    end if

    if (aux_ylo%nc .le. 0) then
       call physbndry_reg_destroy(aux_ylo)
    end if

    if (aux_zlo%nc .le. 0) then
       call physbndry_reg_destroy(aux_zlo)
    end if

    if (aux_xhi%nc .le. 0) then
       call physbndry_reg_destroy(aux_xhi)
    end if

    if (aux_yhi%nc .le. 0) then
       call physbndry_reg_destroy(aux_yhi)
    end if

    if (aux_zhi%nc .le. 0) then
       call physbndry_reg_destroy(aux_zhi)
    end if
  end subroutine nscbc_close

  subroutine nscbc(Q, U, Fdif, rhs, dx)
    type(multifab), intent(in   ) :: Q, U, Fdif
    type(multifab), intent(inout) :: rhs
    double precision, intent(in) :: dx(Q%dim)

    integer :: n, nb, dm, ngq, ngu
    integer :: blo(Q%dim), bhi(Q%dim)
    integer :: dlo(Q%dim), dhi(Q%dim)
    integer ::  lo(Q%dim),  hi(Q%dim)
    integer :: glo(Q%dim), ghi(Q%dim)

    double precision :: proc_Ma2_xlo, proc_Ma2_xhi
    double precision :: proc_Ma2_ylo, proc_Ma2_yhi
    double precision :: proc_Ma2_zlo, proc_Ma2_zhi

    double precision, pointer, dimension(:,:,:,:) :: qp, up, fdp, rhp, auxp

    nb = nboxes(Q)
    dm = Q%dim
    ngq = nghost(Q)
    ngu = nghost(U)

    proc_Ma2_xlo=0.d0 
    proc_Ma2_xhi=0.d0 
    proc_Ma2_ylo=0.d0 
    proc_Ma2_yhi=0.d0 
    proc_Ma2_zlo=0.d0 
    proc_Ma2_zhi=0.d0

    do n=1,nb
       if ( remote(Q,n) ) cycle

       if (isValid(aux_xlo,n)) then
          lo = lwb(get_box(Q,n))
          hi = upb(get_box(Q,n))

          glo = lwb(get_box(aux_xlo%data,n))
          ghi = upb(get_box(aux_xlo%data,n))

          qp => dataptr(Q, n)
          auxp => dataptr(aux_xlo%data,n)

          call compute_aux(lo,hi,ngq,qp,auxp,glo,ghi,proc_Ma2_xlo)
       end if

       if (isValid(aux_xhi,n)) then
          lo = lwb(get_box(Q,n))
          hi = upb(get_box(Q,n))

          glo = lwb(get_box(aux_xhi%data,n))
          ghi = upb(get_box(aux_xhi%data,n))

          qp => dataptr(Q, n)
          auxp => dataptr(aux_xhi%data,n)

          call compute_aux(lo,hi,ngq,qp,auxp,glo,ghi,proc_Ma2_xhi)
       end if

    end do

    if (bcx_lo .eq. OUTLET) then
       call parallel_reduce(Ma2_xlo, proc_Ma2_xlo, MPI_MAX)
    end if

    if (bcx_hi .eq. OUTLET) then
       call parallel_reduce(Ma2_xhi, proc_Ma2_xhi, MPI_MAX)
    end if

    if (bcy_lo .eq. OUTLET) then
       call parallel_reduce(Ma2_ylo, proc_Ma2_ylo, MPI_MAX)
    end if

    if (bcy_hi .eq. OUTLET) then
       call parallel_reduce(Ma2_yhi, proc_Ma2_yhi, MPI_MAX)
    end if

    if (bcz_lo .eq. OUTLET) then
       call parallel_reduce(Ma2_zlo, proc_Ma2_zlo, MPI_MAX)
    end if

    if (bcz_hi .eq. OUTLET) then
       call parallel_reduce(Ma2_zhi, proc_Ma2_zhi, MPI_MAX)
    end if

    do n=1,nb
       if ( remote(Q,n) ) cycle

       if (isValid(aux_xlo,n)) then
          
          qp  => dataptr(Q,n)
          up  => dataptr(U,n)
          rhp => dataptr(rhs,n)
          fdp => dataptr(Fdif, n)
          auxp => dataptr(aux_xlo%data,n)

          lo = lwb(get_box(Q,n))
          hi = upb(get_box(Q,n))

          call get_data_lo_hi(n,dlo,dhi)
          call get_boxbc(n,blo,bhi)

          call nscbc_xlo(lo,hi,ngq,ngu,dx,qp,up,fdp,rhp,auxp(:,lo(1),:,:),dlo,dhi,blo,bhi)

       end if

       if (isValid(aux_xhi,n)) then
          
          qp  => dataptr(Q,n)
          up  => dataptr(U,n)
          rhp => dataptr(rhs,n)
          fdp => dataptr(Fdif, n)
          auxp => dataptr(aux_xhi%data,n)

          lo = lwb(get_box(Q,n))
          hi = upb(get_box(Q,n))

          call get_data_lo_hi(n,dlo,dhi)
          call get_boxbc(n,blo,bhi)

          call nscbc_xhi(lo,hi,ngq,ngu,dx,qp,up,fdp,rhp,auxp(:,hi(1),:,:),dlo,dhi,blo,bhi)

       end if

    end do

  end subroutine nscbc


  subroutine compute_aux(lo,hi,ng,Q,G,glo,ghi,mach2)
    integer, intent(in) :: ng
    integer, dimension(3), intent(in) :: lo, hi, glo, ghi
    double precision, intent(in ) :: Q(-ng+lo(1):hi(1)+ng,-ng+lo(2):hi(2)+ng,-ng+lo(3):hi(3)+ng,nprim)
    double precision, intent(out) :: G(naux,glo(1):hi(1),glo(2):ghi(2),glo(3):ghi(3))
    double precision, intent(inout) :: mach2

    integer :: i, j, k, iwrk
    double precision :: Tt, rwrk, cv, cp, gamma, Wbar, vel2, cs2
    double precision :: Yt(nspecies), wdot(nspecies)

    double precision, parameter :: Ru = 8.31451d7

    do k=glo(3),ghi(3)
       do j=glo(2),ghi(2)
          do i=glo(1),ghi(1)

             Tt = q(i,j,k,qtemp)
             Yt = q(i,j,k,qy1:qy1+nspecies-1)

             call ckcvbs(Tt, Yt, iwrk, rwrk, cv)
             call ckmmwy(Yt, iwrk, rwrk, Wbar)

             cp = cv + Ru/Wbar
             gamma = cp / cv
             cs2 = gamma*q(i,j,k,qpres)/q(i,j,k,qrho)

             vel2 = q(i,j,k,qu)**2 + q(i,j,k,qv)**2 + q(i,j,k,qw)**2
             mach2 = max(mach2, vel2/cs2)

             G(igamma,i,j,k) = gamma
             G(iWbar ,i,j,k) = Wbar
             G(icv   ,i,j,k) = cv
             G(ics   ,i,j,k) = sqrt(cs2)

             call ckwyr(q(i,j,k,qrho), Tt, Yt, iwrk, rwrk, wdot)
             G(iwdot1:,i,j,k) = wdot*molecular_weight

          end do
       end do
    end do

  end subroutine compute_aux


  subroutine nscbc_xlo(lo,hi,ngq,ngu,dx,Q,U,fd,rhs,aux,dlo,dhi,bclo,bchi)
    integer, intent(in) :: lo(3), hi(3), ngq, ngu, dlo(3), dhi(3), bclo(3), bchi(3)
    double precision,intent(in   )::dx(3)
    double precision,intent(in   )::Q  (-ngq+lo(1):hi(1)+ngq,-ngq+lo(2):hi(2)+ngq,-ngq+lo(3):hi(3)+ngq,nprim)
    double precision,intent(in   )::U  (-ngu+lo(1):hi(1)+ngu,-ngu+lo(2):hi(2)+ngu,-ngu+lo(3):hi(3)+ngu,ncons)
    double precision,intent(in   )::fd (lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),ncons)
    double precision,intent(inout)::rhs(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),ncons)
    double precision,intent(in   )::aux(naux,lo(2):hi(2),lo(3):hi(3))

    integer :: i,j,k,n
    double precision :: dxinv(3)
    double precision :: un, dpdn, dudn(3), drhodn, dYdn(nspecies)
    double precision :: L(5+nspecies)
    double precision :: S_p, S_Y(nspecies), d_u, d_v, d_w, d_p, d_Y(nspecies)
    double precision :: scratch, a_mulitD

    double precision, dimension(lo(2):hi(2),lo(3):hi(3)) :: &
         dpdy, dpdz, dudy, dudz, dvdy, dwdz
    
    i = lo(1)

    do n=1,3
       dxinv(n) = 1.0d0 / dx(n)
    end do

    call comp_trans_deriv_x(i,lo,hi,ngq,Q,dxinv,dlo,dhi,dpdy,dudy,dvdy,dpdz,dudz,dwdz)

    ! lo-x boundary
    if (bclo(1) .eq. OUTLET) then

       do k=lo(3),hi(3)
          do j=lo(2),hi(2)
             
             un = q(i,j,k,qu)

             drhodn  = dxinv(1)*first_deriv_rb(q(i:i+3,j,k,qrho))
             dudn(2) = dxinv(1)*first_deriv_rb(q(i:i+3,j,k,qv))
             dudn(3) = dxinv(1)*first_deriv_rb(q(i:i+3,j,k,qw))
             dudn(1) = dxinv(1)*first_deriv_rb(q(i:i+3,j,k,qu))
             dpdn    = dxinv(1)*first_deriv_rb(q(i:i+3,j,k,qpres))
             do n=1,nspecies
                dYdn(n) = dxinv(1)*first_deriv_lb(q(i:i+3,j,k,qy1+n-1))
             end do

             ! Simple 1D LODI 
             L(1) = (un-aux(ics,j,k))*0.5d0*(dpdn-q(i,j,k,qrho)*aux(ics,j,k)*dudn(1))
             L(2) = un*(drhodn-dpdn/aux(ics,j,k)**2)
             L(3) = un*dudn(2)
             L(4) = un*dudn(3)
             L(5) = sigma*aux(ics,j,k)*(1.d0-Ma2_xlo)/(2.d0*Lxdomain)*(q(i,j,k,qpres)-Pinfty)
             L(6:) = un*dYdn

             ! multi-D effects
             a_mulitD = min(0.99d0, (1.0d0-abs(un)/aux(ics,j,k)))
             L(5) = L(5) - a_mulitD*0.5d0*(q(i,j,k,qv)*dpdy(j,k)+q(i,j,k,qw)*dpdz(j,k) &
                  + aux(igamma,j,k)*q(i,j,k,qpres)*(dvdy(j,k)+dwdz(j,k)) &
                  + q(i,j,k,qrho)*aux(ics,j,k)*(q(i,j,k,qv)*dudy(j,k)+q(i,j,k,qw)*dudz(j,k)))
             
             ! viscous and reaction effects
             d_u = fd(i,j,k,imx) / q(i,j,k,qrho)
             d_v = fd(i,j,k,imy) / q(i,j,k,qrho)
             d_w = fd(i,j,k,imz) / q(i,j,k,qrho)

             S_Y = 0.d0  
             S_p = 0.d0
             d_p = 0.d0
             do n=1,nspecies
                scratch = q(i,j,k,qh1+n-1)*(1.d0-aux(igamma,j,k)) + &
                     aux(igamma,j,k)*q(i,j,k,qpres)*aux(iWbar,j,k)/(q(i,j,k,qrho)*molecular_weight(n))
                S_p = S_p + scratch*aux(iwdot1+n-1,j,k)
                d_p = d_p + scratch*fd(i,j,k,iry1+n-1)
                d_Y(n) = fd(i,j,k,iry1+n-1) / q(i,j,k,qrho)
             end do
             d_p = d_p + (aux(igamma,j,k)-1.d0)*(fd(i,j,k,iene) &
                  - d_u*U(i,j,k,imx) - d_v*U(i,j,k,imy) - d_w*U(i,j,k,imz))
             
             L(5) = L(5) + 0.5d0*(S_p + d_p + q(i,j,k,qrho)*aux(ics,j,k)*d_u)

             if (q(i,j,k,qu) > 0.d0) then
                L(2) = -S_p/aux(ics,j,k)**2
                L(3) = 0.d0
                L(4) = 0.d0
                L(5) = 0.5d0*S_p
                L(6:) = S_Y      
             end if

          end do
       end do

    else if (bclo(1).ne.PERIODIC .and. bchi(1).ne.INTERIOR) then
       call bl_error("Unknown boundary")
    end if

  end subroutine nscbc_xlo


  subroutine nscbc_xhi(lo,hi,ngq,ngu,dx,Q,U,fd,rhs,aux,dlo,dhi,bclo,bchi)
    integer, intent(in) :: lo(3), hi(3), ngq, ngu, dlo(3), dhi(3), bclo(3), bchi(3)
    double precision,intent(in   )::dx(3)
    double precision,intent(in   )::Q  (-ngq+lo(1):hi(1)+ngq,-ngq+lo(2):hi(2)+ngq,-ngq+lo(3):hi(3)+ngq,nprim)
    double precision,intent(in   )::U  (-ngu+lo(1):hi(1)+ngu,-ngu+lo(2):hi(2)+ngu,-ngu+lo(3):hi(3)+ngu,ncons)
    double precision,intent(in   )::fd (lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),ncons)
    double precision,intent(inout)::rhs(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),ncons)
    double precision,intent(in   )::aux(naux,lo(2):hi(2),lo(3):hi(3))

    integer :: i,j,k,n
    double precision :: dxinv(3)
    double precision :: un, dpdn, dudn(3), drhodn, dYdn(nspecies)
    double precision :: L(5+nspecies)
    double precision :: S_p, S_Y(nspecies), d_u, d_v, d_w, d_p, d_Y(nspecies)
    double precision :: scratch, a_mulitD

    double precision, dimension(lo(2):hi(2),lo(3):hi(3)) :: &
         dpdy, dpdz, dudy, dudz, dvdy, dwdz
    
    i = hi(1)

    do n=1,3
       dxinv(n) = 1.0d0 / dx(n)
    end do

    call comp_trans_deriv_x(i,lo,hi,ngq,Q,dxinv,dlo,dhi,dpdy,dudy,dvdy,dpdz,dudz,dwdz)

    ! hi-x boundary
    if (bchi(1) .eq. OUTLET) then

       do k=lo(3),hi(3)
          do j=lo(2),hi(2)
             
             un = q(i,j,k,qu)

             drhodn  = dxinv(1)*first_deriv_lb(q(i-3:i,j,k,qrho))
             dudn(2) = dxinv(1)*first_deriv_lb(q(i-3:i,j,k,qv))
             dudn(3) = dxinv(1)*first_deriv_lb(q(i-3:i,j,k,qw))
             dudn(1) = dxinv(1)*first_deriv_lb(q(i-3:i,j,k,qu))
             dpdn    = dxinv(1)*first_deriv_lb(q(i-3:i,j,k,qpres))
             do n=1,nspecies
                dYdn(n) = dxinv(1)*first_deriv_lb(q(i-3:i,j,k,qy1+n-1))
             end do

             ! Simple 1D LODI 
             L(1) = sigma*aux(ics,j,k)*(1.d0-Ma2_xhi)/(2.d0*Lxdomain)*(q(i,j,k,qpres)-Pinfty)
             L(2) = un*(drhodn-dpdn/aux(ics,j,k)**2)
             L(3) = un*dudn(2)
             L(4) = un*dudn(3)
             L(5) = (un+aux(ics,j,k))*0.5d0*(dpdn+q(i,j,k,qrho)*aux(ics,j,k)*dudn(1))
             L(6:) = un*dYdn

             ! multi-D effects
             a_mulitD = min(0.99d0, (1.0d0-abs(un)/aux(ics,j,k)))
             L(1) = L(1) - a_mulitD*0.5d0*(q(i,j,k,qv)*dpdy(j,k)+q(i,j,k,qw)*dpdz(j,k) &
                  + aux(igamma,j,k)*q(i,j,k,qpres)*(dvdy(j,k)+dwdz(j,k)) &
                  - q(i,j,k,qrho)*aux(ics,j,k)*(q(i,j,k,qv)*dudy(j,k)+q(i,j,k,qw)*dudz(j,k)))
             
             ! viscous and reaction effects
             d_u = fd(i,j,k,imx) / q(i,j,k,qrho)
             d_v = fd(i,j,k,imy) / q(i,j,k,qrho)
             d_w = fd(i,j,k,imz) / q(i,j,k,qrho)

             S_Y = 0.d0  
             S_p = 0.d0
             d_p = 0.d0
             do n=1,nspecies
                scratch = q(i,j,k,qh1+n-1)*(1.d0-aux(igamma,j,k)) + &
                     aux(igamma,j,k)*q(i,j,k,qpres)*aux(iWbar,j,k)/(q(i,j,k,qrho)*molecular_weight(n))
                S_p = S_p + scratch*aux(iwdot1+n-1,j,k)
                d_p = d_p + scratch*fd(i,j,k,iry1+n-1)
                d_Y(n) = fd(i,j,k,iry1+n-1) / q(i,j,k,qrho)
             end do
             d_p = d_p + (aux(igamma,j,k)-1.d0)*(fd(i,j,k,iene) &
                  - d_u*U(i,j,k,imx) - d_v*U(i,j,k,imy) - d_w*U(i,j,k,imz))
             
             L(1) = L(1) + 0.5d0*(S_p + d_p - q(i,j,k,qrho)*aux(ics,j,k)*d_u)

             if (q(i,j,k,qu) < 0.d0) then
                L(1) = 0.5d0*S_p
                L(2) = -S_p/aux(ics,j,k)**2
                L(3) = 0.d0
                L(4) = 0.d0
                L(6:) = S_Y      
             end if

          end do
       end do

    else if (bchi(1).ne.PERIODIC .and. bchi(1).ne.INTERIOR) then
       call bl_error("Unknown boundary")
    end if

  end subroutine nscbc_xhi


  subroutine comp_trans_deriv_x(i,lo,hi,ngq,Q,dxinv,dlo,dhi,dpdy,dudy,dvdy,dpdz,dudz,dwdz)
    integer, intent(in) :: i, lo(3), hi(3), ngq,dlo(3),dhi(3)
    double precision, intent(in) :: dxinv(3)
    double precision, intent(in) :: Q(-ngq+lo(1):hi(1)+ngq,-ngq+lo(2):hi(2)+ngq,-ngq+lo(3):hi(3)+ngq,nprim)
    double precision, intent(out):: dpdy(lo(2):hi(2),lo(3):hi(3))
    double precision, intent(out):: dudy(lo(2):hi(2),lo(3):hi(3))
    double precision, intent(out):: dvdy(lo(2):hi(2),lo(3):hi(3))
    double precision, intent(out):: dpdz(lo(2):hi(2),lo(3):hi(3))
    double precision, intent(out):: dudz(lo(2):hi(2),lo(3):hi(3))
    double precision, intent(out):: dwdz(lo(2):hi(2),lo(3):hi(3))

    integer :: j, k, slo(3), shi(3)

    slo = dlo + stencil_ng
    shi = dhi - stencil_ng

    ! d()/dy
    do k=lo(3),hi(3)
       
       do j=slo(2),shi(2)
          dudy(j,k) = dxinv(2)*first_deriv_8(q(i,j-4:j+4,k,qu))
          dvdy(j,k) = dxinv(2)*first_deriv_8(q(i,j-4:j+4,k,qv))
          dpdy(j,k) = dxinv(2)*first_deriv_8(q(i,j-4:j+4,k,qpres))
       end do

       ! lo-y boundary
       if (dlo(2) .eq. lo(2)) then
          j = lo(2)
          ! use completely right-biased stencil
          dudy(j,k) = dxinv(2)*first_deriv_rb(q(i,j:j+3,k,qu))
          dvdy(j,k) = dxinv(2)*first_deriv_rb(q(i,j:j+3,k,qv))
          dpdy(j,k) = dxinv(2)*first_deriv_rb(q(i,j:j+3,k,qpres))          

          j = lo(2)+1
          ! use 3rd-order slightly right-biased stencil
          dudy(j,k) = dxinv(2)*first_deriv_r3(q(i,j-1:j+2,k,qu))
          dvdy(j,k) = dxinv(2)*first_deriv_r3(q(i,j-1:j+2,k,qv))
          dpdy(j,k) = dxinv(2)*first_deriv_r3(q(i,j-1:j+2,k,qpres))          

          j = lo(2)+2
          ! use 4th-order stencil
          dudy(j,k) = dxinv(2)*first_deriv_4(q(i,j-2:j+2,k,qu))
          dvdy(j,k) = dxinv(2)*first_deriv_4(q(i,j-2:j+2,k,qv))
          dpdy(j,k) = dxinv(2)*first_deriv_4(q(i,j-2:j+2,k,qpres))          

          j = lo(2)+3
          ! use 6th-order stencil
          dudy(j,k) = dxinv(2)*first_deriv_6(q(i,j-3:j+3,k,qu))
          dvdy(j,k) = dxinv(2)*first_deriv_6(q(i,j-3:j+3,k,qv))
          dpdy(j,k) = dxinv(2)*first_deriv_6(q(i,j-3:j+3,k,qpres))          
       end if

       ! hi-y bondary
       if (dhi(2) .eq. hi(2)) then
          j = hi(2)-3
          ! use 6th-order stencil
          dudy(j,k) = dxinv(2)*first_deriv_6(q(i,j-3:j+3,k,qu))
          dvdy(j,k) = dxinv(2)*first_deriv_6(q(i,j-3:j+3,k,qv))
          dpdy(j,k) = dxinv(2)*first_deriv_6(q(i,j-3:j+3,k,qpres))          

          j = hi(2)-2
          ! use 4th-order stencil
          dudy(j,k) = dxinv(2)*first_deriv_4(q(i,j-2:j+2,k,qu))
          dvdy(j,k) = dxinv(2)*first_deriv_4(q(i,j-2:j+2,k,qv))
          dpdy(j,k) = dxinv(2)*first_deriv_4(q(i,j-2:j+2,k,qpres))          
          
          j = hi(2)-1
          ! use 3rd-order slightly left-biased stencil
          dudy(j,k) = dxinv(2)*first_deriv_l3(q(i,j-2:j+1,k,qu))
          dvdy(j,k) = dxinv(2)*first_deriv_l3(q(i,j-2:j+1,k,qv))
          dpdy(j,k) = dxinv(2)*first_deriv_l3(q(i,j-2:j+1,k,qpres))          

          j = hi(2)
          ! use completely left-biased stencil
          dudy(j,k) = dxinv(2)*first_deriv_lb(q(i,j-3:j,k,qu))
          dvdy(j,k) = dxinv(2)*first_deriv_lb(q(i,j-3:j,k,qv))
          dpdy(j,k) = dxinv(2)*first_deriv_lb(q(i,j-3:j,k,qpres))  
       end if

    end do

    ! d()/dz
    do k=slo(3),shi(3)
       do j=lo(2),hi(2)
          dudz(j,k) = dxinv(3)*first_deriv_8(q(i,j,k-4:k+4,qu))
          dwdz(j,k) = dxinv(3)*first_deriv_8(q(i,j,k-4:k+4,qw))
          dpdz(j,k) = dxinv(3)*first_deriv_8(q(i,j,k-4:k+4,qpres))
       end do
    end do

    ! lo-z boundary
    if (dlo(3) .eq. lo(3)) then
       k = lo(3)
       ! use completely right-biased stencil
       do j=lo(2),hi(2)
          dudz(j,k) = dxinv(3)*first_deriv_rb(q(i,j,k:k+3,qu))
          dwdz(j,k) = dxinv(3)*first_deriv_rb(q(i,j,k:k+3,qw))
          dpdz(j,k) = dxinv(3)*first_deriv_rb(q(i,j,k:k+3,qpres))
       end do

       k = lo(3)+1
       ! use 3rd-order slightly right-biased stencil
       do j=lo(2),hi(2)
          dudz(j,k) = dxinv(3)*first_deriv_r3(q(i,j,k-1:k+2,qu))
          dwdz(j,k) = dxinv(3)*first_deriv_r3(q(i,j,k-1:k+2,qw))
          dpdz(j,k) = dxinv(3)*first_deriv_r3(q(i,j,k-1:k+2,qpres))
       end do

       k = lo(3)+2
       ! use 4th-order stencil
       do j=lo(2),hi(2)
          dudz(j,k) = dxinv(3)*first_deriv_4(q(i,j,k-2:k+2,qu))
          dwdz(j,k) = dxinv(3)*first_deriv_4(q(i,j,k-2:k+2,qw))
          dpdz(j,k) = dxinv(3)*first_deriv_4(q(i,j,k-2:k+2,qpres))
       end do

       k = lo(3)+3
       ! use 6th-order stencil
       do j=lo(2),hi(2)
          dudz(j,k) = dxinv(3)*first_deriv_6(q(i,j,k-3:k+3,qu))
          dwdz(j,k) = dxinv(3)*first_deriv_6(q(i,j,k-3:k+3,qw))
          dpdz(j,k) = dxinv(3)*first_deriv_6(q(i,j,k-3:k+3,qpres))
       end do
    end if
    
    ! hi-z boundary
    if (dlo(3) .eq. lo(3)) then
       k = hi(3)-3
       ! use 6th-order stencil
       do j=lo(2),hi(2)
          dudz(j,k) = dxinv(3)*first_deriv_6(q(i,j,k-3:k+3,qu))
          dwdz(j,k) = dxinv(3)*first_deriv_6(q(i,j,k-3:k+3,qw))
          dpdz(j,k) = dxinv(3)*first_deriv_6(q(i,j,k-3:k+3,qpres))
       end do

       k = hi(3)-2
       ! use 4th-order stencil
       do j=lo(2),hi(2)
          dudz(j,k) = dxinv(3)*first_deriv_4(q(i,j,k-2:k+2,qu))
          dwdz(j,k) = dxinv(3)*first_deriv_4(q(i,j,k-2:k+2,qw))
          dpdz(j,k) = dxinv(3)*first_deriv_4(q(i,j,k-2:k+2,qpres))
       end do

       k = hi(3)-1
       ! use 3rd-order slightly left-biased stencil
       do j=lo(2),hi(2)
          dudz(j,k) = dxinv(3)*first_deriv_l3(q(i,j,k-2:k+1,qu))
          dwdz(j,k) = dxinv(3)*first_deriv_l3(q(i,j,k-2:k+1,qw))
          dpdz(j,k) = dxinv(3)*first_deriv_l3(q(i,j,k-2:k+1,qpres))
       end do

       k = hi(3)
       ! use completely left-biased stencil
       do j=lo(2),hi(2)
          dudz(j,k) = dxinv(3)*first_deriv_lb(q(i,j,k-3:k,qu))
          dwdz(j,k) = dxinv(3)*first_deriv_lb(q(i,j,k-3:k,qw))
          dpdz(j,k) = dxinv(3)*first_deriv_lb(q(i,j,k-3:k,qpres))
       end do
    end if

  end subroutine comp_trans_deriv_x

end module nscbc_module
