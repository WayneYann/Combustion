module kernels
  use probin, only: nx
  use lmc
  implicit none
  double precision, parameter :: eps = 1.d-6
contains

  !
  ! Calculate cell-centered grad pi from nodal pi
  !
  subroutine calc_grad_pi(gp, press, lo, hi, dx)
    integer,          intent(in   ) :: lo, hi
    double precision, intent(  out) :: gp(lo-1:hi+1)
    double precision, intent(in   ) :: press(lo-1:hi+2), dx

    integer :: i

    do i = lo-1, hi+1
       gp(i) = (press(i+1) - press(i)) / dx
    end do
  end subroutine calc_grad_pi

  !
  ! Calculate edge velocities
  !
  subroutine calc_edge_vel(edgevel, vel, lo, hi, bc)
    integer,          intent(in   ) :: lo, hi, bc(2)
    double precision, intent(  out) :: edgevel(lo:hi+1)
    double precision, intent(in   ) :: vel(lo-2:hi+2)

    integer          :: i
    double precision :: slo, shi, slope(lo-1:hi+1)

    call mkslopes(vel, slope, lo, hi, bc)

    do i = lo, hi+1
       slo = vel(i-1) + 0.5d0 * slope(i-1)
       shi = vel(i  ) - 0.5d0 * slope(i  )
       if ( (slo+shi) .gt. eps) then
          edgevel(i) = slo
       else if ( (slo+shi) .lt. -eps) then
          edgevel(i) = shi
       else ! XXX if ( (abs(slo+shi) .le. eps) .or. (slo .le. 0.d0 .and. shi .ge. 0.d0)) then
          edgevel(i) = 0.d0
       end if
       if (i .eq. lo   .and. bc(1) .eq. 1) edgevel(i) = vel(i-1)
       if (i .eq. hi+1 .and. bc(2) .eq. 2) edgevel(i) = slo
    end do
  end subroutine calc_edge_vel

  !
  ! Project.
  !
  subroutine mac_project(macvel,rho,divu,dx,lo,hi,bc)
    integer,          intent(in   ) :: lo, hi, bc(2)
    double precision, intent(inout) :: macvel(lo:hi+1)
    double precision, intent(in   ) :: rho(lo-2:hi+2), divu(lo-1:hi+1), dx

    integer          :: i
    double precision :: a(nx), b(nx), c(nx), r(nx), u(nx), gam(nx), phi(-1:nx+1)

    do i=lo,hi
       u(i+1) = 0.d0
       r(i+1) = (macvel(i+1)-macvel(i))/dx - divu(i)
       a(i+1) =  (2.d0/dx**2)*(1.d0/(rho(i-1)+rho(i)))
       b(i+1) = -(2.d0/dx**2)*(1.d0/(rho(i-1)+rho(i)) + 1.d0/(rho(i)+rho(i+1)))
       c(i+1) =  (2.d0/dx**2)*(1.d0/(rho(i)+rho(i+1)))

       if (i .eq. lo .and. bc(1) .eq. 1) then
          a(i+1) = 0.d0
          b(i+1) = -(2.d0/dx**2)*(1.d0/(rho(i)+rho(i+1)))
          c(i+1) = -b(i+1)
       end if

       if (i .eq. hi .and. bc(2) .eq. 2) then
          a(i+1) = (1.d0/(3.d0*dx**2))*(1.d0/rho(i)) + (2.d0/dx**2)*(1.d0/(rho(i-1)+rho(i)))
          b(i+1) = -(3.d0/dx**2)*(1.d0/rho(i))       - (2.d0/dx**2)*(1.d0/(rho(i-1)+rho(i)))
          c(i+1) = 0.d0
       end if
    end do

    call tridiag(a,b,c,r,u,gam,hi-lo+1)

    phi(lo:hi) = u(lo+1:hi+1)
    if (bc(1) .eq. 1) phi(lo-1) = phi(lo)
    if (bc(2) .eq. 2) phi(hi+1) = -2.d0*phi(hi) + (1.d0/3.d0)*phi(hi-1)

    do i=lo,hi+1
       macvel(i) = macvel(i) - (2.d0/(rho(i-1)+rho(i)))*(phi(i)-phi(i-1))/dx
    end do
  end subroutine mac_project


  !
  ! Scalar advection.
  !
  subroutine scal_aofs(scal,macvel,aofs,divu,dx,lo,hi,bc)
    integer,          intent(in   ) :: lo, hi, bc(2)
    double precision, intent(in   ) :: dx, scal(lo-2:hi+2,0:nscal), macvel(lo:hi+1), divu(lo-1:hi+1)
    double precision, intent(  out) :: aofs(lo:hi,0:nscal)

    double precision :: sedge(lo:hi+1,0:nscal), slope(lo-1:hi+1), slo, shi, Y(nspec), rwrk, hmix
    integer          :: i, n, iwrk, ispec
    logical          :: compute_comp(0:nscal)

    compute_comp = .true.
    compute_comp(Density) = .false.
    compute_comp(RhoRT)   = .false. ! XXX: MWE: should this be RhoH?
    compute_comp(Temp)    = .false.

    do n = 1,nscal
       if (.not. compute_comp(n)) cycle

       call mkslopes(scal(:,n),slope,lo,hi,bc)

       do i=lo,hi+1
          slo = scal(i-1,n) + 0.5d0 * slope(i-1)
          shi = scal(i  ,n) - 0.5d0 * slope(i  )

          if ( macvel(i) .gt. eps) then
             sedge(i,n) = slo
          else if ( macvel(i) .lt. -eps) then
             sedge(i,n) = shi
          else if ( abs(macvel(i)) .le. eps) then
             sedge(i,n) = 0.5d0 * (slo + shi)
          endif

          if (i .eq. lo   .and. bc(1) .eq. 1) sedge(i,n) = scal(i-1,n)
          if (i .eq. hi+1 .and. bc(2) .eq. 2) sedge(i,n) = slo
       end do
    end do

    do i=lo,hi+1
       sedge(i,Density) = 0.d0
       ! compute Rho on edges as sum of (Rho Y_i) on edges,
       do n = 1,nspec
          ispec = FirstSpec-1+n
          sedge(i,Density) = sedge(i,Density) + sedge(i,ispec)
       end do
       ! compute rho.hmix as sum of (H_i.Rho.Y_i)
       if (.not. compute_comp(RhoH) ) then
          ! XXX: MWE: it seems like we want to do this no matter
          ! what..., otherwise sedge(:,Temp) is bogus down below...
          do n = 1,nspec
             ispec = FirstSpec-1+n
             Y(n) = sedge(i,ispec) / sedge(i,Density)
          end do
          call CKHBMS(sedge(i,Temp),Y,IWRK,RWRK,Hmix)
          sedge(i,RhoH) = Hmix * sedge(i,Density)
       end if
       ! XXX: MWE: what about sedge(:,RhoRT)?
    end do

    do n = 1,nscal
       if (n.eq.Temp) then
          do i=lo,hi
             aofs(i,n) = ( macvel(i+1)*sedge(i+1,n) - macvel(i  )*sedge(i  ,n)) / dx
             aofs(i,n) = aofs(i,n) - &
                  (macvel(i+1)  - macvel(i)) * 0.5d0 * ( sedge(i,n) + sedge(i+1,n)) / dx
          end do
       else
          do i=lo,hi
             aofs(i,n) = ( macvel(i+1)*sedge(i+1,n) - macvel(i  )*sedge(i  ,n)) / dx
          end do
       end if

       ! XXX
       ! make these negative here so we can add as source terms later.
       do i=lo,hi
          aofs(i,n) = -aofs(i,n)
       end do
    end do
  end subroutine scal_aofs


  !
  ! Solve XXX.
  !
  subroutine tridiag(a,b,c,r,u,gam,n)
    integer,          intent(in   ) :: n
    double precision, intent(in   ) :: a(n),b(n),c(n),r(n)
    double precision, intent(inout) :: u(n), gam(n)

    integer          :: j
    double precision :: bet

    if (b(1) .eq. 0) stop "CANT HAVE B(1) = ZERO"

    bet  = b(1)
    u(1) = r(1)/bet

    do j = 2,n
       gam(j) = c(j-1)/bet
       bet    = b(j) - a(j)*gam(j)
       if (bet .eq. 0)  stop 'TRIDIAG FAILED'
       u(j)   = (r(j)-a(j)*u(j-1))/bet
    end do

    do j = n-1,1,-1
       u(j) = u(j) - gam(j+1)*u(j+1)
    end do
  end subroutine tridiag

  !
  ! XXX.
  !
  subroutine mkslopes(scal,slope,lo,hi,bc)
    use lmc, only: unlim
    integer,          intent(in   ) :: lo, hi, bc(2)
    double precision, intent(in   ) :: scal(lo-2:hi+2)
    double precision, intent(  out) :: slope(lo-1:hi+1)

    double precision   :: slo,shi,slim,smid
    integer            :: i

    integer, parameter :: cen=1, lim=2, flag=3, fromm=4

    if (unlim .eq. 0) then
       do i=lo-1,hi+1
          shi = 2.0d0*(scal(i+1) - scal(i  ))
          slo = 2.0d0*(scal(i  ) - scal(i-1))
          smid = 0.5d0*(scal(i+1) - scal(i-1))
          if (i .eq. lo .and. bc(1) .eq. 1) then
             ! inflow: value in ghost cell is value at inflow face
             smid = (scal(1)+3.d0*scal(0)-4.d0*scal(-1))/3.d0
          end if
          slim = min(abs(slo),abs(shi))
          slim = min(abs(smid),slim)
          if (slo*shi .lt. 0.d0) then
             slope(i) = 0.d0
          else
             slope(i) = sign(1.d0,smid)*slim
          end if
       end do
    else
       do i=lo-1,hi+1
          slope(i) = 0.5d0*(scal(i+1) - scal(i-1))
       end do
    end if

    if (bc(1) .eq. 1) slope(lo-1) = 0.d0
    if (bc(2) .eq. 2) slope(hi:)  = 0.d0
  end subroutine mkslopes

end module kernels
