module wchem

  implicit none

  double precision, allocatable, save :: Jac(:,:), A(:,:)
  integer, allocatable, save :: ipvt(:)
  integer, save :: nstep
  
  public :: bechem

contains

  subroutine bechem(rho0, Y0, rho, YT, dt, g, ierr, always_new_J)
    integer, intent(in) :: g
    double precision, intent(in   ) :: rho0, rho, dt
    double precision, intent(in   ) :: Y0(nspecies+1)
    double precision, intent(inout) :: YT(nspecies+1)
    integer, intent(out), optional :: ierr
    logical, intent(in), optional :: always_new_J

    logical :: new_J
    integer :: iwrk, iter, n, info, age
    double precision :: rwrk, rhoinv, cv, rmax, rmin
    double precision, dimension(nspecies) :: uk
    double precision, dimension(nspecies+1) :: YT_init, r, dYTdt
    integer, parameter :: J_int = 5
    logical, save :: compute_new_A

    new_J = .false.
    if (present(always_new_J)) then
       new_J = always_new_J
    end if

    if (g .eq. 1) compute_new_A = .true.

    if (.not. allocated(A)) then
       allocate(Jac(nspecies+1,nspecies+1))
       allocate(A(nspecies+1,nspecies+1))
       allocate(ipvt(nspecies+1))
    end if

    rhoinv = 1.d0/rho

    YT_init = YT
    YT = Y0

    age = 0
    rmin = 1.d50

    do iter = 0, 100

       call vckwyr(1, rho, YT(nspecies+1), YT(1), iwrk, rwrk, dYTdt)
       call ckums(YT(nspecies+1), iwrk, rwrk, uk)
       call ckcvbs(YT(nspecies+1),YT(1),iwrk,rwrk,cv)       
       dYTdt(nspecies+1) = 0.d0
       do n=1,nspecies
          dYTdt(n) = dYTdt(n) * molecular_weight(n) * rhoinv
          dYTdt(nspecies+1) = dYTdt(nspecies+1) - dYTdt(n)*uk(n)
       end do
       dYTdt(nspecies+1) = dYTdt(nspecies+1) / cv

       r = YT - YT_init - dt * dYTdt

       rmax = maxval(abs(r(1:nspecies)))
       if (rmax .le. 1.d-14) then 
          exit 
       endif
       
       if (      new_J          &
            .or. compute_new_A  &
            .or. age.eq.J_int   &
            .or. rmax.ge.rmin ) then
          call LUA(rho, YT, dt)
          age = 0
          compute_new_A = .false.
       else
          age = age + 1
       end if

       rmin = min(rmin,rmax)

       call dgesl(A, nspecies+1, nspecies+1, ipvt, r, 0)

       YT = YT - r
       if ( maxval(YT(1:nspecies)) .gt. 1.d0  .or.  &
            minval(YT(1:nspecies)) .lt. -1.d-8 ) then
          compute_new_A = .true.
       end if

    end do

    if (present(ierr)) then
       if (iter .gt. 100) then
          ierr = iter
          YT = YT_init
       else
          ierr = 0
       end if
    end if

  contains

    subroutine LUA(rho, YT, dt)
      double precision, intent(in) :: rho, YT(nspecies+1), dt
      integer :: i, j, iwrk, info
      double precision :: rwrk
      double precision, dimension(nspecies) :: C
      integer, parameter :: consP = 0

      call ckytcr(rho, YT(nspecies+1), YT(1), iwrk, rwrk, C)
      call DWDOT(Jac, C, YT(nspecies+1), consP)

      do j=1,nspecies
         do i=1,nspecies
            Jac(i,j) = Jac(i,j) * molecular_weight(i) * inv_mwt(j)
         end do
         i=nspecies+1
         Jac(i,j) = Jac(i,j) * inv_mwt(j) * rho
      end do
      !
      j = nspecies+1
      do i=1,nspecies
         Jac(i,j) = Jac(i,j) * molecular_weight(i) * rhoinv
      end do
      
      A = -dt*Jac
      do i=1,nspecies+1
         A(i,i) = 1.d0 + A(i,i)
      end do
      
      call dgefa(A, nspecies+1, nspecies+1, ipvt, info)
      
    end subroutine LUA

  end subroutine beburn

end module wchem
