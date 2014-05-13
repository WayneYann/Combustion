module burner_module

  use chemistry_module, only : nspecies, spec_names, molecular_weight, inv_mwt
  use meth_params_module, only : use_vode

  implicit none

  double precision, allocatable, save :: Jac(:,:), A(:,:)
  integer, allocatable, save :: ipvt(:)
  integer, save :: nstep
  !$omp threadprivate(Jac,A,ipvt,nstep) 

  private

  public :: burn, compute_rhodYdt, splitburn, beburn

contains

  subroutine burn(np, rho, YT, dt, force_new_J, oerr)
    use meth_params_module, only : use_vode
    integer, intent(in) :: np
    double precision, intent(in   ) :: rho(np), dt
    double precision, intent(inout) :: YT(nspecies+1,np)
    logical, intent(in) :: force_new_J
    integer, intent(out), optional :: oerr

    if (use_vode) then
       call burn_vode(np, rho, YT, dt, force_new_J, oerr)
    else
       call burn_bdf(np, rho, YT, dt, force_new_J, oerr)
    end if

  end subroutine burn


  subroutine burn_vode(np, rho, YT, dt, force_new_J, oerr)
    use vode_module, only : verbose, itol, rtol, atol, vode_MF=>MF, always_new_j, &
         voderwork, vodeiwork, lvoderwork, lvodeiwork, voderpar, vodeipar
    integer, intent(in) :: np
    double precision, intent(in   ) :: rho(np), dt
    double precision, intent(inout) :: YT(nspecies+1,np)
    logical, intent(in) :: force_new_J
    integer, intent(out), optional :: oerr

    external f_jac, f_rhs, dvode

    ! vode stuff
    integer, parameter :: itask=1, iopt=1
    integer :: MF, istate, ifail

    double precision :: time
    integer :: g

    if (force_new_J) then
       call setfirst(.true.)
    end if

    do g = 1, np

       voderpar(1) = rho(g)

       istate = 1
       time = 0.d0

       if (always_new_j) call setfirst(.true.)

       MF = vode_MF  ! vode might change its sign!
       call dvode(f_rhs, nspecies+1, YT(:,g), time, dt, itol, rtol, atol, itask, &
            istate, iopt, voderwork, lvoderwork, vodeiwork, lvodeiwork, &
            f_jac, MF, voderpar, vodeipar)

       nstep = vodeiwork(11)

       if (verbose .ge. 1) then
          write(6,*) '......dvode done:'
          write(6,*) ' last successful step size = ',voderwork(11)
          write(6,*) '          next step to try = ',voderwork(12)
          write(6,*) '   integrated time reached = ',voderwork(13)
          write(6,*) '      number of time steps = ',vodeiwork(11)
          write(6,*) '              number of fs = ',vodeiwork(12)
          write(6,*) '              number of Js = ',vodeiwork(13)
          write(6,*) '    method order last used = ',vodeiwork(14)
          write(6,*) '   method order to be used = ',vodeiwork(15)
          write(6,*) '            number of LUDs = ',vodeiwork(19)
          write(6,*) ' number of Newton iterations ',vodeiwork(20)
          write(6,*) ' number of Newton failures = ',vodeiwork(21)
          if (ISTATE.eq.-4 .or. ISTATE.eq.-5) then
             ifail = vodeiwork(16)
             if (ifail .eq. nspecies+1) then
                write(6,*) '   T has the largest error'
             else
                write(6,*) '   spec with largest error = ', trim(spec_names(ifail))
             end if
             call flush(6)
          end if
       end if

       if (istate < 0) then
          print *, 'chemsolv: VODE failed'
          print *, 'istate = ', istate, ' time =', time
          if (present(oerr)) then
             oerr = 1
             return
          else
             call bl_error("ERROR in burn: VODE failed")
          end if
       end if

    end do

    if (present(oerr)) oerr = 0

  end subroutine burn_vode


  subroutine burn_bdf(np, rho_in, YT, dt, force_new_J, oerr)
    use bdf
    use bdf_data, only : ts, reuse_jac
    use feval, only : f_rhs, f_jac, rho
    integer, intent(in) :: np
    double precision, intent(in   ) :: rho_in(np), dt
    double precision, intent(inout) :: YT(nspecies+1,np)
    logical, intent(in) :: force_new_J
    integer, intent(out), optional :: oerr

    double precision :: t0, t1, y1(nspecies+1,np)
    integer :: neq, np_bdf, i, p, ierr
    logical :: reset, reuse_J

    neq = nspecies+1
    np_bdf = ts%npt
    t0 = 0.d0
    t1 = dt

    reset = .true.

    if (mod(np,np_bdf) .ne. 0) then

       print *, "ERROR: Either np or np_bdf has wrong value", np, np_bdf
       call bl_error("ERROR: either np ot np_bdf has wrong value")

    else

       if (force_new_J) then
          reuse_J = .false.
       else
          reuse_J = reuse_jac
       end if

       do i = 1, np, np_bdf

          rho(1:np_bdf) = rho_in(i:i+np_bdf-1)

          call bdf_advance(ts, f_rhs, f_jac, neq, np_bdf, YT(:,i:i+np_bdf-1), t0,  &
               y1(:,i:i+np_bdf-1), t1, dt, reset, reuse_J, ierr)

          nstep = ts%n - 1

          reuse_J = reuse_jac

          if (ierr .ne. 0) then
             print *, 'chemsolv: BDF failed:', errors(ierr)
             print *, 'BDF rtol:', minval(ts%rtol), maxval(ts%rtol)
             print *, 'BDF atol:', minval(ts%atol), maxval(ts%atol)
             print *, 'BDF y:'
             do p = 1, ts%npt
                print *, p, ts%y(:,p)
             end do
             print *, 'BDF y0:'
             print *, YT(:,i:i+np_bdf-1)
             if (present(oerr)) then
                oerr = 1
                return
             else
                call bl_error("ERROR in burn: BDF failed")
             end if
          end if

       end do

    end if

    YT = y1
    
    if (present(oerr)) oerr = 0

  end subroutine burn_bdf


  subroutine compute_rhodYdt(np, rho, T, Y, rdYdt)
    integer, intent(in) :: np
    double precision, intent(in) :: rho(np), T(np), Y(np,nspecies)
    double precision, intent(out) :: rdYdt(np,nspecies)

    integer :: i, n, iwrk
    double precision :: rwrk

    call vckwyr(np, rho, T, Y, iwrk, rwrk, rdYdt)

    do n=1,nspecies
       do i=1,np
          rdYdt(i,n) = rdYdt(i,n) * molecular_weight(n)
       end do
    end do

  end subroutine compute_rhodYdt


  subroutine splitburn(np, rho0, Y0, rho, YT, dt)
    integer, intent(in) :: np
    double precision, intent(in   ) :: rho0, rho(np), dt
    double precision, intent(in   ) :: Y0(nspecies+1)
    double precision, intent(inout) :: YT(nspecies+1,np)

    integer :: i, j, n, iwrk, info
    double precision :: C(nspecies)
    double precision :: Y(np,nspecies), T(np), rdYdt(np,nspecies)
    double precision :: rwrk, dt_step, fac
    integer, parameter :: consP = 0
    
    if (.not. allocated(Jac)) then
       allocate(Jac(nspecies+1,nspecies+1))
       allocate(A(nspecies,nspecies))
       allocate(ipvt(nspecies))
    end if

    ! form A

    call ckytcr(rho0, Y0(nspecies+1), Y0, iwrk, rwrk, C)
    call DWDOT(Jac, C, Y0(nspecies+1), consP)

    dt_step = dt / dble(nstep)

    do j=1,nspecies
       do i=1,nspecies
          Jac(i,j) = Jac(i,j) * molecular_weight(i) * inv_mwt(j)
          A(i,j) = -dt_step*Jac(i,j)
       end do
    end do

    do i=1,nspecies
       A(i,i) = 1.d0 + A(i,i)
    end do
    
    call dgefa(A, nspecies, nspecies, ipvt, info)

    ! form b and store in YT
    
    do i=1,np
       do n=1,nspecies
          Y(i,n) = YT(n,i)
       end do
       T(i) = YT(nspecies+1,i)
    end do

    call vckwyr(np, rho, T, Y, iwrk, rwrk, rdYdt)

    do i=1,np
       fac = dt_step / rho(i)
       do n=1,nspecies
          YT(n,i) = rdYdt(i,n) * molecular_weight(n) * fac
       end do
       
       do j=1,nstep
          call dgesl(A, nspecies, nspecies, ipvt, YT(1:nspecies,i), 0)
       end do
    end do

  end subroutine splitburn


  subroutine beburn(rho0, Y0, rho, YT, dt, g, ierr)
    integer, intent(in) :: g
    double precision, intent(in   ) :: rho0, rho, dt
    double precision, intent(in   ) :: Y0(nspecies+1)
    double precision, intent(inout) :: YT(nspecies+1)
    integer, intent(out) :: ierr

    integer :: iwrk, iter, n, info, age
    double precision :: rwrk, rhoinv, cv, rmax, rmin
    double precision, dimension(nspecies) :: uk
    double precision, dimension(nspecies+1) :: YT_init, r, dYTdt
    integer, parameter :: J_int = 5

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
       
       if ( (iter.eq.0 .and. g.eq.1)    &
            .or. age.eq.J_int           &
            .or. rmax.ge.rmin ) then
          call LUA(rho, YT, dt)
          age = 0
       else
          age = age + 1
       end if

       rmin = min(rmin,rmax)

       call dgesl(A, nspecies+1, nspecies+1, ipvt, r, 0)

       YT = YT - r

    end do

    if (iter .gt. 100) then
       ierr = iter
    else
       ierr = 0
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

end module burner_module
