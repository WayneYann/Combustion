module burner_module

  use chemistry_module, only : nspecies, spec_names, molecular_weight, inv_mwt
  use meth_params_module, only : use_vode

  implicit none

  double precision, allocatable, save :: Jac(:,:), A(:,:)
  integer, allocatable, save :: ipvt(:)
  !$omp threadprivate(Jac,A,ipvt) 

  private

  public :: burn, compute_rhodYdt, init_burn_linear, burn_linear

contains

  subroutine burn(np, rho, YT, dt, force_new_J)
    use meth_params_module, only : use_vode
    use feval, only : f_rhs, rho_feval => rho
    integer, intent(in) :: np
    double precision, intent(in   ) :: rho(np), dt
    double precision, intent(inout) :: YT(nspecies+1,np)
    logical, intent(in) :: force_new_J

    double precision :: YTdot(nspecies+1,np)

    if (dt == 0.d0) then
       rho_feval(1:np) = rho
       call f_rhs(nspecies+1, np, YT, 0.d0, YTdot)
       YT = YTdot
       return
    end if

    if (use_vode) then
       call burn_vode(np, rho, YT, dt, force_new_J)
    else
       call burn_bdf(np, rho, YT, dt, force_new_J)
    end if

  end subroutine burn


  subroutine burn_vode(np, rho, YT, dt, force_new_J)
    use vode_module, only : verbose, itol, rtol, atol, vode_MF=>MF, always_new_j, &
         voderwork, vodeiwork, lvoderwork, lvodeiwork, voderpar, vodeipar
    integer, intent(in) :: np
    double precision, intent(in   ) :: rho(np), dt
    double precision, intent(inout) :: YT(nspecies+1,np)
    logical, intent(in) :: force_new_J

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
          call bl_error("ERROR in burn: VODE failed")
       end if

    end do

  end subroutine burn_vode


  subroutine burn_bdf(np, rho_in, YT, dt, force_new_J)
    use bdf
    use bdf_data, only : ts, reuse_jac
    use feval, only : f_rhs, f_jac, rho
    integer, intent(in) :: np
    double precision, intent(in   ) :: rho_in(np), dt
    double precision, intent(inout) :: YT(nspecies+1,np)
    logical, intent(in) :: force_new_J

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
             call bl_error("ERROR in burn: BDF failed")
          end if

       end do

    end if

    YT = y1

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


  subroutine init_burn_linear(rho0, Y0, dt)
    double precision, intent(in) :: rho0, Y0(nspecies+1), dt

    integer :: i, j, iwrk, info
    double precision :: C(nspecies)
    double precision :: rwrk
    integer, parameter :: consP = 0
    
    if (.not. allocated(Jac)) then
       allocate(Jac(nspecies+1,nspecies+1))
       allocate(A(nspecies,nspecies))
       allocate(ipvt(nspecies))
    end if

    call ckytcr(rho0, Y0(nspecies+1), Y0, iwrk, rwrk, C)
    call DWDOT(Jac, C, Y0(nspecies+1), consP)

    do j=1,nspecies
       do i=1,nspecies
          Jac(i,j) = Jac(i,j) * molecular_weight(i) * inv_mwt(j)
          A(i,j) = -dt*Jac(i,j)
       end do
    end do

    do i=1,nspecies
       A(i,i) = 1.d0 + A(i,i)
    end do
    
    call dgefa(A, nspecies, nspecies, ipvt, info)

  end subroutine init_burn_linear


  subroutine burn_linear(dY)
    double precision, intent(inout) :: dY(nspecies)
    call dgesl(A, nspecies, nspecies, ipvt, dY, 0)
  end subroutine burn_linear

end module burner_module
