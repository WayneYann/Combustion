subroutine chemsolve(lo, hi, &
  uin, il1, il2, ih1, ih2, &
  uou, ol1, ol2, oh1, oh2, &
  stop_time, dt, verbose, use_vode)
  use chemistry_module, only : nspecies
  implicit none
  integer, intent(in) :: lo(2), hi(2)
  integer, intent(in) :: il1, il2, ih1, ih2
  integer, intent(in) :: ol1, ol2, oh1, oh2
  double precision, intent(in ) :: uin(il1:ih1,il2:ih2,nspecies+2)
  double precision, intent(out) :: uou(ol1:oh1,ol2:oh2,nspecies+2)
  double precision, intent(in) :: stop_time, dt
  integer, intent(in) :: verbose, use_vode

  integer :: i, j, n
  double precision :: rhot, Yt(nspecies+1)

  if (use_vode .ne. 0) call setfirst(.true.)

  do j=lo(2),hi(2)
     do i=lo(1),hi(1)

        do n=1,nspecies+1
           Yt(n) = uin(i,j,n)
        end do
        rhot = uin(i,j,nspecies+2)
        
        if (use_vode .ne. 0) then
           call burn_vode(rhot, Yt, stop_time, dt, verbose)
        else
           call burn_bdf(rhot, Yt, stop_time, dt, verbose)           
        end if

        do n=1,nspecies+1
           uou(i,j,n) = Yt(n)
        end do
        uou(i,j,nspecies+2) = rhot

     end do
  end do

end subroutine chemsolve


subroutine burn_vode(rho, YT, stop_time, dt, verbose)
  use vode_module, only : itol, rtol, atol, vode_MF=>MF, always_new_j, &
       voderwork, vodeiwork, lvoderwork, lvodeiwork, voderpar, vodeipar
  use chemistry_module, only : nspecies, spec_names
  integer, intent(in) :: verbose
  double precision, intent(in   ) :: rho, stop_time, dt
  double precision, intent(inout) :: YT(nspecies+1)

  external f_jac, f_rhs, dvode
  
  ! vode stuff
  integer, parameter :: itask=1, iopt=1
  integer :: MF, istate, ifail
  
  double precision :: time
  
  voderpar(1) = rho
  
  istate = 1
  time = 0.d0
  
  if (always_new_j) call setfirst(.true.)
  
  MF = vode_MF  ! vode might change its sign!
  call dvode(f_rhs, nspecies+1, YT, time, stop_time, itol, rtol, atol, itask, &
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
  
end subroutine burn_vode


subroutine burn_bdf(rho_in, YT, stop_time, dt, verbose)
  use chemistry_module, only : nspecies
  use bdf
  use bdf_data, only : ts, reuse_jac
  use feval, only : f_rhs, f_jac, rho
  double precision, intent(in   ) :: rho_in, stop_time, dt
  double precision, intent(inout) :: YT(nspecies+1)
  
  double precision :: t0, t1, y1(nspecies+1)
  integer :: neq, ierr
  logical :: reset
  
  rho = rho_in
  
  neq = nspecies+1
  t0 = 0.d0
  t1 = stop_time

  reset = .true.

  call bdf_advance(ts, f_rhs, f_jac, neq, YT, t0, y1, t1, dt, reset, reuse_jac, ierr)

  if (ierr .ne. 0) then
     print *, 'chemsolv: BDF failed'
     call bl_error("ERROR in burn: BDF failed")       
  end if
  
end subroutine burn_bdf

