	program supernova

        implicit none

        include 'nums.fi'
        include 'sndata.fi'
        include 'probdata.fi'

        integer nx
        integer nsteps
        integer nsteps_taken
        integer at_nstep
        integer plot_int, chk_int
        integer num_init_iters
        integer num_divu_iters
        parameter (num_divu_iters = 3)
        parameter (num_init_iters = 2)
        parameter (nx = 512)

        integer MAX_NSCAL
        parameter (MAX_NSCAL = 7)

        real*8   vel_new(-1:nx  )
        real*8   vel_old(-1:nx  )
        real*8  scal_new(-1:nx  ,MAX_NSCAL)
        real*8  scal_old(-1:nx  ,MAX_NSCAL)
        real*8  scal_hold(-1:nx  ,MAX_NSCAL)
        real*8 press_new(0 :nx  )
        real*8 press_old(0 :nx  )
        real*8  Ydot_new(0 :nx-1,MAX_NSCAL)
        real*8  Ydot_old(0 :nx-1,MAX_NSCAL)
        real*8   rhohalf(0 :nx-1)
        real*8  divu_old(0 :nx-1)
        real*8  divu_new(0 :nx-1)
        real*8  beta_old(0 :nx-1,MAX_NSCAL)
        real*8  beta_new(0 :nx-1,MAX_NSCAL)
        real*8      dsdt(0 :nx-1)
        real*8     intra(0 :nx-1,MAX_NSCAL)


c       Local variables
        real*8 problo,probhi
        real*8 dx
        real*8 time
        real*8 dt_old
        real*8 dt_init,dt_init2
        real*8 init_shrink
        real*8 change_max
        real*8 local_change_max
        real*8 stop_time
        real*8 cfl
        real*8 cfl_used
        real*8 umax
        real*8 dt
        real*8 dt_dummy

        integer i,n,ns,nd,np,ispec
        integer length
        integer lo,hi,lo_lim,hi_lim
        integer do_init

c       New arrays for MISDC.
        real*8    const_src(0 :nx-1,MAX_NSCAL)
        real*8  lin_src_old(0 :nx-1,MAX_NSCAL)
        real*8  lin_src_new(0 :nx-1,MAX_NSCAL)

        character chkfile*(16)
        
        namelist /fortin/ nsteps,stop_time,cfl,
     $                    problo,probhi,chkfile,
     $                    plot_int, chk_int
        open(9,file='probin_sn',form='formatted',status='old')
        read(9,fortin)
        close(unit=9)
        print *,'  '
        print *,' SETTING nsteps    = ',nsteps
        print *,' SETTING stop_time = ',stop_time
        print *,' SETTING cfl       = ',cfl
        print *,' SETTING problo    = ',problo
        print *,' SETTING probhi    = ',probhi
        print *,' SETTING chkfile   = ',chkfile
        print *,'  '

        dx = (probhi-problo)/dble(nx)
        call probinit(problo,probhi)
        nscal = nspec+3
        change_max  = 1.05d0
        init_shrink = 0.1d0

        length = len(chkfile)
        print *,'LENGTH ',length
        print *,'CHKFILE ',chkfile
        if (chkfile(1:1) .eq. 'c') then

          call readcheck(chkfile,nx,vel_new,scal_new,press_new,
     $                   Ydot_new,divu_new,dsdt,intra,
     $                   time,at_nstep,dt_old,cfl_used)
          do_init = 0

          do i = -1,nx
            vel_old(i) =  vel_new(i)
            do n = 1,nscal
              scal_old(i,n) = scal_new(i,n)
            enddo
          enddo
          do i = 0,nx-1
            do n = 1,nspec
              Ydot_old(i,n) =  Ydot_new(i,n)
            enddo
            divu_old(i) = divu_new(i)
          enddo
          do i = 0,nx
            press_old(i) =  press_new(i)
          enddo

          call minmax_vel(nx,vel_new)

        else

          do_init = 1
          time = 0.d0
          at_nstep = 0

          call initdata(nx,vel_new,scal_new,Ydot_new,dx)

c         Define RhoH as function of Temp.
          lo     = 0
          hi     = nx-1
          lo_lim = -1
          hi_lim = nx
          call temp_to_rhoh(lo,hi,lo_lim,hi_lim,scal_new)
          scal_new(nx,RhoH) = scal_new(nx-1,RhoH)

          do i = -1,nx
            vel_old(i) =  vel_new(i)
            do n = 1,nscal
              scal_old(i,n) = scal_new(i,n)
            enddo
          enddo

          call minmax_vel(nx,vel_new)

c         Define density for initial projection.
          do i = 0,nx-1
            rhohalf(i) = scal_old(i,Density)
          enddo

          call calc_diffusivities(nx,scal_new,beta_new)
  
          call calc_divu(nx,scal_new,beta_new,Ydot_new,divu_new,dx)

          call est_dt(nx,vel_new,scal_new,divu_new,dsdt,
     $                cfl,umax,dx,dt)

          print *,'initialVelocityProject: levels = 0  0'
          dt_dummy = -1.d0
          call project(nx,vel_old,vel_new,rhohalf,divu_new,
     $                 press_old,press_new,dx,dt_dummy)

          call est_dt(nx,vel_new,scal_new,divu_new,dsdt,
     $                cfl,umax,dx,dt)

          dt = dt * init_shrink

          dt_init = dt

          do i = -1,nx
             vel_old(i) =  vel_new(i)
          enddo

          do i = 0,nx
            do n = 1,nscal
                scal_hold(i,n) = scal_old(i,n)
            enddo
          enddo

          do i = 0,nx-1
            do n = 1,nscal
                const_src(i,n) = 0.d0
              lin_src_old(i,n) = 0.d0
              lin_src_new(i,n) = 0.d0
            enddo
          enddo

          print *,' '
          print *,'...doing num_divu_iters = ',num_divu_iters 
          print *,' '

          do nd = 1,num_divu_iters
            print *,' '
            write(6,1000) nd
1000        format('...doing ',i2,'th divu_iter ')
            print *,' '
            call strang_chem(nx,scal_old,scal_new,
     $                       const_src,lin_src_old,lin_src_new,
     $                       intra,dt)
            do i = 0,nx-1
              do n = 1,nspec 
                ispec = FirstSpec-1+n
                Ydot_new(i,n) = 
     $           (scal_new(i,ispec)-scal_old(i,ispec))/dt
              enddo
            enddo
            do i = 0,nx
              do n = 1,nscal
                scal_new(i,n) = scal_hold(i,n)
              enddo
            enddo

            call calc_divu(nx,scal_new,beta_new,Ydot_new,divu_new,dx)
            print *,'initialVelocityProject: levels = 0  0'
            dt_dummy = -1.d0
            call project(nx,vel_old,vel_new,rhohalf,divu_new,
     $                   press_old,press_new,dx,dt_dummy)
            print *,' '
            print *,' '
            dt_init = dt
            call est_dt(nx,vel_new,scal_new,divu_new,dsdt,
     $                  cfl,umax,dx,dt_init2)
            dt_init2 = dt_init2 * init_shrink
            dt = min(dt_init,dt_init2)
            print *,' '
          enddo
  
          do i = 0,nx-1
            do n = 1,nspec
              Ydot_old(i,n) = Ydot_new(i,n)
            enddo
             vel_old(i) =  vel_new(i)
            divu_old(i) = divu_new(i)
          enddo
  
          do n = 1,num_init_iters

            call est_dt(nx,vel_new,scal_old,divu_old,dsdt,
     $                  cfl,umax,dx,dt)
            dt = dt * init_shrink
            dt = min(dt,dt_init)
            print *,' '
            print *,'INITIAL PRESSURE ITERATION ',n
            print *,' '
            write(6,1001) time,dt

c           Here we zero out intra before each advance.
            do ns = 1,nscal
             do i = 0,nx-1
                 intra(i,n) = 0.d0
             enddo
            enddo

            call advance(nx,vel_old,vel_new,scal_old,scal_new,
     $                   Ydot_old,Ydot_new,press_old,press_new,
     $                   divu_old,divu_new,dsdt,beta_old,beta_new,
     $                   intra,dx,dt,time)
            call minmax_vel(nx,vel_new)

            do i = 0,nx
                 vel_new(i)  =   vel_old(i)
            enddo

            do i = 0,nx
               do ns = 1,nscal
                 scal_new(i,ns) =  scal_hold(i,ns)
                 scal_old(i,ns) =  scal_hold(i,ns)
               enddo
            enddo
            do i = 0,nx-1
               do ns = 1,nspec
                 Ydot_new(i,ns) =  Ydot_old(i,ns)
               enddo
               divu_new(i) = divu_old(i)
            enddo

            do i = 0,nx
               press_old(i)  = press_new(i)
            enddo

          enddo
1001      format('Advancing level 0 : starting time = ',
     $            e15.9,' with dt = ',e15.9)

          print *,'COMPLETED INITIAL ITERATIONS'
          print *,' '

          cfl_used = cfl * init_shrink

c         Here we zero out intra before the "real" time stepping,
c           but not within advance, so we use a lagged intra each time.
          do n = 1,nscal
           do i = 0,nx-1
               intra(i,n) = 0.d0
           enddo
          enddo

        endif

        nsteps_taken = at_nstep

        do n = at_nstep+1,nsteps

          if (time.ge.stop_time) exit

          call est_dt(nx,vel_new,scal_old,divu_old,dsdt,cfl,umax,dx,dt)
          if (n.eq.1) then
            print *,'MULT BY INIT_SHRINK ',dt,init_shrink
            dt = dt * init_shrink
            print *,'MIN OF ',dt,' AND ',dt_init
            dt = min(dt,dt_init)
            dt_old = dt
          else
            local_change_max = min(change_max,cfl/cfl_used)
            print *,'MIN OF ',dt,' AND ',local_change_max,' * ',dt_old
            dt = min(dt,local_change_max * dt_old)
          endif
          dt = min(dt,stop_time-time)
          print *,' '
          write(6,1001) time,dt
          print *,' '

          call advance(nx,vel_old,vel_new,scal_old,scal_new,
     $                 Ydot_old,Ydot_new,press_old,press_new,
     $                 divu_old,divu_new,dsdt,beta_old,beta_new,
     $                 intra,dx,dt,time)

          call minmax_vel(nx,vel_new)

          do i = 0,nx
            vel_old(i) =  vel_new(i)
            do ns = 1,nscal
              scal_old(i,ns) =  scal_new(i,ns)
            enddo
          enddo

          do i = 0,nx-1
             do ns = 1,nspec
               Ydot_old(i,ns) =  Ydot_new(i,ns)
             enddo
             divu_old(i) = divu_new(i)
          enddo

          do i = 0,nx
             press_old(i)  = press_new(i)
          enddo

          time = time + dt
          dt_old = dt

          cfl_used = umax * dt / dx

          print *,' '
          write(6,1002) n,time,dt
          print *,' '
1002      format('STEP = ',i6,' TIME = ',e16.10,' DT = ',e16.10)

          nsteps_taken = n
          if (time.ge.stop_time) exit 

          np = nsteps_taken / plot_int
          if (nsteps_taken .eq. np*plot_int)
     $    call   writeplt(nx,vel_new,scal_new,press_new,dx,
     $                    nsteps_taken,time)

          np = nsteps_taken / chk_int
          if (nsteps_taken .eq. np*chk_int)
     $      call writecheck(nsteps_taken,nx,vel_new,scal_new,press_new,
     $                      Ydot_new,divu_new,dsdt,intra,dx,time,dt_old,
     $                      cfl_used)

        enddo

        call writecheck(nsteps_taken,nx,vel_new,scal_new,press_new,
     $                  Ydot_new,divu_new,dsdt,intra,
     $                  dx,time,dt_old,cfl_used)
        call   writeplt(nx,vel_new,scal_new,press_new,dx,
     $                  nsteps_taken,time)

        end

        subroutine minmax_scal(nx,scal)

c       Quantities passed in
        integer nx
        real*8  scal(-1:nx,nscal)

c       Local variables
        real*8  smin,smax

        include 'sndata.fi'
        include 'nums.fi'

        print *,' '
        smin = abs(scal(0,Temp))
        smax = smin
        do i = 1,nx-1
          smin = min(smin,abs(scal(i,Temp)))
          smax = max(smax,abs(scal(i,Temp)))
        enddo

        write(6,1001) smin,smax
1001    format(' Min,max temp = ',f13.2,', ',f13.2)

        smin = abs(scal(0,Density))
        smax = smin
        do i = 1,nx-1
          smin = min(smin,abs(scal(i,Density)))
          smax = max(smax,abs(scal(i,Density)))
        enddo

        write(6,1002) smin,smax
1002    format(' Min,max rho  = ',f13.2,', ',f13.2)

        smin = abs(scal(0,RhoH))
        smax = smin
        do i = 1,nx-1
          smin = min(smin,abs(scal(i,RhoH)))
          smax = max(smax,abs(scal(i,RhoH)))
        enddo

        write(6,1003) smin,smax
1003    format(' Min,max rhoh = ',e16.10,', ',e16.10)

        print *,' '

        end

        subroutine minmax_vel(nx,vel)

c       Quantities passed in
        integer nx
        real*8  vel(-1:nx  )

c       Local variables
        real*8  vel_min,vel_max

        vel_min = abs(vel(0))
        vel_max = abs(vel(0))
        do i = 1,nx-1
          vel_min = min(vel_min,abs(vel(i)))
          vel_max = max(vel_max,abs(vel(i)))
        enddo

        write(6,1001) vel_max
1001    format(' LEV = 0 UMAX = ',f21.9)

        end
