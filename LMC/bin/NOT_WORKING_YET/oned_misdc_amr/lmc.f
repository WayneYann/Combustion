      subroutine lmc()

      implicit none

      include 'spec.h'

      integer nsteps
      integer nsteps_taken
      integer at_nstep
      integer plot_int, chk_int
      integer num_init_iters
      integer num_divu_iters
      integer do_initial_projection

      real*8   vel_new(-1:nx  )
      real*8   vel_old(-1:nx  )
      real*8  scal_new(-1:nx  ,maxscal)
      real*8  scal_old(-1:nx  ,maxscal)
      real*8 press_new(0 :nx  )
      real*8 press_old(0 :nx  )
      real*8 I_R     (0:nx-1,0:maxspec)
      real*8 I_R_hold(0:nx-1,0:maxspec)
      real*8   rhohalf( 0:nx-1)
      real*8  divu_old(0 :nx-1)
      real*8  divu_new(0 :nx-1)
      real*8  beta_old(-1 :nx,maxscal)
      real*8  beta_new(-1 :nx,maxscal)
      real*8    mu_new(-1 :nx)
      real*8      dsdt(0 :nx-1)

      real*8 problo,probhi
      real*8 dx
      real*8 time
      real*8 dt_init,dt_new
      real*8 init_shrink
      real*8 change_max
      real*8 local_change_max
      real*8 stop_time
      real*8 cfl
      real*8 cfl_used
      real*8 umax
      real*8 dt
      real*8 fixed_dt
      real*8 dt_dummy

      integer n, nd

c     New arrays for MISDC.
      real*8    const_src(0 :nx-1,maxscal)
      real*8  lin_src_old(0 :nx-1,maxscal)
      real*8  lin_src_new(0 :nx-1,maxscal)
      
      character chkfile*(16)
      real*8 Patm

      namelist /fortin/ nsteps,stop_time,cfl,
     $                  problo,probhi,chkfile,
     $                  plot_int, chk_int, change_max,
     $                  init_shrink, flame_offset,
     $                  dpdt_factor, Patm, coef_avg_harm,
     $                  misdc_iterMAX, predict_temp_for_coeffs,
     $                  do_initial_projection, num_divu_iters, 
     $                  num_init_iters,fixed_dt,
     $                  nochem_hack, use_strang, 
     $                  V_in, lim_rxns,
     $                  LeEQ1, tranfile, TMIN_TRANS, Pr, Sc,
     $                  thickFacTR, thickFacCH, max_vode_subcycles,
     $                  min_vode_timestep

c     Set defaults, change with namelist
      nsteps = 10
      stop_time = 1.e4
      cfl = 0.5
      problo = 0.0
      probhi = 3.5
      chkfile = 'null'
      plot_int = 1
      chk_int = 1
      change_max = 1.05d0
      init_shrink = 0.1d0
      flame_offset = 0.d0
      dpdt_factor = 0.d0
      Patm = 0
      coef_avg_harm = 0
      misdc_iterMAX = 3
      divu_ceiling_flag = 1
      divu_dt_factor    = 0.4d0
      rho_divu_ceiling  = 0.01
      predict_temp_for_coeffs = 1
      do_initial_projection = 1
      num_divu_iters = 3
      num_init_iters = 2
      fixed_dt = -1.d0
      nochem_hack = .false.
      use_strang = .false.
      V_in = 1.d20
      unlim = 0
      lim_rxns = 1

      tranfile = 'tran.asc.grimech30'
      TMIN_TRANS = 0.d0
      Pr = 0.7d0
      Sc = 0.7d0
      LeEQ1 = 0
      thickFacTR = 1.d0
      thickFacCH = 1.d0
      max_vode_subcycles = 15000
      min_vode_timestep = 1.e-19

      divu_old = 0.d0
      press_old = 0.d0

      divu_new = 0.d0
      press_new = 0.d0

      open(9,file='probin',form='formatted',status='old')
      read(9,fortin)
      close(unit=9)
      write(*,fortin)

c     Initialize chem/tran database
      call initchem()

      Pcgs = Patm * P1ATM
      
      dx = (probhi-problo)/DBLE(nx)
      
      call probinit(problo,probhi)
      
      if ( chkfile .ne. 'null') then

         print *,'CHKFILE ',chkfile
         
         call read_check(chkfile,vel_new,scal_new,press_new,
     $                   I_R,divu_new,dsdt,
     $                   time,at_nstep,dt,cfl_used)

         call write_plt(vel_new,scal_new,press_new,divu_new,I_R,
     $                  dx,at_nstep,time)

         at_nstep = at_nstep + 1

         dt_init = dt

         vel_old = vel_new
         scal_old = scal_new
         divu_old = divu_new
         press_old = press_new
                  
      else
         
         time = 0.d0
         at_nstep = 1

C take vals from PMF and fills vel, Y, and Temp
C computes rho and h, fills in rhoH and rhoY
C sets I_R to zero
C fills ghost cells
         call initdata(vel_new,scal_new,I_R,dx)

         call write_plt(vel_new,scal_new,press_new,divu_new,I_R,
     &                  dx,99999,time)

C Note: RhoRT is only a diagnostic
         scal_new(:,RhoRT) = -1.d20
         press_old = 0.d0
         dsdt = 0.d0

         call minmax_vel(nx,vel_new)
         
         call calc_diffusivities(scal_new,beta_new,mu_new,
     $                           dx,time,.true.)

         vel_old = vel_new
         divu_old = divu_new
         scal_old = scal_new
         beta_old = beta_new
c     Define density for initial projection.
         rhohalf(0:nx-1) = scal_old(0:nx-1,Density)
         
         if (do_initial_projection .eq. 1) then

            print *,'initialVelocityProject: '
c     setting dt=-1 ensures we simply project div(u)=S and
c     return zero pressure
            dt_dummy = -1.d0

            call calc_divu(scal_new,beta_new,I_R,divu_new,dx,time)

C     fills vel_new ghost cells
            call project(vel_new,rhohalf,divu_new,
     $                   press_old,press_new,dx,dt_dummy,time)

         end if

         call write_plt(vel_new,scal_new,press_new,divu_new,I_R,
     &                  dx,99998,time)


         if (fixed_dt > 0) then
            dt = fixed_dt
            dt_init = dt
         else
            call est_dt(nx,vel_new,scal_new,divu_new,dsdt,
     $                  cfl,umax,dx,dt)
            
            dt = dt * init_shrink
            dt_init = dt
         endif
         
         vel_old = vel_new

         const_src = 0.d0
         lin_src_old = 0.d0
         lin_src_new = 0.d0

         print *,' '
         print *,' '
         print *,'...doing num_divu_iters = ',num_divu_iters 
         print *,' '
         print *,' '
         do nd = 1,num_divu_iters

            print *,' ...doing divu_iter number',nd,' dt=',dt
            
            call strang_chem(scal_old,scal_new,const_src,lin_src_old,
     $                       lin_src_new,I_R,dt*0.5d0)

            call calc_divu(scal_old,beta_old,I_R,divu_new,dx,time)

            print *,'divu_iters velocity Project: '
c     setting dt=-1 ensures we simply project div(u)=S and
c     return zero pressure
            dt_dummy = -1.d0
            
C vel_old does not get used in proj(), 
C assumes good data in vel_new
            call project(vel_new,rhohalf,divu_new,
     $                   press_old,press_new,dx,dt_dummy,time)

            dt_init = dt
             
            if (fixed_dt > 0) then
               dt = fixed_dt
            else
C CEG:: not sure that this should be scal_new and not scal_old
C   probably doesn't matter that much
               call est_dt(nx,vel_new,scal_new,divu_new,dsdt,
     $                     cfl,umax,dx,dt_new)
               dt_new = dt_new * init_shrink
               dt = min(dt_init,dt_new)
            endif
            print *,' '

            vel_old = vel_new
            divu_old = divu_new

         enddo

         print *,' '
         print *,'...doing num_init_iters = ',num_init_iters 
         print *,' '
         if (num_init_iters .le. 0) then
            initial_iter = 0
         else
            initial_iter = 1
         endif
         do n = 1,num_init_iters

            print *,' '
            print *,'INITIAL PRESSURE ITERATION ',n

            if (fixed_dt > 0) then
               dt = fixed_dt
            else
               call est_dt(nx,vel_new,scal_old,divu_old,dsdt,
     $                     cfl,umax,dx,dt)
               dt = dt * init_shrink
               dt = min(dt,dt_init)
            endif
            write(6,1001) time,dt

            I_R_hold = I_R

            call advance(vel_old,vel_new,scal_old,scal_new,
     $                   I_R,press_old,press_new,
     $                   divu_old,divu_new,dsdt,beta_old,beta_new,
     $                   dx,dt,time)

            call minmax_vel(nx,vel_new)

c     update pressure
            press_old = press_new

c     reset state
            I_R = I_R_hold
            divu_new = divu_old
            vel_new = vel_old
            scal_new = scal_old

            initial_iter = 0          

         enddo

         call write_plt(vel_new,scal_new,press_new,divu_new,I_R,
     &                  dx,0,time)

         cfl_used = cfl * init_shrink

         print *,' '      
         print *,' '      
         print *,'COMPLETED INITIAL ITERATIONS'
         print *,' '      
         print *,'START ADVANCING THE SOLUTION '
         print *,' '            

 1001    format('Advancing: starting time = ',
     $        e15.9,' with dt = ',e15.9)

      endif

C-- Now advance 
      do nsteps_taken = at_nstep, nsteps

         if (time.ge.stop_time) exit

         if (time > 0.d0) then 
            if (fixed_dt > 0) then
               dt = fixed_dt
            else
               call est_dt(nx,vel_new,scal_old,divu_old,dsdt,
     $                     cfl,umax,dx,dt_new)
               if (nsteps_taken .eq. 1) then
                  dt = dt*init_shrink
                  dt = min(dt,dt_init)
               else
                  local_change_max = min(change_max,cfl/cfl_used)
                  dt = dt * local_change_max 
                  dt = min(dt,dt_new)
               endif            
            endif
         endif

          dt = min(dt,stop_time-time)

         write(6,*)
         write(6,1001 )time,dt
         write(6,*)'STEP = ',nsteps_taken
         
         call advance(vel_old,vel_new,scal_old,scal_new,
     $                I_R,press_old,press_new,
     $                divu_old,divu_new,dsdt,beta_old,beta_new,
     $                dx,dt,time)

         call minmax_vel(nx,vel_new)
c     update state, time
         vel_old = vel_new
         scal_old = scal_new
         divu_old = divu_new
         press_old = press_new

         time = time + dt

         cfl_used = umax*dt/dx 

         if (MOD(nsteps_taken,plot_int).eq.0 .OR. 
     &        nsteps_taken.eq.nsteps) then 
            call write_plt(vel_new,scal_new,press_new,divu_new,I_R,
     $           dx,nsteps_taken,time)
         endif
         if (MOD(nsteps_taken,chk_int).eq.0 .OR.
     &        nsteps_taken.eq.nsteps) then 
            call write_check(nsteps_taken,vel_new,scal_new,press_new,
     $           I_R,divu_new,dsdt,dx,time,dt,cfl_used)
         endif
      enddo



      print *,' '      
      print *,'COMPLETED SUCCESSFULLY'
      print *,' '      

      end

      subroutine minmax_vel(nx,vel)
      integer nx
      real*8  vel(-1:nx  )

      real*8  vel_min,vel_max

      vel_min = abs(vel(0))
      vel_max = abs(vel(0))
      do i = 1,nx-1
         vel_min = min(vel_min,abs(vel(i)))
         vel_max = max(vel_max,abs(vel(i)))
      enddo
      
      write(6,1001) vel_min
      write(6,1002) vel_max
 1001 format(' UMAX = ',f21.9)
 1002 format(' UMAX = ',f21.9)
      
      end
      
