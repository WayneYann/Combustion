      subroutine lmc()

      implicit none

      include 'spec.h'

      integer nsteps
      integer nsteps_taken
      integer at_nstep
      integer plot_int
      integer chk_int
      integer do_initial_projection
      integer num_divu_iters
      integer num_init_iters

!     cell-centered, 2 ghost cells
      real*8, allocatable ::   vel_new(:)
      real*8, allocatable ::   vel_old(:)
      real*8, allocatable ::  scal_new(:,:)
      real*8, allocatable ::  scal_old(:,:)
      real*8, allocatable :: scal_hold(:,:)

!     cell-centered, 1 ghost cell
      real*8, allocatable ::      I_R(:,:)
      real*8, allocatable :: beta_old(:,:)
      real*8, allocatable :: beta_new(:,:)
      real*8, allocatable :: mu_dummy(:)

!     cell-centered, no ghost cells
      real*8, allocatable ::    divu_old(:)
      real*8, allocatable ::    divu_new(:)
      real*8, allocatable ::        dsdt(:)
      real*8, allocatable ::   const_src(:,:)
      real*8, allocatable :: lin_src_old(:,:)
      real*8, allocatable :: lin_src_new(:,:)

!     nodal, 1 ghost cell
      real*8, allocatable :: press_new(:)
      real*8, allocatable :: press_old(:)

      real*8 problo,probhi
      real*8 dx,time
      real*8 dt_init,dt_new
      real*8 init_shrink
      real*8 change_max
      real*8 local_change_max
      real*8 stop_time
      real*8 cfl
      real*8 cfl_used
      real*8 umax
      real*8 dt,fixed_dt
      real*8 Patm

      integer divu_iter,init_iter

      character chkfile*(16)

      namelist /fortin/ nx,nsteps,stop_time,cfl,
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
     $                  max_vode_subcycles,
     $                  min_vode_timestep, divu_ceiling_flag,
     $                  divu_dt_factor, rho_divu_ceiling, unlim

c     Set defaults, change with namelist
      nx = 256
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
      Patm = 1.d0
      coef_avg_harm = 0
      misdc_iterMAX = 3
      predict_temp_for_coeffs = 1
      do_initial_projection = 1
      num_divu_iters = 3
      num_init_iters = 2
      fixed_dt = -1.d0
      nochem_hack = .false.
      use_strang = .false.
      V_in = 1.d20
      lim_rxns = 1
      LeEQ1 = 0
      tranfile = 'tran.asc.grimech30'
      TMIN_TRANS = 0.d0
      Pr = 0.7d0
      Sc = 0.7d0
      max_vode_subcycles = 15000
      min_vode_timestep = 1.e-19
      divu_ceiling_flag = 1
      divu_dt_factor    = 0.4d0
      rho_divu_ceiling  = 0.01
      unlim = 0

      open(9,file='probin',form='formatted',status='old')
      read(9,fortin)
      close(unit=9)
      write(*,fortin)

c     Initialize chem/tran database
      call initchem()

      Pcgs = Patm * P1ATM

      call probinit(problo,probhi)

!     cell-centered, 2 ghost cells
      allocate(   vel_new(-2:nx+1))
      allocate(   vel_old(-2:nx+1))
      allocate(  scal_new(-2:nx+1,nscal))
      allocate(  scal_old(-2:nx+1,nscal))
      allocate( scal_hold(-2:nx+1,nscal))

!     cell-centered, 1 ghost cell
      allocate(      I_R(-1:nx,0:Nspec))
      allocate( beta_old(-1:nx,nscal))
      allocate( beta_new(-1:nx,nscal))
      allocate( mu_dummy(-1:nx))

!     cell-centered, no ghost cells
      allocate(    divu_old(0:nx-1))
      allocate(    divu_new(0:nx-1))
      allocate(        dsdt(0:nx-1))
      allocate(   const_src(0:nx-1,nscal))
      allocate( lin_src_old(0:nx-1,nscal))
      allocate( lin_src_new(0:nx-1,nscal))

!     nodal, 1 ghost cell
      allocate( press_new(-1:nx+1))
      allocate( press_old(-1:nx+1))

      divu_old = 0.d0
      divu_new = 0.d0

      press_old = 0.d0
      press_new = 0.d0

      dsdt = 0.d0
      
      dx = (probhi-problo)/DBLE(nx)

      lo(0) = 0
      hi(0) = nx-1

!     0=interior; 1=inflow; 2=outflow
      bc_lo(0) = 1
      bc_hi(0) = 2
      
      if ( chkfile .ne. 'null') then

         print *,'CHKFILE ',chkfile
         
         call read_check(chkfile,vel_old,scal_old,press_old,
     $                   I_R,divu_old,dsdt,
     $                   time,at_nstep,dt,cfl_used)

         call write_plt(vel_old,scal_old,press_old,divu_old,I_R,
     $                  dx,at_nstep,time)

         at_nstep = at_nstep + 1

         dt_init = dt

c     needed for seed to EOS after first strang_chem call
         scal_new(:,Temp) = scal_old(:,Temp)
                  
      else
         
         time = 0.d0
         at_nstep = 1

C take vals from PMF and fills vel, Y, and Temp
C computes rho and h, fills in rhoH and rhoY
C sets I_R to zero
         call initdata(vel_old,scal_old,I_R,dx)

c     needed for seed to EOS after first strang_chem call
         scal_new(:,Temp) = scal_old(:,Temp)

         call write_plt(vel_old,scal_old,press_old,divu_old,I_R,
     &                  dx,99999,time)

         call calc_diffusivities(scal_old,beta_old,mu_dummy,dx,time)
         
         if (do_initial_projection .eq. 1) then

            print *,'initialVelocityProject: '
            call calc_divu(scal_old,beta_old,I_R,divu_old,dx)

c     passing in dt=-1 ensures we simply project div(u)=S and
c     return zero pressure
            call project(vel_old,scal_old(0:,Density),divu_old,
     $                   press_old,press_new,dx,-1.d0,time)

         end if

         call write_plt(vel_old,scal_old,press_old,divu_old,I_R,
     &                  dx,99998,time)

         if (fixed_dt > 0) then
            dt = fixed_dt
            dt_init = dt
         else
            call est_dt(vel_old,scal_old,divu_old,dsdt,
     $                  cfl,umax,dx,dt)
            
            dt = dt * init_shrink
            dt_init = dt
         endif

         const_src = 0.d0
         lin_src_old = 0.d0
         lin_src_new = 0.d0

         print *,' '
         print *,' '
         print *,'...doing num_divu_iters = ',num_divu_iters 
         print *,' '
         print *,' '
         do divu_iter=1,num_divu_iters

            print *,' ...doing divu_iter number',divu_iter,' dt=',dt
            
            call strang_chem(scal_old,scal_new,const_src,lin_src_old,
     $                       lin_src_new,I_R,dt*0.5d0,dx,time)

c     reset temperature just in case strang_chem call is not well poased
            scal_new(:,Temp) = scal_old(:,Temp)

            call calc_divu(scal_old,beta_old,I_R,divu_old,dx)

            print *,'divu_iters velocity Project: '
            
c     passing in dt=-1 ensures we simply project div(u)=S and
c     return zero pressure
            call project(vel_old,scal_old(0:,Density),divu_old,
     $                   press_old,press_new,dx,-1.d0,time)

            dt_init = dt
             
            if (fixed_dt > 0) then
               dt = fixed_dt
            else
               call est_dt(vel_old,scal_old,divu_old,dsdt,
     $                     cfl,umax,dx,dt_new)
               dt_new = dt_new * init_shrink
               dt = min(dt_init,dt_new)
            endif
            print *,' '

         enddo

         call write_plt(vel_old,scal_old,press_old,divu_old,I_R,
     &                  dx,99997,time)

         print *,' '
         print *,'...doing num_init_iters = ',num_init_iters 
         print *,' '
         if (num_init_iters .le. 0) then
            is_first_initial_iter = 0
         else
            is_first_initial_iter = 1
         endif
         do init_iter=1,num_init_iters

            print *,' '
            print *,'INITIAL PRESSURE ITERATION ',init_iter

            if (fixed_dt > 0) then
               dt = fixed_dt
            else
               call est_dt(vel_old,scal_old,divu_old,dsdt,
     $                     cfl,umax,dx,dt)
               dt = dt * init_shrink
               dt = min(dt,dt_init)
            endif
            write(6,1001) time,dt

c     strang split overwrites scal_old so we preserve it
            if (use_strang) then
               scal_hold = scal_old
            end if

            call advance(vel_old,vel_new,scal_old,scal_new,
     $                   I_R,press_old,press_new,
     $                   divu_old,divu_new,dsdt,beta_old,beta_new,
     $                   dx,dt,time)

c     restore scal_old
            if (use_strang) then
               scal_old = scal_hold
            end if

c     update pressure and I_R
            press_old = press_new

            is_first_initial_iter = 0          

         enddo

         call write_plt(vel_old,scal_old,press_old,divu_old,I_R,
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
               call est_dt(vel_new,scal_old,divu_old,dsdt,
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
