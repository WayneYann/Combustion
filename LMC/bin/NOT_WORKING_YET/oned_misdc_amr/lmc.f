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
      real*8  scal_hold(-1:nx ,maxscal)
      real*8 press_new(0 :nx  )
      real*8 press_old(0 :nx  )
      real*8 I_R_new(0:nx-1,0:maxspec)
      real*8 I_R_old(0:nx-1,0:maxspec)
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

      integer i, n, nd, ns

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
     $                   I_R_new,divu_new,dsdt,
     $                   time,at_nstep,dt,cfl_used)

         call write_plt(vel_new,scal_new,press_new,divu_new,I_R_new,
     $                  dx,at_nstep,time)

         at_nstep = at_nstep + 1

         dt_init = dt
            
         do i = -1,nx
            vel_old(i) =  vel_new(i)
            do n = 1,nscal
               scal_old(i,n) = scal_new(i,n)
            enddo
         enddo
         do i = 0,nx-1
            do n = 0,Nspec
               I_R_old(i,n) =  I_R_new(i,n)
            enddo
            divu_old(i) = divu_new(i)
         enddo
         do i = 0,nx
            press_old(i) =  press_new(i)
         enddo
                  
      else
         
         time = 0.d0
         at_nstep = 1

C take vals from PMF and fills vel, spec (rhoY), Temp
C computes rho, rhoH
C sets I_R to zero
C Does NOT fill ghost cells
         call initdata(vel_new,scal_new,I_R_new,dx)

         call write_plt(vel_new,scal_new,press_new,divu_new,I_R_new,
     &                  dx,99999,time)

C FIXME?
C I don't think scal(RhoRT) ever actually gets used for anything,
C  But scal_aofs still wants to compute an advection term for it,
C  so initialize here to a riduculous number for now
         do i = 0,nx-1
            press_old(i) =  0.d0
            dsdt(i) =  0.d0
            dsdt(i) =  0.d0
            scal_new(i,RhoRT) = -1.d20
         enddo
         scal_new(-1,RhoRT) = -1.d20
         scal_new(nx,RhoRT) = -1.d20
         press_old(nx) =  0.d0
  
C Fills in ghost cells for rho, Y, Temp, rhoH, but not RhoRT 
         call set_bc_s(scal_new,dx,time)

         call minmax_vel(nx,vel_new)
         
         call calc_diffusivities(scal_new,beta_new,mu_new,
     $                           dx,time,.true.)

         do i=0,nx-1
            vel_old(i) = vel_new(i)
            divu_old(i) = divu_new(i)
            do n = 1,nscal
               scal_old(i,n) = scal_new(i,n)
               beta_old(i,n) = beta_new(i,n) 
            enddo
c     Define density for initial projection.
            rhohalf(i) = scal_old(i,Density)
         enddo
         
         if (do_initial_projection .eq. 1) then

            print *,'initialVelocityProject: '
            dt_dummy = -1.d0

            call calc_divu(scal_new,beta_new,I_R_new,divu_new,dx,time)

C     fills vel_new ghost cells, but not for vel_old 
            call project(vel_new,rhohalf,divu_new,
     $                   press_old,press_new,dx,dt_dummy,time)

         end if

         call write_plt(vel_new,scal_new,press_new,divu_new,I_R_new,
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
C  CEG:: needed for strang chemistry
         do i = 0,nx-1
            do n = 1,nscal
               scal_hold(i,n) = scal_old(i,n)
            enddo
         enddo
         
         do i = -1,nx
            vel_old(i) =  vel_new(i)
         enddo
         
         do i = 0,nx-1
            do n=1,nscal
               const_src(i,n) = 0.d0
               lin_src_old(i,n) = 0.d0
               lin_src_new(i,n) = 0.d0
            enddo
         enddo

         print *,' '
         print *,' '
         print *,'...doing num_divu_iters = ',num_divu_iters 
         print *,' '
         print *,' '
         do nd = 1,num_divu_iters

            print *,' ...doing divu_iter number',nd,' dt=',dt
            
C            if (use_strang) then
               call strang_chem(scal_old,scal_new,
     $              const_src,lin_src_old,lin_src_new,
     $              I_R_new,dt*0.5d0)
C            else 
C Using strang vs SDC seems to have little effect in the long run
C maybe sdc needs a better estimate of IR here
C increasing sdc iters did not help
C               sdc_iter = misdc_iterMAX
C               misdc_iterMAX = 10
C               call advance(vel_old,vel_new,scal_old,scal_new,
C     $                   I_R_new,press_old,press_new,
C     $                   divu_old,divu_new,dsdt,beta_old,beta_new,
C     $                   dx,0.5d0*dt,time)
C               do i = 0,nx-1
C                  vel_new(i) =  vel_old(i)
C               enddo
C               misdc_iterMAX = sdc_iter
C            endif
           
            call calc_divu(scal_old,beta_old,I_R_new,
     &                     divu_new,dx,time)

            print *,'divu_iters velocity Project: '
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

            do i = 0,nx-1
               do n = 0,Nspec
                  I_R_old(i,n) = I_R_new(i,n)
               enddo
               vel_old(i) =  vel_new(i)
               divu_old(i) = divu_new(i)
            enddo

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

            call advance(vel_old,vel_new,scal_old,scal_new,
     $                   I_R_new,press_old,press_new,
     $                   divu_old,divu_new,dsdt,beta_old,beta_new,
     $                   dx,dt,time)

            call minmax_vel(nx,vel_new)

            do i = 0,nx
               press_old(i)  = press_new(i)
            enddo

c     Reset state, I_R
            do i = 0,nx-1
               vel_new(i)  =   vel_old(i)
               do ns = 1,nscal
                  scal_new(i,ns) =  scal_hold(i,ns)
                  scal_old(i,ns) =  scal_hold(i,ns)
               enddo
               do ns = 0,Nspec
                  I_R_new(i,ns) =  I_R_old(i,ns)
               enddo
C ceg:: don't think this is needed.  advance overwrites what's in divu_new
               divu_new(i) = divu_old(i)
            enddo
            initial_iter = 0          
         enddo


CCCCCCCCCCCCCCCCCC
C CEG doesn't seem to make any real difference where i do this
C         do i = 0,nx-1 
C            do ns = 0,Nspec
C john's code set I_R = 0 here so time lagged I_R always gets used
C  this performs worse with SDC iters
C               I_R_new(i,ns) =  I_R_old(i,ns)
C            enddo
C         enddo
CCCCCCCCCCCCCCCCCCCc

         call write_plt(vel_new,scal_new,press_new,divu_new,I_R_new,
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
     $                I_R_new,press_old,press_new,
     $                divu_old,divu_new,dsdt,beta_old,beta_new,
     $                dx,dt,time)

         call minmax_vel(nx,vel_new)
c     update state, I_R, time
         do i = 0,nx-1
            vel_old(i)  =   vel_new(i)
            do ns = 1,nscal
               scal_old(i,ns) = scal_new(i,ns)   
            enddo
            do ns = 0,Nspec
               I_R_old(i,ns) = I_R_new(i,ns) 
            enddo
            divu_old(i) = divu_new(i)
         enddo
         do i = 0,nx
            press_old(i)  = press_new(i)
         enddo
         time = time + dt

         cfl_used = umax*dt/dx 

         if (MOD(nsteps_taken,plot_int).eq.0 .OR. 
     &        nsteps_taken.eq.nsteps) then 
            call write_plt(vel_new,scal_new,press_new,divu_new,I_R_new,
     $           dx,nsteps_taken,time)
         endif
         if (MOD(nsteps_taken,chk_int).eq.0 .OR.
     &        nsteps_taken.eq.nsteps) then 
            call write_check(nsteps_taken,vel_new,scal_new,press_new,
     $           I_R_new,divu_new,dsdt,dx,time,dt,cfl_used)
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
      
      write(6,1001) vel_max
 1001 format(' UMAX = ',f21.9)
      
      end
      
