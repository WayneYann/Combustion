      subroutine lmc()
      implicit none
      include 'spec.h'
      integer nsteps
      integer nsteps_taken
      integer at_nstep
      integer plot_int, chk_int
      integer num_init_iters
      integer num_divu_iters
      
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
      real*8    mu_old(-1 :nx)
      real*8    mu_new(-1 :nx)
      real*8      dsdt(0 :nx-1)
      real*8    tforce(0 :nx-1,maxscal)


c     Local variables
      real*8 problo,probhi
      real*8 dx
      real*8 time
      real*8 dt_old
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

      real*8 divu_max

      integer do_init, is, i, n, nd, ns

C CEG debgging REMOVE ME
      integer hi, lo, ncomp,j
      real*8 ptherm(-1:nx)
      real*8 Y(maxspec)
      
c     New arrays for MISDC.
      real*8    const_src(0 :nx-1,maxscal)
      real*8  lin_src_old(0 :nx-1,maxscal)
      real*8  lin_src_new(0 :nx-1,maxscal)
      
      character chkfile*(16)
      real*8 Patm

      namelist /fortin/ nsteps,stop_time,cfl,
     $                  problo,probhi,chkfile,
     $                  plot_int, chk_int, change_max,
     $                  init_shrink, probtype,flame_offset,
     $                  dpdt_factor, Patm, coef_avg_harm,
     $                  misdc_iterMAX, predict_temp_for_coeffs,
     $                  num_divu_iters, num_init_iters,fixed_dt


c     Initialize chem/tran database
      call initchem()

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
      probtype = 1
      flame_offset = 0.d0
      dpdt_factor = 0.d0
      Patm = 0
      coef_avg_harm = 0
      misdc_iterMAX = 3
      divu_ceiling_flag = 1
      divu_dt_factor    = 0.4d0
      rho_divu_ceiling  = 0.01
      predict_temp_for_coeffs = 1
      num_divu_iters = 3
      num_init_iters = 2
      fixed_dt = -1.d0

      open(9,file='probin',form='formatted',status='old')
      read(9,fortin)
      close(unit=9)
c      write(*,fortin)

      nochem_hack = .true.

      Pcgs = Patm * P1ATM
      
      dx = (probhi-problo)/DBLE(nx)
      
      call probinit(problo,probhi)
      
      if ( chkfile .ne. 'null') then

         print *,'CHKFILE ',chkfile
         
         call read_check(chkfile,nx,vel_new,scal_new,press_new,
     $                   I_R_new,divu_new,dsdt,
     $                   time,at_nstep,dt_old,cfl_used)
         do_init = 0
            
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
         
c     call minmax_vel(nx,vel_new)
         
      else
         
         do_init = 1
         time = 0.d0
         at_nstep = 0

C take vals from PMF and fills vel, spec (rhoY), Temp
C                              computes rho, rhoH, I_R
C Does NOT fill ghost cells
C vel and press have uninitialized (ie grabage in) ghost cells coming out
         call initdata(vel_new,scal_new,I_R_new(0,0),dx)
C CEG:: not sure where/if this get initialized if not done here
C FIXME
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
  
C Fills in ghost cells for rho, Y, Temp, rhoH 
         call set_bc_grow_s(scal_new,dx,time)
         do i = -1,nx
            do n = 1,nscal
               scal_old(i,n) = scal_new(i,n)
            enddo
         enddo

         call minmax_vel(nx,vel_new)
         
         call calc_diffusivities(scal_new,beta_new,mu_new)

         call calc_divu(scal_new,beta_new,I_R_new,divu_new,dx,time)
         call write_plt(vel_new,scal_new,press_new,divu_new,
     &                  dx,dt,99999,time)
        
         print *,'initialVelocityProject: '
         dt_dummy = -1.d0

         do i=0,nx-1
            vel_old(i) = vel_new(i)
C CEG:: for the event that divu_iters = 0
            divu_old(i) = divu_new(i)
         enddo

c     Define density for initial projection.
         do i = 0,nx-1
            rhohalf(i) = scal_old(i,Density)
         enddo

C fills vel_new ghost cells, but not for vel_old 
C press_new gets vals everywhere
         call project(vel_old,vel_new,rhohalf,divu_new,
     $                press_old,press_new,dx,dt_dummy,time)

C debugging FIXME
C 1006 FORMAT((I5,1X),11(E22.15,1X))      
C         hi = 255
C         lo = 0
C         ncomp = 11
C         call compute_pthermo(scal_new,ptherm)
C         open(UNIT=11, FILE='after.dat', STATUS = 'REPLACE')
C         write(11,*)'# ', hi-lo, ncomp 
C         do j=lo,hi
C            do n = 1,Nspec
C               Y(n) = scal_new(j,FirstSpec+n-1)*1.d3
C            enddo
C            write(11,1006) j, vel_new(j)*1.d-2, 
C     &                     scal_new(j,Density)*1.d3,
C     &                     (Y(n),n=1,Nspec),
C     $                     scal_new(j,RhoH)*1.d-1,
C     $                     scal_new(j,Temp),
C     $                     ptherm(j)*1.d-1
C         enddo
C         close(11)
C         stop
CCCCCCCCCCCCC

         if (fixed_dt > 0) then
            dt = fixed_dt
         else
            call est_dt(nx,vel_new,scal_new,divu_new,dsdt,
     $                  cfl,umax,dx,dt)
            
            dt = dt * init_shrink
            dt_init = dt
         endif
C CEG:: why no ghost cells here?-- everything is set up to fill ghost cells
C  right before they are used
C  I don't know if this is neccessary.  I don't think scal_old gets changed
         do i = 0,nx-1
            do n = 1,nscal
               scal_hold(i,n) = scal_old(i,n)
            enddo
         enddo
         
C CEG:: do i really need to change the ghosts here???
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
            
            call strang_chem(scal_old,scal_new,
     $                       const_src,lin_src_old,lin_src_new,
     $                       I_R_new,dt)
           
            call calc_divu(scal_new,beta_new,I_R_new,
     &                     divu_new,dx,time)

            print *,'divu_iters velocity Project: '
            dt_dummy = -1.d0
            
            call project(vel_old,vel_new,rhohalf,divu_new,
     $                   press_old,press_new,dx,dt_dummy,time)

            dt_init = dt
             
            if (fixed_dt > 0) then
               dt = fixed_dt
            else
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

C debugging FIXME
C$$$         hi = 255
C$$$         lo = 0
C$$$         ncomp = 11
C$$$         call compute_pthermo(scal_new,ptherm)
C$$$         open(UNIT=11, FILE='before.dat', STATUS = 'REPLACE')
C$$$         write(11,*)'# ', hi-lo, ncomp 
C$$$         do j=lo,hi
C$$$            do n = 1,Nspec
C$$$               Y(n) = scal_new(i,FirstSpec+n-1)/scal_new(i,Density)
C$$$            enddo
C$$$            write(11,1006) j, vel_new(i)*1.d-2, 
C$$$     &                     scal_new(i,Density)*1.d3,
C$$$     &                     (Y(n),n=1,Nspec),
C$$$     $                     scal_new(i,RhoH)*1.d-4,
C$$$     $                     scal_new(i,Temp),
C$$$     $                     ptherm(i)*1.d-1
C$$$         enddo
C$$$         close(11)
CCCCCCCCCCCCC

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
     $                   I_R_new,I_R_new,press_old,press_new,
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
               do ns = 1,Nspec
                  I_R_new(i,ns) =  I_R_old(i,ns)
               enddo
               divu_new(i) = divu_old(i)
            enddo
                        
         enddo

         call write_plt(vel_new,scal_new,press_new,divu_new,
     &                  dx,dt,0,time)

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
C-- Done with initialization--------------------

C-- Now advance 
      do nsteps_taken = 1, nsteps

         if (time > 0.d0) then 
            if (fixed_dt > 0) then
               dt = fixed_dt
            else
               call est_dt(nx,vel_new,scal_old,divu_old,dsdt,
     $                     cfl,umax,dx,dt_new)
               dt = dt * change_max 
               dt = min(dt,dt_new)            
            endif
         endif

         write(6,*)
         write(6,1001 )time,dt
         write(6,*)'STEP = ',nsteps_taken
         
         call advance(vel_old,vel_new,scal_old,scal_new,
     $                I_R_new,I_R_new,press_old,press_new,
     $                divu_old,divu_new,dsdt,beta_old,beta_new,
     $                dx,dt,time)

c     update state, I_R, time
         do i = 0,nx-1
            vel_old(i)  =   vel_new(i)
            do ns = 1,nscal
               scal_old(i,ns) = scal_new(i,ns)   
            enddo
            do ns = 1,Nspec
               I_R_old(i,ns) = I_R_new(i,ns) 
            enddo
            divu_old(i) = divu_new(i)
         enddo
         do i = 0,nx
            press_old(i)  = press_new(i)
         enddo
         time = time + dt

         call write_plt(vel_new,scal_new,press_new,divu_new,dx,dt,
     &                  nsteps_taken,time)
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
      
