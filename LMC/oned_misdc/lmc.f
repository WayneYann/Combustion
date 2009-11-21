      subroutine lmc()
      implicit none
      include 'spec.h'
      integer nx
      integer nsteps
      integer nsteps_taken
      integer at_nstep
      integer plot_int, chk_int
      integer num_init_iters
      integer num_divu_iters
      parameter (nx = 64)
      

c     Room for rho, rhoH, Temp, RhoRT + species (rho.Y)
      integer MAX_NSCAL
      parameter (MAX_NSCAL = maxspec + 4)

      real*8   vel_new(-1:nx  )
      real*8   vel_old(-1:nx  )
      real*8  scal_new(-1:nx  ,MAX_NSCAL)
      real*8  scal_old(-1:nx  ,MAX_NSCAL)
      real*8  scal_hold(-1:nx ,MAX_NSCAL)
      real*8 press_new(0 :nx  )
      real*8 press_old(0 :nx  )
      real*8  Ydot_new(0 :nx-1,maxspec)
      real*8  Ydot_old(0 :nx-1,maxspec)
      real*8   rhohalf(-1:nx  )
      real*8  divu_old(0 :nx-1)
      real*8  divu_new(0 :nx-1)
      real*8  beta_old(-1 :nx,MAX_NSCAL)
      real*8  beta_new(-1 :nx,MAX_NSCAL)
      real*8    mu_old(-1 :nx)
      real*8    mu_new(-1 :nx)
      real*8      dsdt(0 :nx-1)
      real*8     intra(0 :nx-1,MAX_NSCAL)


c     Local variables
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

      integer do_init, is, i, n, nd, ns
      
c     New arrays for MISDC.
      real*8    const_src(0 :nx-1,MAX_NSCAL)
      real*8  lin_src_old(0 :nx-1,MAX_NSCAL)
      real*8  lin_src_new(0 :nx-1,MAX_NSCAL)
      
      character chkfile*(16)
      real*8 Patm

      namelist /fortin/ nsteps,stop_time,cfl,
     $                  problo,probhi,chkfile,
     $                  plot_int, chk_int, change_max,
     $                  init_shrink, probtype,flame_offset,
     $                  dpdt_factor, Patm, coef_avg_harm,
     $                  misdc_iterMAX, predict_temp_for_coeffs,
     $                  num_divu_iters, num_init_iters

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
      coef_avg_harm = 1
      misdc_iterMAX = 3
      divu_ceiling_flag = 1
      divu_dt_factor    = 0.4d0
      rho_divu_ceiling  = 0.01
      predict_temp_for_coeffs = 1
      num_divu_iters = 3
      num_init_iters = 2

      open(9,file='probin',form='formatted',status='old')
      read(9,fortin)
      close(unit=9)
c      write(*,fortin)

      Pcgs = Patm * P1ATM
      
      dx = (probhi-problo)/DBLE(nx)
      nscal = Nspec+4

      call probinit(problo,probhi)
      
      if ( chkfile .ne. 'null') then

         print *,'CHKFILE ',chkfile
         
         call read_check(chkfile,nx,vel_new,scal_new,press_new,
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
         
c         call minmax_vel(nx,vel_new)
         
      else

          do_init = 1
          time = 0.d0
          at_nstep = 0

          call initdata(nx,vel_new,scal_new,Ydot_new,dx)
          call set_bc_grow_s(nx,scal_new,dx,time)
          do i = -1,nx
            do n = 1,nscal
              scal_old(i,n) = scal_new(i,n)
            enddo
          enddo

          call minmax_vel(nx,vel_new)

c         Define density for initial projection.
          do i = -1,nx
            rhohalf(i) = scal_old(i,Density)
          enddo

c     FIXME: Currently computes beta on bc using grow vals (see above)
c              Should rather use bc vals, and modify stencil for laplacians
          call calc_diffusivities(nx,scal_new,beta_new,mu_new)

          call calc_divu(nx,scal_new,beta_new,Ydot_new,divu_new,dx,time)

          call est_dt(nx,vel_new,scal_new,divu_new,dsdt,
     $                cfl,umax,dx,dt)

          print *,'initialVelocityProject: '
          dt_dummy = -1.d0

          do i=0,nx-1
             vel_old(i) = vel_new(i)
          enddo
          call project(nx,vel_old,vel_new,rhohalf,divu_new,
     $                 press_old,press_new,dx,dt_dummy,time)

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
             print *,' ...doing divu_iter number',nd,'dt=',dt

             call strang_chem(nx,scal_old,scal_new,
     $                        const_src,lin_src_old,lin_src_new,
     $                        intra,dt)

             do i = 0,nx-1
                do n = 1,nspec 
                   is = FirstSpec-1+n
                   Ydot_new(i,n) = 
     $                  (scal_new(i,is)-scal_old(i,is))/dt
                enddo
             enddo
             do i = 0,nx
                do n = 1,nscal
                   scal_new(i,n) = scal_hold(i,n)
                enddo
             enddo

             call calc_divu(nx,scal_new,beta_new,Ydot_new,
     &            divu_new,dx,time)

             print *,'initialVelocityProject: '
             dt_dummy = -1.d0

             call project(nx,vel_old,vel_new,rhohalf,divu_new,
     $            press_old,press_new,dx,dt_dummy,time)

             print *,' '
             print *,' '
             dt_init = dt

             call est_dt(nx,vel_new,scal_new,divu_new,dsdt,
     $                   cfl,umax,dx,dt_init2)
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
     $                   cfl,umax,dx,dt)
             dt = dt * init_shrink
             dt = min(dt,dt_init)
             print *,' '
             print *,'INITIAL PRESSURE ITERATION ',n
             print *,' '
             write(6,1001) time,dt

c     Here we zero out intra before each advance. FIXME: WHY?
             do ns = 1,nscal
                do i = 0,nx-1
                   intra(i,n) = 0.d0
                enddo
             enddo
             
             call advance(nx,vel_old,vel_new,scal_old,scal_new,
     $                    Ydot_old,Ydot_new,press_old,press_new,
     $                    divu_old,divu_new,dsdt,beta_old,beta_new,
     $                    intra,dx,dt,time)
             call minmax_vel(nx,vel_new)

             do i = 0,nx-1
                vel_new(i)  =   vel_old(i)
             enddo

             do i = 0,nx-1
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
 1001    format('Advancing: starting time = ',
     $        e15.9,' with dt = ',e15.9)
         
         print *,'COMPLETED INITIAL ITERATIONS'
         print *,' '
         
         cfl_used = cfl * init_shrink
         
c     Here we zero out intra before the "real" time stepping,
c     but not within advance, so we use a lagged intra each time.
c     FIXME: Why start on a bad note?
         do n = 1,nscal
            do i = 0,nx-1
               intra(i,n) = 0.d0
            enddo
         enddo
      endif
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
      
