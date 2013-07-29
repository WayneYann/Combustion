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
      real*8, allocatable ::   vel_new(:,:)
      real*8, allocatable ::   vel_old(:,:)
      real*8, allocatable ::  scal_new(:,:,:)
      real*8, allocatable ::  scal_old(:,:,:)
      real*8, allocatable :: scal_hold(:,:,:)

!     cell-centered, 1 ghost cell
      real*8, allocatable ::      I_R(:,:,:)
      real*8, allocatable :: beta_old(:,:,:)
      real*8, allocatable :: beta_new(:,:,:)
      real*8, allocatable :: beta_for_Y_old(:,:,:)
      real*8, allocatable :: beta_for_Y_new(:,:,:)
      real*8, allocatable :: beta_for_Wbar_old(:,:,:)
      real*8, allocatable :: beta_for_Wbar_new(:,:,:)
      real*8, allocatable :: mu_dummy(:,:)
      real*8, allocatable :: divu_old(:,:)
      real*8, allocatable :: divu_new(:,:)

!     cell-centered, no ghost cells
      real*8, allocatable ::        dSdt(:,:)
      real*8, allocatable ::   delta_chi(:,:)
      real*8, allocatable ::   const_src(:,:,:)
      real*8, allocatable :: lin_src_old(:,:,:)
      real*8, allocatable :: lin_src_new(:,:,:)

!     nodal, 1 ghost cell
      real*8, allocatable :: press_new(:,:)
      real*8, allocatable :: press_old(:,:)

      integer, allocatable :: lo(:), hi(:), bc(:,:)

      real*8, allocatable :: dx(:), dt(:)

      real*8 problo,probhi
      real*8 time
      real*8 init_shrink
      real*8 stop_time
      real*8 fixed_dt
      real*8 Patm

      integer l,divu_iter,init_iter

      character chkfile*(16)

      namelist /fortin/ nx,nlevs,rr,subcycling,nsteps,stop_time,
     $                  problo,probhi,chkfile,
     $                  plot_int, chk_int,
     $                  init_shrink, flame_offset,
     $                  dpdt_factor, 
     $                  Patm, coef_avg_harm,
     $                  probtype,
     $                  misdc_iterMAX,
     $                  do_initial_projection, num_divu_iters, 
     $                  num_init_iters,fixed_dt,
     $                  V_in, lim_rxns,
     $                  LeEQ1, tranfile, TMIN_TRANS, Pr, Sc,
     $                  max_vode_subcycles,
     $                  min_vode_timestep, divu_ceiling_flag,
     $                  divu_dt_factor, rho_divu_ceiling, unlim

c     Set defaults, change with namelist
      nx = 256
      nlevs = 1
      rr = 2
      subcycling = .false.
      nsteps = 10
      stop_time = 1.e4
      problo = 0.0
      probhi = 3.5
      chkfile = 'null'
      plot_int = 1
      chk_int = 1
      init_shrink = 0.1d0
      flame_offset = 0.d0
      dpdt_factor = 0.d0
      Patm = 1.d0
      coef_avg_harm = 0
      misdc_iterMAX = 3
      do_initial_projection = 1
      num_divu_iters = 3
      num_init_iters = 2
      fixed_dt = -1.d0
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
      probtype = 1

      open(9,file='probin',form='formatted',status='old')
      read(9,fortin)
      close(unit=9)

      if (probtype.ne.1 .and. probtype.ne.2) then
         print *,'Unknown probtype:',probtype,'  Must be 1 or 2' 
         stop
      endif

      write(*,fortin)

c     number of cells at finest level
c     assumes rr is the same between all levels
      nfine = nx * rr**(nlevs-1)

c     Initialize chem/tran database and nspec
      call initchem()

      Pcgs = Patm * P1ATM

c     defines Density, Temp, RhoH, RhoRT, FirstSpec, LastSpec, nscal,
c     u_bc, T_bc, Y_bc, h_bc, and rho_bc
      call probinit(problo,probhi)

!     cell-centered, 2 ghost cells
      allocate(  vel_new(0:nlevs-1,-2:nfine+1))
      allocate(  vel_old(0:nlevs-1,-2:nfine+1))
      allocate( scal_new(0:nlevs-1,-2:nfine+1,nscal))
      allocate( scal_old(0:nlevs-1,-2:nfine+1,nscal))
      allocate(scal_hold(0:nlevs-1,-2:nfine+1,nscal))

!     cell-centered, 1 ghost cell
      allocate(     I_R(0:nlevs-1,-1:nfine,0:Nspec))
      allocate(beta_old(0:nlevs-1,-1:nfine,nscal))
      allocate(beta_new(0:nlevs-1,-1:nfine,nscal))
      allocate(beta_for_Y_old(0:nlevs-1,-1:nfine,nscal))
      allocate(beta_for_Y_new(0:nlevs-1,-1:nfine,nscal))
      allocate(beta_for_Wbar_old(0:nlevs-1,-1:nfine,nscal))
      allocate(beta_for_Wbar_new(0:nlevs-1,-1:nfine,nscal))
      allocate(mu_dummy(0:nlevs-1,-1:nfine))
      allocate(divu_old(0:nlevs-1,-1:nfine))
      allocate(divu_new(0:nlevs-1,-1:nfine))

!     cell-centered, no ghost cells
      allocate(       dSdt(0:nlevs-1,0:nfine-1))
      allocate(  delta_chi(0:nlevs-1,0:nfine-1))
      allocate(  const_src(0:nlevs-1,0:nfine-1,nscal))
      allocate(lin_src_old(0:nlevs-1,0:nfine-1,nscal))
      allocate(lin_src_new(0:nlevs-1,0:nfine-1,nscal))

!     nodal, 1 ghost cell
      allocate(press_new(0:nlevs-1,-1:nfine+1))
      allocate(press_old(0:nlevs-1,-1:nfine+1))

      allocate(lo(0:nlevs-1))
      allocate(hi(0:nlevs-1))
      allocate(bc(0:nlevs-1,2))

      allocate(dx(0:nlevs-1))
      allocate(dt(0:nlevs-1))

!     only need to zero these so plotfile has sensible data
      divu_old = 0.d0
      divu_new = 0.d0

!     must zero this or else RHS in mac project could be undefined
      dSdt = 0.d0
      delta_chi = 0.d0
      
!     initialize dx
      dx(0) = (probhi-problo)/DBLE(nx)
      do l=1,nlevs-1
         dx(l) = dx(l-1) / dble(rr)
      end do

!     initialize dt
      if (fixed_dt .le. 0.d0) then
         print*,'Error: must specify fixed_dt'
         stop
      else
         dt(0) = fixed_dt
         do l=1,nlevs-1
            if (subcycling) then
               dt(l) = dt(l-1) / dble(rr)
            else
               dt(l) = dt(l-1)
            end if
         end do
      end if

!     initialize lo and hi at each level
      lo(0) = 0
      hi(0) = nx-1
!     for now, the fine grid covers central 50% of domain and does not move
      if (nlevs .gt. 1) then
         lo(1) = (nx/4)*rr
         hi(1) = (3*nx/4)*rr - 1
      end if
      if (nlevs .gt. 2) then
         print*,'Error: grids only specified for nlevs = 2'
         stop
      end if

!     initialize boundary conditions
!     0=interior; 1=inflow; 2=outflow
      bc(0,1) = 1
      bc(0,2) = 2
      do l=1,nlevs-1
         bc(l,1:2) = 0
      end do
      
      if ( chkfile .ne. 'null') then

         print *,'CHKFILE ',chkfile
         
         call read_check(chkfile,vel_old,scal_old,press_old,
     $                   I_R,divu_old,dSdt,
     $                   time,at_nstep,dt,lo,hi)

         call write_plt(vel_old,scal_old,press_old,divu_old,I_R,
     $                  dx,at_nstep,time,lo,hi,bc)

         at_nstep = at_nstep + 1

c     needed for seed to EOS after first strang_chem call
         scal_new(:,:,Temp) = scal_old(:,:,Temp)
                  
      else
         
         time = 0.d0
         at_nstep = 1

C take vals from PMF and fills vel, Y, and Temp
C computes rho and h, fills in rhoH and rhoY
C sets I_R to zero

         call initdata(vel_old,scal_old,I_R,dx,lo,hi,bc)

c     needed for seed to EOS after first strang_chem call
         scal_new(:,:,Temp) = scal_old(:,:,Temp)

         press_old = 0.d0

         call write_plt(vel_old,scal_old,press_old,divu_old,I_R,
     &                  dx,99999,time,lo,hi,bc)

         do l=0,nlevs-1
            call calc_diffusivities(scal_old(l,:,:),beta_old(l,:,:),
     &                              beta_for_Y_old(l,:,:),
     &                              beta_for_Wbar_old(l,:,:),
     &                              mu_dummy(l,:),lo(l),hi(l))
         end do
         
         if (do_initial_projection .eq. 1) then

            print *,'initialVelocityProject: '
            do l=0,nlevs-1
               call calc_divu(scal_old(l,:,:),beta_old(l,:,:),
     &                        I_R(l,:,:),divu_old(l,:),dx(l),
     &                        lo(l),hi(l))
            end do

c     passing in dt=-1 ensures we simply project div(u)=S and
c     return zero pressure
            call project_level(vel_old(0,:),scal_old(:,0:,Density),
     $                         divu_old(0,:),press_old(0,:),
     $                         press_new(0,:),dx(0),-1.d0,
     $                         lo(0),hi(0),bc(0,:))

         end if

         call write_plt(vel_old,scal_old,press_old,divu_old,I_R,
     &                  dx,99998,time,lo,hi,bc)

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
            
            do l=0,nlevs-1
               call strang_chem(scal_old(l,:,:),scal_new(l,:,:),
     $                          const_src(l,:,:),lin_src_old(l,:,:),
     $                          lin_src_new(l,:,:),I_R(l,:,:),
     $                          0.5d0*dt(l),lo(l),hi(l),bc(l,:))
            end do

c     reset temperature just in case strang_chem call is not well poased
            scal_new(:,:,Temp) = scal_old(:,:,Temp)

            do l=0,nlevs-1
               call calc_divu(scal_old(l,:,:),beta_old(l,:,:),
     &                        I_R(l,:,:),divu_old(l,:),dx(l),
     &                        lo(l),hi(l))
            end do

            print *,'divu_iters velocity Project: '
            
c     passing in dt=-1 ensures we simply project div(u)=S and
c     return zero pressure
            call project_level(vel_old(0,:),scal_old(:,0:,Density),
     $                         divu_old(0,:),press_old(0,:),
     $                         press_new(0,:),dx(0),-1.d0,
     $                         lo(0),hi(0),bc(0,:))

         enddo

         call write_plt(vel_old,scal_old,press_old,divu_old,I_R,
     &                  dx,99997,time,lo,hi,bc)

         print *,' '
         print *,'...doing num_init_iters = ',num_init_iters 
         print *,' '
         if (num_init_iters .le. 0) then
            is_first_initial_iter = 0
         else
            is_first_initial_iter = 1
         endif
         do init_iter=1,num_init_iters

            doing_init_iters = 1

            print *,' '
            print *,'INITIAL PRESSURE ITERATION ',init_iter

            call advance(vel_old,vel_new,scal_old,scal_new,
     $                   I_R,press_old,press_new,
     $                   divu_old,divu_new,dSdt,beta_old,beta_new,
     $                   beta_for_Y_old,beta_for_Y_new,
     $                   beta_for_Wbar_old,beta_for_Wbar_new,
     $                   dx,dt,lo,hi,bc,delta_chi,-init_iter)

c     update pressure and I_R
            press_old = press_new

            is_first_initial_iter = 0          

         enddo

         doing_init_iters = 0

         call write_plt(vel_old,scal_old,press_old,divu_old,I_R,
     &                  dx,0,time,lo,hi,bc)

         print *,' '      
         print *,' '      
         print *,'COMPLETED INITIAL ITERATIONS'
         print *,' '      
         print *,'START ADVANCING THE SOLUTION '
         print *,' '            

 1001    format('Advancing: starting time = ',
     $        e15.9,' with dt = ',e15.9)

      endif

      call write_plt(vel_new,scal_new,press_new,divu_new,I_R,
     $     dx,at_nstep,time,lo,hi,bc)
      call write_check(at_nstep,vel_new,scal_new,press_new,
     $     I_R,divu_new,dSdt,dx,time,dt,lo,hi)

C-- Now advance 
      do nsteps_taken = at_nstep, nsteps

         if (time.ge.stop_time) exit

         dt = min(dt,stop_time-time)

         write(6,*)
         write(6,1001 )time,dt
         write(6,*)'STEP = ',nsteps_taken
         
         call advance(vel_old,vel_new,scal_old,scal_new,
     $                I_R,press_old,press_new,
     $                divu_old,divu_new,dSdt,beta_old,beta_new,
     $                beta_for_Y_old,beta_for_Y_new,
     $                beta_for_Wbar_old,beta_for_Wbar_new,
     $                dx,dt,lo,hi,bc,delta_chi,nsteps_taken)

c     update state, time
         vel_old = vel_new
         scal_old = scal_new
         divu_old = divu_new
         press_old = press_new

         time = time + dt(0)

         if (MOD(nsteps_taken,plot_int).eq.0 .OR. 
     &        nsteps_taken.eq.nsteps) then 
            call write_plt(vel_new,scal_new,press_new,divu_new,I_R,
     $                     dx,nsteps_taken,time,lo,hi,bc)
         endif
         if (MOD(nsteps_taken,chk_int).eq.0 .OR.
     &        nsteps_taken.eq.nsteps) then 
            call write_check(nsteps_taken,vel_new,scal_new,press_new,
     $                       I_R,divu_new,dSdt,dx,time,dt,lo,hi)
         endif
      enddo

      print *,' '      
      print *,'COMPLETED SUCCESSFULLY'
      print *,' '      

      end
