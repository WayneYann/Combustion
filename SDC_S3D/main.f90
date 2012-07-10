program main

  use boxlib
  use pfasst
  use parallel
  use multifab_module
  use bl_IO_module
  use layout_module
  use init_data_module
  use write_plotfile_module
  use advance_module
  use bl_prof_module
  use sdcquad_module
  use mpi
  use feval

  implicit none

  !
  ! We only support 3-D.
  !
  integer, parameter :: DM = 3
  !
  ! We need four grow cells.
  !
  integer, parameter :: NG  = 4
  !
  ! We have five components (no species).
  !
  integer, parameter :: NC  = 5

  integer            :: nsteps, plot_int, n_cell, max_grid_size
  integer            :: un, farg, narg
  logical            :: need_inputs_file, found_inputs_file
  character(len=128) :: inputs_file_name
  integer            :: i, lo(DM), hi(DM), istep, method
  double precision   :: prob_lo(DM), prob_hi(DM), cfl, eta, alam
  double precision   :: dx(DM), dt, time, tfinal, start_time, end_time
  logical            :: is_periodic(DM)
  type(box)          :: bx
  type(boxarray)     :: ba
  type(layout)       :: la
  type(multifab)     :: U
  type(feval_ctx_t), target :: ctx
  type(sdcquad)      :: sdc
  type(pf_pfasst_t)  :: pf
  type(pf_comm_t)    :: tcomm
  integer            :: nvars(2), nnodes(2)

  type(bl_prof_timer), save :: bpt, bpt_init_data

  !
  ! What's settable via an inputs file.
  !
  namelist /probin/ tfinal, nsteps, dt, cfl, plot_int, &
       n_cell, max_grid_size, eta, alam, &
       method


  !
  ! Inititalize BoxLib.
  !
  call parallel_initialize(MPI_COMM_SELF)
  ! call boxlib_initialize()
  call bl_prof_initialize(on = .true.)

  call build(bpt, "bpt_main")

  start_time = parallel_wtime()


  !
  ! Namelist default values -- overwritable via inputs file.
  !
  tfinal        = 2.5d-7        ! Final time
  nsteps        = 100           ! Maximum number of time steps
  dt            = tfinal/nsteps ! Time step size (ignored if CFL > 0)
  cfl           = 0.5d0         ! Desired CFL number (use fixed size steps if CFL < 0)
  plot_int      = 10            ! Plot interval (time steps)
  n_cell        = 32            ! Number of grid cells per dimension
  max_grid_size = 32
  eta           = 1.8d-4        ! Diffusion coefficient
  alam          = 1.5d2         ! Diffusion coefficient
  method        = 1             ! Time-stepping scheme (1=RK3, 2=SDC, 3=PFASST)


  !
  ! Read inputs file and overwrite any default values.
  !
  narg = command_argument_count()
  need_inputs_file = .true.
  farg = 1
  if ( need_inputs_file .AND. narg >= 1 ) then
     call get_command_argument(farg, value = inputs_file_name)
     inquire(file = inputs_file_name, exist = found_inputs_file )
     if ( found_inputs_file ) then
        farg = farg + 1
        un = unit_new()
        open(unit=un, file = inputs_file_name, status = 'old', action = 'read')
        read(unit=un, nml = probin)
        close(unit=un)
        need_inputs_file = .false.
     end if
  end if

  if ( parallel_IOProcessor() ) then
     write(6,probin)
  end if


  !
  ! Physical problem is a box on (-1,-1) to (1,1), periodic on all sides.
  !
  prob_lo     = -0.1d0
  prob_hi     =  0.1d0
  is_periodic = .true.


  !
  ! Create a box from (0,0) to (n_cell-1,n_cell-1).
  !
  lo = 0
  hi = n_cell-1
  bx = make_box(lo,hi)

  do i = 1,DM
     dx(i) = (prob_hi(i)-prob_lo(i)) / n_cell
  end do

  call boxarray_build_bx(ba,bx)
  call boxarray_maxsize(ba,max_grid_size)
  call layout_build_ba(la,ba,boxarray_bbox(ba),pmask=is_periodic)
  call destroy(ba)


  !
  ! Compute initial condition.
  !
  call multifab_build(U,la,NC,NG)
  
  call build(bpt_init_data, "bpt_init_data")
  call init_data(U,dx,prob_lo,prob_hi)
  call destroy(bpt_init_data)

  istep  = 0
  time   = 0.d0

  if (plot_int > 0) then
     call write_plotfile(U,istep,dx,time,prob_lo,prob_hi)
  end if


  !
  ! Create SDC or PFASST context
  !
  if (method == 2) then 
     open(unit=un, file=inputs_file_name, status='old', action='read')
     call build(sdc, un)
     close(unit=un)
     call mk_imex_smats(sdc)
  end if

  if (method == 3) then
     if (.not. parallel_nprocs() > 1) then
        stop "ERROR: Not enough processors for PFASST run."
     end if

     ! XXX: take nvars from volume of initial condition...
     nvars  = NC * n_cell**DM
     nnodes = [ 5, 3 ] 

     call create(tcomm, MPI_COMM_WORLD)
     call create(pf, tcomm, 2, nvars, nnodes)

     pf%niters = 12
     pf%qtype  = 1
     pf%levels(2)%nsweeps = 2

     pf%dt   = dt
     pf%tend = pf%dt * tcomm%nproc

     ! XXX: check dt?

     call setup(tcomm, pf)
     call setup(pf)

     pf%levels(1)%ctx = c_loc(ctx)
     pf%levels(2)%ctx = c_loc(ctx)
  end if


  !
  ! Create feval context
  !
  ctx%la = la
  ctx%nc = nc
  ctx%ng = ng
  ctx%dx = dx

  ctx%eta  = eta
  ctx%alam = alam


  !
  ! Main time loop.
  !
  if (method == 1 .or. method == 2) then
     
     do istep=1,nsteps

        if (parallel_IOProcessor()) then
           print*,'Advancing time step',istep,'time = ',time
        end if

        call advance(U,dt,dx,cfl,time,tfinal,method,ctx,sdc)

        time = time + dt

        if (plot_int > 0) then
           if ( mod(istep,plot_int) .eq. 0 &
                .or. istep .eq. nsteps &
                .or. time >= tfinal ) then
              call write_plotfile(U,istep,dx,time,prob_lo,prob_hi)
           end if
        end if

        if (time >= tfinal) then
           exit
        end if

     end do

  else if (method == 3) then

     ! XXX: set PFASST initial condition

     print *, 'PFASST NOT IMPLEMENTED YET'

     ! call run(pf)

  end if


  !
  ! Destroy/finalize everything.
  !
  call destroy(U)
  call destroy(la)

  end_time = parallel_wtime()

  if ( parallel_IOProcessor() ) then
     print*,"Run time (s) =",end_time-start_time
  end if

  if (method == 2) then
     call destroy(sdc)
  else if (method == 3) then
     call destroy(pf)
     call destroy(tcomm)
  end if

  call destroy(bpt)
  call bl_prof_glean("bl_prof_report")
  call bl_prof_finalize()

  call boxlib_finalize()

end program main
