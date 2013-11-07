program sdclmc
  use sdclib
  use iso_c_binding
  use encap
  use feval
  use lmc
  use probin

  implicit none

  type(lmc_encap), target :: U0
  type(sdc_encap)         :: enc
  type(sdc_imex)          :: imex
  type(sdc_nodes)         :: nodes

  character(256)   :: probin_fname
  integer          :: k, nstep, ierr
  double precision :: t

  ! init lmc

  if (command_argument_count() == 1) then
     call get_command_argument(1, value=probin_fname)
  else
     probin_fname = "probin.nml"
  end if
  call read_probin(probin_fname)
  call init_chem()
  call init_prob(prob_lo)

  ! init sdc
  call lmc_encap_build(enc, nx, nscal)
  call sdc_nodes_build(nodes, 3, SDC_GAUSS_LOBATTO, ierr)
  call sdc_imex_build(imex, nodes, c_funloc(f1eval), c_funloc(f2eval), c_funloc(f2comp), ierr)
  call sdc_imex_setup(imex, enc, c_null_ptr, ierr)
  call sdc_imex_allocate(imex, ierr)
  call lmc_encap_create_simple(U0, nx, nscal)

  ! initial condition
  call init_data(U0%vel, U0%scal, dx, lo, hi, bc)

  ! XXX: initial projection, initial divu and pressure


  call sdc_imex_set_q0(imex, c_loc(U0))
  call sdc_imex_spread(imex, t)

  ! time loop
  do nstep = 1, nsteps

     print *, "Current step: time:", t, ", dt:", dt, ", step:", nstep

     if (nsteps > 1) call sdc_imex_spread_qend(imex)

     ! sdc iters
     do k = 1, 8
        call sdc_imex_sweep(imex, t, dt, 0)
     end do
     call sdc_imex_get_qend(imex, c_loc(U0))

     t = t + dt

     ! if (MOD(nsteps_taken,plot_int).eq.0 .OR.
     !    &        nsteps_taken.eq.nsteps) then
     !    call write_plt(vel_new,scal_new,press_new,divu_new,I_R,
     !    $                     dx,nsteps_taken,time,lo,hi,bc)
     ! endif
     ! if (MOD(nsteps_taken,chk_int).eq.0 .OR.
     !    &        nsteps_taken.eq.nsteps) then
     !    call write_check(nsteps_taken,vel_new,scal_new,press_new,
     !    $                       I_R,divu_new,dx,time,dt,lo,hi)
     ! endif
  enddo

end program sdclmc
