subroutine PROBINIT (init,name,namlen,problo,probhi)

  use probdata_module
  implicit none

  integer init, namlen
  integer name(namlen)
  double precision problo(3), probhi(3)

  integer untin,i

  namelist /fortin/ pamb, v_cf, T_cf, v_jet, T_jet, X_H2_jet, r_jet, &
       max_tracerr_lev, tracerr, max_vorterr_lev, vorterr, max_tempgrad_lev, tempgrad

  integer maxlen
  parameter (maxlen=256)
  character probin*(maxlen)
  
  if (namlen .gt. maxlen) then
     write(6,*) 'probin file name too long'
     stop
  end if

  do i = 1, namlen
     probin(i:i) = char(name(i))
  end do

  ! defaults
  pamb = 1.01325d6  ! 1 Patm
  v_cf = 5.5d3
  T_cf = 750.d0
  v_jet = 2.5d4
  T_jet = 420d0
  X_H2_jet = 0.7d0
  r_jet = 0.05d0

  ! Read namelists
  untin = 9
  open(untin,file=probin(1:namlen),form='formatted',status='old')
  read(untin,fortin)
  close(unit=untin)

end subroutine PROBINIT


subroutine rns_initdata(level,time,lo,hi,nscal, &
     state,state_l1,state_l2,state_l3,state_h1,state_h2,state_h3, &
     delta,xlo,xhi)

  use probdata_module, only : probdata_initialized, init_probdata, state_cf
  use meth_params_module, only : NVAR

  implicit none

  integer level, nscal
  integer lo(3), hi(3)
  integer state_l1,state_l2,state_l3,state_h1,state_h2,state_h3
  double precision xlo(3), xhi(3), time, delta(3)
  double precision state(state_l1:state_h1,state_l2:state_h2,state_l3:state_h3,NVAR)
  
  ! local variables
  integer :: i, j, k

  if (.not.probdata_initialized) call init_probdata()

  !$omp parallel do private(i,j,k)
  do k = lo(3), hi(3)
     do j = lo(2), hi(2)
        do i = lo(1), hi(1)
           state(i,j,k,:) = state_cf
        end do
     end do
  end do

end subroutine rns_initdata
