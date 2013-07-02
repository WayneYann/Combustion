
subroutine PROBINIT (init,name,namlen,problo,probhi)

  use probdata_module
  implicit none

  integer init, namlen
  integer name(namlen)
  double precision problo(2), probhi(2)

  xshock0_hi = xshock0_lo + (probhi(2)-problo(2))/tan(thetashock)
  vshock_x = vshock/sin(thetashock)

end subroutine PROBINIT

! ::: -----------------------------------------------------------
! ::: This routine is called at problem setup time and is used
! ::: to initialize data on each grid.  
! ::: 
! ::: NOTE:  all arrays have one cell of ghost zones surrounding
! :::        the grid interior.  Values in these cells need not
! :::        be set here.
! ::: 
! ::: INPUTS/OUTPUTS:
! ::: 
! ::: level     => amr level of grid
! ::: time      => time at which to init data             
! ::: lo,hi     => index limits of grid interior (cell centered)
! ::: nstate    => number of state components.  You should know
! :::		   this already!
! ::: state     <=  Scalar array
! ::: delta     => cell size
! ::: xlo,xhi   => physical locations of lower left and upper
! :::              right hand corner of grid.  (does not include
! :::		   ghost region).
! ::: -----------------------------------------------------------
subroutine rns_initdata(level,time,lo,hi,nscal, &
     state,state_l1,state_l2,state_h1,state_h2, &
     delta,xlo,xhi)

  use eos_module, only : gamma_const
  use probdata_module
  use meth_params_module, only : NVAR, URHO, UMX, UMY, UEDEN, UTEMP, UFS, NSPEC
  use chemistry_module, only : Patm, nspecies

  implicit none

  integer level, nscal
  integer lo(2), hi(2)
  integer state_l1,state_l2,state_h1,state_h2
  double precision xlo(2), xhi(2), time, delta(2)
  double precision state(state_l1:state_h1,state_l2:state_h2,NVAR)
  
  ! local variables
  integer :: i, j
  double precision :: xcen, ycen, ei0, ei1

  logical, save :: first_call = .true.

  if (first_call) then

     first_call = .false.

     ei0 = p0/((gamma_const-1.d0)*rho0)
     ei1 = p1/((gamma_const-1.d0)*rho1)

     allocate(state0(NVAR))
     allocate(state1(NVAR))

     state0(URHO)  = rho0
     state0(UMX)   =  rho0*v0*sin(thetashock)
     state0(UMY)   = -rho0*v0*cos(thetashock)
     state0(UEDEN) = rho0*(ei0 + 0.5d0*v0**2)
     state0(UTEMP) = 0.d0
     state0(UFS)   = 0.25d0 * rho0
     state0(UFS+1) = 0.75d0 * rho0

     state1(URHO)  = rho1
     state1(UMX)   =  rho1*v1*sin(thetashock)
     state1(UMY)   = -rho1*v1*cos(thetashock)
     state1(UEDEN) = rho1*(ei1 + 0.5d0*v1**2)
     state1(UTEMP) = 0.d0
     state1(UFS)   = 0.75d0 * rho1
     state1(UFS+1) = 0.25d0 * rho1

  end if

  if (NSPEC.ne. 2) then
     write(6,*)"nspec .ne. 2", NSPEC
     stop
  end if

  do j = state_l2, state_h2
     ycen = xlo(2) + delta(2)*(dble(j-lo(2)) + 0.5d0)

     do i = state_l1, state_h1
        xcen = xlo(1) + delta(1)*(dble(i-lo(1)) + 0.5d0)

        if (ycen .gt. tan(thetashock)*(xcen-xshock0_lo)) then
           state(i,j,:) = state1
        else
           state(i,j,:) = state0
        end if

     end do
  end do

end subroutine rns_initdata


! ::: -----------------------------------------------------------

     subroutine rns_hypfill(adv,adv_l1,adv_l2,adv_h1,adv_h2, &
                           domlo,domhi,delta,xlo,time,bc)
 
     use meth_params_module, only : NVAR
     use probdata_module

     implicit none
     include 'bc_types.fi'
     integer adv_l1,adv_l2,adv_h1,adv_h2
     integer bc(2,2,*)
     integer domlo(2), domhi(2)
     double precision delta(2), xlo(2), time
     double precision adv(adv_l1:adv_h1,adv_l2:adv_h2,NVAR)

     integer i, j, n
     double precision :: xcen, xshock

     do n = 1,NVAR
        call filcc(adv(adv_l1,adv_l2,n), &
             adv_l1,adv_l2,adv_h1,adv_h2, &
             domlo,domhi,delta,xlo,bc(1,1,n))
     enddo
     
     ! XLO
     if (adv_l1.lt.domlo(1)) then
        if (bc(1,1,1).eq.EXT_DIR) then
           do n=1,NVAR
              do j = adv_l2,adv_h2  ! fill the corners too
                 do i = adv_l1, domlo(1)-1
                    adv(i,j,n) = state1(n)
                 end do
              end do
           end do
        else
           print *,'SHOULD NEVER GET HERE bc(1,1,1) .ne. EXT_DIR) '
           stop
        end if
     end if
     
     ! YLO
     if (adv_l2 .lt. domlo(2)) then
        do n=1,NVAR
           do j = adv_l2, domlo(2)-1
              do i = adv_l1, adv_h1
                 xcen = delta(1)*(i + 0.5d0)
                 if (xcen < xshock0_lo) then
                    adv(i,j,n) = state1(n)
                 end if
              end do
           end do
        end do
     end if
     
     ! YHI
     if (adv_h2 .gt. domhi(2)) then
        xshock = xshock0_hi + vshock_x*time
        do n=1,NVAR
           do j = domhi(2)+1, adv_h2
              do i = adv_l1, adv_h1
                 xcen = delta(1)*(i + 0.5d0)
                 if (xcen < xshock) then
                    adv(i,j,n) = state1(n)
                 else
                    adv(i,j,n) = state0(n)
                 end if
              end do
           end do
        end do
     end if
     
   end subroutine rns_hypfill

! ::: -----------------------------------------------------------

      subroutine rns_denfill(adv,adv_l1,adv_l2,adv_h1,adv_h2, &
           domlo,domhi,delta,xlo,time,bc)

      implicit none
      include 'bc_types.fi'
      integer adv_l1,adv_l2,adv_h1,adv_h2
      integer bc(2,2,*)
      integer domlo(2), domhi(2)
      double precision delta(2), xlo(2), time
      double precision adv(adv_l1:adv_h1,adv_l2:adv_h2)

      print *, 'rns_denfill: SHOULD NEVER GET HERE'
      stop

!     Note: this function should not be needed, technically, but is provided
!     to filpatch because there are many times in the algorithm when just
!     the density is needed.  We try to rig up the filling so that the same
!     function is called here and in hypfill where all the states are filled.

      call filcc(adv,adv_l1,adv_l2,adv_h1,adv_h2,domlo,domhi,delta,xlo,bc)

!     XLO
      if ( bc(1,1,1).eq.EXT_DIR .and. adv_l1.lt.domlo(1)) then
         print *,'SHOULD NEVER GET HERE bc(1,1,1) .eq. EXT_DIR) '
         stop
      end if

!     XHI
      if ( bc(1,2,1).eq.EXT_DIR .and. adv_l1.lt.domlo(1)) then
         print *,'SHOULD NEVER GET HERE bc(1,2,1) .eq. EXT_DIR) '
         stop
      end if

!     YLO
      if ( bc(2,1,1).eq.EXT_DIR .and. adv_l2.lt.domlo(2)) then
         print *,'SHOULD NEVER GET HERE bc(2,1,1) .eq. EXT_DIR) '
         stop
      end if

!     YHI
      if ( bc(2,2,1).eq.EXT_DIR .and. adv_l2.lt.domlo(2)) then
         print *,'SHOULD NEVER GET HERE bc(2,2,1) .eq. EXT_DIR) '
         stop
      end if

      end subroutine rns_denfill
