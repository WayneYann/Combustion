
subroutine PROBINIT (init,name,namlen,problo,probhi)

  use probdata_module
  implicit none

  integer init, namlen
  integer name(namlen)
  double precision problo(2), probhi(2)

  integer untin,i

  namelist /fortin/ prob_type, pertmag, rfire, uinit, vinit, winit

!
!     Build "probin" filename -- the name of file containing fortin namelist.
!     
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
         
! set namelist defaults
  prob_type = 1

! problem type 1
  pertmag = 0.d0
  rfire   = 0.15d0
  uinit   = 0.d0
  vinit   = 0.d0
  winit   = 0.d0

!     Read namelists
  untin = 9
  open(untin,file=probin(1:namlen),form='formatted',status='old')
  read(untin,fortin)
  close(unit=untin)

  center(1) = 0.5d0*(problo(1) + probhi(1))
  center(2) = 0.5d0*(problo(2) + probhi(2))

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

  use eos_module, only : eos_given_PTY
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
  integer :: i, j, n, iwrk
  double precision :: xcen, ycen, r, rfront
  double precision :: pmf_vals(NSPEC+3), Xt(nspec), Yt(nspec)
  double precision :: rhot, et, Pt, Tt, u1t, u2t, rwrk

  if (nspecies .ne. NSPEC) then
     write(6,*)"nspecies, nspec ", nspecies, NSPEC
     stop
  end if

  do j = lo(2), hi(2)
     ycen = xlo(2) + delta(2)*(dble(j-lo(2)) + 0.5d0)


     do i = state_l1, state_h1
        xcen = xlo(1) + delta(1)*(dble(i-lo(1)) + 0.5d0)

        r = sqrt((xcen-center(1))**2 + (ycen-center(2))**2)
        rfront = rfire - r + 3.011d0 ! 3.011d0 is roughly the sufrace of fire for pmf.

        call pmf(rfront,rfront,pmf_vals,n)
     
        if (n .ne. nspec+3) then
           write(6,*)"n, nspec ", n, nspec
           stop
        end if
        
        Xt = pmf_vals(4:)

        Pt  = patm
        Tt  = pmf_vals(1)
        u1t = uinit
        u2t = vinit
     
        call eos_given_PTY(rhot, et, Yt, Pt, Tt, Xt)

        state(i,j,URHO ) = rhot
        state(i,j,UMX  ) = rhot*u1t
        state(i,j,UMX  ) = rhot*u2t
        state(i,j,UEDEN) = rhot*(et + 0.5d0*(u1t**2+u2t**2))
        state(i,j,UTEMP) = Tt

        do n=1,nspec
           state(i,j,UFS+n-1) = Yt(n)*rhot
        end do
     end do
  end do

end subroutine rns_initdata


! ::: -----------------------------------------------------------

     subroutine rns_hypfill(adv,adv_l1,adv_l2,adv_h1,adv_h2, &
                           domlo,domhi,delta,xlo,time,bc)
 
     use meth_params_module, only : NVAR

     implicit none
     include 'bc_types.fi'
     integer adv_l1,adv_l2,adv_h1,adv_h2
     integer bc(2,2,*)
     integer domlo(2), domhi(2)
     double precision delta(2), xlo(2), time
     double precision adv(adv_l1:adv_h1,adv_l2:adv_h2,NVAR)

     integer n

      do n = 1,NVAR
         call filcc(adv(adv_l1,adv_l2,n), &
              adv_l1,adv_l2,adv_h1,adv_h2, &
              domlo,domhi,delta,xlo,bc(1,1,n))
      enddo

      do n = 1,NVAR

!        XLO
         if ( bc(1,1,n).eq.EXT_DIR .and. adv_l1.lt.domlo(1)) then
            print *,'SHOULD NEVER GET HERE bc(1,1,n) .eq. EXT_DIR) '
            stop
         end if

!        XHI
         if ( bc(1,2,n).eq.EXT_DIR .and. adv_l1.lt.domlo(1)) then
            print *,'SHOULD NEVER GET HERE bc(1,2,n) .eq. EXT_DIR) '
            stop
         end if

!        YLO
         if ( bc(2,1,n).eq.EXT_DIR .and. adv_l2.lt.domlo(2)) then
            print *,'SHOULD NEVER GET HERE bc(2,1,n) .eq. EXT_DIR) '
            stop
         end if

!        YHI
         if ( bc(2,2,n).eq.EXT_DIR .and. adv_l2.lt.domlo(2)) then
            print *,'SHOULD NEVER GET HERE bc(2,2,n) .eq. EXT_DIR) '
            stop
         end if

      end do

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
