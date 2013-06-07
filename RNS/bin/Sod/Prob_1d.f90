
subroutine PROBINIT (init,name,namlen,problo,probhi)

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
subroutine cns_initdata(level,time,lo,hi,nscal, &
     state,state_l1,state_h1,delta,xlo,xhi)

end subroutine cns_initdata


! ::: 
! ::: -----------------------------------------------------------
! :::

     subroutine ca_hypfill(adv,adv_l1,adv_h1, &
                           domlo,domhi,delta,xlo,time,bc)

     use meth_params_module, only : NVAR
     implicit none
     include 'bc_types.fi'

     integer          :: adv_l1,adv_h1
     integer          :: bc(1,2,*)
     integer          :: domlo(1), domhi(1)
     double precision :: delta(1), xlo(1), time
     double precision :: adv(adv_l1:adv_h1,NVAR)

     integer n

     do n = 1,NVAR
        call filcc(adv(adv_l1,n),adv_l1,adv_h1, &
                   domlo,domhi,delta,xlo,bc(1,1,n))
     enddo

     do n = 1, NVAR

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

     end do

     end subroutine ca_hypfill

! ::: 
! ::: -----------------------------------------------------------
! :::

      subroutine ca_denfill(adv,adv_l1,adv_h1,domlo,domhi,delta,xlo,time,bc)

      implicit none
      include 'bc_types.fi'

      integer          :: adv_l1,adv_h1
      integer          :: bc(1,2,*)
      integer          :: domlo(1), domhi(1)
      double precision :: delta(1), xlo(1), time
      double precision :: adv(adv_l1:adv_h1)

!     Note: this function should not be needed, technically, but is provided
!     to filpatch because there are many times in the algorithm when just
!     the density is needed.  We try to rig up the filling so that the same
!     function is called here and in hypfill where all the states are filled.

      call filcc(adv,adv_l1,adv_h1,domlo,domhi,delta,xlo,bc)

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

      end subroutine ca_denfill
