
subroutine PROBINIT (init,name,namlen,problo,probhi)

  use probdata_module
  implicit none

  integer init, namlen
  integer name(namlen)
  double precision problo(2), probhi(2)

  integer untin,i

  namelist /fortin/ Minfty, Rvortex, Cvortex

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
  Minfty  = 0.05d0
  Rvortex = 0.1d0
  Cvortex = -0.0025d0
  
!     Read namelists
  untin = 9
  open(untin,file=probin(1:namlen),form='formatted',status='old')
  read(untin,fortin)
  close(unit=untin)

  Lx = probhi(1) - problo(1)

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
  use chemistry_module, only : Patm, Ru, nspecies

  implicit none

  integer level, nscal
  integer lo(2), hi(2)
  integer state_l1,state_l2,state_h1,state_h2
  double precision xlo(2), xhi(2), time, delta(2)
  double precision state(state_l1:state_h1,state_l2:state_h2,NVAR)
  
  ! local variables
  integer :: i, j, ii, jj, n, iwrk
  double precision :: xcen, ycen, xg, yg
  double precision :: pmf_vals(NSPEC+3), Xt(nspec), Yt(nspec)
  double precision :: rhot, et, Pt, Tt, u1t, u2t, rwrk, Cv, Cp, gamma, cs, exptmp, &
       Rc, Cc, uinfty, Wbar

  double precision, parameter :: gp(2) = (/ -1.d0/sqrt(3.d0), 1.d0/sqrt(3.d0) /)

  if (nspecies .ne. NSPEC) then
     write(6,*)"nspecies, nspec ", nspecies, NSPEC
     stop
  end if
  
  call pmf(-1.d0,-1.d0,pmf_vals,n)

  if (n .ne. nspec+3) then
     write(6,*)"n, nspec ", n, nspec
     stop
  end if
  
  Tt = pmf_vals(1)

  do n = 1,nspecies
     Xt(n) = pmf_vals(3+n)
  end do
  CALL CKXTY (Xt, IWRK, RWRK, Yt)
  CALL CKRHOY(patm,Tt,Yt,IWRK,RWRK,rhot)
  ! we now have rho
  
  call CKCVBL(Tt, Xt, iwrk, rwrk, Cv)
  Cp = Cv + Ru
  gamma = Cp / Cv
  cs = sqrt(gamma*patm/rhot)
  
  uinfty = Minfty*cs
  Rc = Rvortex*Lx
  Cc = Cvortex*cs*Lx
    
  do j = state_l2, state_h2
     ycen = xlo(2) + delta(2)*(dble(j-lo(2)) + 0.5d0)
     
     do i = state_l1, state_h1
        xcen = xlo(1) + delta(1)*(dble(i-lo(1)) + 0.5d0)
        
        do n=1,NVAR
           state(i,j,n) = 0.d0
        end do
        
        do jj = 1, 2
           yg = ycen + 0.5d0*delta(2)*gp(jj)
           do ii = 1, 2
              xg = xcen + 0.5d0*delta(1)*gp(ii)
              
              exptmp = exp(-(xg**2+yg**2)/(2.d0*Rc**2))
              
              u1t = uinfty - Cc*exptmp*yg/Rc**2
              u2t = Cc*exptmp*xg/Rc**2
              
              Pt = patm - rhot*(Cc/Rc)**2*exptmp
              
              call CKMMWX(Xt, iwrk, rwrk, Wbar)
              Tt = Pt*Wbar / (rhot*Ru)
              
              call CKUBMS(Tt,Yt,IWRK,RWRK,et)
              
              state(i,j,URHO ) = state(i,j,URHO ) + 0.25d0*rhot
              state(i,j,UMX  ) = state(i,j,UMX  ) + 0.25d0*rhot*u1t
              state(i,j,UMY  ) = state(i,j,UMY  ) + 0.25d0*rhot*u2t
              state(i,j,UEDEN) = state(i,j,UEDEN) + 0.25d0*rhot*(et + 0.5d0*(u1t**2+u2t**2))
              state(i,j,UTEMP) = state(i,j,UTEMP) + 0.25d0*Tt
              do n=1, NSPEC
                 state(i,j,UFS+n-1) = state(i,j,UFS+n-1) + 0.25d0*rhot*Yt(n)
              end do
              
              if (.not. sinftysaved) then
                 Tt = pmf_vals(1)
                 call CKUBMS(Tt,Yt,IWRK,RWRK,et)
                 
                 allocate(stateinfty(NVAR))
                 
                 stateinfty(URHO ) = rhot
                 stateinfty(UMX  ) = rhot*uinfty
                 stateinfty(UMY  ) = 0.d0
                 stateinfty(UEDEN) = rhot*(et + 0.5d0*uinfty**2)
                 stateinfty(UTEMP) = Tt
                 do n=1,NSPEC
                    stateinfty(UFS+n-1) = rhot*Yt(n)
                 end do
                 
                 sinftysaved = .true.
              end if
              
           end do
        end do

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

      do n = 1,NVAR
         call filcc(adv(adv_l1,adv_l2,n), &
              adv_l1,adv_l2,adv_h1,adv_h2, &
              domlo,domhi,delta,xlo,bc(1,1,n))
      enddo

!        XLO
      if (adv_l1.lt.domlo(1)) then
         if (bc(1,1,1).eq.EXT_DIR) then
            do n=1,NVAR
               do j = adv_l2,adv_h2  ! fill the corners too
                  do i = adv_l1, domlo(1)-1
                     adv(i,j,n) = stateinfty(n)
                  end do
               end do
            end do
         else
            print *,'SHOULD NEVER GET HERE bc(1,1,1) .ne. EXT_DIR) '
            stop
         end if
      end if
      
      do n = 1, NVAR
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

      print *, 'rns_denfill: SHOULD NEVER GET HERE'
      stop

      call filcc(adv,adv_l1,adv_l2,adv_h1,adv_h2,domlo,domhi,delta,xlo,bc)

      !     XLO
      if ( bc(1,1,1).eq.EXT_DIR .and. adv_l1.lt.domlo(1)) then
         print *,'SHOULD NEVER GET HERE bc(1,1,1) .eq. EXT_DIR) '
         stop
      end if
      
      !     XHI
      if ( bc(1,2,1).eq.EXT_DIR .and. adv_h1.gt.domhi(1)) then
         print *,'SHOULD NEVER GET HERE bc(1,2,1) .eq. EXT_DIR) '
         stop
      end if
      
      !     YLO
      if ( bc(2,1,1).eq.EXT_DIR .and. adv_l2.lt.domlo(2)) then
         print *,'SHOULD NEVER GET HERE bc(2,1,1) .eq. EXT_DIR) '
         stop
      end if
      
      !     YHI
      if ( bc(2,2,1).eq.EXT_DIR .and. adv_h2.gt.domhi(2)) then
         print *,'SHOULD NEVER GET HERE bc(2,2,1) .eq. EXT_DIR) '
         stop
      end if

      end subroutine rns_denfill
