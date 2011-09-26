
      subroutine PROBINIT (init,name,namlen,problo,probhi)

      use probdata_module
      implicit none

      integer init, namlen
      integer name(namlen)
      double precision problo(3), probhi(3)

      double precision vctr
      integer untin,i


      namelist /fortin/ zstandoff, pertmag, pAmb, idir, probtype, &
           denerr,  dengrad,  max_denerr_lev,  max_dengrad_lev, &
           velgrad,  max_velgrad_lev, &
           presserr,pressgrad,max_presserr_lev,max_pressgrad_lev

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


      idir = 1                ! direction across which to jump

      denerr = 1.d20
      dengrad = 1.d20
      max_denerr_lev = -1
      max_dengrad_lev = -1

      presserr = 1.d20
      pressgrad = 1.d20
      max_presserr_lev = -1
      max_pressgrad_lev = -1

      velgrad = 1.d20
      max_velgrad_lev = -1

      pAmb = 101325.0

!     Read namelists
      untin = 9
      open(untin,file=probin(1:namlen),form='formatted',status='old')
      read(untin,fortin)
      close(unit=untin)

! compute the internal energy (erg/cc) for the left and right state

      
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
     subroutine ca_initdata(level,time,lo,hi,nscal, &
        state,state_l1,state_l2,state_l3,state_h1,state_h2,state_h3, &
        delta,xlo,xhi)

     use cdwrk_module
     use probdata_module
     use meth_params_module, only : NVAR, URHO, UMX, UMY, UMZ, UEDEN, UEINT, UTEMP, UFS
     use prob_params_module
     implicit none

     integer level, nscal
     integer lo(3), hi(3)
     integer state_l1,state_l2,state_l3,state_h1,state_h2,state_h3
     double precision xlo(3), xhi(3), time, delta(3)
     double precision state(state_l1:state_h1,state_l2:state_h2, &
                            state_l3:state_h3,NVAR)

     double precision xcen,ycen,zcen

     double precision Xt(maxspec), Yt(maxspec)
     double precision z1,z2,Lx,Ly,pert
     double precision pt,rhot,u1t,u2t,u3t,Wavg,Tt,blend,et,Cvt
     double precision pmf_vals(maxspec+3)
     double precision pi
     integer i,j,k,n

      pi = 4.d0*atan2(1.d0,1.d0)
      do k = lo(3), hi(3)
         zcen = xlo(3) + delta(3)*(float(k-lo(3)) + 0.5d0)

         do j = lo(2), hi(2)
            ycen = xlo(2) + delta(2)*(float(j-lo(2)) + 0.5d0)

            do i = lo(1), hi(1)
               xcen = xlo(1) + delta(1)*(float(i-lo(1)) + 0.5d0)

               if(probtype .eq. 1 )then

                pert = 0.d0
                if (pertmag .gt. 0.d0) then
                   Lx = phys_prob_hi(1) - phys_prob_lo(1)
                   Ly = phys_prob_hi(2) - phys_prob_lo(2)
                   pert = pertmag*(1.000 * sin(2*Pi*4*xcen/Lx)             * sin(2*Pi*5*ycen/Ly)   &
                                    + 1.023 * sin(2*Pi*2*(xcen-.004598)/Lx)   * sin(2*Pi*4*(ycen-.0053765)/Ly)   &
                                    + 0.945 * sin(2*Pi*3*(xcen-.00712435)/Lx) * sin(2*Pi*3*(ycen-.02137)/Ly)   &
                                    + 1.017 * sin(2*Pi*5*(xcen-.0033)/Lx)     * sin(2*Pi*6*(ycen-.018)/Ly)   &
                                                + .982 * sin(2*Pi*5*(xcen-.014234)/Lx) )
                  endif

                  z1 = (zcen-0.d50*delta(3)-zstandoff+pert)*100.d0
                  z2 = (zcen+0.d50*delta(3)-zstandoff+pert)*100.d0

                  call pmf(z1,z2,pmf_vals,n)

                  if (n.ne.Nspec+3) then
                     write(6,*)"n,Nspec",n,Nspec
                     call bl_abort('INITDATA: n .ne. Nspec+3')
                  endif

                  Tt = pmf_vals(1)
                  do n = 1,Nspec
                     Xt(n) = pmf_vals(3+n)
                  end do
                  u1t = 0.d0
                  u2t = 0.d0
                  u3t = pmf_vals(2)*1.d-2
                  CALL CKXTY (Xt, IWRK, RWRK, Yt)
                  CALL CKRHOY(pAmb*1.d1,Tt,Yt,IWRK(ckbi),RWRK(ckbr),rhot)
                  rhot = rhot*1.e3
                  call CKUBMS(Tt,Yt,IWRK,RWRK,et)
                  et = et*1.d-4

                  state(i,j,k,URHO) = rhot
                  state(i,j,k,UTEMP) = Tt
                  state(i,j,k,UMX) = rhot*u1t
                  state(i,j,k,UMY) = rhot*u2t
                  state(i,j,k,UMZ) = rhot*u3t

                  state(i,j,k,UEDEN) = rhot*(et + 0.5d0*(u1t**2 + u2t**2 + u3t**2))

                  state(i,j,k,UEINT) = rhot*et

                  do n=1,Nspec
                     state(i,j,k,UFS+n-1) = Yt(n)*rhot
                  end do

              else

                  call bl_abort('unknown probtype in INITDATA')

              endif

            enddo
         enddo
      enddo



      end subroutine ca_initdata


! ::: -----------------------------------------------------------
      subroutine ca_hypfill(adv,adv_l1,adv_l2,adv_l3,adv_h1,adv_h2, &
           adv_h3,domlo,domhi,delta,xlo,time,bc)

     use meth_params_module, only : NVAR
     implicit none

     include 'bc_types.fi'
     integer adv_l1,adv_l2,adv_l3,adv_h1,adv_h2,adv_h3
     integer bc(3,2,*)
     integer domlo(3), domhi(3)
     double precision delta(3), xlo(3), time
     double precision adv(adv_l1:adv_h1,adv_l2:adv_h2,adv_l3:adv_h3,NVAR)

      integer n

      do n = 1,NVAR
         call filcc(adv(adv_l1,adv_l2,adv_l3,n), &
              adv_l1,adv_l2,adv_l3,adv_h1,adv_h2,adv_h3, &
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

!        YLO
         if ( bc(3,1,n).eq.EXT_DIR .and. adv_l3.lt.domlo(3)) then
            print *,'SHOULD NEVER GET HERE bc(3,1,n) .eq. EXT_DIR) '
            stop
         end if

!        YHI
         if ( bc(3,2,n).eq.EXT_DIR .and. adv_l3.lt.domlo(3)) then
            print *,'SHOULD NEVER GET HERE bc(3,2,n) .eq. EXT_DIR) '
            stop
         end if

      end do

      end subroutine ca_hypfill


      subroutine ca_denfill(adv,adv_l1,adv_l2,adv_l3,adv_h1,adv_h2, &
           adv_h3,domlo,domhi,delta,xlo,time,bc)

      implicit none

      include 'bc_types.fi'
      integer adv_l1,adv_l2,adv_l3,adv_h1,adv_h2,adv_h3
      integer bc(3,2,*)
      integer domlo(3), domhi(3)
      double precision delta(3), xlo(3), time
      double precision adv(adv_l1:adv_h1,adv_l2:adv_h2,adv_l3:adv_h3)

!     Note: this function should not be needed, technically, but is provided
!     to filpatch because there are many times in the algorithm when just
!     the density is needed.  We try to rig up the filling so that the same
!     function is called here and in hypfill where all the states are filled.

      call filcc(adv,adv_l1,adv_l2,adv_l3,adv_h1,adv_h2,adv_h3, &
                 domlo,domhi,delta,xlo,bc)

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

!     ZLO
      if ( bc(3,1,1).eq.EXT_DIR .and. adv_l3.lt.domlo(3)) then
         print *,'SHOULD NEVER GET HERE bc(3,1,1) .eq. EXT_DIR) '
         stop
      end if

!     ZHI
      if ( bc(3,2,1).eq.EXT_DIR .and. adv_l3.lt.domlo(3)) then
         print *,'SHOULD NEVER GET HERE bc(3,2,1) .eq. EXT_DIR) '
         stop
      end if

      end subroutine ca_denfill


! ::: -----------------------------------------------------------
 
      subroutine ca_reactfill(react,react_l1,react_l2,react_l3,&
           react_h1,react_h2,react_h3, &
           domlo,domhi,delta,xlo,time,bc)
 
      use probdata_module
      implicit none
      include 'bc_types.fi'
 
      integer :: react_l1,react_l2,react_l3,react_h1,react_h2,react_h3
      integer :: bc(3,2,*)
      integer :: domlo(3), domhi(3)
      double precision delta(3), xlo(3), time
      double precision react(react_l1:react_h1,react_l2:react_h2,react_l3:react_h3)
 
      call filcc(react,react_l1,react_l2,react_l3,react_h1,react_h2,react_h3,domlo,domhi,delta,xlo,bc)
 
      end subroutine ca_reactfill

