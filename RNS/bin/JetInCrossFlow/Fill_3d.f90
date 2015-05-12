! Fill the entire state
subroutine rns_grpfill(adv,adv_l1,adv_l2,adv_l3,adv_h1,adv_h2,adv_h3, &
     domlo,domhi,delta,xlo,time,bc_in)
 
  use probdata_module, only : state_cf, state_jet, r_jet, T_cf
  use meth_params_module, only : NVAR, URHO, UMX, UMY, UMZ, UEDEN, UTEMP, UFS, NSPEC
  use sdc_boundary_module, only : isFEval
  
  implicit none
  include 'bc_types.fi'
  integer adv_l1,adv_l2,adv_l3,adv_h1,adv_h2,adv_h3
  integer bc_in(3,2,*)
  integer domlo(3), domhi(3)
  double precision delta(3), xlo(3), time
  double precision adv(adv_l1:adv_h1,adv_l2:adv_h2,adv_l3:adv_h3,NVAR)
  
  integer n, i, j, k, ii, kk, iwrk
  integer bc(3,2,NVAR)
  double precision :: x, z, rhoInv, Yt(nspec), rwrk, et
  double precision, parameter :: gp(2) = (/ -1.d0/sqrt(3.d0), 1.d0/sqrt(3.d0) /)
  double precision, parameter :: wgt = 0.25d0

  bc = bc_in(:,:,1:NVAR)
  if (isFEval) then
     where (bc .eq. EXT_DIR) bc = FOEXTRAP
  end if

  do n = 1,NVAR
     call filcc(adv(adv_l1,adv_l2,adv_l3,n), &
          adv_l1,adv_l2,adv_l3,adv_h1,adv_h2,adv_h3, &
          domlo,domhi,delta,xlo,bc(1,1,n))
  enddo

  if (isFEval) return

  ! XLO: inflow
  if (adv_l1.lt.domlo(1)) then
     if (bc(1,1,1).eq.EXT_DIR) then
        do k = adv_l3, adv_h3
           do j = adv_l2, adv_h2
              do i = adv_l1, domlo(1)-1
                 adv(i,j,k,:) = state_cf
              end do
           end do
        end do
     else
        print *,'grpfill: SHOULD NEVER GET HERE bc(1,1,1) .ne. EXT_DIR'
        stop        
     end if
  end if

  ! YLO: Wall
  if (adv_l2.lt.domlo(2)) then
     do k = adv_l3, adv_h3
        z = (DBLE(k-adv_l3)+.5d0)*delta(3)+xlo(3)
        do j = adv_l2, domlo(2)-1
           do i = adv_l1, adv_h1
              x = (DBLE(i-adv_l1)+.5d0)*delta(1)+xlo(1)
              if (x**2+z**2 .le. r_jet**2) then
                 adv(i,j,k,:) = state_jet
              else
                 ! fixed T
                 adv(i,j,k,UTEMP) = T_cf
                 rhoInv = 1.d0/adv(i,j,k,URHO)
                 Yt = adv(i,j,k,UFS:UFS+NSPEC-1)*rhoInv
                 call CKUBMS(T_cf, Yt, iwrk, rwrk, et)
                 adv(i,j,k,UEDEN) = adv(i,j,k,URHO)*et  &
                      + 0.5d0*rhoInv*(adv(i,j,k,UMX)**2+adv(i,j,k,UMY)**2+adv(i,j,k,UMZ)**2)
              end if
           end do
        end do
     end do
  end if

  do n = 1, NVAR

     !     XLO
!     if ( bc(1,1,n).eq.EXT_DIR .and. adv_l1.lt.domlo(1)) then
!        print *,'grpfill: SHOULD NEVER GET HERE bc(1,1,n) .eq. EXT_DIR) '
!        stop
!     end if
     
     !     XHI
     if ( bc(1,2,n).eq.EXT_DIR .and. adv_h1.gt.domhi(1)) then
        print *,'grpfill: SHOULD NEVER GET HERE bc(1,2,n) .eq. EXT_DIR) '
        stop
     end if
     
     !     YLO
!     if ( bc(2,1,n).eq.EXT_DIR .and. adv_l2.lt.domlo(2)) then
!        print *,'grpfill: SHOULD NEVER GET HERE bc(2,1,n) .eq. EXT_DIR) '
!        stop
!     end if
     
     !     YHI
     if ( bc(2,2,n).eq.EXT_DIR .and. adv_h2.gt.domhi(2)) then
        print *,'grpfill: SHOULD NEVER GET HERE bc(2,2,n) .eq. EXT_DIR) '
        stop
     end if

     !     ZLO
     if ( bc(3,1,n).eq.EXT_DIR .and. adv_l3.lt.domlo(3)) then
        print *,'grpfill: SHOULD NEVER GET HERE bc(3,1,n) .eq. EXT_DIR) '
        stop
     end if
     
     !     ZHI
     if ( bc(3,2,n).eq.EXT_DIR .and. adv_h3.gt.domhi(3)) then
        print *,'grpfill: SHOULD NEVER GET HERE bc(3,2,n) .eq. EXT_DIR) '
        stop
     end if
     
  end do
  
end subroutine rns_grpfill

! Fill one variable
subroutine rns_regfill(adv,adv_l1,adv_l2,adv_l3,adv_h1,adv_h2,adv_h3, &
     domlo,domhi,delta,xlo,time,bc)
  
  implicit none
  include 'bc_types.fi'
  integer adv_l1,adv_l2,adv_l3,adv_h1,adv_h2,adv_h3
  integer bc(3,2,*)
  integer domlo(3), domhi(3)
  double precision delta(3), xlo(3), time
  double precision adv(adv_l1:adv_h1,adv_l2:adv_h2,adv_l3:adv_h3)

  call filcc(adv,adv_l1,adv_l2,adv_l3,adv_h1,adv_h2,adv_h3, &
       domlo,domhi,delta,xlo,bc)
  
  !     XLO
  if ( bc(1,1,1).eq.EXT_DIR .and. adv_l1.lt.domlo(1)) then
     print *,'regfill: SHOULD NEVER GET HERE bc(1,1,1) .eq. EXT_DIR) '
     stop
  end if
      
  !     XHI
  if ( bc(1,2,1).eq.EXT_DIR .and. adv_h1.gt.domhi(1)) then
     print *,'regfill: SHOULD NEVER GET HERE bc(1,2,1) .eq. EXT_DIR) '
     stop
  end if
  
  !     YLO
  if ( bc(2,1,1).eq.EXT_DIR .and. adv_l2.lt.domlo(2)) then
     print *,'regfill: SHOULD NEVER GET HERE bc(2,1,1) .eq. EXT_DIR) '
     stop
  end if
  
  !     YHI
  if ( bc(2,2,1).eq.EXT_DIR .and. adv_h2.gt.domhi(2)) then
     print *,'regfill: SHOULD NEVER GET HERE bc(2,2,1) .eq. EXT_DIR) '
     stop
  end if

  !     ZLO
  if ( bc(3,1,1).eq.EXT_DIR .and. adv_l3.lt.domlo(3)) then
     print *,'regfill: SHOULD NEVER GET HERE bc(3,1,1) .eq. EXT_DIR) '
     stop
  end if
  
  !     ZHI
  if ( bc(3,2,1).eq.EXT_DIR .and. adv_h3.gt.domhi(3)) then
     print *,'regfill: SHOULD NEVER GET HERE bc(3,2,1) .eq. EXT_DIR) '
     stop
  end if
  
end subroutine rns_regfill


! Fill density
subroutine rns_denfill(adv,adv_l1,adv_l2,adv_l3,adv_h1,adv_h2,adv_h3, &
     domlo,domhi,delta,xlo,time,bc)
  
  implicit none
  include 'bc_types.fi'
  integer adv_l1,adv_l2,adv_l3,adv_h1,adv_h2,adv_h3
  integer bc(3,2,*)
  integer domlo(3), domhi(3)
  double precision delta(3), xlo(3), time
  double precision adv(adv_l1:adv_h1,adv_l2:adv_h2,adv_l3:adv_h3)

  call filcc(adv,adv_l1,adv_l2,adv_l3,adv_h1,adv_h2,adv_h3, &
       domlo,domhi,delta,xlo,bc)
  
  !     XLO
  if ( bc(1,1,1).eq.EXT_DIR .and. adv_l1.lt.domlo(1)) then
     print *,'denfill: SHOULD NEVER GET HERE bc(1,1,1) .eq. EXT_DIR) '
     stop
  end if
  
  !     XHI
  if ( bc(1,2,1).eq.EXT_DIR .and. adv_h1.gt.domhi(1)) then
     print *,'denfill: SHOULD NEVER GET HERE bc(1,2,1) .eq. EXT_DIR) '
     stop
  end if
  
  !     YLO
  if ( bc(2,1,1).eq.EXT_DIR .and. adv_l2.lt.domlo(2)) then
     print *,'denfill: SHOULD NEVER GET HERE bc(2,1,1) .eq. EXT_DIR) '
     stop
  end if
  
  !     YHI
  if ( bc(2,2,1).eq.EXT_DIR .and. adv_h2.gt.domhi(2)) then
     print *,'denfill: SHOULD NEVER GET HERE bc(2,2,1) .eq. EXT_DIR) '
     stop
  end if

  !     ZLO
  if ( bc(3,1,1).eq.EXT_DIR .and. adv_l3.lt.domlo(3)) then
     print *,'denfill: SHOULD NEVER GET HERE bc(3,1,1) .eq. EXT_DIR) '
     stop
  end if
  
  !     ZHI
  if ( bc(3,2,1).eq.EXT_DIR .and. adv_h3.gt.domhi(3)) then
     print *,'denfill: SHOULD NEVER GET HERE bc(3,2,1) .eq. EXT_DIR) '
     stop
  end if
  
end subroutine rns_denfill


! Fill x-momentum
subroutine rns_mxfill(adv,adv_l1,adv_l2,adv_l3,adv_h1,adv_h2,adv_h3, &
     domlo,domhi,delta,xlo,time,bc)
  
  implicit none
  include 'bc_types.fi'
  integer adv_l1,adv_l2,adv_l3,adv_h1,adv_h2,adv_h3
  integer bc(3,2,*)
  integer domlo(3), domhi(3)
  double precision delta(3), xlo(3), time
  double precision adv(adv_l1:adv_h1,adv_l2:adv_h2,adv_l3:adv_h3)
  
  call filcc(adv,adv_l1,adv_l2,adv_l3,adv_h1,adv_h2,adv_h3, &
       domlo,domhi,delta,xlo,bc)
  
  !     XLO
  if ( bc(1,1,1).eq.EXT_DIR .and. adv_l1.lt.domlo(1)) then
     print *,'mxfill: SHOULD NEVER GET HERE bc(1,1,1) .eq. EXT_DIR) '
     stop
  end if
  
  !     XHI
  if ( bc(1,2,1).eq.EXT_DIR .and. adv_h1.gt.domhi(1)) then
     print *,'mxfill: SHOULD NEVER GET HERE bc(1,2,1) .eq. EXT_DIR) '
     stop
  end if
  
  !     YLO
  if ( bc(2,1,1).eq.EXT_DIR .and. adv_l2.lt.domlo(2)) then
     print *,'mxfill: SHOULD NEVER GET HERE bc(2,1,1) .eq. EXT_DIR) '
     stop
  end if
  
  !     YHI
  if ( bc(2,2,1).eq.EXT_DIR .and. adv_h2.gt.domhi(2)) then
     print *,'mxfill: SHOULD NEVER GET HERE bc(2,2,1) .eq. EXT_DIR) '
     stop
  end if

  !     ZLO
  if ( bc(3,1,1).eq.EXT_DIR .and. adv_l3.lt.domlo(3)) then
     print *,'mxfill: SHOULD NEVER GET HERE bc(3,1,1) .eq. EXT_DIR) '
     stop
  end if
  
  !     ZHI
  if ( bc(3,2,1).eq.EXT_DIR .and. adv_h3.gt.domhi(3)) then
     print *,'mxfill: SHOULD NEVER GET HERE bc(3,2,1) .eq. EXT_DIR) '
     stop
  end if
  
end subroutine rns_mxfill


! Fill y-momentum
subroutine rns_myfill(adv,adv_l1,adv_l2,adv_l3,adv_h1,adv_h2,adv_h3, &
     domlo,domhi,delta,xlo,time,bc)
  
  implicit none
  include 'bc_types.fi'
  integer adv_l1,adv_l2,adv_l3,adv_h1,adv_h2,adv_h3
  integer bc(3,2,*)
  integer domlo(3), domhi(3)
  double precision delta(3), xlo(3), time
  double precision adv(adv_l1:adv_h1,adv_l2:adv_h2,adv_l3:adv_h3)

  call filcc(adv,adv_l1,adv_l2,adv_l3,adv_h1,adv_h2,adv_h3, &
       domlo,domhi,delta,xlo,bc)
  
  !     XLO
  if ( bc(1,1,1).eq.EXT_DIR .and. adv_l1.lt.domlo(1)) then
     print *,'myfill: SHOULD NEVER GET HERE bc(1,1,1) .eq. EXT_DIR) '
     stop
  end if
  
  !     XHI
  if ( bc(1,2,1).eq.EXT_DIR .and. adv_h1.gt.domhi(1)) then
     print *,'myfill: SHOULD NEVER GET HERE bc(1,2,1) .eq. EXT_DIR) '
     stop
  end if
  
  !     YLO
  if ( bc(2,1,1).eq.EXT_DIR .and. adv_l2.lt.domlo(2)) then
     print *,'myfill: SHOULD NEVER GET HERE bc(2,1,1) .eq. EXT_DIR) '
     stop
  end if
  
  !     YHI
  if ( bc(2,2,1).eq.EXT_DIR .and. adv_h2.gt.domhi(2)) then
     print *,'myfill: SHOULD NEVER GET HERE bc(2,2,1) .eq. EXT_DIR) '
     stop
  end if

  !     ZLO
  if ( bc(3,1,1).eq.EXT_DIR .and. adv_l3.lt.domlo(3)) then
     print *,'myfill: SHOULD NEVER GET HERE bc(3,1,1) .eq. EXT_DIR) '
     stop
  end if
  
  !     ZHI
  if ( bc(3,2,1).eq.EXT_DIR .and. adv_h3.gt.domhi(3)) then
     print *,'myfill: SHOULD NEVER GET HERE bc(3,2,1) .eq. EXT_DIR) '
     stop
  end if
  
end subroutine rns_myfill


! Fill z-momentum
subroutine rns_mzfill(adv,adv_l1,adv_l2,adv_l3,adv_h1,adv_h2,adv_h3, &
     domlo,domhi,delta,xlo,time,bc)
  
  implicit none
  include 'bc_types.fi'
  integer adv_l1,adv_l2,adv_l3,adv_h1,adv_h2,adv_h3
  integer bc(3,2,*)
  integer domlo(3), domhi(3)
  double precision delta(3), xlo(3), time
  double precision adv(adv_l1:adv_h1,adv_l2:adv_h2,adv_l3:adv_h3)

  call filcc(adv,adv_l1,adv_l2,adv_l3,adv_h1,adv_h2,adv_h3, &
       domlo,domhi,delta,xlo,bc)
  
  !     XLO
  if ( bc(1,1,1).eq.EXT_DIR .and. adv_l1.lt.domlo(1)) then
     print *,'mzfill: SHOULD NEVER GET HERE bc(1,1,1) .eq. EXT_DIR) '
     stop
  end if
  
  !     XHI
  if ( bc(1,2,1).eq.EXT_DIR .and. adv_h1.gt.domhi(1)) then
     print *,'mzfill: SHOULD NEVER GET HERE bc(1,2,1) .eq. EXT_DIR) '
     stop
  end if
  
  !     YLO
  if ( bc(2,1,1).eq.EXT_DIR .and. adv_l2.lt.domlo(2)) then
     print *,'mzfill: SHOULD NEVER GET HERE bc(2,1,1) .eq. EXT_DIR) '
     stop
  end if
  
  !     YHI
  if ( bc(2,2,1).eq.EXT_DIR .and. adv_h2.gt.domhi(2)) then
     print *,'mzfill: SHOULD NEVER GET HERE bc(2,2,1) .eq. EXT_DIR) '
     stop
  end if

  !     ZLO
  if ( bc(3,1,1).eq.EXT_DIR .and. adv_l3.lt.domlo(3)) then
     print *,'mzfill: SHOULD NEVER GET HERE bc(3,1,1) .eq. EXT_DIR) '
     stop
  end if
  
  !     ZHI
  if ( bc(3,2,1).eq.EXT_DIR .and. adv_h3.gt.domhi(3)) then
     print *,'mzfill: SHOULD NEVER GET HERE bc(3,2,1) .eq. EXT_DIR) '
     stop
  end if
  
end subroutine rns_mzfill


! Fill temperature
subroutine rns_tempfill(adv,adv_l1,adv_l2,adv_l3,adv_h1,adv_h2,adv_h3, &
     domlo,domhi,delta,xlo,time,bc)
  
  implicit none
  include 'bc_types.fi'
  integer adv_l1,adv_l2,adv_l3,adv_h1,adv_h2,adv_h3
  integer bc(3,2,*)
  integer domlo(3), domhi(3)
  double precision delta(3), xlo(3), time
  double precision adv(adv_l1:adv_h1,adv_l2:adv_h2,adv_l3:adv_h3)

  call filcc(adv,adv_l1,adv_l2,adv_l3,adv_h1,adv_h2,adv_h3, &
       domlo,domhi,delta,xlo,bc)
  
  !     XLO
  if ( bc(1,1,1).eq.EXT_DIR .and. adv_l1.lt.domlo(1)) then
     print *,'tempfill: SHOULD NEVER GET HERE bc(1,1,1) .eq. EXT_DIR) '
     stop
  end if
  
  !     XHI
  if ( bc(1,2,1).eq.EXT_DIR .and. adv_h1.gt.domhi(1)) then
     print *,'tempfill: SHOULD NEVER GET HERE bc(1,2,1) .eq. EXT_DIR) '
     stop
  end if
  
  !     YLO
  if ( bc(2,1,1).eq.EXT_DIR .and. adv_l2.lt.domlo(2)) then
     print *,'tempfill: SHOULD NEVER GET HERE bc(2,1,1) .eq. EXT_DIR) '
     stop
  end if
  
  !     YHI
  if ( bc(2,2,1).eq.EXT_DIR .and. adv_h2.gt.domhi(2)) then
     print *,'tempfill: SHOULD NEVER GET HERE bc(2,2,1) .eq. EXT_DIR) '
     stop
  end if

  !     ZLO
  if ( bc(3,1,1).eq.EXT_DIR .and. adv_l3.lt.domlo(3)) then
     print *,'tempfill: SHOULD NEVER GET HERE bc(3,1,1) .eq. EXT_DIR) '
     stop
  end if
  
  !     ZHI
  if ( bc(3,2,1).eq.EXT_DIR .and. adv_h3.gt.domhi(3)) then
     print *,'tempfill: SHOULD NEVER GET HERE bc(3,2,1) .eq. EXT_DIR) '
     stop
  end if
  
end subroutine rns_tempfill
