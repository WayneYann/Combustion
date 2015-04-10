! Fill the entire state
subroutine rns_grpfill(adv,adv_l1,adv_h1, &
     domlo,domhi,delta,xlo,time,bc_in)
  
  use eos_module, only : gamma_const, eos_get_T
  use probdata_module
  use meth_params_module, only : NVAR, URHO, UMX, UEDEN, UTEMP, UFS, NSPEC
  use sdc_boundary_module, only : isFEval
  implicit none
  include 'bc_types.fi'
  
  integer          :: adv_l1,adv_h1
  integer          :: bc_in(1,2,*)
  integer          :: domlo(1), domhi(1)
  double precision :: delta(1), xlo(1), time
  double precision :: adv(adv_l1:adv_h1,NVAR)
  
  integer n
  integer :: bc(1,2,NVAR)

  double precision xg, xcen, rho, u, p, Yt(2), Tt, et
  integer i, ii

  double precision, parameter :: gp(2) = (/ -1.d0/sqrt(3.d0), 1.d0/sqrt(3.d0) /)

  bc = bc_in(:,:,1:NVAR)
  if (isFEval) then
     do n=1,NVAR
        if (bc(1,1,n) .eq. EXT_DIR) bc(1,1,n) = FOEXTRAP
        if (bc(1,2,n) .eq. EXT_DIR) bc(1,2,n) = FOEXTRAP
     end do
  end if

  do n = 1,NVAR
     call filcc(adv(adv_l1,n),adv_l1,adv_h1, &
          domlo,domhi,delta,xlo,bc(1,1,n))
  enddo
  
  if (isFEval) return

  if (adv_l1 .lt. domlo(1)) then
     if (bc(1,1,1) .eq. EXT_DIR) then
        do i = adv_l1, domlo(1)-1
           
           adv(i,:) = 0.d0

           rho = 3.857143d0
           u   = 2.629369d0
           p   = 10.33333d0

           et = p/((gamma_const-1.d0)*rho)
           Yt = 0.5d0
           call eos_get_T(Tt, et, Yt)
           
           adv(i,URHO ) = adv(i,URHO ) + rho
           adv(i,UMX  ) = adv(i,UMX  ) + rho*u
           adv(i,UEDEN) = adv(i,UEDEN) + rho*(et+0.5d0*u*u)
           adv(i,UTEMP) = adv(i,UTEMP) + Tt
           do n=1, NSPEC
              adv(i,UFS+n-1) = adv(i,UFS+n-1) + rho*Yt(n)
           end do
        end do        
     end if
  end if

  if (adv_h1 .gt. domhi(1)) then
     if (bc(1,2,1) .eq. EXT_DIR) then
        
        do i = domhi(1)+1, adv_h1
           xcen = xlo(1) + delta(1)*(DBLE(i-adv_l1)+.5d0)
           
           adv(i,:) = 0.d0
           
           do ii = 1, 2
              xg = xcen + 0.5d0*delta(1)*gp(ii)
              
              rho = 1.d0+0.2d0*sin(5.d0*xg)
              u = 0.d0
              p = 1.d0
              
              et = p/((gamma_const-1.d0)*rho)
              Yt = 0.5d0
              call eos_get_T(Tt, et, Yt)
              
              adv(i,URHO ) = adv(i,URHO ) + 0.5d0*rho
              adv(i,UMX  ) = adv(i,UMX  ) + 0.5d0*rho*u
              adv(i,UEDEN) = adv(i,UEDEN) + 0.5d0*rho*(et+0.5d0*u*u)
              adv(i,UTEMP) = adv(i,UTEMP) + 0.5d0*Tt
              do n=1, NSPEC
                 adv(i,UFS+n-1) = adv(i,UFS+n-1) + 0.5d0*rho*Yt(n)
              end do
           end do
        end do
        
     end if
  end if

  do n = 1, NVAR
     
     !        XLO
     ! if ( bc(1,1,n).eq.EXT_DIR .and. adv_l1.lt.domlo(1)) then
     !    print *,'grpfill: SHOULD NEVER GET HERE bc(1,1,n) .eq. EXT_DIR) '
     !    stop
     ! end if
     
     !        XHI
     ! if ( bc(1,2,n).eq.EXT_DIR .and. adv_h1.gt.domhi(1)) then
     !    print *,'grpfill: SHOULD NEVER GET HERE bc(1,2,n) .eq. EXT_DIR) '
     !    stop
     ! end if
     
  end do
  
end subroutine rns_grpfill

! Fill one variable
subroutine rns_regfill(adv,adv_l1,adv_h1,domlo,domhi,delta,xlo,time,bc)
  
  implicit none
  include 'bc_types.fi'
  
  integer          :: adv_l1,adv_h1
  integer          :: bc(1,2,*)
  integer          :: domlo(1), domhi(1)
  double precision :: delta(1), xlo(1), time
  double precision :: adv(adv_l1:adv_h1)
  
  call filcc(adv,adv_l1,adv_h1,domlo,domhi,delta,xlo,bc)
  
  ! !     XLO
  ! if ( bc(1,1,1).eq.EXT_DIR .and. adv_l1.lt.domlo(1)) then
  !    print *,'regfill: SHOULD NEVER GET HERE bc(1,1,1) .eq. EXT_DIR) '
  !    stop
  ! end if
  
  !     XHI
  ! if ( bc(1,2,1).eq.EXT_DIR .and. adv_h1.gt.domhi(1)) then
  !    print *,'regfill: SHOULD NEVER GET HERE bc(1,2,1) .eq. EXT_DIR) '
  !    stop
  ! end if
  
end subroutine rns_regfill

! Fill density
subroutine rns_denfill(adv,adv_l1,adv_h1,domlo,domhi,delta,xlo,time,bc)
  
  implicit none
  include 'bc_types.fi'
  
  integer          :: adv_l1,adv_h1
  integer          :: bc(1,2,*)
  integer          :: domlo(1), domhi(1)
  double precision :: delta(1), xlo(1), time
  double precision :: adv(adv_l1:adv_h1)
  
  call filcc(adv,adv_l1,adv_h1,domlo,domhi,delta,xlo,bc)
  
  ! !     XLO
  ! if ( bc(1,1,1).eq.EXT_DIR .and. adv_l1.lt.domlo(1)) then
  !    print *,'denfill: SHOULD NEVER GET HERE bc(1,1,1) .eq. EXT_DIR) '
  !    stop
  ! end if
  
  !     XHI
  ! if ( bc(1,2,1).eq.EXT_DIR .and. adv_h1.gt.domhi(1)) then
  !    print *,'denfill: SHOULD NEVER GET HERE bc(1,2,1) .eq. EXT_DIR) '
  !    stop
  ! end if
  
end subroutine rns_denfill

! Fill x-momentum
subroutine rns_mxfill(adv,adv_l1,adv_h1,domlo,domhi,delta,xlo,time,bc)
  
  implicit none
  include 'bc_types.fi'
  
  integer          :: adv_l1,adv_h1
  integer          :: bc(1,2,*)
  integer          :: domlo(1), domhi(1)
  double precision :: delta(1), xlo(1), time
  double precision :: adv(adv_l1:adv_h1)
  
  call filcc(adv,adv_l1,adv_h1,domlo,domhi,delta,xlo,bc)
  
  ! !     XLO
  ! if ( bc(1,1,1).eq.EXT_DIR .and. adv_l1.lt.domlo(1)) then
  !    print *,'mxfill: SHOULD NEVER GET HERE bc(1,1,1) .eq. EXT_DIR) '
  !    stop
  ! end if
  
  !     XHI
  ! if ( bc(1,2,1).eq.EXT_DIR .and. adv_h1.gt.domhi(1)) then
  !    print *,'mxfill: SHOULD NEVER GET HERE bc(1,2,1) .eq. EXT_DIR) '
  !    stop
  ! end if
  
end subroutine rns_mxfill

! Fill temperature
subroutine rns_tempfill(adv,adv_l1,adv_h1,domlo,domhi,delta,xlo,time,bc)
  
  implicit none
  include 'bc_types.fi'
  
  integer          :: adv_l1,adv_h1
  integer          :: bc(1,2,*)
  integer          :: domlo(1), domhi(1)
  double precision :: delta(1), xlo(1), time
  double precision :: adv(adv_l1:adv_h1)
  
  call filcc(adv,adv_l1,adv_h1,domlo,domhi,delta,xlo,bc)
  
  ! !     XLO
  ! if ( bc(1,1,1).eq.EXT_DIR .and. adv_l1.lt.domlo(1)) then
  !    print *,'tempfill: SHOULD NEVER GET HERE bc(1,1,1) .eq. EXT_DIR) '
  !    stop
  ! end if
  
  !     XHI
  ! if ( bc(1,2,1).eq.EXT_DIR .and. adv_h1.gt.domhi(1)) then
  !    print *,'tempfill: SHOULD NEVER GET HERE bc(1,2,1) .eq. EXT_DIR) '
  !    stop
  ! end if
  
end subroutine rns_tempfill
