! Fill the entire state
subroutine rns_grpfill(adv,adv_l1,adv_h1, &
     domlo,domhi,delta,xlo,time,bc_in)
  
  use meth_params_module, only : NVAR
  use sdc_boundary_module, only : isFEval
  use inflow_module
  implicit none
  include 'bc_types.fi'
  
  integer          :: adv_l1,adv_h1
  integer          :: bc_in(1,2,*)
  integer          :: domlo(1), domhi(1)
  double precision :: delta(1), xlo(1), time
  double precision :: adv(adv_l1:adv_h1,NVAR)
  
  integer n, i
  integer bc(1,2,NVAR)

  bc = bc_in(:,:,1:NVAR)
  if (isFEval) then
     where (bc .eq. EXT_DIR) bc = FOEXTRAP
  end if

  do n = 1,NVAR
     call filcc(adv(adv_l1,n),adv_l1,adv_h1, &
          domlo,domhi,delta,xlo,bc(1,1,n))
  enddo

  if (isFEval) return

  !        XLO
  if (adv_l1.lt.domlo(1)) then
     if ( bc(1,1,1).eq.EXT_DIR) then
        do n = 1, NVAR
           do i = adv_l1, domlo(1)-1 
              adv(i,n) = inflow_state(n)
           end do
        end do
     else
!zzz        print *,'grpfill: SHOULD NEVER GET HERE bc(1,1,n) .eq. EXT_DIR) '
!zzz        stop
     end if
  end if
     
  !        XHI
  if ( bc(1,2,1).eq.EXT_DIR .and. adv_h1.gt.domhi(1)) then
     print *,'grpfill: SHOULD NEVER GET HERE bc(1,2,n) .eq. EXT_DIR) '
     stop
  end if
  
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
  
end subroutine rns_regfill

! Fill density
subroutine rns_denfill(adv,adv_l1,adv_h1,domlo,domhi,delta,xlo,time,bc)
  
  use meth_params_module, only : URHO
  use inflow_module
  implicit none
  include 'bc_types.fi'
  
  integer          :: adv_l1,adv_h1
  integer          :: bc(1,2,*)
  integer          :: domlo(1), domhi(1)
  double precision :: delta(1), xlo(1), time
  double precision :: adv(adv_l1:adv_h1)

  call filcc(adv,adv_l1,adv_h1,domlo,domhi,delta,xlo,bc)
  
  !     XLO
  if (adv_l1.lt.domlo(1)) then
     if ( bc(1,1,1).eq.EXT_DIR) then
        adv(adv_l1:domlo(1)-1) = inflow_state(URHO)
     else
!zzz        print *,'denfill: SHOULD NEVER GET HERE bc(1,1,1) .eq. EXT_DIR) '
!zzz        stop
     end if
  end if
  
  !     XHI
  if ( bc(1,2,1).eq.EXT_DIR .and. adv_h1.gt.domhi(1)) then
     print *,'denfill: SHOULD NEVER GET HERE bc(1,2,1) .eq. EXT_DIR) '
     stop
  end if
  
end subroutine rns_denfill

! Fill x-momentum
subroutine rns_mxfill(adv,adv_l1,adv_h1,domlo,domhi,delta,xlo,time,bc)
  
  use meth_params_module, only : UMX
  use inflow_module
  implicit none
  include 'bc_types.fi'
  
  integer          :: adv_l1,adv_h1
  integer          :: bc(1,2,*)
  integer          :: domlo(1), domhi(1)
  double precision :: delta(1), xlo(1), time
  double precision :: adv(adv_l1:adv_h1)
  
  call filcc(adv,adv_l1,adv_h1,domlo,domhi,delta,xlo,bc)
  
  !     XLO
  if (adv_l1.lt.domlo(1)) then
     if ( bc(1,1,1).eq.EXT_DIR) then
        adv(adv_l1:domlo(1)-1) = inflow_state(UMX)
     else
!zzz        print *,'mxfill: SHOULD NEVER GET HERE bc(1,1,1) .eq. EXT_DIR) '
!zzz        stop
     end if
  end if
  
  !     XHI
  if ( bc(1,2,1).eq.EXT_DIR .and. adv_h1.gt.domhi(1)) then
     print *,'mxfill: SHOULD NEVER GET HERE bc(1,2,1) .eq. EXT_DIR) '
     stop
  end if
  
end subroutine rns_mxfill

! Fill temperature
subroutine rns_tempfill(adv,adv_l1,adv_h1,domlo,domhi,delta,xlo,time,bc)
  
  use meth_params_module, only : UTEMP
  use inflow_module
  implicit none
  include 'bc_types.fi'
  
  integer          :: adv_l1,adv_h1
  integer          :: bc(1,2,*)
  integer          :: domlo(1), domhi(1)
  double precision :: delta(1), xlo(1), time
  double precision :: adv(adv_l1:adv_h1)
  
  call filcc(adv,adv_l1,adv_h1,domlo,domhi,delta,xlo,bc)
  
  !     XLO
  if (adv_l1.lt.domlo(1)) then
     if ( bc(1,1,1).eq.EXT_DIR) then
        adv(adv_l1:domlo(1)-1) = inflow_state(UTEMP)
     else
!zzz        print *,'tempfill: SHOULD NEVER GET HERE bc(1,1,1) .eq. EXT_DIR) '
!zzz        stop
     endif
  end if
  
  !     XHI
  if ( bc(1,2,1).eq.EXT_DIR .and. adv_h1.gt.domhi(1)) then
     print *,'tempfill: SHOULD NEVER GET HERE bc(1,2,1) .eq. EXT_DIR) '
     stop
  end if
  
end subroutine rns_tempfill
