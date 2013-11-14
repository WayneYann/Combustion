! Fill the entire state
subroutine rns_grpfill(adv,adv_l1,adv_h1, &
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
        print *,'grpfill: SHOULD NEVER GET HERE bc(1,1,n) .eq. EXT_DIR) '
        stop
     end if
     
     !        XHI
     if ( bc(1,2,n).eq.EXT_DIR .and. adv_h1.gt.domhi(1)) then
        print *,'grpfill: SHOULD NEVER GET HERE bc(1,2,n) .eq. EXT_DIR) '
        stop
     end if
     
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
     print *,'denfill: SHOULD NEVER GET HERE bc(1,1,1) .eq. EXT_DIR) '
     stop
  end if
  
  !     XHI
  if ( bc(1,2,1).eq.EXT_DIR .and. adv_h1.gt.domhi(1)) then
     print *,'denfill: SHOULD NEVER GET HERE bc(1,2,1) .eq. EXT_DIR) '
     stop
  end if
  
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
  
end subroutine rns_mxfill
