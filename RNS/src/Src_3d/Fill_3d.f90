! Fill the entire state
subroutine rns_grpfill(adv,adv_l1,adv_l2,adv_l3,adv_h1,adv_h2,adv_h3, &
     domlo,domhi,delta,xlo,time,bc)
 
  use meth_params_module, only : NVAR
  
  implicit none
  include 'bc_types.fi'
  integer adv_l1,adv_l2,adv_l3,adv_h1,adv_h2,adv_h3
  integer bc(3,2,*)
  integer domlo(3), domhi(3)
  double precision delta(3), xlo(3), time
  double precision adv(adv_l1:adv_h1,adv_l2:adv_h2,adv_l3:adv_h3,NVAR)
  
  integer n

  !$omp parallel do private(n)
  do n = 1,NVAR
     call filcc(adv(adv_l1,adv_l2,adv_l3,n), &
          adv_l1,adv_l2,adv_l3,adv_h1,adv_h2,adv_h3, &
          domlo,domhi,delta,xlo,bc(1,1,n))
  enddo
  !$omp end parallel do

  do n = 1, NVAR

     !     XLO
     if ( bc(1,1,n).eq.EXT_DIR .and. adv_l1.lt.domlo(1)) then
        print *,'grpfill: SHOULD NEVER GET HERE bc(1,1,n) .eq. EXT_DIR) '
        stop
     end if
     
     !     XHI
     if ( bc(1,2,n).eq.EXT_DIR .and. adv_h1.gt.domhi(1)) then
        print *,'grpfill: SHOULD NEVER GET HERE bc(1,2,n) .eq. EXT_DIR) '
        stop
     end if
     
     !     YLO
     if ( bc(2,1,n).eq.EXT_DIR .and. adv_l2.lt.domlo(2)) then
        print *,'grpfill: SHOULD NEVER GET HERE bc(2,1,n) .eq. EXT_DIR) '
        stop
     end if
     
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
