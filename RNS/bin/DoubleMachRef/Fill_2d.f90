! Fill the entire state
subroutine rns_grpfill(adv,adv_l1,adv_l2,adv_h1,adv_h2, &
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
     do n=1,NVAR
        if (bc(1,1,n).eq.EXT_DIR) then
           do j = adv_l2,adv_h2  ! fill the corners too
              do i = adv_l1, domlo(1)-1
                 adv(i,j,n) = state1(n)
              end do
           end do
        else
           print *,'grpfill: SHOULD NEVER GET HERE bc(1,1,n) .ne. EXT_DIR) '
           stop
        end if
     end do
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
  
end subroutine rns_grpfill

! Fill one variable
subroutine rns_regfill(adv,adv_l1,adv_l2,adv_h1,adv_h2, &
     domlo,domhi,delta,xlo,time,bc)
  
  implicit none
  include 'bc_types.fi'
  integer adv_l1,adv_l2,adv_h1,adv_h2
  integer bc(2,2,*)
  integer domlo(2), domhi(2)
  double precision delta(2), xlo(2), time
  double precision adv(adv_l1:adv_h1,adv_l2:adv_h2)
  
  call filcc(adv,adv_l1,adv_l2,adv_h1,adv_h2,domlo,domhi,delta,xlo,bc)
  
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
  
end subroutine rns_regfill


! Fill density
subroutine rns_denfill(adv,adv_l1,adv_l2,adv_h1,adv_h2, &
     domlo,domhi,delta,xlo,time,bc)
  
  implicit none
  include 'bc_types.fi'
  integer adv_l1,adv_l2,adv_h1,adv_h2
  integer bc(2,2,*)
  integer domlo(2), domhi(2)
  double precision delta(2), xlo(2), time
  double precision adv(adv_l1:adv_h1,adv_l2:adv_h2)

  
  call filcc(adv,adv_l1,adv_l2,adv_h1,adv_h2,domlo,domhi,delta,xlo,bc)
  
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
  
end subroutine rns_denfill


! Fill x-momentum
subroutine rns_mxfill(adv,adv_l1,adv_l2,adv_h1,adv_h2, &
     domlo,domhi,delta,xlo,time,bc)
  
  implicit none
  include 'bc_types.fi'
  integer adv_l1,adv_l2,adv_h1,adv_h2
  integer bc(2,2,*)
  integer domlo(2), domhi(2)
  double precision delta(2), xlo(2), time
  double precision adv(adv_l1:adv_h1,adv_l2:adv_h2)

  
  call filcc(adv,adv_l1,adv_l2,adv_h1,adv_h2,domlo,domhi,delta,xlo,bc)
  
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
  
end subroutine rns_mxfill


! Fill y-momentum
subroutine rns_myfill(adv,adv_l1,adv_l2,adv_h1,adv_h2, &
     domlo,domhi,delta,xlo,time,bc)
  
  implicit none
  include 'bc_types.fi'
  integer adv_l1,adv_l2,adv_h1,adv_h2
  integer bc(2,2,*)
  integer domlo(2), domhi(2)
  double precision delta(2), xlo(2), time
  double precision adv(adv_l1:adv_h1,adv_l2:adv_h2)

  
  call filcc(adv,adv_l1,adv_l2,adv_h1,adv_h2,domlo,domhi,delta,xlo,bc)
  
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
  
end subroutine rns_myfill
