! Fill the entire state
subroutine rns_grpfill(adv,adv_l1,adv_h1, &
     domlo,domhi,delta,xlo,time,bc)
  
  use eos_module, only : gamma_const, eos_get_T
  use probdata_module
  use meth_params_module, only : NVAR, URHO, UMX, UEDEN, UTEMP, UFS, NSPEC
  implicit none
  include 'bc_types.fi'
  
  integer          :: adv_l1,adv_h1
  integer          :: bc(1,2,*)
  integer          :: domlo(1), domhi(1)
  double precision :: delta(1), xlo(1), time
  double precision :: adv(adv_l1:adv_h1,NVAR)

  double precision, parameter :: gp(2) = (/ -1.d0/sqrt(3.d0), 1.d0/sqrt(3.d0) /)
  double precision, parameter :: Pi = 3.1415926535897932d0
  
  integer :: i, ii, n
  double precision :: xcen, xg
  double precision :: rhot, Pt, et, Tt, Yt(2), ekt
  
  do n = 1,NVAR
     call filcc(adv(adv_l1,n),adv_l1,adv_h1, &
          domlo,domhi,delta,xlo,bc(1,1,n))
  enddo
  
!  do n = 1, NVAR
     
     !        XLO
     if ( bc(1,1,1).eq.EXT_DIR .and. adv_l1.lt.domlo(1)) then
!        print *,'grpfill: SHOULD NEVER GET HERE bc(1,1,n) .eq. EXT_DIR) '
!        stop

        do i=adv_l1,domlo(1)-1

           xcen = xlo(1) + delta(1)*(dble(i-adv_l1) + 0.5d0) - u0*time

           adv(i,:) = 0.d0
           
           do ii = 1, 2
              xg = xcen + 0.5d0*delta(1)*gp(ii)
              
              rhot = rho0 + drho0*sin(2.d0*Pi/length(1)*xg)
              
              Pt = p0
              
              et = Pt/((gamma_const-1.d0)*rhot)
              
              ekt = 0.5d0*u0*u0
        
              Yt(1) = 1.d0
              Yt(2) = 0.d0
              
              call eos_get_T(Tt, et, Yt)
        
              adv(i,URHO ) = adv(i,URHO ) + 0.5d0*rhot
              adv(i,UMX  ) = adv(i,UMX  ) + 0.5d0*rhot*u0
              adv(i,UEDEN) = adv(i,UEDEN) + 0.5d0*rhot*(et+ekt)
              adv(i,UTEMP) = adv(i,UTEMP) + 0.5d0*Tt
              do n=1, NSPEC
                 adv(i,UFS+n-1) = adv(i,UFS+n-1) + 0.5d0*rhot*Yt(n)
              end do
              
           end do
           
        end do

     end if
     
     !        XHI
     if ( bc(1,2,1).eq.EXT_DIR .and. adv_h1.gt.domhi(1)) then
!        print *,'grpfill: SHOULD NEVER GET HERE bc(1,2,n) .eq. EXT_DIR) '
!        stop

        do i=domhi(1)+1,adv_h1

           xcen = xlo(1) + delta(1)*(dble(i-adv_l1) + 0.5d0) - u0*time

           adv(i,:) = 0.d0
           
           do ii = 1, 2
              xg = xcen + 0.5d0*delta(1)*gp(ii)
              
              rhot = rho0 + drho0*sin(2.d0*Pi/length(1)*xg)
              
              Pt = p0
              
              et = Pt/((gamma_const-1.d0)*rhot)
              
              ekt = 0.5d0*u0*u0
        
              Yt(1) = 1.d0
              Yt(2) = 0.d0
              
              call eos_get_T(Tt, et, Yt)
        
              adv(i,URHO ) = adv(i,URHO ) + 0.5d0*rhot
              adv(i,UMX  ) = adv(i,UMX  ) + 0.5d0*rhot*u0
              adv(i,UEDEN) = adv(i,UEDEN) + 0.5d0*rhot*(et+ekt)
              adv(i,UTEMP) = adv(i,UTEMP) + 0.5d0*Tt
              do n=1, NSPEC
                 adv(i,UFS+n-1) = adv(i,UFS+n-1) + 0.5d0*rhot*Yt(n)
              end do
              
           end do
           
        end do

     end if
     
!  end do
  
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
  
end subroutine rns_tempfill
