! Fill the entire state
subroutine rns_grpfill(adv,adv_l1,adv_l2,adv_l3,adv_h1,adv_h2,adv_h3, &
     domlo,domhi,delta,xlo,time,bc_in)
 
  use meth_params_module
  use probdata_module
  use sdc_boundary_module, only : isFEval
  use turbinflow_module
  
  implicit none
  include 'bc_types.fi'
  integer adv_l1,adv_l2,adv_l3,adv_h1,adv_h2,adv_h3
  integer bc_in(3,2,*)
  integer domlo(3), domhi(3)
  double precision delta(3), xlo(3), time
  double precision adv(adv_l1:adv_h1,adv_l2:adv_h2,adv_l3:adv_h3,NVAR)
  
  integer i,j,k,n,ii,jj,iwrk
  double precision x,y,xg,yg,r,eta,Yt(nspec),Tt,u1t,u2t,u3t,rhot,et,rwrk
  double precision, parameter :: gp(2) = (/ -1.d0/sqrt(3.d0), 1.d0/sqrt(3.d0) /)
  double precision, parameter :: wgt = 0.25d0
  integer bc(3,2,NVAR)
  double precision, allocatable :: vturb(:,:,:), xturb(:), yturb(:), etaturb(:,:)
  double precision :: zturb, Ek, Eknew, rhoeta

  if (.not. jet_initialized) call init_jet()

  if (.not. turbinflow_initialized) call init_turbinflow(turbfile)

  bc = bc_in(:,:,1:NVAR)
  if (isFEval) bc(3,1,:) = FOEXTRAP

  !$omp parallel do private(n)
  do n = 1,NVAR
     call filcc(adv(adv_l1,adv_l2,adv_l3,n), &
          adv_l1,adv_l2,adv_l3,adv_h1,adv_h2,adv_h3, &
          domlo,domhi,delta,xlo,bc(1,1,n))
  enddo
  !$omp end parallel do

  if (isFEval) return

  ! zlo
  if (adv_l3.lt.domlo(3)) then
     if (bc(3,1,1) .eq. EXT_DIR) then

        !$omp parallel do private(i,j,k,n,ii,jj,iwrk,x,y,xg,yg,r,eta) &
        !$omp private(Yt,Tt,u1t,u2t,u3t,rhot,et,rwrk) &
        !$omp collapse(2)
        do k = adv_l3, domlo(3)-1
        do j = adv_l2, adv_h2
        do i = adv_l1, adv_h1

           x = (DBLE(i-adv_l1)+.5d0)*delta(1)+xlo(1)
           y = (DBLE(j-adv_l2)+.5d0)*delta(2)+xlo(2)

           adv(i,j,k,:) = 0.d0

           do jj=1,2
              yg = y + 0.5d0*delta(2)*gp(jj)
              do ii=1,2
                 xg = x + 0.5d0*delta(1)*gp(ii)
                    
                 r = sqrt(xg*xg+yg*yg)
                 
                 eta = 0.5d0*(1.d0-tanh((r -splitx)/xfrontw))

                 do n=1,nspec
                    Yt(n) = eta*fuel_Y(n) + (1.d0-eta)*air_Y(n)
                 end do
                 Tt  = eta * T_in + (1.d0-eta) * T_co
                 u1t = 0.d0
                 u2t = 0.d0
                 u3t = eta * vn_in + (1.d0-eta) * vn_co                    
                    

                 CALL CKRHOY(pamb,Tt,Yt,IWRK,RWRK,rhot)
                 call CKUBMS(Tt,Yt,IWRK,RWRK,et)

                 adv(i,j,k,URHO ) = adv(i,j,k,URHO ) + wgt*rhot
                 adv(i,j,k,UMX  ) = adv(i,j,k,UMX  ) + wgt*rhot*u1t
                 adv(i,j,k,UMY  ) = adv(i,j,k,UMY  ) + wgt*rhot*u2t
                 adv(i,j,k,UMZ  ) = adv(i,j,k,UMZ  ) + wgt*rhot*u3t
                 adv(i,j,k,UEDEN) = adv(i,j,k,UEDEN) + wgt*rhot*(et + 0.5d0*(u1t**2+u2t**2+u3t**2))
                 adv(i,j,k,UTEMP) = adv(i,j,k,UTEMP) + wgt*Tt
                 do n=1, NSPEC
                    adv(i,j,k,UFS+n-1) = adv(i,j,k,UFS+n-1) + wgt*rhot*Yt(n)
                 end do

              end do
           end do

        end do
        end do
        end do
        !$omp end parallel do

        ! add turbuence
        allocate(vturb  (adv_l1:adv_h1,adv_l2:adv_h2,3))
        allocate(xturb  (adv_l1:adv_h1))
        allocate(yturb  (adv_l2:adv_h2))
        allocate(etaturb(adv_l1:adv_h1,adv_l2:adv_h2))

        do i = adv_l1, adv_h1
           xturb(i) = (DBLE(i-adv_l1)+.5d0)*delta(1)+xlo(1)
        end do

        do j = adv_l2, adv_h2
           yturb(j) = (DBLE(j-adv_l2)+.5d0)*delta(2)+xlo(2)
        end do

        do j = adv_l2, adv_h2
           do i = adv_l1, adv_h1
              r = sqrt(xturb(i)**2+yturb(j)**2)
              etaturb(i,j) = 0.5d0*(1.d0-tanh((r-splitx)/xfrontw))
           end do
        end do

        do k = adv_l3, domlo(3)-1
           zturb = (DBLE(k-adv_l3)+.5d0)*delta(3)+xlo(3) + vn_in*time

           call get_turbvelocity(adv_l1,adv_l2,adv_h1,adv_h2,xturb,yturb,zturb,vturb)

           do j = adv_l2, adv_h2
              do i = adv_l1, adv_h1
                 Ek = adv(i,j,k,UMX)**2+adv(i,j,k,UMY)**2+adv(i,j,k,UMZ)**2
                 rhoeta = adv(i,j,k,URHO)*etaturb(i,j)*turb_boost_factor
                 adv(i,j,k,UMX) = adv(i,j,k,UMX) + vturb(i,j,1)*rhoeta
                 adv(i,j,k,UMY) = adv(i,j,k,UMY) + vturb(i,j,2)*rhoeta
                 adv(i,j,k,UMZ) = adv(i,j,k,UMZ) + vturb(i,j,3)*rhoeta
                 Eknew = adv(i,j,k,UMX)**2+adv(i,j,k,UMY)**2+adv(i,j,k,UMZ)**2
                 adv(i,j,k,UEDEN) = adv(i,j,k,UEDEN) + &
                      0.5d0*(Eknew-Ek)/adv(i,j,k,URHO)
              end do
           end do
        end do

        deallocate(vturb,xturb,yturb,etaturb)

     else
        print *,'grpfill: SHOULD NEVER GET HERE bc(3,1,1) .ne. EXT_DIR) '
        stop
     end if
  end if

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
     ! if ( bc(2,1,n).eq.EXT_DIR .and. adv_l2.lt.domlo(2)) then
     !    print *,'grpfill: SHOULD NEVER GET HERE bc(2,1,n) .eq. EXT_DIR) '
     !    stop
     ! end if
     
     !     YHI
     if ( bc(2,2,n).eq.EXT_DIR .and. adv_h2.gt.domhi(2)) then
        print *,'grpfill: SHOULD NEVER GET HERE bc(2,2,n) .eq. EXT_DIR) '
        stop
     end if

     !     ZLO
     ! if ( bc(3,1,n).eq.EXT_DIR .and. adv_l3.lt.domlo(3)) then
     !    print *,'grpfill: SHOULD NEVER GET HERE bc(3,1,n) .eq. EXT_DIR) '
     !    stop
     ! end if
     
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
  
  use meth_params_module
  use probdata_module

  implicit none
  include 'bc_types.fi'
  integer adv_l1,adv_l2,adv_l3,adv_h1,adv_h2,adv_h3
  integer bc(3,2,*)
  integer domlo(3), domhi(3)
  double precision delta(3), xlo(3), time
  double precision adv(adv_l1:adv_h1,adv_l2:adv_h2,adv_l3:adv_h3)

  integer i,j,k,n,ii,jj,iwrk
  double precision x,y,xg,yg,r,eta,Yt(nspec),Tt,rhot,rwrk
  double precision, parameter :: gp(2) = (/ -1.d0/sqrt(3.d0), 1.d0/sqrt(3.d0) /)
  double precision, parameter :: wgt = 0.25d0

  if (.not. jet_initialized) then
     call init_jet()
  end if

  call filcc(adv,adv_l1,adv_l2,adv_l3,adv_h1,adv_h2,adv_h3, &
       domlo,domhi,delta,xlo,bc)

  ! zlo
  if (adv_l3.lt.domlo(3)) then
     if (bc(3,1,1) .eq. EXT_DIR) then

        !$omp parallel do private(i,j,k,n,ii,jj,iwrk,x,y,xg,yg,r,eta) &
        !$omp private(Yt,Tt,rhot,rwrk) &
        !$omp collapse(2)
        do k = adv_l3, domlo(3)-1
        do j = adv_l2, adv_h2
        do i = adv_l1, adv_h1

           x = (DBLE(i-adv_l1)+.5d0)*delta(1)+xlo(1)
           y = (DBLE(j-adv_l2)+.5d0)*delta(2)+xlo(2)

           adv(i,j,k) = 0.d0

           do jj=1,2
              yg = y + 0.5d0*delta(2)*gp(jj)
              do ii=1,2
                 xg = x + 0.5d0*delta(1)*gp(ii)
                    
                 r = sqrt(xg*xg+yg*yg)

                 eta = 0.5d0*(1.d0-tanh((r -splitx)/xfrontw))

                 do n=1,nspec
                    Yt(n) = eta*fuel_Y(n) + (1.d0-eta)*air_Y(n)
                 end do
                 Tt  = eta * T_in + (1.d0-eta) * T_co

                 CALL CKRHOY(pamb,Tt,Yt,IWRK,RWRK,rhot)
                 adv(i,j,k) = adv(i,j,k) + wgt*rhot

              end do
           end do

        end do
        end do
        end do
        !$omp end parallel do

     else
        print *,'denfill: SHOULD NEVER GET HERE bc(3,1,1) .ne. EXT_DIR) '
        stop
     end if
  end if
  
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

  ! !     ZLO
  ! if ( bc(3,1,1).eq.EXT_DIR .and. adv_l3.lt.domlo(3)) then
  !    print *,'denfill: SHOULD NEVER GET HERE bc(3,1,1) .eq. EXT_DIR) '
  !    stop
  ! end if
  
  !     ZHI
  if ( bc(3,2,1).eq.EXT_DIR .and. adv_h3.gt.domhi(3)) then
     print *,'denfill: SHOULD NEVER GET HERE bc(3,2,1) .eq. EXT_DIR) '
     stop
  end if
  
end subroutine rns_denfill


! Fill x-momentum
subroutine rns_mxfill(adv,adv_l1,adv_l2,adv_l3,adv_h1,adv_h2,adv_h3, &
     domlo,domhi,delta,xlo,time,bc)

  use meth_params_module
  use probdata_module
  
  implicit none
  include 'bc_types.fi'
  integer adv_l1,adv_l2,adv_l3,adv_h1,adv_h2,adv_h3
  integer bc(3,2,*)
  integer domlo(3), domhi(3)
  double precision delta(3), xlo(3), time
  double precision adv(adv_l1:adv_h1,adv_l2:adv_h2,adv_l3:adv_h3)

  integer :: i,j,k
  
  call filcc(adv,adv_l1,adv_l2,adv_l3,adv_h1,adv_h2,adv_h3, &
       domlo,domhi,delta,xlo,bc)

  ! zlo
  if (adv_l3.lt.domlo(3)) then
     if (bc(3,1,1) .eq. EXT_DIR) then
        do k = adv_l3, domlo(3)-1
           do j = adv_l2, adv_h2
              do i = adv_l1, adv_h1
                 adv(i,j,k) = 0.d0
              end do
           end do
        end do
     end if
  end if
  
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
  ! if ( bc(3,1,1).eq.EXT_DIR .and. adv_l3.lt.domlo(3)) then
  !    print *,'mxfill: SHOULD NEVER GET HERE bc(3,1,1) .eq. EXT_DIR) '
  !    stop
  ! end if
  
  !     ZHI
  if ( bc(3,2,1).eq.EXT_DIR .and. adv_h3.gt.domhi(3)) then
     print *,'mxfill: SHOULD NEVER GET HERE bc(3,2,1) .eq. EXT_DIR) '
     stop
  end if
  
end subroutine rns_mxfill


! Fill y-momentum
subroutine rns_myfill(adv,adv_l1,adv_l2,adv_l3,adv_h1,adv_h2,adv_h3, &
     domlo,domhi,delta,xlo,time,bc)

  use meth_params_module
  use probdata_module
  
  implicit none
  include 'bc_types.fi'
  integer adv_l1,adv_l2,adv_l3,adv_h1,adv_h2,adv_h3
  integer bc(3,2,*)
  integer domlo(3), domhi(3)
  double precision delta(3), xlo(3), time
  double precision adv(adv_l1:adv_h1,adv_l2:adv_h2,adv_l3:adv_h3)

  integer :: i,j,k

  call filcc(adv,adv_l1,adv_l2,adv_l3,adv_h1,adv_h2,adv_h3, &
       domlo,domhi,delta,xlo,bc)

  ! zlo
  if (adv_l3.lt.domlo(3)) then
     if (bc(3,1,1) .eq. EXT_DIR) then
        do k = adv_l3, domlo(3)-1
           do j = adv_l2, adv_h2
              do i = adv_l1, adv_h1
                 adv(i,j,k) = 0.d0
              end do
           end do
        end do
     end if
  end if
  
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
  ! if ( bc(3,1,1).eq.EXT_DIR .and. adv_l3.lt.domlo(3)) then
  !    print *,'myfill: SHOULD NEVER GET HERE bc(3,1,1) .eq. EXT_DIR) '
  !    stop
  ! end if
  
  !     ZHI
  if ( bc(3,2,1).eq.EXT_DIR .and. adv_h3.gt.domhi(3)) then
     print *,'myfill: SHOULD NEVER GET HERE bc(3,2,1) .eq. EXT_DIR) '
     stop
  end if
  
end subroutine rns_myfill


! Fill z-momentum
subroutine rns_mzfill(adv,adv_l1,adv_l2,adv_l3,adv_h1,adv_h2,adv_h3, &
     domlo,domhi,delta,xlo,time,bc)

  use meth_params_module
  use probdata_module
  
  implicit none
  include 'bc_types.fi'
  integer adv_l1,adv_l2,adv_l3,adv_h1,adv_h2,adv_h3
  integer bc(3,2,*)
  integer domlo(3), domhi(3)
  double precision delta(3), xlo(3), time
  double precision adv(adv_l1:adv_h1,adv_l2:adv_h2,adv_l3:adv_h3)

  integer i,j,k,n,ii,jj,iwrk
  double precision x,y,xg,yg,r,eta,Yt(nspec),Tt,u3t,rhot,rwrk
  double precision, parameter :: gp(2) = (/ -1.d0/sqrt(3.d0), 1.d0/sqrt(3.d0) /)
  double precision, parameter :: wgt = 0.25d0

  if (.not. jet_initialized) then
     call init_jet()
  end if

  call filcc(adv,adv_l1,adv_l2,adv_l3,adv_h1,adv_h2,adv_h3, &
       domlo,domhi,delta,xlo,bc)

  ! zlo
  if (adv_l3.lt.domlo(3)) then
     if (bc(3,1,1) .eq. EXT_DIR) then

        !$omp parallel do private(i,j,k,n,ii,jj,iwrk,x,y,xg,yg,r,eta) &
        !$omp private(Yt,Tt,u3t,rhot,rwrk) &
        !$omp collapse(2)
        do k = adv_l3, domlo(3)-1
        do j = adv_l2, adv_h2
        do i = adv_l1, adv_h1

           x = (DBLE(i-adv_l1)+.5d0)*delta(1)+xlo(1)
           y = (DBLE(j-adv_l2)+.5d0)*delta(2)+xlo(2)

           adv(i,j,k) = 0.d0

           do jj=1,2
              yg = y + 0.5d0*delta(2)*gp(jj)
              do ii=1,2
                 xg = x + 0.5d0*delta(1)*gp(ii)
                    
                 r = sqrt(xg*xg+yg*yg)

                 eta = 0.5d0*(1.d0-tanh((r -splitx)/xfrontw))

                 do n=1,nspec
                    Yt(n) = eta*fuel_Y(n) + (1.d0-eta)*air_Y(n)
                 end do
                 Tt  = eta * T_in + (1.d0-eta) * T_co
                 u3t = eta * vn_in + (1.d0-eta) * vn_co                    
                    
                 CALL CKRHOY(pamb,Tt,Yt,IWRK,RWRK,rhot)

                 adv(i,j,k) = adv(i,j,k) + wgt*rhot*u3t

              end do
           end do

        end do
        end do
        end do
        !$omp end parallel do

     else
        print *,'mzfill: SHOULD NEVER GET HERE bc(3,1,1) .ne. EXT_DIR) '
        stop
     end if
  end if
  
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
  ! if ( bc(3,1,1).eq.EXT_DIR .and. adv_l3.lt.domlo(3)) then
  !    print *,'mzfill: SHOULD NEVER GET HERE bc(3,1,1) .eq. EXT_DIR) '
  !    stop
  ! end if
  
  !     ZHI
  if ( bc(3,2,1).eq.EXT_DIR .and. adv_h3.gt.domhi(3)) then
     print *,'mzfill: SHOULD NEVER GET HERE bc(3,2,1) .eq. EXT_DIR) '
     stop
  end if
  
end subroutine rns_mzfill


! Fill temperature
subroutine rns_tempfill(adv,adv_l1,adv_l2,adv_l3,adv_h1,adv_h2,adv_h3, &
     domlo,domhi,delta,xlo,time,bc)
  
  use meth_params_module
  use probdata_module

  implicit none
  include 'bc_types.fi'
  integer adv_l1,adv_l2,adv_l3,adv_h1,adv_h2,adv_h3
  integer bc(3,2,*)
  integer domlo(3), domhi(3)
  double precision delta(3), xlo(3), time
  double precision adv(adv_l1:adv_h1,adv_l2:adv_h2,adv_l3:adv_h3)

  integer i,j,k,ii,jj
  double precision x,y,xg,yg,r,eta,Tt
  double precision, parameter :: gp(2) = (/ -1.d0/sqrt(3.d0), 1.d0/sqrt(3.d0) /)
  double precision, parameter :: wgt = 0.25d0

  if (.not. jet_initialized) then
     call init_jet()
  end if

  call filcc(adv,adv_l1,adv_l2,adv_l3,adv_h1,adv_h2,adv_h3, &
       domlo,domhi,delta,xlo,bc)
  
  ! zlo
  if (adv_l3.lt.domlo(3)) then
     if (bc(3,1,1) .eq. EXT_DIR) then

        !$omp parallel do private(i,j,k,ii,jj,x,y,xg,yg,r,eta,Tt) &
        !$omp collapse(2)
        do k = adv_l3, domlo(3)-1
        do j = adv_l2, adv_h2
        do i = adv_l1, adv_h1

           x = (DBLE(i-adv_l1)+.5d0)*delta(1)+xlo(1)
           y = (DBLE(j-adv_l2)+.5d0)*delta(2)+xlo(2)

           adv(i,j,k) = 0.d0

           do jj=1,2
              yg = y + 0.5d0*delta(2)*gp(jj)
              do ii=1,2
                 xg = x + 0.5d0*delta(1)*gp(ii)
                    
                 r = sqrt(xg*xg+yg*yg)
                 
                 eta = 0.5d0*(1.d0-tanh((r -splitx)/xfrontw))

                 Tt  = eta * T_in + (1.d0-eta) * T_co

                 adv(i,j,k) = adv(i,j,k) + wgt*Tt

              end do
           end do

        end do
        end do
        end do
        !$omp end parallel do

     else
        print *,'grpfill: SHOULD NEVER GET HERE bc(3,1,1) .ne. EXT_DIR) '
        stop
     end if
  end if

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
  ! if ( bc(3,1,1).eq.EXT_DIR .and. adv_l3.lt.domlo(3)) then
  !    print *,'tempfill: SHOULD NEVER GET HERE bc(3,1,1) .eq. EXT_DIR) '
  !    stop
  ! end if
  
  !     ZHI
  if ( bc(3,2,1).eq.EXT_DIR .and. adv_h3.gt.domhi(3)) then
     print *,'tempfill: SHOULD NEVER GET HERE bc(3,2,1) .eq. EXT_DIR) '
     stop
  end if
  
end subroutine rns_tempfill
