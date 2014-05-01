! Fill the entire state
subroutine rns_grpfill(adv,adv_l1,adv_l2,adv_h1,adv_h2, &
     domlo,domhi,delta,xlo,time,bc_in)
 
  use meth_params_module
  use probdata_module
  use sdc_boundary_module, only : isFEval
  
  implicit none
  include 'bc_types.fi'
  integer adv_l1,adv_l2,adv_h1,adv_h2
  integer bc_in(2,2,*)
  integer domlo(2), domhi(2)
  double precision delta(2), xlo(2), time
  double precision adv(adv_l1:adv_h1,adv_l2:adv_h2,NVAR)
  
  integer i, j, n, iwrk, ii
  double precision :: x, xg, facx, fact, sigma, eta, Pi, eta1
  double precision rhot,u1t,u2t,Tt,et,Yt(NSPEC),rwrk
  double precision, parameter :: gp(2) = (/ -1.d0/sqrt(3.d0), 1.d0/sqrt(3.d0) /)
  double precision, parameter :: wgt = 0.5d0

  integer bc(2,2,NVAR)

  if (.not. dmejet_initialized) then
     call init_DME_jet()
  end if

  bc = bc_in(:,:,1:NVAR)
  if (isFEval) bc(2,1,:) = FOEXTRAP

  !$omp parallel do private(n)
  do n = 1,NVAR
     call filcc(adv(adv_l1,adv_l2,n), &
          adv_l1,adv_l2,adv_h1,adv_h2, &
          domlo,domhi,delta,xlo,bc(1,1,n))
  enddo
  !$omp end parallel do

  if (isFEval) return

!        YLO
  if (adv_l2.lt.domlo(2)) then
     if (bc(2,1,1).eq.EXT_DIR) then

        Pi = 4.d0*atan(1.d0)
        facx = 2.d0*Pi/((domhi(1)-domlo(1)+1)*delta(1))
        fact = sin(2.d0*Pi*time/inflow_period)

        sigma = 2.5d0*xfrontw*splitx

        ! fill the corners too
        !$omp parallel do private(i,j,n,iwrk,ii,x,xg,eta,eta1,rhot,u1t,u2t,Tt,et,Yt,rwrk) &
        !$omp collapse(2)
        do j = adv_l2, domlo(2)-1 
           do i = adv_l1,adv_h1
        
              x = (DBLE(i-adv_l1)+.5d0)*delta(1)+xlo(1)
              
              adv(i,j,:) = 0.d0

              do ii=1,2
                 xg = x + 0.5d0*delta(1)*gp(ii)

                 if (prob_type .eq. 0) then

                    eta = 0.5d0 * (tanh((xg + splitx)/sigma)   &
                         &       - tanh((xg - splitx)/sigma))

                    do n=1,nspec
                       Yt(n) = eta*fuel_Y(n) + (1.d0-eta)*air_Y(n)
                    end do
                    Tt  = eta * T_in + (1.d0-eta) * T_co
                    u1t = 0.d0
                    u2t = eta *vn_in + (1.d0-eta) *vn_co &
                         + inflow_vnmag*eta*sin(xg*facx)*fact

                 else if (prob_type .eq. 1) then

                    eta = 0.5d0 * (tanh((xg + splitx)/Tfrontw)  &
                         &       - tanh((xg - splitx)/Tfrontw))
                    eta1 = 0.5d0 * (tanh((xg + blobr)/xfrontw)  &
                         &        - tanh((xg - blobr)/xfrontw))

                    do n=1,nspec
                       Yt(n) = eta*fuel_Y(n) + (1.d0-eta)*air_Y(n)
                    end do
                    Tt  = eta * T_in + (1.d0-eta) * T_co
                    u1t = 0.d0
                    u2t = eta1 * vn_in + (1.d0-eta1) * vn_co

                 else 

                    eta = 0.5d0 * (tanh((xg + splitx)/xfrontw)  &
                         &       - tanh((xg - splitx)/xfrontw))

                    do n=1,nspec
                       Yt(n) = eta*fuel_Y(n) + (1.d0-eta)*air_Y(n)
                    end do
                    Tt  = eta * T_in + (1.d0-eta) * T_co
                    u1t = 0.d0
                    u2t = eta * vn_in + (1.d0-eta) * vn_co                    
                    
                 end if
       
                 CALL CKRHOY(pamb,Tt,Yt,IWRK,RWRK,rhot)
                 call CKUBMS(Tt,Yt,IWRK,RWRK,et)

                 adv(i,j,URHO ) = adv(i,j,URHO ) + wgt*rhot
                 adv(i,j,UMX  ) = adv(i,j,UMX  ) + wgt*rhot*u1t
                 adv(i,j,UMY  ) = adv(i,j,UMY  ) + wgt*rhot*u2t
                 adv(i,j,UEDEN) = adv(i,j,UEDEN) + wgt*rhot*(et + 0.5d0*(u1t**2+u2t**2))
                 adv(i,j,UTEMP) = adv(i,j,UTEMP) + wgt*Tt
                 do n=1, NSPEC
                    adv(i,j,UFS+n-1) = adv(i,j,UFS+n-1) + wgt*rhot*Yt(n)
                 end do
              
              end do

           end do
        end do
        !$omp end parallel do
     else
        print *,'grpfill: SHOULD NEVER GET HERE bc(2,1,1) .ne. EXT_DIR) '
        stop
     end if
  end if
  
  do n = 1, NVAR

     !     XLO
     if ( bc(1,1,n).eq.EXT_DIR .and. adv_l1.gt.domlo(1)) then
        print *,'grpfill: SHOULD NEVER GET HERE bc(1,1,n) .eq. EXT_DIR) '
        stop
     end if

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
     
  end do
  
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
  
  use meth_params_module
  use probdata_module

  implicit none
  include 'bc_types.fi'
  integer adv_l1,adv_l2,adv_h1,adv_h2
  integer bc(2,2,*)
  integer domlo(2), domhi(2)
  double precision delta(2), xlo(2), time
  double precision adv(adv_l1:adv_h1,adv_l2:adv_h2)

  integer :: i, j, n, iwrk, ii
  double precision :: x, xg, sigma, eta
  double precision rhot,Tt,Yt(NSPEC),rwrk
  double precision, parameter :: gp(2) = (/ -1.d0/sqrt(3.d0), 1.d0/sqrt(3.d0) /)
  double precision, parameter :: wgt = 0.5d0

  if (.not. dmejet_initialized) then
     call init_DME_jet()
  end if

  call filcc(adv,adv_l1,adv_l2,adv_h1,adv_h2,domlo,domhi,delta,xlo,bc)

  !        YLO
  if (adv_l2.lt.domlo(2)) then
     if (bc(2,1,1).eq.EXT_DIR) then

        sigma = 2.5d0*xfrontw*splitx

        ! fill the corners too
        !$omp parallel do private(i,j,n,iwrk,ii,x,xg,eta,rhot,Tt,Yt,rwrk) &
        !$omp collapse(2)
        do j = adv_l2, domlo(2)-1 
           do i = adv_l1,adv_h1

              x = (DBLE(i-adv_l1)+.5d0)*delta(1)+xlo(1)
              
              adv(i,j) = 0.d0

              do ii=1,2

                 xg = x + 0.5d0*delta(1)*gp(ii)

                 if (prob_type .eq. 0) then
                    eta = 0.5d0 * (tanh((xg + splitx)/sigma)   &
                         &       - tanh((xg - splitx)/sigma))
                 else if (prob_type .eq. 1) then
                    eta = 0.5d0 * (tanh((xg + splitx)/Tfrontw)  &
                         &       - tanh((xg - splitx)/Tfrontw))
                 else if (prob_type .eq. 2) then
                    eta = 0.5d0 * (tanh((xg + splitx)/xfrontw)  &
                         &       - tanh((xg - splitx)/xfrontw))
                 end if

                 do n=1,nspec
                    Yt(n) = eta*fuel_Y(n) + (1.d0-eta)*air_Y(n)
                 end do
                 Tt  = eta * T_in + (1.d0-eta) * T_co
       
                 CALL CKRHOY(pamb,Tt,Yt,IWRK,RWRK,rhot)
                 adv(i,j) = adv(i,j) + wgt*rhot
              
              end do

           end do
        end do
        !$omp end parallel do
     else
        print *,'SHOULD NEVER GET HERE bc(2,1,1) .ne. EXT_DIR) '
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
!  if ( bc(2,1,1).eq.EXT_DIR .and. adv_l2.lt.domlo(2)) then
!     print *,'denfill: SHOULD NEVER GET HERE bc(2,1,1) .eq. EXT_DIR) '
!     stop
!  end if
  
  !     YHI
  if ( bc(2,2,1).eq.EXT_DIR .and. adv_h2.gt.domhi(2)) then
     print *,'denfill: SHOULD NEVER GET HERE bc(2,2,1) .eq. EXT_DIR) '
     stop
  end if
  
end subroutine rns_denfill


! Fill x-momentum
subroutine rns_mxfill(adv,adv_l1,adv_l2,adv_h1,adv_h2, &
     domlo,domhi,delta,xlo,time,bc)

  use meth_params_module
  use probdata_module
  
  implicit none
  include 'bc_types.fi'
  integer adv_l1,adv_l2,adv_h1,adv_h2
  integer bc(2,2,*)
  integer domlo(2), domhi(2)
  double precision delta(2), xlo(2), time
  double precision adv(adv_l1:adv_h1,adv_l2:adv_h2)

  integer :: i, j

  if (.not. dmejet_initialized) then
     call init_DME_jet()
  end if
  
  call filcc(adv,adv_l1,adv_l2,adv_h1,adv_h2,domlo,domhi,delta,xlo,bc)

  !        YLO
  if (adv_l2.lt.domlo(2)) then
     if (bc(2,1,1).eq.EXT_DIR) then

        ! fill the corners too
        do j = adv_l2, domlo(2)-1 
           do i = adv_l1,adv_h1

              adv(i,j) = 0.d0
              
           end do
        end do
     else
        print *,'SHOULD NEVER GET HERE bc(1,1,1) .ne. EXT_DIR) '
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
     print *,'mxfill: SHOULD NEVER GET HERE bc(1,2,1) .eq. EXT_DIR) '
     stop
  end if
  
!  !     YLO
!  if ( bc(2,1,1).eq.EXT_DIR .and. adv_l2.lt.domlo(2)) then
!     print *,'mxfill: SHOULD NEVER GET HERE bc(2,1,1) .eq. EXT_DIR) '
!     stop
!  end if
  
  !     YHI
  if ( bc(2,2,1).eq.EXT_DIR .and. adv_h2.gt.domhi(2)) then
     print *,'mxfill: SHOULD NEVER GET HERE bc(2,2,1) .eq. EXT_DIR) '
     stop
  end if
  
end subroutine rns_mxfill


! Fill y-momentum
subroutine rns_myfill(adv,adv_l1,adv_l2,adv_h1,adv_h2, &
     domlo,domhi,delta,xlo,time,bc)

  use meth_params_module
  use probdata_module
  
  implicit none
  include 'bc_types.fi'
  integer adv_l1,adv_l2,adv_h1,adv_h2
  integer bc(2,2,*)
  integer domlo(2), domhi(2)
  double precision delta(2), xlo(2), time
  double precision adv(adv_l1:adv_h1,adv_l2:adv_h2)
  
  integer :: i, j, n, iwrk, ii
  double precision :: x, xg, facx, fact, sigma, eta, Pi, eta1
  double precision rhot,u2t,Tt,Yt(NSPEC),rwrk
  double precision, parameter :: gp(2) = (/ -1.d0/sqrt(3.d0), 1.d0/sqrt(3.d0) /)
  double precision, parameter :: wgt = 0.5d0

  if (.not. dmejet_initialized) then
     call init_DME_jet()
  end if

  call filcc(adv,adv_l1,adv_l2,adv_h1,adv_h2,domlo,domhi,delta,xlo,bc)

!        YLO
  if (adv_l2.lt.domlo(2)) then
     if (bc(2,1,1).eq.EXT_DIR) then

        Pi = 4.d0*atan(1.d0)
        facx = 2.d0*Pi/((domhi(1)-domlo(1)+1)*delta(1))
        fact = sin(2.d0*Pi*time/inflow_period)

        sigma = 2.5d0*xfrontw*splitx

        ! fill the corners too
        !$omp parallel do private(i,j,n,iwrk,ii,x,xg,eta,eta1,rhot,u2t,Tt,Yt,rwrk) &
        !$omp collapse(2)
        do j = adv_l2, domlo(2)-1 
           do i = adv_l1,adv_h1

              x = (DBLE(i-adv_l1)+.5d0)*delta(1)+xlo(1)
              
              adv(i,j) = 0.d0

              do ii=1,2
                 xg = x + 0.5d0*delta(1)*gp(ii)
                 
                 if (prob_type .eq. 0) then
                    eta = 0.5d0 * (tanh((xg + splitx)/sigma)   &
                         &       - tanh((xg - splitx)/sigma))
                 else if (prob_type .eq. 1) then
                    eta = 0.5d0 * (tanh((xg + splitx)/Tfrontw)  &
                         &       - tanh((xg - splitx)/Tfrontw))
                    eta1 = 0.5d0 * (tanh((xg + blobr)/xfrontw)  &
                         &        - tanh((xg - blobr)/xfrontw))
                 else if (prob_type .eq. 2) then
                    eta = 0.5d0 * (tanh((xg + splitx)/xfrontw)  &
                         &       - tanh((xg - splitx)/xfrontw))
                 end if

                 do n=1,nspec
                    Yt(n) = eta*fuel_Y(n) + (1.d0-eta)*air_Y(n)
                 end do
                 Tt  = eta * T_in + (1.d0-eta) * T_co
                 if (prob_type .eq. 0) then 
                    u2t = eta *vn_in + (1.d0-eta) *vn_co &
                         + inflow_vnmag*eta*sin(xg*facx)*fact
                 else if (prob_type .eq. 1) then
                    u2t = eta1 * vn_in + (1.d0-eta1) * vn_co
                 else if (prob_type .eq. 2) then
                    u2t = eta * vn_in + (1.d0-eta) * vn_co
                 end if
       
                 CALL CKRHOY(pamb,Tt,Yt,IWRK,RWRK,rhot)

                 adv(i,j) = adv(i,j) + wgt*rhot*u2t
              
              end do

           end do
        end do
        !$omp end parallel do
     else
        print *,'SHOULD NEVER GET HERE bc(1,1,1) .ne. EXT_DIR) '
        stop
     end if
  end if

  !     XLO
  if ( bc(1,1,1).eq.EXT_DIR .and. adv_l1.lt.domlo(1)) then
     print *,'regfill: SHOULD NEVER GET HERE bc(1,1,1) .eq. EXT_DIR) '
     stop
  end if
  
  !     XHI
  if ( bc(1,2,1).eq.EXT_DIR .and. adv_h1.gt.domhi(1)) then
     print *,'myfill: SHOULD NEVER GET HERE bc(1,2,1) .eq. EXT_DIR) '
     stop
  end if
  
  !     YLO
!  if ( bc(2,1,1).eq.EXT_DIR .and. adv_l2.lt.domlo(2)) then
!     print *,'myfill: SHOULD NEVER GET HERE bc(2,1,1) .eq. EXT_DIR) '
!     stop
!  end if
  
  !     YHI
  if ( bc(2,2,1).eq.EXT_DIR .and. adv_h2.gt.domhi(2)) then
     print *,'myfill: SHOULD NEVER GET HERE bc(2,2,1) .eq. EXT_DIR) '
     stop
  end if
  
end subroutine rns_myfill


! Fill temperature
subroutine rns_tempfill(adv,adv_l1,adv_l2,adv_h1,adv_h2, &
     domlo,domhi,delta,xlo,time,bc)

  use probdata_module
  
  implicit none
  include 'bc_types.fi'
  integer adv_l1,adv_l2,adv_h1,adv_h2
  integer bc(2,2,*)
  integer domlo(2), domhi(2)
  double precision delta(2), xlo(2), time
  double precision adv(adv_l1:adv_h1,adv_l2:adv_h2)
  
  integer :: i, j, ii
  double precision :: x, xg, sigma, eta
  double precision :: Tt
  double precision, parameter :: gp(2) = (/ -1.d0/sqrt(3.d0), 1.d0/sqrt(3.d0) /)
  double precision, parameter :: wgt = 0.5d0

  if (.not. dmejet_initialized) then
     call init_DME_jet()
  end if

  call filcc(adv,adv_l1,adv_l2,adv_h1,adv_h2,domlo,domhi,delta,xlo,bc)
  
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
   if (adv_l2.lt.domlo(2)) then
     if (bc(2,1,1).eq.EXT_DIR) then

        sigma = 2.5d0*xfrontw*splitx

        ! fill the corners too
        !$omp parallel do private(i,j,ii,x,xg,eta,Tt) collapse(2)
        do j = adv_l2, domlo(2)-1 
           do i = adv_l1,adv_h1

              x = (DBLE(i-adv_l1)+.5d0)*delta(1)+xlo(1)
              
              adv(i,j) = 0.d0

              do ii=1,2

                 xg = x + 0.5d0*delta(1)*gp(ii)

                 if (prob_type .eq. 0) then
                    eta = 0.5d0 * (tanh((xg + splitx)/sigma)   &
                         &       - tanh((xg - splitx)/sigma))
                 else if (prob_type .eq. 1) then
                    eta = 0.5d0 * (tanh((xg + splitx)/Tfrontw)  &
                         &       - tanh((xg - splitx)/Tfrontw))
                 else if (prob_type .eq. 2) then
                    eta = 0.5d0 * (tanh((xg + splitx)/xfrontw)  &
                         &       - tanh((xg - splitx)/xfrontw))
                 end if

                 Tt  = eta * T_in + (1.d0-eta) * T_co

                 adv(i,j) = adv(i,j) + wgt*Tt
              
              end do

           end do
        end do
        !$omp end parallel do
     else
        print *,'SHOULD NEVER GET HERE bc(2,1,1) .ne. EXT_DIR) '
        stop
     end if
  end if
  
  !     YHI
  if ( bc(2,2,1).eq.EXT_DIR .and. adv_h2.gt.domhi(2)) then
     print *,'tempfill: SHOULD NEVER GET HERE bc(2,2,1) .eq. EXT_DIR) '
     stop
  end if
  
end subroutine rns_tempfill
