cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     NOTE: The comments and equation references correspond to the final
c           published version available online at:
c
c     http://www.tandfonline.com/doi/full/10.1080/13647830.2012.701019
c     WILL'S EDIT TO THE CODE
c     Created: 06/03/2015
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      
      function quadrature_1(f0, f1, f2, dt) result(q)
         implicit none
         
         real*8 f0, f1, f2, dt
         real*8 q
         
         q = (5.0*f0/24.0 + f1/3.0 - f2/24.0)*dt
      end function quadrature_1


      function quadrature_2(f0, f1, f2, dt) result(q)
         implicit none
         
         real*8 f0, f1, f2, dt
         real*8 q
         
         q = (-f0/24.0 + f1/3.0 + 5*f2/24.0)*dt
      end function quadrature_2
      
      subroutine advance(vel_old,vel_new,scal_old,scal_new,
     $                   I_R,press_old,press_new,
     $                   divu_old,divu_new,dSdt,beta_n,beta_2,
     $                   dx,dt,lo,hi,bc,delta_chi,istep)

      implicit none

      include 'spec.h'

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     state variables held in scal_old and scal_new contain:
c
c     1           = Density
c     2           = Temp
c     3           = RhoH
c     4           = RhoRT (i.e., "ptherm")
c     5:5+Nspec-1 = FirstSpec:LastSpec (rho*Y_k)
c
c     For the CHEMH mechanism this code defaults to, Nspec=9 and nscal=13
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c     cell-centered, 2 ghost cells
      real*8    vel_old(0:nlevs-1,-2:nfine+1)
      real*8    vel_new(0:nlevs-1,-2:nfine+1)
      real*8   scal_old(0:nlevs-1,-2:nfine+1,nscal)
      real*8   scal_new(0:nlevs-1,-2:nfine+1,nscal)
      
      real*8   scal_1_prev(0:nlevs-1,-2:nfine+1,nscal)
      real*8   scal_2_prev(0:nlevs-1,-2:nfine+1,nscal)
      real*8   scal_1_AD(0:nlevs-1,-2:nfine+1,nscal)
      real*8   scal_2_AD(0:nlevs-1,-2:nfine+1,nscal)
      
      real*8   scal_1_curr(0:nlevs-1,-2:nfine+1,nscal)
      real*8   scal_2_curr(0:nlevs-1,-2:nfine+1,nscal)
      
c     cell-centered, 1 ghost cell
      real*8        I_R(0:nlevs-1,-1:nfine  ,0:Nspec)
      real*8   beta_n(0:nlevs-1,-1:nfine  ,nscal)
      real*8   beta_1(0:nlevs-1,-1:nfine  ,nscal)
      real*8   beta_2(0:nlevs-1,-1:nfine  ,nscal)
      real*8   divu_old(0:nlevs-1,-1:nfine)
      real*8   divu_new(0:nlevs-1,-1:nfine)

c     cell-centered, no ghost cells
      real*8       dSdt(0:nlevs-1, 0:nfine-1)
      real*8  delta_chi(0:nlevs-1, 0:nfine-1)
      
      real*8  delta_chi_0(0:nlevs-1, 0:nfine-1)
      real*8  delta_chi_1(0:nlevs-1, 0:nfine-1)
      real*8  delta_chi_2(0:nlevs-1, 0:nfine-1)
      real*8  delta_chi_3(0:nlevs-1, 0:nfine-1)

c     nodal, 1 ghost cell
      real*8  press_old(0:nlevs-1,-1:nfine+1)
      real*8  press_new(0:nlevs-1,-1:nfine+1)

      integer lo(0:nlevs-1)
      integer hi(0:nlevs-1)
      integer bc(0:nlevs-1,2)
      real*8  dx(0:nlevs-1)
      real*8  dt(0:nlevs-1)
      integer istep

c     local variables

c     cell-centered, 1 ghost cell
      real*8       mu_old(0:nlevs-1,-1:nfine)
      real*8       mu_new(0:nlevs-1,-1:nfine)
      real*8     mu_dummy(0:nlevs-1,-1:nfine)
      real*8           gp(0:nlevs-1,-1:nfine)
      real*8         visc(0:nlevs-1,-1:nfine)
      real*8     I_R_divu(0:nlevs-1,-1:nfine,  0:Nspec)
      real*8     I_R_temp(0:nlevs-1,-1:nfine,  0:Nspec)
      real*8     diff_n(0:nlevs-1,-1:nfine,  nscal)
      real*8     diff_1(0:nlevs-1,-1:nfine,  nscal)
      real*8     diff_1_curr(0:nlevs-1,-1:nfine,  nscal)
      real*8     diff_2(0:nlevs-1,-1:nfine,  nscal)
      real*8     diff_hat_1(0:nlevs-1,-1:nfine,  nscal)
      real*8     diff_hat_2(0:nlevs-1,-1:nfine,  nscal)
      real*8     diff_tmp(0:nlevs-1,-1:nfine,  nscal)
      real*8 diffdiff_n(0:nlevs-1,-1:nfine)
      real*8 diffdiff_1(0:nlevs-1,-1:nfine)
      real*8 diffdiff_2(0:nlevs-1,-1:nfine)
      real*8 diffdiff_1_curr(0:nlevs-1,-1:nfine)
      real*8  divu_effect(0:nlevs-1,-1:nfine)
      
      real*8  divu_1(0:nlevs-1,-1:nfine)      
      real*8  divu_2(0:nlevs-1,-1:nfine)

c     cell-centered, no ghost cells
      real*8      rhohalf(0:nlevs-1, 0:nfine-1)
      real*8        alpha(0:nlevs-1, 0:nfine-1)
      real*8      vel_Rhs(0:nlevs-1, 0:nfine-1)
      
c     WILL: store old and new values for A in order to pass 
c           piecewise linear to the ODE solver
c      real*8         aofs(0:nlevs-1, 0:nfine-1,nscal)
      real*8         aofs_n(0:nlevs-1, 0:nfine-1, nscal)
      real*8         aofs_1(0:nlevs-1, 0:nfine-1, nscal)
      real*8         aofs_2(0:nlevs-1, 0:nfine-1, nscal)
      real*8         aofs_avg(0:nlevs-1, 0:nfine-1, nscal)
      real*8         aofs_curr(0:nlevs-1, 0:nfine-1, nscal)
      
c     WILLL: right-hand side used in backward euler solve
      real*8       berhs(0:nfine-1,nscal)

      real*8 gamma_lo(0:nlevs-1, 0:nfine-1,Nspec)
      real*8 gamma_hi(0:nlevs-1, 0:nfine-1,Nspec)
      real*8    const_src(0:nlevs-1, 0:nfine-1,nscal)
      real*8          Rhs(0:nlevs-1, 0:nfine-1,nscal)
      real*8         dRhs(0:nlevs-1, 0:nfine-1,0:Nspec)

c     nodal, no ghost cells
c     WILL: changed this to store new and old velocities
c      real*8       macvel(0:nlevs-1, 0:nfine  )
      real*8       macvel_n(0:nlevs-1, 0:nfine  )
      real*8       macvel_1(0:nlevs-1, 0:nfine  )
      real*8       macvel_2(0:nlevs-1, 0:nfine  )
      real*8       macvel_avg(0:nlevs-1, 0:nfine)
      real*8      veledge(0:nlevs-1, 0:nfine    )

      real*8 Y(Nspec),WDOTK(Nspec),C(Nspec),RWRK
      real*8 wdot_n(0:nfine-1, Nspec)
      real*8 wdot_1(0:nfine-1, Nspec)
      real*8 wdot_2(0:nfine-1, Nspec)
      real*8 cpmix,rhocp,vel_theta,be_cn_theta
      
      integer i,is,misdc,n,rho_flag,IWRK
      
      real*8 dt1, dt2
      logical const_quad, do_be, provide_wdot
      
      real*8 :: quadrature_1, quadrature_2
      
      dt1 = dt(0)/2.0
      dt2 = dt(0)/2.0
      
      provide_wdot = .true.

c     "diffdiff" means "differential diffusion", which corresponds to
c     sum_m div [ h_m (rho D_m - lambda/cp) grad Y_m ]
c     in equation (3)
      
      diffdiff_n(0,:) = 0.d0
      diffdiff_1(0,:) = 0.d0
      diffdiff_2(0,:) = 0.d0

      print *,'advance: at start of time step'

ccccccccccccccccccccccccccccccccccccccccccc
c     Step 1: Compute advection velocities
ccccccccccccccccccccccccccccccccccccccccccc
 
c     compute cell-centered grad pi from nodal pi
      do i=lo(0)-1,hi(0)+1
         gp(0,i) = (press_old(0,i+1) - press_old(0,i)) / dx(0)
      enddo

      print *,'... predict edge velocities'

c     compute U^{ADV,*}
      call pre_mac_predict(vel_old(0,:),scal_old(0,:,:),gp(0,:),
     $                     macvel_n(0,:),dx(0),dt(0),lo(0),hi(0),bc(0,:))

c     compute ptherm = p(rho,T,Y)
c     this is needed for any dpdt-based correction scheme
      call compute_pthermo(scal_old(0,:,:),lo(0),hi(0),bc(0,:))

c     reset delta_chi
      delta_chi = 0.d0
      delta_chi_0 = 0.d0
      delta_chi_1 = 0.d0
      delta_chi_2 = 0.d0
      delta_chi_3 = 0.d0

c     delta_chi = delta_chi + (peos-p0)/(dt*peos) + (1/peos) u dot grad peos
      call add_dpdt(scal_old(0,:,:),scal_old(0,:,RhoRT),
     $              delta_chi_0(0,:),macvel_n(0,:),dx(0),dt(0),
     $              lo(0),hi(0),bc(0,:))

c     S_hat^{n+1/2} = S^{n+1/2} + delta_chi
      do i=lo(0),hi(0)
         divu_effect(0,i) = divu_old(0,i) + delta_chi_0(0,i)
      end do

c     mac projection
c     macvel_n will now satisfy div(umac) = S_hat^{n+1/2}
      call macproj(macvel_n(0,:),scal_old(0,:,Density),
     &             divu_effect(0,:),dx,lo(0),hi(0),bc(0,:))

c     WILL: compute A^n
      call scal_aofs(scal_old(0,:,:),macvel_n(0,:),aofs_n(0,:,:),
     $                  dx(0),dt(0),lo(0),hi(0),bc(0,:))

ccccccccccccccccccccccccccccccccccccccccccc
c     Step 2: Advance thermodynamic variables
ccccccccccccccccccccccccccccccccccccccccccc

c     WILL: removed the strang algorithm (for simplicity)
      if (use_strang) then
         print *,'REMOVED STRANG ALGORITHM'
         stop
      else

ccccccccccccccccccccccccccccccccccccccccccc
c     Step 2: Advance thermodynamic variables (SDC algorithm)
ccccccccccccccccccccccccccccccccccccccccccc

c     diffusion solves in predictor are regular Crank-Nicolson
         be_cn_theta = 0.5d0

c     compute transport coefficients
c        rho D_m     (for species)
c        lambda / cp (for enthalpy)
c        lambda      (for temperature)
         call calc_diffusivities(scal_old(0,:,:),beta_n(0,:,:),
     &                           mu_old(0,:),lo(0),hi(0))

c     compute diffusion terms at time n
         print *,'... creating the diffusive terms with old data'

c     compute conservatively corrected div gamma_m 
c     also save gamma_m for computing diffdiff terms later
         call get_spec_visc_terms(scal_old(0,:,:),beta_n(0,:,:),
     &                            diff_n(0,:,FirstSpec:),
     &                            gamma_lo(0,:,:),gamma_hi(0,:,:),
     &                            dx(0),lo(0),hi(0))
c     compute div lambda/cp grad h (no differential diffusion)
         call get_rhoh_visc_terms(scal_old(0,:,:),beta_n(0,:,:),
     &                            diff_n(0,:,RhoH),dx(0),lo(0),hi(0))

         if (LeEQ1 .eq. 0) then
c     calculate differential diffusion "diffdiff" terms, i.e.,
c     sum_m div [ h_m (rho D_m - lambda/cp) grad Y_m ]
c     we pass in conservative gamma_m via gamma
c     we take lambda / cp from beta
c     we compute h_m using T from the first argument
c     we compute grad Y_m using Y_m from the second argument
            call get_diffdiff_terms(scal_old(0,:,:),scal_old(0,:,:),
     $                              gamma_lo(0,:,:),
     $                              gamma_hi(0,:,:),beta_n(0,:,:),
     $                              diffdiff_n(0,:),dx(0),lo(0),hi(0))
         end if

         do i=lo(0),hi(0)
            do n=1,Nspec
               C(n) = scal_old(0,i,FirstSpec+n-1)*invmwt(n)
            end do
            call CKWC(scal_old(0,i,Temp),C,IWRK,RWRK,wdot_n(i,:))
            do n=1,Nspec
               I_R(0,i,n) = wdot_n(i,n)*mwt(n)
            end do
         end do
         
C----------------------------------------------------------------
c     Begin initial predictor
C----------------------------------------------------------------

         if (fancy_predictor .eq. 1) then
            print *,'REMOVED FANCY PREDICTOR'
            stop
         end if 
c        non-fancy predictor that simply sets scal_new = scal_old
         
         scal_1_prev = scal_old
         scal_2_prev = scal_old
         scal_1_curr = scal_old
         scal_2_curr = scal_old
         
         wdot_1 = wdot_n
         wdot_2 = wdot_n
         
         aofs_1 = aofs_n
         aofs_2 = aofs_n
         aofs_curr = aofs_n
         
         diff_1 = diff_n
         diff_2 = diff_n
         
         macvel_1 = macvel_n
         macvel_2 = macvel_n
         
         divu_1   = divu_old
         divu_2   = divu_old
         divu_new = divu_old

C----------------------------------------------------------------
c     End initial predictor
C----------------------------------------------------------------

C----------------------------------------------------------------
c     Begin MISDC iterations
C----------------------------------------------------------------

c     diffusion solves in SDC iterations are iterative corrections
c     that have a backward Euler character
         be_cn_theta = 1.d0

         do misdc = 1, misdc_iterMAX
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!! START SUBSTEP ONE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            do i=lo(0),hi(0)
               do n = 1,Nspec
                  Y(n) = scal_1_prev(0,i,FirstSpec+n-1)/scal_1_prev(0,i,Density)
               enddo
               ! compute the production rates from the previous iterate
               call CKWYR(scal_1_prev(0,i,Density), scal_1_prev(0,i,Temp), Y, iwrk, rwrk, wdot_1(i,:))
               
               do n=1,Nspec
                  I_R_divu(0,i,n) = wdot_1(i,n)*mwt(n)
               end do
            end do
c     compute transport coefficients
c        rho D_m     (for species)
c        lambda / cp (for enthalpy)
c        lambda      (for temperature)
           call calc_diffusivities(scal_1_prev(0,:,:),beta_1(0,:,:),
     &                              mu_dummy(0,:),lo(0),hi(0))
     
c     compute a conservative div gamma_m
c     save gamma_m for differential diffusion computation
            call get_spec_visc_terms(scal_1_prev(0,:,:),beta_1(0,:,:),
     &                               diff_1(0,:,FirstSpec:),
     &                               gamma_lo(0,:,:),
     &                               gamma_hi(0,:,:),
     &                               dx(0),lo(0),hi(0))
c     compute div lambda/cp grad h (no differential diffusion)
            call get_rhoh_visc_terms(scal_1_prev(0,:,:),beta_1(0,:,:),
     &                               diff_1(0,:,RhoH),dx(0),lo(0),hi(0))

            if (LeEQ1 .eq. 0) then
c     calculate differential diffusion "diffdiff" terms, i.e.,
c     sum_m div [ h_m (rho D_m - lambda/cp) grad Y_m ]
c     we pass in conservative gamma_m via gamma
c     we take lambda / cp from beta
c     we compute h_m using T from the first argument
c     we compute grad Y_m using Y_m from the second argument
               call get_diffdiff_terms(scal_1_prev(0,:,:),scal_1_prev(0,:,:),
     $                                 gamma_lo(0,:,:),
     $                                 gamma_hi(0,:,:),beta_1(0,:,:),
     $                                 diffdiff_1(0,:),dx(0),lo(0),hi(0))
            end if

c     compute advective flux divergence
            call scal_aofs(scal_1_prev(0,:,:),macvel_1(0,:),aofs_1(0,:,:),
     $                     dx(0),dt1,lo(0),hi(0),bc(0,:))

            
            
c     either the mac velocities have changed, or this was not called earlier
c     because we are not using the fancy predictor

            scal_1_AD = scal_1_prev
            do i=lo(0),hi(0)
               do n=1,nscal
                  aofs_avg(0,i,n) = quadrature_1(aofs_n(0,i,n),
     $                                           aofs_1(0,i,n),
     $                                           aofs_2(0,i,n), dt)/dt1
               enddo
            enddo
            
            print *,'... update rho'
            call update_rho(scal_old(0,:,:),scal_1_AD(0,:,:),aofs_avg(0,:,:),
     &                      dt1,lo(0),hi(0),bc(0,:))
            aofs_avg = 0.5d0*(aofs_n + aofs_1)
            
c     update rhoY_m with advection terms and set up RHS for equation (47) C-N solve
            print *,'... do correction diffusion solve for species'
            do i=lo(0),hi(0)
               do n=1,Nspec
                  is = FirstSpec + n - 1
c     includes deferred correction term for species
                  dRhs(0,i,n) = -dt1*diff_1(0,i,is)
     &                 + quadrature_1(wdot_n(i,n)*mwt(n)+diff_n(0,i,is)+aofs_n(0,i,is),
     &                                wdot_1(i,n)*mwt(n)+diff_1(0,i,is)+aofs_1(0,i,is),
     &                                wdot_2(i,n)*mwt(n)+diff_2(0,i,is)+aofs_2(0,i,is),
     &                                dt(0))
               enddo
c     includes deferred correction term for enthalpy
c     differential diffusion will be added later
               dRhs(0,i,0) = -dt1*diff_1(0,i,RhoH)
     &              + quadrature_1(diff_n(0,i,RhoH)+aofs_n(0,i,RhoH),
     &                             diff_1(0,i,RhoH)+aofs_1(0,i,RhoH),
     &                             diff_2(0,i,RhoH)+aofs_2(0,i,RhoH),
     &                             dt(0))
            enddo
c     WILL: TODO: make sure this should be aofs_avg
            call update_spec(scal_old(0,:,:),scal_1_prev(0,:,:),
     &                       alpha(0,:),beta_n(0,:,:),
     &                       dRhs(0,0:,1:),Rhs(0,0:,FirstSpec:),
     &                       dx(0),dt1,be_cn_theta,lo(0),hi(0),bc(0,:))

c     Solve C-N system in equation (47) for \tilde{Y}_{m,AD}^{(k+1)}
            rho_flag = 2
            do n=1,Nspec
               is = FirstSpec + n - 1
               call cn_solve(scal_1_AD(0,:,:),alpha(0,:),beta_1(0,:,:),
     $                       Rhs(0,:,is),dx(0),dt1,is,be_cn_theta,
     $                       rho_flag,.false.,lo(0),hi(0),bc(0,:))
            enddo
            
            
            if (LeEQ1 .eq. 1) then
c     extract div gamma^n
               do i=lo(0),hi(0)
                  do n=1,Nspec
                     is = FirstSpec + n - 1
                     diff_hat_1(0,i,is) = (scal_1_AD(0,i,is)-scal_old(0,i,is))/dt1 
     $                    - dRhs(0,i,n)/dt1
                  enddo
               enddo

            else

c     compute conservatively corrected div gamma_m 
c     also save gamma_m for computing diffdiff terms later
               call get_spec_visc_terms(scal_1_AD(0,:,:),beta_1(0,:,:),
     &                                  diff_hat_1(0,:,FirstSpec:),
     &                                  gamma_lo(0,:,:),
     &                                  gamma_hi(0,:,:),
     &                                  dx(0),lo(0),hi(0))

c     add differential diffusion to forcing for enthalpy solve in equation (49)
               do i=lo(0),hi(0)
                  dRhs(0,i,0) = dRhs(0,i,0) 
     $                 + quadrature_1(diffdiff_n(0,i),
     $                                diffdiff_1(0,i),
     $                                diffdiff_2(0,i), dt(0))
               end do
            end if
            
            print *,'... do correction diffusion solve for rhoh'

c     WILL: TODO: double check this
c     update rhoh with advection terms and set up RHS for equation (49) C-N solve
            call update_rhoh(scal_old(0,:,:),scal_1_prev(0,:,:),
     &                       alpha(0,:),beta_n(0,:,:),
     &                       dRhs(0,:,0),Rhs(0,:,RhoH),dx(0),dt1,
     &                       be_cn_theta,lo(0),hi(0),bc(0,:))

c     Solve C-N system in equation (49) for h_{AD}^{(k+1)}
            rho_flag = 2
            call cn_solve(scal_1_AD(0,:,:),alpha(0,:),beta_1(0,:,:),
     $                    Rhs(0,:,RhoH),dx(0),dt1,RhoH,be_cn_theta,
     $                    rho_flag,.false.,lo(0),hi(0),bc(0,:))
            
c     extract D for RhoH
c     WILL: TODO: make sure this is right
            do i=lo(0),hi(0)
               diff_hat_1(0,i,RhoH) = (scal_1_AD(0,i,RhoH)-scal_old(0,i,RhoH))/dt1 
     $              - dRhs(0,i,0)/dt1
            enddo
            
            print *,'... react with const sources'
            
c     compute A+D source terms for reaction integration
            do n = 1,nscal
               do i=lo(0),hi(0)
                  const_src(0,i,n) = diff_hat_1(0,i,n) - diff_1(0,i,n)
     $               + quadrature_1(aofs_n(0,i,n) + diff_n(0,i,n),
     $                              aofs_1(0,i,n) + diff_1(0,i,n),
     $                              aofs_2(0,i,n) + diff_2(0,i,n), dt(0))/dt1
               enddo
            enddo
c     add differential diffusion
            do i=lo(0),hi(0)
               const_src(0,i,RhoH) = const_src(0,i,RhoH)
     $                 + quadrature_1(diffdiff_n(0,i),
     $                                diffdiff_1(0,i),
     $                                diffdiff_2(0,i), dt(0))/dt1
               do n=1,Nspec
                  const_src(0,i,FirstSpec+n-1) = const_src(0,i,FirstSpec+n-1) 
     $              - wdot_1(i,n)*mwt(n) 
     $              + quadrature_1(wdot_n(i,n), wdot_1(i,n), wdot_2(i,n), dt(0))*mwt(n)/dt1
               end do
            end do

c     solve equations (50), (51) and (52)
            call strang_chem(scal_old(0,:,:), scal_1_curr(0,:,:),
     $                       provide_wdot, const_src(0,:,:),
     $                       I_R(0,:,:),dt1,lo(0),hi(0),bc(0,:))

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
cccccccccccccccccccccccccccccccccccc
c     compute S^{n+1}
cccccccccccccccccccccccccccccccccccc

               print *,'... compute S^{n+1}'
               print *,'    old and new'

c     instantaneous omegadot for divu calc
               do i=lo(0),hi(0)
                  do n = 1,Nspec
                     Y(n) = scal_1_curr(0,i,FirstSpec+n-1)/scal_1_curr(0,i,Density)
                  enddo
                  ! compute the production rates from the current iterate
                  call CKWYR(scal_1_curr(0,i,Density), scal_1_curr(0,i,Temp), Y, iwrk, rwrk, WDOTK)
                  do n=1,Nspec
                     I_R_divu(0,i,n) = wdotk(n)*mwt(n)
                  end do
               end do
               
               call calc_diffusivities(scal_1_curr(0,:,:),beta_1(0,:,:),
     &                              mu_dummy(0,:),lo(0),hi(0))
c     divu
               call calc_divu(scal_1_curr(0,:,:),beta_1(0,:,:),I_R_divu(0,:,:),
     &                           divu_new(0,:),dx(0),lo(0),hi(0))
               
               divu_1 = divu_new
cccccccccccccccccccccccccccccccccccc
c     update delta_chi
cccccccccccccccccccccccccccccccccccc

               print *,'... updating S^{n+1} and macvel_1'
               print *,'    using fancy delta_chi'

c     compute ptherm = p(rho,T,Y)
c     this is needed for any dpdt-based correction scheme
               call compute_pthermo(scal_1_curr(0,:,:),lo(0),hi(0),bc(0,:))
               
c     delta_chi = delta_chi + (peos-p0)/(dt*peos) + (1/peos) u dot grad peos
c     completely reevaluate delta_chi_2
               delta_chi_2 = 0.d0
               call add_dpdt(scal_1_curr(0,:,:),scal_1_curr(0,:,RhoRT),
     $                       delta_chi_2(0,:),macvel_1(0,:),dx(0),dt1,
     $                       lo(0),hi(0),bc(0,:))
c     and use that quantity to increment delta_chi_1

               delta_chi_1 = delta_chi_1 + delta_chi_2
               
c     S_hat^{n+1} = S^{n+1} + delta_chi
               do i=lo(0),hi(0)
                  divu_effect(0,i) = divu_new(0,i) + delta_chi_0(0,i) + delta_chi_1(0,i)
               end do

c     macvel_1 will now satisfy div(umac) = S_hat^{n+1}
               call macproj(macvel_1(0,:),scal_1_curr(0,:,Density),
     &                      divu_effect(0,:),dx,lo(0),hi(0),bc(0,:))


               call scal_aofs(scal_1_curr(0,:,:),macvel_1(0,:),aofs_curr(0,:,:),
     $                     dx(0),dt1,lo(0),hi(0),bc(0,:))
            
c     compute conservatively corrected div gamma_m 
c     also save gamma_m for computing diffdiff terms later
              call get_spec_visc_terms(scal_1_curr(0,:,:),beta_1(0,:,:),
     &                                 diff_1_curr(0,:,FirstSpec:),
     &                                 gamma_lo(0,:,:),gamma_hi(0,:,:),
     &                                 dx(0),lo(0),hi(0))
              if (LeEQ1 .eq. 0) then
c     calculate differential diffusion "diffdiff" terms, i.e.,
c     sum_m div [ h_m (rho D_m - lambda/cp) grad Y_m ]
c     we pass in conservative gamma_m via gamma
c     we take lambda / cp from beta
c     we compute h_m using T from the first argument
c     we compute grad Y_m using Y_m from the second argument
                 call get_diffdiff_terms(scal_1_curr(0,:,:),scal_1_curr(0,:,:),
     $                                   gamma_lo(0,:,:),
     $                                   gamma_hi(0,:,:),beta_1(0,:,:),
     $                                   diffdiff_1_curr(0,:),dx(0),lo(0),hi(0))
              end if
            
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!! END FIRST SUBSTEP !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!! START SECOND SUBSTEP !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            
            do i=lo(0),hi(0)
               do n = 1,Nspec
                  Y(n) = scal_2_prev(0,i,FirstSpec+n-1)/scal_2_prev(0,i,Density)
               enddo
               ! compute the production rates from the previous iterate
               call CKWYR(scal_2_prev(0,i,Density), scal_2_prev(0,i,Temp), Y, iwrk, rwrk, wdot_2(i,:))
            end do
            
            print *,'... doing SDC iter ',misdc
            
            print *,'... compute lagged diff_new, D(U^{n+1,k-1})'

c     compute transport coefficients
c        rho D_m     (for species)
c        lambda / cp (for enthalpy)
c        lambda      (for temperature)
            call calc_diffusivities(scal_2_prev(0,:,:),beta_2(0,:,:),
     &                              mu_dummy(0,:),lo(0),hi(0))
c     compute a conservative div gamma_m
c     save gamma_m for differential diffusion computation
            call get_spec_visc_terms(scal_2_prev(0,:,:),beta_2(0,:,:),
     &                               diff_2(0,:,FirstSpec:),
     &                               gamma_lo(0,:,:),
     &                               gamma_hi(0,:,:),
     &                               dx(0),lo(0),hi(0))
c     compute div lambda/cp grad h (no differential diffusion)
            call get_rhoh_visc_terms(scal_2_prev(0,:,:),beta_2(0,:,:),
     &                               diff_2(0,:,RhoH),dx(0),lo(0),hi(0))

            if (LeEQ1 .eq. 0) then
c     calculate differential diffusion "diffdiff" terms, i.e.,
c     sum_m div [ h_m (rho D_m - lambda/cp) grad Y_m ]
c     we pass in conservative gamma_m via gamma
c     we take lambda / cp from beta
c     we compute h_m using T from the first argument
c     we compute grad Y_m using Y_m from the second argument
               call get_diffdiff_terms(scal_2_prev(0,:,:),scal_2_prev(0,:,:),
     $                                 gamma_lo(0,:,:),
     $                                 gamma_hi(0,:,:),beta_2(0,:,:),
     $                                 diffdiff_2(0,:),dx(0),lo(0),hi(0))
            end if

cccccccccccccccccccccccccccccccccccc
c     compute S^{n+1}
cccccccccccccccccccccccccccccccccccc

               print *,'... compute S^{n+1}'
               print *,'    old and new'

c     instantaneous omegadot for divu calc
               do i=lo(0),hi(0)
                  do n=1,Nspec
                     I_R_divu(0,i,n) = wdot_2(i,n)*mwt(n)
                  end do
               end do

c     divu
               call calc_divu(scal_2_prev(0,:,:),beta_2(0,:,:),I_R_divu(0,:,:),
     &                           divu_2(0,:),dx(0),lo(0),hi(0))
            
cccccccccccccccccccccccccccccccccccc
c     update delta_chi and project
cccccccccccccccccccccccccccccccccccc

            print *,'... updating S^{n+1} and macvel_2'
            print *,'    using fancy delta_chi'

c     S_hat^{n+1} = S^{n+1} + delta_chi
            do i=lo(0),hi(0)
               divu_effect(0,i) = divu_new(0,i) + delta_chi_2(0,i) + delta_chi_3(0,i)
            end do

c     macvel_2 will now satisfy div(umac) = S_hat^{n+1}
c     WILL: changed scal_old to scal_1
            call macproj(macvel_2(0,:),scal_2_prev(0,:,Density),
     &                      divu_effect(0,:),dx,lo(0),hi(0),bc(0,:))       
            
c     compute advective flux divergence
            call scal_aofs(scal_2_prev(0,:,:),macvel_2(0,:),aofs_2(0,:,:),
     $                     dx(0),dt2,lo(0),hi(0),bc(0,:))

            do i=lo(0),hi(0)
               do n=1,nscal
                  aofs_avg(0,i,n) = aofs_curr(0,i,n) - aofs_1(0,i,n)
     $                            + quadrature_2(aofs_n(0,i,n),
     $                                           aofs_1(0,i,n),
     $                                           aofs_2(0,i,n), dt)/dt2
               enddo
            enddo
            
            scal_2_AD = scal_2_prev
            print *,'... update rho'
            call update_rho(scal_1_prev(0,:,:),scal_2_AD(0,:,:),aofs_avg(0,:,:),
     &                      dt2,lo(0),hi(0),bc(0,:))
            aofs_avg = 0.5d0*(aofs_1 + aofs_2)
            
c     update rhoY_m with advection terms and set up RHS for equation (47) C-N solve
            print *,'... do correction diffusion solve for species'
            do i=lo(0),hi(0)
               do n=1,Nspec
                  is = FirstSpec + n - 1
                  dRhs(0,i,n) =   dt2*(aofs_curr(0,i,is)-aofs_1(0,i,is)-diff_2(0,i,is))
     &                 + quadrature_2(wdot_n(i,n)*mwt(n)+diff_n(0,i,is)+aofs_n(0,i,is),
     &                                wdot_1(i,n)*mwt(n)+diff_1(0,i,is)+aofs_1(0,i,is),
     &                                wdot_2(i,n)*mwt(n)+diff_2(0,i,is)+aofs_2(0,i,is),
     &                                dt(0))
               enddo
c     includes deferred correction term for enthalpy
c     differential diffusion will be added later
               dRhs(0,i,0) = dt1*(aofs_curr(0,i,RhoH)-aofs_1(0,i,RhoH)-diff_2(0,i,RhoH))
     &              + quadrature_2(diff_n(0,i,RhoH)+aofs_n(0,i,RhoH),
     &                             diff_1(0,i,RhoH)+aofs_1(0,i,RhoH),
     &                             diff_2(0,i,RhoH)+aofs_2(0,i,RhoH),
     &                             dt(0))
            enddo
            
c     WILL: TODO: make sure this should be aofs_avg
            call update_spec(scal_1_AD(0,:,:),scal_2_AD(0,:,:),
     &                       alpha(0,:),beta_1(0,:,:),
     &                       dRhs(0,0:,1:),Rhs(0,0:,FirstSpec:),
     &                       dx(0),dt2,be_cn_theta,lo(0),hi(0),bc(0,:))

c     Solve C-N system in equation (47) for \tilde{Y}_{m,AD}^{(k+1)}
            rho_flag = 2
            do n=1,Nspec
               is = FirstSpec + n - 1
               call cn_solve(scal_2_AD(0,:,:),alpha(0,:),beta_2(0,:,:),
     $                       Rhs(0,:,is),dx(0),dt2,is,be_cn_theta,
     $                       rho_flag,.false.,lo(0),hi(0),bc(0,:))
            enddo
            
            if (LeEQ1 .eq. 1) then
c     extract div gamma^n
               do i=lo(0),hi(0)
                  do n=1,Nspec
                     is = FirstSpec + n - 1
                     diff_hat_2(0,i,is) = (scal_2_AD(0,i,is)-scal_1_AD(0,i,is)
     $                                   - dRhs(0,i,n))/dt2
                  enddo
               enddo

            else

c     compute conservatively corrected div gamma_m 
c     also save gamma_m for computing diffdiff terms later
            call get_spec_visc_terms(scal_2_AD(0,:,:),beta_2(0,:,:),
     &                               diff_hat_2(0,:,FirstSpec:),
     &                               gamma_lo(0,:,:),
     &                               gamma_hi(0,:,:),
     &                               dx(0),lo(0),hi(0))

c     add differential diffusion to forcing for enthalpy solve in equation (49)
c     should this use the Gauss-Lobatto quadrature rule??
               do i=lo(0),hi(0)
                  dRhs(0,i,0) = dRhs(0,i,0) + dt2*(diffdiff_1_curr(0,i) - diffdiff_1(0,i))
     $                 + quadrature_2(diffdiff_n(0,i),
     $                                diffdiff_1(0,i),
     $                                diffdiff_2(0,i), dt(0))
               end do

            end if
            
            print *,'... do correction diffusion solve for rhoh'

c     update rhoh with advection terms and set up RHS for equation (49) C-N solve
            call update_rhoh(scal_1_AD(0,:,:),scal_2_AD(0,:,:),
     &                       alpha(0,:),beta_1(0,:,:),
     &                       dRhs(0,:,0),Rhs(0,:,RhoH),dx(0),dt2,
     &                       be_cn_theta,lo(0),hi(0),bc(0,:))

c     Solve C-N system in equation (49) for h_{AD}^{(k+1)}
            rho_flag = 2
            call cn_solve(scal_2_AD(0,:,:),alpha(0,:),beta_2(0,:,:),
     $                    Rhs(0,:,RhoH),dx(0),dt2,RhoH,be_cn_theta,
     $                    rho_flag,.false.,lo(0),hi(0),bc(0,:))
            
c     extract D for RhoH
c     WILL: TODO: make sure this is right
            do i=lo(0),hi(0)
               diff_hat_2(0,i,RhoH) = (scal_2_AD(0,i,RhoH)-scal_1_AD(0,i,RhoH)
     $                               - dRhs(0,i,0))/dt2
            enddo
            
            print *,'... react with const sources'
            
c     compute A+D source terms for reaction integration
            do n = 1,nscal
               do i=lo(0),hi(0)
                  const_src(0,i,n) = aofs_curr(0,i,n)  - aofs_1(0,i,n)
     $                             + diff_hat_2(0,i,n) - diff_2(0,i,n)
     $               + quadrature_2(aofs_n(0,i,n) + diff_n(0,i,n),
     $                              aofs_1(0,i,n) + diff_1(0,i,n),
     $                              aofs_2(0,i,n) + diff_2(0,i,n), dt(0))/dt2
               enddo
            enddo
c     add differential diffusion
            do i=lo(0),hi(0)
               const_src(0,i,RhoH) = const_src(0,i,RhoH)
     $                 + diffdiff_1_curr(0,i) - diffdiff_1(0,i)
     $                 + quadrature_2(diffdiff_n(0,i),
     $                                diffdiff_1(0,i),
     $                                diffdiff_2(0,i), dt(0))/dt2
c     add the reaction forcing term
               do n=1,Nspec
                  const_src(0,i,FirstSpec+n-1) = const_src(0,i,FirstSpec+n-1) 
     $              - wdot_2(i,n)*mwt(n)
     $              + quadrature_2(wdot_n(i,n), wdot_1(i,n), wdot_2(i,n), dt(0))*mwt(n)/dt2
               end do
            end do

c     solve equations (50), (51) and (52)
            call strang_chem(scal_1_curr(0,:,:), scal_2_curr(0,:,:),
     $                       provide_wdot, const_src(0,:,:),
     $                       I_R(0,:,:),dt2,lo(0),hi(0),bc(0,:))

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

cccccccccccccccccccccccccccccccccccc
c     compute S^{n+1}
cccccccccccccccccccccccccccccccccccc

               print *,'... compute S^{n+1}'
               print *,'    old and new'

c     instantaneous omegadot for divu calc
               do i=lo(0),hi(0)
                  do n = 1,Nspec
                     Y(n) = scal_2_curr(0,i,FirstSpec+n-1)/scal_2_curr(0,i,Density)
                  enddo
                  ! compute the production rates from the current iterate
                  call CKWYR(scal_2_curr(0,i,Density), scal_2_curr(0,i,Temp), Y, iwrk, rwrk, WDOTK)
                  do n=1,Nspec
                     I_R_divu(0,i,n) = wdotk(n)*mwt(n)
                  end do
               end do
               call calc_diffusivities(scal_2_curr(0,:,:),beta_2(0,:,:),
     &                              mu_dummy(0,:),lo(0),hi(0))
c     divu
               call calc_divu(scal_2_curr(0,:,:),beta_2(0,:,:),I_R_divu(0,:,:),
     &                           divu_new(0,:),dx(0),lo(0),hi(0))

cccccccccccccccccccccccccccccccccccc
c     update delta_chi and project
cccccccccccccccccccccccccccccccccccc

               print *,'... updating S^{n+1} and macvel_2'
               print *,'    using fancy delta_chi'

c     compute ptherm = p(rho,T,Y)
c     this is needed for any dpdt-based correction scheme
               call compute_pthermo(scal_2_curr(0,:,:),lo(0),hi(0),bc(0,:))
               
c     delta_chi = delta_chi + (peos-p0)/(dt*peos) + (1/peos) u dot grad peos
c     increment delta_chi_3
               call add_dpdt(scal_2_curr(0,:,:),scal_2_curr(0,:,RhoRT),
     $                       delta_chi_3(0,:),macvel_2(0,:),dx(0),dt2,
     $                       lo(0),hi(0),bc(0,:))

            scal_1_prev = scal_1_curr
            scal_2_prev = scal_2_curr

!         scal_1_prev = scal_1_curr
!         scal_new = scal_1_curr

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!! END SECOND SUBSTEP !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

C----------------------------------------------------------------
c     End MISDC iterations
C----------------------------------------------------------------

         enddo

      end if
      
      scal_new = scal_2_curr
C----------------------------------------------------------------
c     Step 3: Advance the velocity
C----------------------------------------------------------------

c     omegadot for divu_new computation is instantaneous
c     value of omegadot at t^{n+1}
      do i=lo(0),hi(0)
         do n=1,Nspec
            C(n) = scal_new(0,i,FirstSpec+n-1)*invmwt(n)
         end do
         call CKWC(scal_new(0,i,Temp),C,IWRK,RWRK,WDOTK)
         
         do n=1,Nspec
            I_R_divu(0,i,n) = WDOTK(n)*mwt(n)
         end do
      end do
      
c     compute transport coefficients
c        rho D_m     (for species)
c        lambda / cp (for enthalpy)
c        lambda      (for temperature)       
      call calc_diffusivities(scal_new(0,:,:),beta_2(0,:,:),
     &                        mu_new(0,:),lo(0),hi(0))
      
c     calculate S
      call calc_divu(scal_new(0,:,:),beta_2(0,:,:),I_R_divu(0,:,:),
     &               divu_new(0,:),dx(0),lo(0),hi(0))
      
c     calculate dSdt
      do i=lo(0),hi(0)
         dSdt(0,i) = (divu_new(0,i) - divu_old(0,i)) / dt(0)
      enddo

      print *,'... update velocities'

      vel_theta = 0.5d0

c     get velocity visc terms to use as a forcing term for advection
      call get_vel_visc_terms(vel_old(0,:),mu_old(0,:),visc(0,:),dx(0),
     $                        lo(0),hi(0))

      do i=lo(0),hi(0)
         visc(0,i) = visc(0,i)/scal_old(0,i,Density)
      enddo

c     compute velocity edge states
      call vel_edge_states(vel_old(0,:),scal_old(0,:,Density),gp(0,:),
     $                     macvel_n(0,:),veledge(0,:),dx(0),dt(0),
     $                     visc(0,:),lo(0),hi(0),bc(0,:))

c     alculate rhohalf
      do i=lo(0),hi(0)
         !rhohalf(0,i) = 0.5d0*(scal_old(0,i,Density)+scal_new(0,i,Density))
         rhohalf(0,i) = scal_1_curr(0,i,Density)
      enddo      
      
      !macvel_avg = 0.5d0*(macvel_n + macvel_2)
      macvel_avg = macvel_1
      !macvel_avg = macvel_n/6.0 + macvel_1/3.0 + macvel_2/6.0

c     update velocity and set up RHS for C-N diffusion solve
      call update_vel(vel_old(0,:),vel_new(0,:),gp(0,:),rhohalf(0,:),
     &                macvel_avg(0,:),veledge(0,:),alpha(0,:),mu_old(0,:),
     &                vel_Rhs(0,:),dx(0),dt(0),vel_theta,
     &                lo(0),hi(0),bc(0,:))

      if (is_first_initial_iter .eq. 1) then

c     during the first pressure initialization step, use an
c     explicit update for diffusion
         call get_vel_visc_terms(vel_old(0,:),mu_old(0,:),visc(0,:),
     $                           dx(0),lo(0),hi(0))
         do i=lo(0),hi(0)
            vel_new(0,i) = vel_new(0,i) + visc(0,i)*dt(0)/rhohalf(0,i)
         enddo

      else

c     crank-nicolson viscous solve
         rho_flag = 1
         call cn_solve(vel_new(0,:),alpha(0,:),mu_new(0,:),
     $                 vel_Rhs(0,:),dx(0),dt(0),1,vel_theta,rho_flag,
     $                 .true.,lo(0),hi(0),bc(0,:))

      endif

c     compute ptherm = p(rho,T,Y)
c     this is needed for any dpdt-based correction scheme
      call compute_pthermo(scal_new(0,:,:),lo(0),hi(0),bc(0,:))

c     S_hat^{n+1} = S^{n+1} + dpdt_factor*(ptherm-p0)/(gamma*dt*p0)
c                           + dpdt_factor*(u dot grad p)/(gamma*p0)
      call add_dpdt_nodal(scal_new(0,:,:),scal_new(0,:,RhoRT),
     &                    divu_new(0,:),vel_new(0,:),dx(0),dt(0),
     &                    lo(0),hi(0),bc(0,:))

c     project cell-centered velocities
      print *,'...nodal projection...'
      call project_level(vel_new(0,:),rhohalf(0,:),divu_new(0,:),
     &                   press_old(0,:),press_new(0,:),dx(0),dt(0),
     &                   lo(0),hi(0),bc(0,:))

      end
