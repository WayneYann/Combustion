cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     NOTE: The comments and equation references correspond to the final
c           published version available online at:
c
c     http://www.tandfonline.com/doi/full/10.1080/13647830.2012.701019
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine advance(vel_old,vel_new,scal_old,scal_new,
     $                   I_R,press_old,press_new,
     $                   divu_old,divu_new,beta_old,beta_new,
     $                   beta_for_Y_old,beta_for_Y_new,
     $                   beta_for_Wbar_old,beta_for_Wbar_new,
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
      real*8   scal_new(0:nlevs-1,-2:nfine+1,nscal)
      real*8   scal_old(0:nlevs-1,-2:nfine+1,nscal)

c     cell-centered, 1 ghost cell
      real*8        I_R(0:nlevs-1,-1:nfine  ,0:Nspec)
      real*8   beta_old(0:nlevs-1,-1:nfine  ,nscal)
      real*8   beta_new(0:nlevs-1,-1:nfine  ,nscal)
      real*8   beta_for_Y_old(0:nlevs-1,-1:nfine  ,nscal)
      real*8   beta_for_Y_new(0:nlevs-1,-1:nfine  ,nscal)
      real*8   beta_for_Wbar_old(0:nlevs-1,-1:nfine  ,nscal)
      real*8   beta_for_Wbar_new(0:nlevs-1,-1:nfine  ,nscal)
      real*8   divu_old(0:nlevs-1,-1:nfine)
      real*8   divu_new(0:nlevs-1,-1:nfine)

c     cell-centered, no ghost cells
      real*8  delta_chi(0:nlevs-1,0:nfine-1)
      real*8     deltaT(0:nlevs-1,0:nfine-1)

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
      real*8           gp(0:nlevs-1,-1:nfine)
      real*8         visc(0:nlevs-1,-1:nfine)
      real*8  I_R_instant(0:nlevs-1,-1:nfine,  0:Nspec)
      real*8     diff_old(0:nlevs-1,-1:nfine,  nscal)
      real*8     diff_new(0:nlevs-1,-1:nfine,  nscal)
      real*8     diff_hat(0:nlevs-1,-1:nfine,  nscal)
      real*8     diff_tmp(0:nlevs-1,-1:nfine,  nscal)
      real*8 diffdiff_old(0:nlevs-1,-1:nfine)
      real*8 diffdiff_new(0:nlevs-1,-1:nfine)
      real*8 diffdiff_hat(0:nlevs-1,-1:nfine)
      real*8  divu_effect(0:nlevs-1,-1:nfine)

c     cell-centered, no ghost cells
      real*8       rhohalf(0:nlevs-1, 0:nfine-1)
      real*8         alpha(0:nlevs-1, 0:nfine-1)
      real*8        rho_cp(0:nlevs-1, 0:nfine-1)
      real*8       vel_Rhs(0:nlevs-1, 0:nfine-1)
      real*8      aofs_old(0:nlevs-1, 0:nfine-1,nscal)
      real*8      aofs_new(0:nlevs-1, 0:nfine-1,nscal)
      real*8      gamma_lo(0:nlevs-1, 0:nfine-1,Nspec)
      real*8      gamma_hi(0:nlevs-1, 0:nfine-1,Nspec)
      real*8 gamma_Wbar_lo(0:nlevs-1, 0:nfine-1,Nspec)
      real*8 gamma_Wbar_hi(0:nlevs-1, 0:nfine-1,Nspec)
      real*8     const_src(0:nlevs-1, 0:nfine-1,nscal)
      real*8   lin_src_old(0:nlevs-1, 0:nfine-1,nscal)
      real*8   lin_src_new(0:nlevs-1, 0:nfine-1,nscal)
      real*8           Rhs(0:nlevs-1, 0:nfine-1,nscal)
      real*8          dRhs(0:nlevs-1, 0:nfine-1,0:Nspec)

c     nodal, no ghost cells
      real*8       macvel(0:nlevs-1, 0:nfine  )
      real*8      veledge(0:nlevs-1, 0:nfine  )

      real*8 Y(Nspec),WDOTK(Nspec),C(Nspec),RWRK
      
      integer i,is,misdc,n,rho_flag,IWRK,l,j

      real*8   scal_tmp(0:nlevs-1,-2:nfine+1,nscal)
      real*8 norm(Nspec),deltaTsum

      print *,'advance: at start of time step',istep

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
     $                     macvel(0,:),dx(0),dt(0),lo(0),hi(0),bc(0,:))

c     reset delta_chi
      delta_chi = 0.d0

c     compute ptherm = p(rho,T,Y)
c     this is needed for any dpdt-based correction scheme
      call compute_pthermo(scal_old(0,:,:),lo(0),hi(0),bc(0,:))

c     delta_chi = delta_chi + dpdt_factor*(peos-p0)/(dt*peos)
      call add_dpdt(scal_old(0,:,:),scal_old(0,:,RhoRT),
     $              delta_chi(0,:),macvel(0,:),dx(0),dt(0),
     $              lo(0),hi(0),bc(0,:))

c     S_hat^n = S^n + delta_chi
      do i=lo(0),hi(0)
         divu_effect(0,i) = divu_old(0,i) + delta_chi(0,i)
      end do

c     mac projection
c     macvel will now satisfy div(umac) = S_hat
      call macproj(macvel(0,:),scal_old(0,:,Density),
     &             divu_effect(0,:),dx,lo(0),hi(0),bc(0,:))

c     compute A^n
      print *,'... creating the advective terms with old data'
      call scal_aofs(scal_old(0,:,:),macvel(0,:),aofs_old(0,:,:),
     $               divu_effect(0,:),dx(0),dt(0),
     $               lo(0),hi(0),bc(0,:))

ccccccccccccccccccccccccccccccccccccccccccc
c     Step 2: Advance thermodynamic variables
ccccccccccccccccccccccccccccccccccccccccccc

c     compute transport coefficients at t^n
c        rho D_m     (for species)
c        lambda / cp (for enthalpy)
c        lambda      (for temperature)
      call calc_diffusivities(scal_old(0,:,:),beta_old(0,:,:),
     &                        beta_for_Y_old(0,:,:),
     &                        beta_for_Wbar_old(0,:,:),
     &                        mu_old(0,:),lo(0),hi(0))

c     compute diffusion terms at t^n
      print *,'... creating the diffusive terms with old data'

c     compute div lambda grad T
      diff_old(0,:,Temp) = 0.d0
      call addDivLambdaGradT(scal_old(0,:,:),beta_old(0,:,:),
     &                       diff_old(0,:,Temp),dx(0),lo(0),hi(0))

c     compute conservatively corrected div gamma_m 
c     also save Gamma_m for computing diffdiff = div h_m Gamma_m later
      call get_spec_visc_terms(scal_old(0,:,:),beta_old(0,:,:),
     &                         diff_old(0,:,FirstSpec:),
     &                         gamma_lo(0,:,:),gamma_hi(0,:,:),
     &                         dx(0),lo(0),hi(0))

      if (LeEQ1 .eq. 0) then
c     compute div h_m Gamma_m
c     we pass in conservative Gamma_m via gamma
c     we compute h_m using T from the scalar argument
         call get_diffdiff_terms(scal_old(0,:,:),
     $                           gamma_lo(0,:,:),gamma_hi(0,:,:),
     $                           diffdiff_old(0,:),dx(0),lo(0),hi(0))
      else
         diffdiff_old = 0.d0
         diffdiff_new = 0.d0
      end if

c     If istep > 1, I_R is instantaneous value at t^n
c     Otherwise,    I_R is I_R^kmax from previous pressure iteration
      if (istep .gt. 1) then
         do i=lo(0),hi(0)
            do n=1,Nspec
               C(n) = scal_old(0,i,FirstSpec+n-1)*invmwt(n)
            end do
            call CKWC(scal_old(0,i,Temp),C,IWRK,RWRK,WDOTK)
            do n=1,Nspec
               I_R(0,i,n) = WDOTK(n)*mwt(n)
            end do
         end do
      end if

c     non-fancy predictor that simply sets scal_new = scal_old
      scal_new = scal_old
      divu_new = divu_old
      beta_new = beta_old
      beta_for_Y_new = beta_for_Y_old
      beta_for_Wbar_new = beta_for_Wbar_old
      diff_new = diff_old
      aofs_new = aofs_old
      diffdiff_new = diffdiff_old

C----------------------------------------------------------------
c     Begin MISDC iterations
C----------------------------------------------------------------

      do misdc = 1, misdc_iterMAX

         print *,'... doing SDC iter ',misdc

         if (misdc .gt. 1) then

            print *,'... compute lagged diff_new, D^{n+1,(k-1)}'

c     compute transport coefficients
c        rho D_m     (for species)
c        lambda / cp (for enthalpy)
c        lambda      (for temperature)
            call calc_diffusivities(scal_new(0,:,:),beta_new(0,:,:),
     &                              beta_for_Y_new(0,:,:),
     &                              beta_for_Wbar_new(0,:,:),
     &                              mu_new(0,:),lo(0),hi(0))

c     compute div lambda grad T
            diff_new(0,:,Temp) = 0.d0
            call addDivLambdaGradT(scal_new(0,:,:),beta_new(0,:,:),
     &                             diff_new(0,:,Temp),dx(0),lo(0),hi(0))

c     compute a conservative div gamma_m
c     save gamma_m for differential diffusion computation
            call get_spec_visc_terms(scal_new(0,:,:),beta_new(0,:,:),
     &                               diff_new(0,:,FirstSpec:),
     &                               gamma_lo(0,:,:),gamma_hi(0,:,:),
     &                               dx(0),lo(0),hi(0))

            if (LeEQ1 .eq. 0) then
c     compute div h_m Gamma_m
c     we pass in conservative gamma_m via gamma
c     we compute h_m using T from the scalar argument
               call get_diffdiff_terms(scal_new(0,:,:),
     $                                 gamma_lo(0,:,:),gamma_hi(0,:,:),
     $                                 diffdiff_new(0,:),dx(0),lo(0),hi(0))

            end if

cccccccccccccccccccccccccccccccccccc
c     re-compute S^{n+1/2} by averaging old and new
cccccccccccccccccccccccccccccccccccc

            print *,'... recompute S^{n+1/2} by averaging'
            print *,'    old and new'

c     instantaneous omegadot for divu calc
            do i=lo(0),hi(0)
               do n=1,Nspec
                  C(n) = scal_new(0,i,FirstSpec+n-1)*invmwt(n)
               end do
               call CKWC(scal_new(0,i,Temp),C,IWRK,RWRK,WDOTK)
               do n=1,Nspec
                  I_R_instant(0,i,n) = WDOTK(n)*mwt(n)
               end do
            end do

c     divu
            call calc_divu(scal_new(0,:,:),beta_new(0,:,:),I_R_instant(0,:,:),
     &                     divu_new(0,:),dx(0),lo(0),hi(0))

cccccccccccccccccccccccccccccccccccc
c     update delta_chi and project
cccccccccccccccccccccccccccccccccccc

c     compute ptherm = p(rho,T,Y)
c     this is needed for any dpdt-based correction scheme
            call compute_pthermo(scal_new(0,:,:),lo(0),hi(0),bc(0,:))
               
c     delta_chi = delta_chi + dpdt_factor*(peos-p0)/(dt*peos)
            call add_dpdt(scal_new(0,:,:),scal_new(0,:,RhoRT),
     $                    delta_chi(0,:),macvel(0,:),dx(0),dt(0),
     $                    lo(0),hi(0),bc(0,:))

c     S_hat^{n+1,(k)} = S^{n+1,(k)} + delta_chi
            do i=lo(0),hi(0)
               divu_new(0,i) = divu_new(0,i) + delta_chi(0,i)
            end do

c     macvel will now satisfy div(umac) = S_hat^{n+1,(k)}
            call macproj(macvel(0,:),scal_new(0,:,Density),
     &                   divu_new(0,:),dx,lo(0),hi(0),bc(0,:))

            print *,'... creating the advective terms with new data'

c     compute A^{n+1,(k)}
            call scal_aofs(scal_new(0,:,:),macvel(0,:),aofs_new(0,:,:),
     $                     divu_new(0,:),dx(0),dt(0),
     $                     lo(0),hi(0),bc(0,:))

         end if

cccccccccccccccccccccccccccccccccccc
c     update delta_chi and project
cccccccccccccccccccccccccccccccccccc

         print *,'... updating S^{n+1/2} and macvel'
         print *,'    using fancy delta_chi'

c     update density
         print *,'... update rho'
         call update_rho(scal_old(0,:,:),scal_new(0,:,:),aofs_old(0,:,:),aofs_new(0,:,:),
     &                   dt(0),lo(0),hi(0),bc(0,:))

c     this is used as the alpha coefficient for species and velocity solver
         do i=lo(0),hi(0)
            alpha(0,i) = scal_new(0,i,Density)
         end do

c     compute deferred correcion terms
         do i=lo(0),hi(0)
            do n=1,Nspec
               is = FirstSpec + n - 1
c     includes deferred correction term for species
c     dRhs for species now holds dt*(I_R + (1/2) div Gamma_m^n + (1/2) div Gamma_m^{(k)} )
               dRhs(0,i,n) = dt(0)*(I_R(0,i,n) 
     &              + 0.5d0*(diff_old(0,i,is) - diff_new(0,i,is)))
            enddo
c     includes deferred correction term for alternate enthalpy formulation
c     this is the lambda grad T part and the h_m Gamma_m part
c     dRhs for enthalpy now holds :
c        (dt/2) div (lambda^n grad T^n - lambda^(k) grad T^(k))
c       +(dt/2) div (h_m^n gamma_m^n - h_m^(k) gamma_m^(k))
            dRhs(0,i,0) = dt(0)*(
     &           + 0.5d0*(diff_old(0,i,Temp) - diff_new(0,i,Temp))
     &           + 0.5d0*(diffdiff_old(0,i) - diffdiff_new(0,i)))
         enddo

c     new iterative coupled species/enthalpy diffusion algorithm
         do l=1,Wbar_iter

            print*,'Wbar iter',l

c     compute div ( beta_for_Wbar^{(k)} grad Wbar_{AD}^{(k+1),l} )
c     also need to save the fluxes themselves for constructing Gamma_m later
            call get_spec_visc_terms_Wbar(scal_new(0,:,:),beta_for_Wbar_new(0,:,:),
     &                                    diff_tmp(0,:,FirstSpec:),
     &                                    gamma_Wbar_lo(0,:,:),
     &                                    gamma_Wbar_hi(0,:,:),
     &                                    dx(0),lo(0),hi(0))

c     construct Rhs for implicit system
            do i=lo(0),hi(0)
               do n=1,Nspec
                  is = FirstSpec + n - 1
                  Rhs(0,i,is) = scal_old(0,i,is) 
     &                 + 0.5d0*dt(0)*(aofs_old(0,i,is) +aofs_new(0,i,is))
     &                 + dRhs(0,i,n) + dt(0)*diff_tmp(0,i,is)
               end do
            end do

c     Solve implicit system
            rho_flag = 2
            do n=1,Nspec
               is = FirstSpec + n - 1
               call cn_solve(scal_new(0,:,:),alpha(0,:),beta_for_Y_new(0,:,:),
     $                       Rhs(0,:,is),dx(0),dt(0),is,1.d0,
     $                       rho_flag,.false.,lo(0),hi(0),bc(0,:))
            enddo

            call set_bc_s(scal_new(0,:,:),lo(0),hi(0),bc(0,:))

c     compute conservatively corrected version of div gamma_m
c     where gamma_m =   beta_for_Y^{(k)} grad \tilde Y_{m,AD}^{(k+1),l+1} 
c                     + beta_for_Wbar^{(k)} grad Wbar_{AD}^{(k+1),l}
c     the latter term is already available from the get_spec_visc_terms_Wbar call above
c     we save the total fluxes for calculating diffdiff terms for the enthalpy solve
            call get_spec_visc_terms_Y_and_Wbar(scal_new(0,:,:),
     &                                          beta_for_Y_new(0,:,:),
     &                                          diff_hat(0,:,FirstSpec:),
     &                                          gamma_Wbar_lo(0,:,:),
     &                                          gamma_Wbar_hi(0,:,:),
     &                                          gamma_lo(0,:,:),
     &                                          gamma_hi(0,:,:),
     &                                          dx(0),lo(0),hi(0))

c     compute rho^{(k+1)}*Y_{m,AD}^{(k+1),l+1}
            do i=lo(0),hi(0)
               do n=1,Nspec
                  is = FirstSpec + n -1
                  scal_new(0,i,is) = scal_old(0,i,is) 
     &                 + 0.5d0*dt(0)*(aofs_old(0,i,is)+aofs_new(0,i,is))
     &                 + dRhs(0,i,n) + dt(0)*diff_hat(0,i,is)
               end do
            end do

c     diagnostic stuff
            if (l .gt. 1) then
               norm = 0.d0
               do i=lo(0),hi(0)
                  do n=1,Nspec
                     is = FirstSpec + n - 1
                     norm(n) = norm(n) + abs(scal_new(0,i,is)-scal_tmp(0,i,is))
                  end do
               end do               
               print*,'change in rhoY relative to previous iter'
               write(*,1000) (norm(1:Nspec))
 1000          format (1000E11.3)
            end if
            scal_tmp = scal_new

            call set_bc_s(scal_new(0,:,:),lo(0),hi(0),bc(0,:))

         end do

c     set Rhs(RhoH) to (rhoh)^n + dt*A + 
c        (dt/2) div (lambda^n grad T^n - lambda^(k) grad T^(k))
c       +(dt/2) div (h_m^n gamma_m^n - h_m^(k) gamma_m^(k)
         do i=lo(0),hi(0)
            Rhs(0,i,RhoH) = scal_old(0,i,RhoH) 
     &           + 0.5d0*dt(0)*(aofs_old(0,i,RhoH)+aofs_new(0,i,RhoH)) 
     &           + dRhs(0,i,0)
         end do

         if (LeEQ1 .eq. 0) then

c     compute div h_m^{(k)} Gamma_{m,AD}^{(k+1)}
c     we pass in conservative gamma_m via gamma
c     we compute h_m using T from the scalar argument
            call get_diffdiff_terms(scal_new(0,:,:),
     $                              gamma_lo(0,:,:),gamma_hi(0,:,:),
     $                              diffdiff_hat(0,:),dx(0),lo(0),hi(0))

         end if

         do j=1,deltaT_iter

            print*,'deltaT iter',j

c     compute rho^{(k+1)}*cp_{AD}^{(k+1),l}
            do i=lo(0),hi(0)
               do n = 1,Nspec
                  Y(n) = scal_new(0,i,FirstSpec+n-1) / scal_new(0,i,Density)
               enddo
               call CKCPBS(scal_new(0,i,Temp),Y,IWRK,RWRK,rho_cp(0,i))
               rho_cp(0,i) = rho_cp(0,i)*scal_new(0,i,Density)
            end do

c     compute div lambda^{(k)} grad T_{AD}^{(k+1),l}
            diff_hat(:,:,Temp) = 0.d0
            call addDivLambdaGradT(scal_new(0,:,:),beta_new(0,:,:),
     $                             diff_hat(0,:,Temp),dx(0),lo(0),hi(0))

c     build rhs for delta T solve and store it in Rhs(Temp)
c     Rhs(RhoH) already holds (rhoh)^n + dt*A
c       + (dt/2) div (lambda^n grad T^n - lambda^(k) grad T^(k))
c       + (dt/2) div (h_m^n gamma_m^n - h_m^(k) gamma_m^(k))
c     make a copy of Rhs(RhoH)
            Rhs(:,:,Temp) = Rhs(:,:,RhoH)

c     need to subtract rho^(k+1) h_AD^{(k+1),l} from Rhs(Temp)
            do i=lo(0),hi(0)
               Rhs(0,i,Temp) = Rhs(0,i,Temp) - scal_new(0,i,RhoH)
            end do

c     need to add dt*div lambda^{(k)} grad T_AD^{(k+1),l} to Rhs(Temp)
            do i=lo(0),hi(0)
               Rhs(0,i,Temp) = Rhs(0,i,Temp) + dt(0)*diff_hat(0,i,Temp)
            end do

c     add dt*div h_m^{(k)} Gamma_{m,AD}^{(k+1)}
            do i=lo(0),hi(0)
               Rhs(0,i,Temp) = Rhs(0,i,Temp) + dt(0)*diffdiff_hat(0,i)
            end do

c     Solve C-N system for delta T
            deltaT = 0.d0
            call cn_solve_deltaT(deltaT(0,:),rho_cp(0,:),
     $                           beta_new(0,:,Temp),
     $                           Rhs(0,:,Temp),dx(0),dt(0),
     $                           1.d0,lo(0),hi(0),bc(0,:))

c     diagnostic stuff
            deltaTsum = 0.d0
            do i=lo(0),hi(0)
               deltaTsum = deltaTsum + abs(deltaT(0,i))
            end do
            print*,'deltaTsum',deltaTsum

            do i=lo(0),hi(0)
c     update temperature and use EOS to get enthalpy
               scal_new(0,i,Temp) = scal_new(0,i,Temp) + deltaT(0,i)
               do n = 1,Nspec
                  Y(n) = scal_new(0,i,FirstSpec+n-1) / scal_new(0,i,Density)
               enddo
               call CKHBMS(scal_new(0,i,Temp),Y,IWRK,RWRK,scal_new(0,i,RhoH))
               scal_new(0,i,RhoH) = scal_new(0,i,RhoH) * scal_new(0,i,Density)
            end do

c     update enthalpy and use EOS to get temperature
c            do i=lo(0),hi(0)
c               scal_new(0,i,RhoH) = scal_new(0,i,RhoH) + rho_cp(0,i)*deltaT(0,i)
c            end do
c            call rhoh_to_temp(scal_new(0,:,:),lo(0),hi(0))

            call set_bc_s(scal_new(0,:,:),lo(0),hi(0),bc(0,:))

c     end loop over m
         end do

c     dRhs for was holding :
c        (dt/2) div (lambda^n grad T^n - lambda^(k) grad T^(k))
c       +(dt/2) div (h_m^n gamma_m^n - h_m^(k) gamma_m^(k))

         print *,'... react with const sources'

c     compute A+D source terms for reaction integration
c     do this in alternate enthalpy formulation for diff term
         do n = FirstSpec,LastSpec
            do i=lo(0),hi(0)
               const_src(0,i,n) = diff_hat(0,i,n) - diff_new(0,i,n)
               lin_src_old(0,i,n) = aofs_old(0,i,n) + diff_old(0,i,n)
               lin_src_new(0,i,n) = aofs_new(0,i,n) + diff_new(0,i,n)
            enddo
         enddo
         do i=lo(0),hi(0)
            const_src(0,i,RhoH) = diff_hat(0,i,Temp) - diff_new(0,i,Temp)
     &           + diffdiff_hat(0,i) - diffdiff_new(0,i)
            lin_src_old(0,i,RhoH) = aofs_old(0,i,RhoH) + diff_old(0,i,Temp) + diffdiff_old(0,i)
            lin_src_new(0,i,RhoH) = aofs_new(0,i,RhoH) + diff_new(0,i,Temp) + diffdiff_new(0,i)
         enddo
         
c     solve equations (50), (51) and (52)
         call strang_chem(scal_old(0,:,:),scal_new(0,:,:),
     $                    const_src(0,:,:),lin_src_old(0,:,:),
     $                    lin_src_new(0,:,:),
     $                    I_R(0,:,:),dt(0),lo(0),hi(0),bc(0,:))
            
C----------------------------------------------------------------
c     End MISDC iterations
C----------------------------------------------------------------

      enddo

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
            I_R_instant(0,i,n) = WDOTK(n)*mwt(n)
         end do
      end do

c     compute transport coefficients
c        rho D_m     (for species)
c        lambda / cp (for enthalpy)
c        lambda      (for temperature)       
      call calc_diffusivities(scal_new(0,:,:),beta_new(0,:,:),
     &                        beta_for_Y_new(0,:,:),
     &                        beta_for_Wbar_new(0,:,:),
     &                        mu_new(0,:),lo(0),hi(0))

c     calculate S
      call calc_divu(scal_new(0,:,:),beta_new(0,:,:),I_R_instant(0,:,:),
     &               divu_new(0,:),dx(0),lo(0),hi(0))

      print *,'... update velocities'

c     get velocity visc terms to use as a forcing term for advection
      call get_vel_visc_terms(vel_old(0,:),mu_old(0,:),visc(0,:),dx(0),
     $                        lo(0),hi(0))

      do i=lo(0),hi(0)
         visc(0,i) = visc(0,i)/scal_old(0,i,Density)
      enddo

c     compute velocity edge states
      call vel_edge_states(vel_old(0,:),scal_old(0,:,Density),gp(0,:),
     $                     macvel(0,:),veledge(0,:),dx(0),dt(0),
     $                     visc(0,:),lo(0),hi(0),bc(0,:))

c     calculate rhohalf
      do i=lo(0),hi(0)
         rhohalf(0,i) = 0.5d0*(scal_old(0,i,Density)+scal_new(0,i,Density))
      enddo      

c     update velocity and set up RHS for C-N diffusion solve
      call update_vel(vel_old(0,:),vel_new(0,:),gp(0,:),rhohalf(0,:),
     &                macvel(0,:),veledge(0,:),alpha(0,:),mu_old(0,:),
     &                vel_Rhs(0,:),dx(0),dt(0),0.5d0,
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
     $                 vel_Rhs(0,:),dx(0),dt(0),1,0.5d0,rho_flag,
     $                 .true.,lo(0),hi(0),bc(0,:))

      endif

c     compute ptherm = p(rho,T,Y)
c     this is needed for any dpdt-based correction scheme
      call compute_pthermo(scal_new(0,:,:),lo(0),hi(0),bc(0,:))

c     reset delta_chi for new-time projection
      delta_chi = 0.d0

c     S_hat^{n+1} = S^{n+1} + dpdt_factor*(ptherm-p0)/(gamma*dt*p0)
c                           + dpdt_factor*(u dot grad p)/(gamma*p0)
      call add_dpdt_nodal(scal_new(0,:,:),scal_new(0,:,RhoRT),
     &                    delta_chi(0,:),vel_new(0,:),dx(0),dt(0),
     &                    lo(0),hi(0),bc(0,:))

c     use divu_effect as a temporary holding place for divu_new + delta_chi
      do i=lo(0),hi(0)
         divu_effect(0,i) = divu_new(0,i) + delta_chi(0,i)
      end do

c     project cell-centered velocities
      print *,'...nodal projection...'
      call project_level(vel_new(0,:),rhohalf(0,:),divu_effect(0,:),
     &                   press_old(0,:),press_new(0,:),dx(0),dt(0),
     &                   lo(0),hi(0),bc(0,:))

      end
