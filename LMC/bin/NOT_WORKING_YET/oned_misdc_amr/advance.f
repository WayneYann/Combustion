cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     NOTE: The comments and equation references correspond to the final
c           published version available online at:
c
c     http://www.tandfonline.com/doi/full/10.1080/13647830.2012.701019
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine advance(vel_old,vel_new,scal_old,scal_new,
     $                   I_R,press_old,press_new,
     $                   divu_old,divu_new,dSdt,beta_old,beta_new,
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
      real*8   divu_old(0:nlevs-1,-1:nfine)
      real*8   divu_new(0:nlevs-1,-1:nfine)

c     cell-centered, no ghost cells
      real*8       dSdt(0:nlevs-1, 0:nfine-1)
      real*8  delta_chi(0:nlevs-1, 0:nfine-1)

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
      real*8     diff_old(0:nlevs-1,-1:nfine,  nscal)
      real*8     diff_new(0:nlevs-1,-1:nfine,  nscal)
      real*8     diff_hat(0:nlevs-1,-1:nfine,  nscal)
      real*8     diff_tmp(0:nlevs-1,-1:nfine,  nscal)
      real*8       tforce(0:nlevs-1,-1:nfine,  nscal)
      real*8 diffdiff_old(0:nlevs-1,-1:nfine)
      real*8 diffdiff_new(0:nlevs-1,-1:nfine)
      real*8  divu_extrap(0:nlevs-1,-1:nfine)
      real*8  divu_effect(0:nlevs-1,-1:nfine)

c     cell-centered, no ghost cells
      real*8      rhohalf(0:nlevs-1, 0:nfine-1)
      real*8        alpha(0:nlevs-1, 0:nfine-1)
      real*8      vel_Rhs(0:nlevs-1, 0:nfine-1)
      real*8         aofs(0:nlevs-1, 0:nfine-1,nscal)
      real*8 spec_flux_lo(0:nlevs-1, 0:nfine-1,Nspec)
      real*8 spec_flux_hi(0:nlevs-1, 0:nfine-1,Nspec)
      real*8    const_src(0:nlevs-1, 0:nfine-1,nscal)
      real*8  lin_src_old(0:nlevs-1, 0:nfine-1,nscal)
      real*8  lin_src_new(0:nlevs-1, 0:nfine-1,nscal)
      real*8          Rhs(0:nlevs-1, 0:nfine-1,nscal)
      real*8         dRhs(0:nlevs-1, 0:nfine-1,0:Nspec)

c     nodal, no ghost cells
      real*8       macvel(0:nlevs-1, 0:nfine  )
      real*8      veledge(0:nlevs-1, 0:nfine  )

      real*8 Y(Nspec),WDOTK(Nspec),C(Nspec),RWRK
      real*8 cpmix,rhocp,vel_theta,be_cn_theta
      
      integer i,is,misdc,n,rho_flag,IWRK

c     "diffdiff" means "differential diffusion", which corresponds to
c     sum_m div [ h_m (rho D_m - lambda/cp) grad Y_m ]
c     in equation (3)
      diffdiff_old(0,:) = 0.d0
      diffdiff_new(0,:) = 0.d0

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
     $                     macvel(0,:),dx(0),dt(0),lo(0),hi(0),bc(0,:))

c     extrapolate S^{n-1} and S^n to get S^{n+1/2}
      do i=lo(0),hi(0)
         divu_extrap(0,i) = divu_old(0,i) + 0.5d0*dt(0)*dSdt(0,i)
      end do

c     compute ptherm = p(rho,T,Y)
c     this is needed for any dpdt-based correction scheme
      call compute_pthermo(scal_old(0,:,:),lo(0),hi(0),bc(0,:))

      if (fancy_dpdt_fix .eq. 1 .and. fancy_predictor .eq. 1) then

ccccccccccccccccccccccccccccccccccc
c     new fancy delta chi algorithm
ccccccccccccccccccccccccccccccccccc

         do i=lo(0),hi(0)

c     delta_chi = delta_chi + (ptherm-p0)/(dt*p0)
            delta_chi(0,i) = delta_chi(0,i) 
     $           + dpdt_factor*(scal_old(0,i,RhoRT)-pcgs)/(dt(0)*pcgs)

c     S_hat^{n+1/2} = S^{n+1/2} + delta_chi
            divu_effect(0,i) = divu_extrap(0,i) + delta_chi(0,i)

         end do

      else

ccccccccccccccccccccccccccccccccccc
c     original dpdt_factor algorithm
ccccccccccccccccccccccccccccccccccc

c     S_hat^{n+1/2} = S^{n+1/2} + dpdt_factor*(ptherm-p0)/(gamma*dt*p0)
c                               + dpdt_factor*(u dot grad p)/(gamma*p0)
         do i=lo(0),hi(0)
            divu_effect(0,i) = divu_extrap(0,i)
         end do
         call add_dpdt(scal_old(0,:,:),scal_old(0,:,RhoRT),
     $                 divu_effect(0,:),macvel(0,:),dx(0),dt(0),
     $                 lo(0),hi(0),bc(0,:))

      end if

c     mac projection
c     macvel will now satisfy div(umac) = S_hat^{n+1/2}
      call macproj(macvel(0,:),scal_old(0,:,Density),
     &             divu_effect(0,:),dx,lo(0),hi(0),bc(0,:))


ccccccccccccccccccccccccccccccccccccccccccc
c     Step 2: Advance thermodynamic variables
ccccccccccccccccccccccccccccccccccccccccccc

      if (use_strang) then

ccccccccccccccccccccccccccccccccccccccccccc
c     Step 2: Advance thermodynamic variables (Strang split algorithm)
ccccccccccccccccccccccccccccccccccccccccccc

         be_cn_theta = 0.5d0

ccccccccccccccccccccccccccccccccccccccccccc
c     Strang Step 2A: First reaction step
ccccccccccccccccccccccccccccccccccccccccccc

         print *,'... react for dt/2'

c     react for dt/2
         do n = 1,nscal
            do i=lo(0),hi(0)
               const_src(0,i,n) = 0.d0
               lin_src_old(0,i,n) = 0.d0
               lin_src_new(0,i,n) = 0.d0
            enddo
         enddo
         call strang_chem(scal_old(0,:,:),scal_new(0,:,:),
     $                    const_src(0,:,:),lin_src_old(0,:,:),
     $                    lin_src_new(0,:,:),
     $                    I_R(0,:,:),dt(0)/2.d0,lo(0),hi(0),bc(0,:))
         
c     set scal_old to be the results after reactions
c     note that rho and rhoh remain constant, only species changed
         do n = FirstSpec,LastSpec
            scal_old(0,:,n) = scal_new(0,:,n)
         enddo

c     call the EOS to get consistent temperature
         call rhoh_to_temp(scal_old(0,:,:),lo(0),hi(0))

ccccccccccccccccccccccccccccccccccccccccccc
c     Strang Step 2B: Advection-diffusion predictor
ccccccccccccccccccccccccccccccccccccccccccc

         print *,'... creating the diffusive terms with old data'

c     compute transport coefficients at t^{<1>}
c        rho D_m     (for species)
c        lambda / cp (for enthalpy)
c        lambda      (for temperature)
         call calc_diffusivities(scal_old(0,:,:),beta_old(0,:,:),
     &                           mu_dummy(0,:),lo(0),hi(0))

c     compute div lambda grad T + gamma_m dot grad h_m, where
c     gamma_m has been conservatively corrected
         call get_temp_visc_terms(scal_old(0,:,:),beta_old(0,:,:),
     &                            diff_old(0,:,Temp),dx(0),lo(0),hi(0))
c     compute conservatively corrected div gamma_m 
c     also save gamma_m for computing diffdiff terms later
         call get_spec_visc_terms(scal_old(0,:,:),beta_old(0,:,:),
     &                            diff_old(0,:,FirstSpec:),
     &                            spec_flux_lo(0,:,:),
     &                            spec_flux_hi(0,:,:),
     &                            dx(0),lo(0),hi(0))
c     compute div lambda/cp grad h (no differential diffusion)
         call get_rhoh_visc_terms(scal_old(0,:,:),beta_old(0,:,:),
     &                            diff_old(0,:,RhoH),dx(0),lo(0),hi(0))

         if (LeEQ1 .eq. 0) then
c     calculate differential diffusion "diffdiff" terms, i.e.,
c     sum_m div [ h_m (rho D_m - lambda/cp) grad Y_m ]
c     we pass in conservative gamma_m via spec_flux
c     we take lambda / cp from beta
c     we compute h_m using T from the first argument
c     we compute grad Y_m using Y_m from the second argument
            call get_diffdiff_terms(scal_old(0,:,:),scal_old(0,:,:),
     $                              spec_flux_lo(0,:,:),
     $                              spec_flux_hi(0,:,:),beta_old(0,:,:),
     $                              diffdiff_old(0,:),dx(0),lo(0),hi(0))
         end if
            
         print *,'... computing aofs with explicit D forcing'

         do i=lo(0),hi(0)
            do n = 1,Nspec
               is = FirstSpec + n - 1
               tforce(0,i,is) = diff_old(0,i,is)
            enddo
            tforce(0,i,RhoH) = diff_old(0,i,RhoH) + diffdiff_old(0,i)
         enddo

c     compute advective flux divergence
         call scal_aofs(scal_old(0,:,:),macvel(0,:),aofs(0,:,:),
     $                  divu_effect(0,:),tforce(0,:,:),dx(0),dt(0),
     $                  lo(0),hi(0),bc(0,:))

         print *,'... update rho'

c     update density via equation (55)
         call update_rho(scal_old(0,:,:),scal_new(0,:,:),aofs(0,:,:),
     &                   dt(0),lo(0),hi(0),bc(0,:))

c     create forcing term for temperature edge state prediction
         do i=lo(0),hi(0)
            do n = 1,Nspec
               Y(n) = scal_old(0,i,FirstSpec+n-1) / scal_old(0,i,Density)
            enddo
            call CKCPBS(scal_old(0,i,Temp),Y,IWRK,RWRK,cpmix)
            rhocp = cpmix * 
     &           (scal_old(0,i,Density) + scal_new(0,i,Density)) / 2.d0
            tforce(0,i,Temp) = diff_old(0,i,Temp)/rhocp
         end do

         print *,'... updating new temperature with lagged coeffs'

c     update T with advection term and set up RHS for equation (56) C-N solve
         call update_temp(scal_old(0,:,:),scal_new(0,:,:),
     $                    aofs(0,:,:),alpha(0,:),
     $                    beta_old(0,:,:),beta_old(0,:,:),
     $                    Rhs(0,:,Temp),dx(0),dt(0),be_cn_theta,
     $                    lo(0),hi(0),bc(0,:))

c     Solve C-N system in equation (56) for \tilde{T}_{pred}^{<2>}
         rho_flag = 1
         call cn_solve(scal_new(0,:,:),alpha(0,:),beta_old(0,:,:),
     $                 Rhs(0,:,Temp),dx(0),dt(0),Temp,be_cn_theta,
     $                 rho_flag,.false.,lo(0),hi(0),bc(0,:))

         print *,'... compute new coeffs after temperature update'

c     compute transport coefficients
c        rho D_m     (for species)
c        lambda / cp (for enthalpy)
c        lambda      (for temperature)
         call calc_diffusivities(scal_new(0,:,:),beta_new(0,:,:),
     &                           mu_dummy(0,:),lo(0),hi(0))

         print *,'... do predictor for species'

c     update rhoY_m with advection terms and set up RHS for equation (57) C-N solve
         do i=lo(0),hi(0)
            dRhs(0,i,0) = 0.0d0
            do n=1,Nspec
               dRhs(0,i,n) = 0.d0
            enddo
         enddo
         call update_spec(scal_old(0,:,:),scal_new(0,:,:),aofs(0,:,:),
     &                    alpha(0,:),beta_old(0,:,:),
     &                    dRhs(0,0:,1:),Rhs(0,0:,FirstSpec:),dx(0),dt(0),
     &                    be_cn_theta,lo(0),hi(0),bc(0,:))

c     Solve C-N system in equation (57) for \tilde{Y}_{m,pred}^{<2>}
         rho_flag = 2
         do n=1,Nspec
            is = FirstSpec + n - 1
            call cn_solve(scal_new(0,:,:),alpha(0,:),beta_new(0,:,:),
     $                    Rhs(0,:,is),dx(0),dt(0),is,be_cn_theta,
     $                    rho_flag,.false.,lo(0),hi(0),bc(0,:))
         enddo
      
         if (LeEQ1 .eq. 0) then

c     compute conservatively corrected div gamma_m 
c     also save gamma_m for computing diffdiff terms later
            call get_spec_visc_terms(scal_new(0,:,:),beta_new(0,:,:),
     &                               diff_tmp(0,:,FirstSpec:),
     &                               spec_flux_lo(0,:,:),
     &                               spec_flux_hi(0,:,:),
     &                               dx(0),lo(0),hi(0))

c     update species with conservative diffusion fluxes using equation (58)
            do i=lo(0),hi(0)
               do n=1,Nspec
                  is = FirstSpec + n - 1
                  scal_new(0,i,is) = scal_old(0,i,is) + 
     $                 dt(0)*(aofs(0,i,is)
     $                 + 0.5d0*diff_old(0,i,is) + 0.5d0*diff_tmp(0,i,is))
               end do
            end do
            call set_bc_s(scal_new(0,:,:),lo(0),hi(0),bc(0,:))
            
         end if

c     compute transport coefficients
c        rho D_m     (for species)
c        lambda / cp (for enthalpy)
c        lambda      (for temperature)
         call calc_diffusivities(scal_new(0,:,:),beta_new(0,:,:),
     &                           mu_dummy(0,:),lo(0),hi(0))

         if (LeEQ1 .eq. 0) then

c     calculate differential diffusion "diffdiff" terms, i.e.,
c     sum_m div [ h_m (rho D_m - lambda/cp) grad Y_m ]
c     we pass in conservative gamma_m via spec_flux
c     we take lambda / cp from beta
c     we compute h_m using T from the first argument
c     we compute grad Y_m using Y_m from the second argument
            call get_diffdiff_terms(scal_new(0,:,:),scal_new(0,:,:),
     $                              spec_flux_lo(0,:,:),
     $                              spec_flux_hi(0,:,:),beta_new(0,:,:),
     $                              diffdiff_new(0,:),dx(0),lo(0),hi(0))
         
            do i=lo(0),hi(0)
               dRhs(0,i,0) = dRhs(0,i,0)
     $              + 0.5d0*dt(0)*(diffdiff_old(0,i) + diffdiff_new(0,i))
            end do
            
         end if

         print *,'... do predictor for rhoh'

c     update rhoh with advection terms and set up RHS for equation (59) C-N solve
         call update_rhoh(scal_old(0,:,:),scal_new(0,:,:),aofs(0,:,:),
     &                    alpha(0,:),beta_old(0,:,:),
     &                    dRhs(0,:,0),Rhs(0,:,RhoH),dx(0),dt(0),
     &                    be_cn_theta,lo(0),hi(0),bc(0,:))

c     Solve C-N system in equation (59) for h_{pred}^{<2>}
         rho_flag = 2
         call cn_solve(scal_new(0,:,:),alpha(0,:),beta_new(0,:,:),
     $                 Rhs(0,:,RhoH),dx(0),dt(0),RhoH,be_cn_theta,
     $                 rho_flag,.false.,lo(0),hi(0),bc(0,:))

c     call the EOS to get consistent temperature
         call rhoh_to_temp(scal_new(0,:,:),lo(0),hi(0))

C----------------------------------------------------------------
C     Corrector

         print *,'... compute new coeffs after predictor'

c     compute transport coefficients
c        rho D_m     (for species)
c        lambda / cp (for enthalpy)
c        lambda      (for temperature)
         call calc_diffusivities(scal_new(0,:,:),beta_new(0,:,:),
     &                           mu_dummy(0,:),lo(0),hi(0))

         print *,'... do corrector for species'

c     update rhoY_m with advection terms and set up RHS for equation (60) C-N solve
         do i=lo(0),hi(0)
            dRhs(0,i,0) = 0.0d0
            do n=1,Nspec
               dRhs(0,i,n) = 0.d0
            enddo
         enddo
         call update_spec(scal_old(0,:,:),scal_new(0,:,:),aofs(0,:,:),
     &                    alpha(0,:),beta_old(0,:,:),
     &                    dRhs(0,0:,1:),Rhs(0,0:,FirstSpec:),dx(0),dt(0),
     &                    be_cn_theta,lo(0),hi(0),bc(0,:))

c     Solve C-N system in equation (60) for \tilde{Y}_m^{<2>}
         rho_flag = 2
         do n=1,Nspec
            is = FirstSpec + n - 1
            call cn_solve(scal_new(0,:,:),alpha(0,:),beta_new(0,:,:),
     $                    Rhs(0,:,is),dx(0),dt(0),is,be_cn_theta,
     $                    rho_flag,.false.,lo(0),hi(0),bc(0,:))
         enddo

         if (LeEQ1 .eq. 0) then

c     compute conservatively corrected div gamma_m 
c     also save gamma_m for computing diffdiff terms later
            call get_spec_visc_terms(scal_new(0,:,:),beta_new(0,:,:),
     &                               diff_tmp(0,:,FirstSpec:),
     &                               spec_flux_lo(0,:,:),
     &                               spec_flux_hi(0,:,:),
     &                               dx(0),lo(0),hi(0))

c     update species with conservative diffusion fluxes using equation (61)
            do i=lo(0),hi(0)
               do n=1,Nspec
                  is = FirstSpec + n - 1
                  scal_new(0,i,is) = scal_old(0,i,is) + 
     $                 dt(0)*(aofs(0,i,is)
     $                 + 0.5d0*diff_old(0,i,is) + 0.5d0*diff_tmp(0,i,is))
               end do
            end do
            call set_bc_s(scal_new(0,:,:),lo(0),hi(0),bc(0,:))

c     calculate differential diffusion "diffdiff" terms, i.e.,
c     sum_m div [ h_m (rho D_m - lambda/cp) grad Y_m ]
c     we pass in conservative gamma_m via spec_flux
c     we take lambda / cp from beta
c     we compute h_m using T from the first argument
c     we compute grad Y_m using Y_m from the second argument
            call get_diffdiff_terms(scal_new(0,:,:),scal_new(0,:,:),
     $                              spec_flux_lo(0,:,:),
     $                              spec_flux_hi(0,:,:),beta_new(0,:,:),
     $                              diffdiff_new(0,:),dx(0),lo(0),hi(0))

            do i=lo(0),hi(0)
               dRhs(0,i,0) = dRhs(0,i,0)
     $              + 0.5d0*dt(0)*(diffdiff_old(0,i) + diffdiff_new(0,i))
            end do
            
         end if
         
         print *,'... do corrector for rhoh'

c     update rhoh with advection terms and set up RHS for equation (62) C-N solve
         call update_rhoh(scal_old(0,:,:),scal_new(0,:,:),aofs(0,:,:),
     &                    alpha(0,:),beta_old(0,:,:),
     &                    dRhs(0,:,0),Rhs(0,:,RhoH),dx(0),dt(0),
     &                    be_cn_theta,lo(0),hi(0),bc(0,:))

c     Solve C-N system in equation (62) for h_^{<2>}
         rho_flag = 2
         call cn_solve(scal_new(0,:,:),alpha(0,:),beta_new(0,:,:),
     $                 Rhs(0,:,RhoH),dx(0),dt(0),RhoH,be_cn_theta,
     $                 rho_flag,.false.,lo(0),hi(0),bc(0,:))

c     call the EOS to get consistent temperature
         call rhoh_to_temp(scal_new(0,:,:),lo(0),hi(0))

         print *,'... react for dt/2'

         do i=lo(0),hi(0)
            do n = FirstSpec,LastSpec
               scal_old(0,i,n) = scal_new(0,i,n)
            enddo
            scal_old(0,i,Temp) = scal_new(0,i,Temp)
            scal_old(0,i,Density) = scal_new(0,i,Density)
         enddo
         call strang_chem(scal_old(0,:,:),scal_new(0,:,:),
     $                    const_src(0,:,:),lin_src_old(0,:,:),
     $                    lin_src_new(0,:,:),
     $                    I_R_temp(0,:,:),dt(0)/2.d0,
     $                    lo(0),hi(0),bc(0,:))

c     call the EOS to get consistent temperature
         call rhoh_to_temp(scal_new(0,:,:),lo(0),hi(0))
         
         I_R(0,:,:) = I_R(0,:,:) + I_R_temp(0,:,:)
         I_R(0,:,:) = I_R(0,:,:) / 2.d0

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
      call calc_diffusivities(scal_old(0,:,:),beta_old(0,:,:),
     &                        mu_old(0,:),lo(0),hi(0))

c     compute diffusion terms at time n
         print *,'... creating the diffusive terms with old data'

c     compute conservatively corrected div gamma_m 
c     also save gamma_m for computing diffdiff terms later
         call get_spec_visc_terms(scal_old(0,:,:),beta_old(0,:,:),
     &                            diff_old(0,:,FirstSpec:),
     &                            spec_flux_lo(0,:,:),spec_flux_hi(0,:,:),
     &                            dx(0),lo(0),hi(0))
c     compute div lambda/cp grad h (no differential diffusion)
         call get_rhoh_visc_terms(scal_old(0,:,:),beta_old(0,:,:),
     &                            diff_old(0,:,RhoH),dx(0),lo(0),hi(0))

         if (LeEQ1 .eq. 0) then
c     calculate differential diffusion "diffdiff" terms, i.e.,
c     sum_m div [ h_m (rho D_m - lambda/cp) grad Y_m ]
c     we pass in conservative gamma_m via spec_flux
c     we take lambda / cp from beta
c     we compute h_m using T from the first argument
c     we compute grad Y_m using Y_m from the second argument
            call get_diffdiff_terms(scal_old(0,:,:),scal_old(0,:,:),
     $                              spec_flux_lo(0,:,:),
     $                              spec_flux_hi(0,:,:),beta_old(0,:,:),
     $                              diffdiff_old(0,:),dx(0),lo(0),hi(0))
         end if

c     If .true., use I_R in predictor is instantaneous value at t^n
c     If .false., use I_R^lagged = I_R^kmax from previous time step
         if (.false.) then
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

         print *,'... computing A forcing = D^n + I_R^{n-1,kmax}'

c     compute advective forcing term
         do i=lo(0),hi(0)
            do n = 1,Nspec
               is = FirstSpec + n - 1
               tforce(0,i,is) = diff_old(0,i,is) + I_R(0,i,n)
            enddo
            tforce(0,i,RhoH) = diff_old(0,i,RhoH) + diffdiff_old(0,i)
         enddo

         if (fancy_predictor .eq. 1) then

c     compute advective flux divergence
            call scal_aofs(scal_old(0,:,:),macvel(0,:),aofs(0,:,:),
     $                     divu_effect(0,:),tforce(0,:,:),dx(0),dt(0),
     $                     lo(0),hi(0),bc(0,:))

            print *,'... update rho'

c     update density via equation (40)
            call update_rho(scal_old(0,:,:),scal_new(0,:,:),aofs(0,:,:),
     &                      dt(0),lo(0),hi(0),bc(0,:))

c     update rhoY_m with advection terms and set up RHS for equation (41) C-N solve
            do i=lo(0),hi(0)
               dRhs(0,i,0) = 0.0d0
               do n=1,Nspec
                  dRhs(0,i,n) = dt(0)*I_R(0,i,n)
               enddo
            enddo
            call update_spec(scal_old(0,:,:),scal_new(0,:,:),aofs(0,:,:),
     &                       alpha(0,:),beta_old(0,:,:),
     &                       dRhs(0,0:,1:),Rhs(0,0:,FirstSpec:),dx(0),dt(0),
     &                       be_cn_theta,lo(0),hi(0),bc(0,:))

            print *,'... do initial diffusion solve for species'

c     Solve C-N system in equation (41) for \tilde{Y}_{m,AD}^{(0)}
            rho_flag = 2
            do n=1,Nspec
               is = FirstSpec + n - 1
               call cn_solve(scal_new(0,:,:),alpha(0,:),beta_old(0,:,:),
     $                       Rhs(0,:,is),dx(0),dt(0),is,be_cn_theta,rho_flag,
     $                       .false.,lo(0),hi(0),bc(0,:))
            enddo

            if (LeEQ1 .eq. 1) then
            
c     extract div gamma^n
               do i=lo(0),hi(0)
                  do n=1,Nspec
                     is = FirstSpec + n - 1
                     diff_hat(0,i,is) = 2.d0*((scal_new(0,i,is)-scal_old(0,i,is))/dt(0)
     $                    - aofs(0,i,is) - I_R(0,i,n) - 0.5d0*diff_old(0,i,is))
                  enddo
               end do
               
            else
            
c     compute conservatively corrected div gamma_m 
c     also save gamma_m for computing diffdiff terms later
               call get_spec_visc_terms(scal_new(0,:,:),beta_old(0,:,:),
     &                                  diff_hat(0,:,FirstSpec:),
     &                                  spec_flux_lo(0,:,:),
     &                                  spec_flux_hi(0,:,:),
     &                                  dx(0),lo(0),hi(0))

c     update species with conservative diffusion fluxes using equation (42)
               do i=lo(0),hi(0)
                  do n=1,Nspec
                     is = FirstSpec + n - 1
                     scal_new(0,i,is) = scal_old(0,i,is) + 
     $                    dt(0)*(aofs(0,i,is) + I_R(0,i,n)
     $                    + 0.5d0*diff_old(0,i,is) + 0.5d0*diff_hat(0,i,is))
                  end do
               end do
               call set_bc_s(scal_new(0,:,:),lo(0),hi(0),bc(0,:))
         
c     calculate differential diffusion "diffdiff" terms, i.e.,
c     sum_m div [ h_m (rho D_m - lambda/cp) grad Y_m ]
c     we pass in conservative gamma_m via spec_flux
c     we take lambda / cp from beta
c     we compute h_m using T from the first argument
c     we compute grad Y_m using Y_m from the second argument
               call get_diffdiff_terms(scal_old(0,:,:),scal_new(0,:,:),
     $                                 spec_flux_lo(0,:,:),
     $                                 spec_flux_hi(0,:,:),beta_old(0,:,:),
     $                                 diffdiff_new(0,:),dx(0),lo(0),hi(0))

c     add differential diffusion to forcing for enthalpy solve
               do i=lo(0),hi(0)
                  dRhs(0,i,0) = dRhs(0,i,0) 
     $                 + 0.5d0*dt(0)*(diffdiff_old(0,i) + diffdiff_new(0,i))
               end do
            
            end if

c     update rhoh with advection terms and set up RHS for equation (43) C-N solve
            call update_rhoh(scal_old(0,:,:),scal_new(0,:,:),aofs(0,:,:),
     &                       alpha(0,:),beta_old(0,:,:),
     &                       dRhs(0,:,0),Rhs(0,:,RhoH),dx(0),dt(0),
     &                       be_cn_theta,lo(0),hi(0),bc(0,:))

            print *,'... do initial diffusion solve for rhoh'

c     Solve C-N system in equation (43) for h_{AD}^{(0)}
            rho_flag = 2
            call cn_solve(scal_new(0,:,:),alpha(0,:),beta_old(0,:,:),
     $                    Rhs(0,:,RhoH),dx(0),dt(0),RhoH,be_cn_theta,rho_flag,
     $                    .false.,lo(0),hi(0),bc(0,:))

c     extract D for RhoH
            do i=lo(0),hi(0)
               diff_hat(0,i,RhoH) = 2.d0*((scal_new(0,i,RhoH)-scal_old(0,i,RhoH))/dt(0) 
     $              - aofs(0,i,RhoH) - dRhs(0,i,0)/dt(0) - 0.5d0*diff_old(0,i,RhoH) )
            enddo

            print *,'... react with constant sources'

c     compute A+D source terms for reaction integration
            do n = 1,nscal
               do i=lo(0),hi(0)
                  const_src(0,i,n) = aofs(0,i,n) 
     $                 + 0.5d0*diff_hat(0,i,n) + 0.5d0*diff_old(0,i,n)
                  lin_src_old(0,i,n) = 0.d0
                  lin_src_new(0,i,n) = 0.d0
               enddo
            enddo
c     add differential diffusion
            do i=lo(0),hi(0)
               const_src(0,i,RhoH) = const_src(0,i,RhoH)
     $              + 0.5d0*(diffdiff_old(0,i)+diffdiff_new(0,i))
            end do
            
c     solve equations (44), (45) and (46)
            call strang_chem(scal_old(0,:,:),scal_new(0,:,:),
     $                       const_src(0,:,:),lin_src_old(0,:,:),
     $                       lin_src_new(0,:,:),
     $                       I_R(0,:,:),dt(0),lo(0),hi(0),bc(0,:))
            
         else

c     non-fancy predictor that simply sets scal_new = scal_old

            scal_new(0,:,:) = scal_old(0,:,:)

         end if

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
            print *,'... doing SDC iter ',misdc
            
            if (fancy_dpdt_fix .eq. 1) then

cccccccccccccccccccccccccccccccccccc
c     new fancy delta chi algorithm
cccccccccccccccccccccccccccccccccccc
               
c     compute ptherm = p(rho,T,Y)
c     this is needed for any dpdt-based correction scheme
               call compute_pthermo(scal_new(0,:,:),lo(0),hi(0),bc(0,:))
               
               do i=lo(0),hi(0)

c     delta_chi = delta_chi + (ptherm-p0)/(dt*p0)
                  delta_chi(0,i) = delta_chi(0,i) 
     $                 + dpdt_factor*(scal_new(0,i,RhoRT)-pcgs)/(dt(0)*pcgs)
                  
c     S_hat^{n+1/2} = S^{n+1/2} + delta_chi
                  divu_effect(0,i) = divu_extrap(0,i) + delta_chi(0,i)

               end do

c     mac projection
c     macvel will now satisfy div(umac) = S_hat^{n+1/2}
               call macproj(macvel(0,:),scal_old(0,:,Density),
     &                      divu_effect(0,:),dx,lo(0),hi(0),bc(0,:))
               
            end if
            
            print *,'... compute lagged diff_new, D(U^{n+1,k-1})'

c     compute transport coefficients
c        rho D_m     (for species)
c        lambda / cp (for enthalpy)
c        lambda      (for temperature)
            call calc_diffusivities(scal_new(0,:,:),beta_new(0,:,:),
     &                              mu_dummy(0,:),lo(0),hi(0))
c     compute a conservative div gamma_m
c     save gamma_m for differential diffusion computation
            call get_spec_visc_terms(scal_new(0,:,:),beta_new(0,:,:),
     &                               diff_new(0,:,FirstSpec:),
     &                               spec_flux_lo(0,:,:),
     &                               spec_flux_hi(0,:,:),
     &                               dx(0),lo(0),hi(0))
c     compute div lambda/cp grad h (no differential diffusion)
            call get_rhoh_visc_terms(scal_new(0,:,:),beta_new(0,:,:),
     &                               diff_new(0,:,RhoH),dx(0),lo(0),hi(0))

            if (LeEQ1 .eq. 0) then
c     calculate differential diffusion "diffdiff" terms, i.e.,
c     sum_m div [ h_m (rho D_m - lambda/cp) grad Y_m ]
c     we pass in conservative gamma_m via spec_flux
c     we take lambda / cp from beta
c     we compute h_m using T from the first argument
c     we compute grad Y_m using Y_m from the second argument
               call get_diffdiff_terms(scal_new(0,:,:),scal_new(0,:,:),
     $                                 spec_flux_lo(0,:,:),
     $                                 spec_flux_hi(0,:,:),beta_new(0,:,:),
     $                                 diffdiff_new(0,:),dx(0),lo(0),hi(0))
            end if

            print *,'... computing A forcing term = D^n + I_R^{k-1}'

c     note: no need to recompute advective forcing for RhoH since it
c     doesn't change from the predictor
            do i=lo(0),hi(0)
               do n = 1,Nspec
                  is = FirstSpec + n - 1
                  tforce(0,i,is) = diff_old(0,i,is) + I_R(0,i,n)
               enddo
            enddo
            
c     compute advective flux divergence
            call scal_aofs(scal_old(0,:,:),macvel(0,:),aofs(0,:,:),
     $                     divu_effect(0,:),tforce(0,:,:),dx(0),dt(0),
     $                     lo(0),hi(0),bc(0,:))


c     update rhoY_m with advection terms and set up RHS for equation (47) C-N solve
            print *,'... do correction diffusion solve for species'
            do i=lo(0),hi(0)
               do n=1,Nspec
                  is = FirstSpec + n - 1
c     includes deferred correction term for species
                  dRhs(0,i,n) = dt(0)*(I_R(0,i,n) 
     &                 + 0.5d0*(diff_old(0,i,is) - diff_new(0,i,is)))
               enddo
c     includes deferred correction term for enthalpy
c     differential diffusion will be added later
               dRhs(0,i,0) = dt(0)*(
     &              + 0.5d0*(diff_old(0,i,RhoH) - diff_new(0,i,RhoH)))
            enddo
            call update_spec(scal_old(0,:,:),scal_new(0,:,:),aofs(0,:,:),
     &                       alpha(0,:),beta_old(0,:,:),
     &                       dRhs(0,0:,1:),Rhs(0,0:,FirstSpec:),
     &                       dx(0),dt(0),be_cn_theta,lo(0),hi(0),bc(0,:))

c     Solve C-N system in equation (47) for \tilde{Y}_{m,AD}^{(k+1)}
            rho_flag = 2
            do n=1,Nspec
               is = FirstSpec + n - 1
               call cn_solve(scal_new(0,:,:),alpha(0,:),beta_new(0,:,:),
     $                       Rhs(0,:,is),dx(0),dt(0),is,be_cn_theta,
     $                       rho_flag,.false.,lo(0),hi(0),bc(0,:))
            enddo
            
            if (LeEQ1 .eq. 1) then

c     extract div gamma^n
               do i=lo(0),hi(0)
                  do n=1,Nspec
                     is = FirstSpec + n - 1
                     diff_hat(0,i,is) = (scal_new(0,i,is)-scal_old(0,i,is))/dt(0) 
     $                    - aofs(0,i,is) - dRhs(0,i,n)/dt(0)
                  enddo
               enddo

            else

c     compute conservatively corrected div gamma_m 
c     also save gamma_m for computing diffdiff terms later
               call get_spec_visc_terms(scal_new(0,:,:),beta_new(0,:,:),
     &                                  diff_hat(0,:,FirstSpec:),
     &                                  spec_flux_lo(0,:,:),
     &                                  spec_flux_hi(0,:,:),
     &                                  dx(0),lo(0),hi(0))

c     add differential diffusion to forcing for enthalpy solve in equation (49)
               do i=lo(0),hi(0)
                  dRhs(0,i,0) = dRhs(0,i,0) 
     $                 + 0.5d0*dt(0)*(diffdiff_old(0,i) + diffdiff_new(0,i))
               end do

            end if
            
            print *,'... do correction diffusion solve for rhoh'

c     update rhoh with advection terms and set up RHS for equation (49) C-N solve
            call update_rhoh(scal_old(0,:,:),scal_new(0,:,:),aofs(0,:,:),
     &                       alpha(0,:),beta_old(0,:,:),
     &                       dRhs(0,:,0),Rhs(0,:,RhoH),dx(0),dt(0),
     &                       be_cn_theta,lo(0),hi(0),bc(0,:))

c     Solve C-N system in equation (49) for h_{AD}^{(k+1)}
            rho_flag = 2
            call cn_solve(scal_new(0,:,:),alpha(0,:),beta_new(0,:,:),
     $                    Rhs(0,:,RhoH),dx(0),dt(0),RhoH,be_cn_theta,
     $                    rho_flag,.false.,lo(0),hi(0),bc(0,:))
            
c     extract D for RhoH
            do i=lo(0),hi(0)
               diff_hat(0,i,RhoH) = (scal_new(0,i,RhoH)-scal_old(0,i,RhoH))/dt(0) 
     $              - aofs(0,i,RhoH) - dRhs(0,i,0)/dt(0)
            enddo
            
            print *,'... react with const sources'

c     compute A+D source terms for reaction integration
            do n = 1,nscal
               do i=lo(0),hi(0)
                  const_src(0,i,n) = aofs(0,i,n)
     $                 + 0.5d0*(diff_old(0,i,n)+diff_new(0,i,n))
     $                 + diff_hat(0,i,n) - diff_new(0,i,n)
                  lin_src_old(0,i,n) = 0.d0
                  lin_src_new(0,i,n) = 0.d0
               enddo
            enddo
c     add differential diffusion
            do i=lo(0),hi(0)
               const_src(0,i,RhoH) = const_src(0,i,RhoH)
     $              + 0.5d0*(diffdiff_old(0,i)+diffdiff_new(0,i))
            end do

c     solve equations (50), (51) and (52)
            call strang_chem(scal_old(0,:,:),scal_new(0,:,:),
     $                       const_src(0,:,:),lin_src_old(0,:,:),
     $                       lin_src_new(0,:,:),
     $                       I_R(0,:,:),dt(0),lo(0),hi(0),bc(0,:))
            
C----------------------------------------------------------------
c     End MISDC iterations
C----------------------------------------------------------------

         enddo

      end if

C----------------------------------------------------------------
c     Step 3: Advance the velocity
C----------------------------------------------------------------

      if (use_strang) then

c     omegadot for divu_new computation is average omegadot
c     from both dt/2 strang calls
         I_R_divu = I_R

      else

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

      end if

c     compute transport coefficients
c        rho D_m     (for species)
c        lambda / cp (for enthalpy)
c        lambda      (for temperature)       
      call calc_diffusivities(scal_new(0,:,:),beta_new(0,:,:),
     &                        mu_new(0,:),lo(0),hi(0))

c     calculate S
      call calc_divu(scal_new(0,:,:),beta_new(0,:,:),I_R_divu(0,:,:),
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
     $                     macvel(0,:),veledge(0,:),dx(0),dt(0),
     $                     visc(0,:),lo(0),hi(0),bc(0,:))

c     calculate rhohalf
      do i=lo(0),hi(0)
         rhohalf(0,i) = 0.5d0*(scal_old(0,i,Density)+scal_new(0,i,Density))
      enddo      

c     update velocity and set up RHS for C-N diffusion solve
      call update_vel(vel_old(0,:),vel_new(0,:),gp(0,:),rhohalf(0,:),
     &                macvel(0,:),veledge(0,:),alpha(0,:),mu_old(0,:),
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
