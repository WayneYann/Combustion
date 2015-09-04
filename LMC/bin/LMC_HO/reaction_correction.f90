      subroutine reaction_correction(scal_mp1_avg, scal_m_avg, &
                                     scal_k_avg, misdc_k, advection_kp1, &
                                     advection_k, diffusion_kp1, diffusion_k, &
                                     wdot_k, wdot_kp1, I_k_avg, I_k_cc, dtm)
         use wchem_module
         use ghost_cells_module
         use cell_conversions_module
         implicit none
         include 'spec.h'
         double precision, intent(inout) ::  scal_mp1_avg(-2:nx+1,nscal)
         double precision, intent(in   ) ::    scal_m_avg(-2:nx+1,nscal)
         double precision, intent(in   ) ::    scal_k_avg(-2:nx+1,nscal)
         integer,          intent(in   ) ::       misdc_k
         double precision, intent(in   ) :: advection_kp1( 0:nx-1,nscal)
         double precision, intent(in   ) ::   advection_k( 0:nx-1,nscal)
         double precision, intent(in   ) :: diffusion_kp1( 0:nx-1,nscal)
         double precision, intent(in   ) ::   diffusion_k( 0:nx-1,nscal)
         double precision, intent(in   ) ::        wdot_k( 0:nx-1,Nspec)
         double precision, intent(out  ) ::      wdot_kp1( 0:nx-1,Nspec)
         double precision, intent(in   ) ::       I_k_avg( 0:nx-1,nscal)
         double precision, intent(in   ) ::        I_k_cc( 0:nx-1,nscal)
         double precision, intent(in   ) ::           dtm
         
         
         double precision ::   wdot_k_avg(0:nx-1,Nspec)
         double precision :: wdot_kp1_avg(0:nx-1,Nspec)
         
         double precision :: rho_cc(-2:nx+1)
         double precision ::   h_cc(-2:nx+1)
         
         double precision :: avg_term(0:nx-1)
         double precision ::  cc_term(0:nx-1, Nspec+1)
         
         double precision :: rhs(Nspec)
         double precision :: guess(Nspec)
         double precision :: rhohguess
         
         double precision :: rhoYold(Nspec)
         double precision :: rhohold
         
         integer :: i, is, n
         
         !integer :: iwrk
         !double precision :: rwrk
         !double precision :: hmix
         !double precision :: T
         double precision :: Y(Nspec)
         
         
         !integer NiterMax, Niter, ifail
         !parameter (NiterMax = 30)
         !double precision res(NiterMax), errMax
         ! chemsolve stuff
         integer :: FuncCount, do_diag, ifail
         double precision :: diag(Nreac)
         
         ! turn off diagnostics
         do_diag = 0
         
         ! we perform the reaction corrections at cell-centers
         ! but advection and diffusion are passed in as cell-averaged 
         ! quantities. we therefore convert them to cell-averages
         do n = 1,Nspec
            is = FirstSpec+n-1
            do i=0,nx-1
               avg_term(i) = scal_m_avg(i, is) + dtm*(advection_kp1(i,is) - advection_k(i,is) &
                          + diffusion_kp1(i,is) - diffusion_k(i,is))
            end do
            call extrapolate_avg_to_cc(cc_term(:,n), avg_term)
            cc_term(:,n) = cc_term(:,n) - dtm*wdot_k(:, n) + I_k_cc(:,is)
         end do
         
         do i=0,nx-1
            avg_term(i) = scal_m_avg(i, RhoH) + dtm*(advection_kp1(i,RhoH) - advection_k(i,RhoH) &
                       + diffusion_kp1(i,RhoH) - diffusion_k(i,RhoH))
         end do
         call extrapolate_avg_to_cc(cc_term(:,Nspec+1), avg_term)
         cc_term(:,Nspec+1) = cc_term(:,Nspec+1) + I_k_cc(:, RhoH)
         
         call avg_to_cc(rho_cc, scal_mp1_avg(:,Density), rho_bc(on_lo))
         call avg_to_cc(h_cc, scal_mp1_avg(:,RhoH), h_bc(on_lo)*rho_bc(on_lo))
         do i=-2,nx+1
            h_cc(i) = h_cc(i)/rho_cc(i)
         end do
         
         do i=0,nx-1
            ! set up the parameters used for the VODE and BE solve
            do n = 1,Nspec
               is = FirstSpec+n-1
               
               rhs(n) = cc_term(i,n)/rho_cc(i)
               
               c_0(n) = rhs(n)
               c_1(n) = 0.d0
               
               rhoYold(n) = scal_m_avg(i, is)
            end do
            
            c_0(0) = cc_term(i,Nspec+1)
            c_1(0) = 0.d0
            
            rhohold = scal_m_avg(i, RhoH)
            
            rhoh_INIT = scal_m_avg(i, RhoH)
            T_INIT    = scal_m_avg(i, Temp)
            
            if(misdc_k .eq. 1) then
               ! call VODE to solve the ODE, to get a guess for the BE Newton solve
               call chemsolve(guess, rhohguess, rhoYold, rhohold, FuncCount, &
                              dtm, diag, do_diag, ifail, i)
            
               ! use the result from VODE as the intial guess for the Newton solve
               do n=1,Nspec
                  guess(n) = guess(n)/rho_cc(i)
               end do
            else
               do n=1,Nspec
                  guess(n) = scal_k_avg(i, FirstSpec+n-1)/scal_k_avg(i, Density)
               end do
            end if
            
            ! call the nonlinear backward Euler solver
            call bechem(Y, guess, rho_cc(i), h_cc(i), rhs, dtm)
            
            do n=1,Nspec
               wdot_kp1(i,n) = rho_cc(i)*(Y(n) - rhs(n))/dtm
            end do
         enddo
         
         do n=1,Nspec
            call extrapolate_cc_to_avg(wdot_k_avg(:,n),   wdot_k(:,n))
            call extrapolate_cc_to_avg(wdot_kp1_avg(:,n), wdot_kp1(:,n))
         end do
         
         do i=0,nx-1
            do n=1,Nspec
               is = FirstSpec+n-1
               scal_mp1_avg(i,is) = scal_m_avg(i,is) &
                  + dtm*(advection_kp1(i,is) - advection_k(i,is) &
                       + diffusion_kp1(i,is) - diffusion_k(i,is) &
                       +  wdot_kp1_avg(i,n ) -  wdot_k_avg(i,n)) + I_k_avg(i,is)
            end do
         end do
         
         call fill_scal_avg_ghost_cells(scal_mp1_avg)
      end subroutine reaction_correction