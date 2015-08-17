      subroutine reaction_correction(scal_mp1_cc, scal_m_cc, advection_kp1, &
                                     advection_k, diffusion_kp1, diffusion_k, &
                                     wdot_k, I_k, dtm)
      use wchem_module
      use ghost_cells_module
      use cell_conversions_module
      implicit none
      include 'spec.h'
      double precision, intent(inout) ::   scal_mp1_cc(-2:nx+1,nscal)
      double precision, intent(in   ) ::     scal_m_cc(-2:nx+1,nscal)
      double precision, intent(in   ) :: advection_kp1( 0:nx-1,nscal)
      double precision, intent(in   ) ::   advection_k( 0:nx-1,nscal)
      double precision, intent(in   ) :: diffusion_kp1( 0:nx-1,nscal)
      double precision, intent(in   ) ::   diffusion_k( 0:nx-1,nscal)
      double precision, intent(in   ) ::        wdot_k( 0:nx-1,Nspec)
      double precision, intent(in   ) ::           I_k( 0:nx-1,nscal)
      double precision, intent(in   ) ::           dtm
      
      double precision :: avg_term(0:nx-1)
      double precision ::  cc_term(0:nx-1, Nspec+1)
      
      double precision :: rhs(Nspec+1)
      double precision :: guess(Nspec+1)
      double precision :: rhohguess
      
      double precision :: rhoYold(Nspec)
      double precision :: rhohold
      
      integer :: i, is, n
      
      double precision :: hmix
      double precision :: Y(Nspec)
      
      ! chemsolve stuff
      integer NiterMax, Niter, ifail
      parameter (NiterMax = 30)
      double precision res(NiterMax), errMax
      integer :: FuncCount, do_diag
      double precision :: diag(Nreac)
      
      ! turn off diagnostics
      do_diag = 0
      errMax = hmix_Typ*1.e-20
      
      
      do n = 1,Nspec
         is = FirstSpec+n-1
         avg_term = dtm*(advection_kp1(:,is) - advection_k(:,is) &
                       + diffusion_kp1(:,is) - diffusion_k(:,is)) + I_k(:,is)
         
         call extrapolate_avg_to_cc(cc_term(:,n), avg_term)
      end do
      
      avg_term = dtm*(advection_kp1(:,RhoH) - advection_k(:,RhoH) &
                       + diffusion_kp1(:,RhoH) - diffusion_k(:,RhoH)) + I_k(:,RhoH)
      call extrapolate_avg_to_cc(cc_term(:,Nspec+1), avg_term)
      
      do i=0,nx-1
         ! need to convert advection and diffusion to cell-centered 
         ! quantities here
         ! set up the parameters used for the VODE and BE solve
         do n = 1,Nspec
            is = FirstSpec+n-1
            
            rhs(n) = scal_m_cc(i, is) + cc_term(i,n) - dtm*wdot_k(i, n)
            
            c_0(n) = rhs(n)
            c_1(n) = 0.d0
            
            rhoYold(n) = scal_m_cc(i, is)
         end do
         
         rhs(Nspec+1) = scal_m_cc(i, RhoH) + cc_term(i,Nspec+1)
         c_0(0) = rhs(Nspec+1)
         c_1(0) = 0.d0
         
         rhohold = scal_m_cc(i, RhoH)
         
         rhoh_INIT = scal_m_cc(i, RhoH)
         T_INIT    = scal_m_cc(i, Temp)
         
         ! call VODE to solve the ODE, to get a guess for the BE Newton solve
         call chemsolve(guess, rhohguess, rhoYold, rhohold, FuncCount, &
                        dtm, diag, do_diag, ifail, i)
         
         ! use the result from VODE as the intial guess for the Newton solve
         guess(Nspec+1) = rhohguess
         
         ! call the nonlinear backward Euler solver
         call bechem(guess, scal_mp1_cc(i, Density), rhs, dtm)
         
         ! set the result as the new value
         do n = 1,Nspec
            scal_mp1_cc(i,FirstSpec+n-1) = rhs(n)
         enddo
         scal_mp1_cc(i, RhoH) = rhs(Nspec + 1)
         
         ! compute the temperature from h and Y
         do n = 1,Nspec
            Y(n) = scal_mp1_cc(i,FirstSpec+n-1)/scal_mp1_cc(i,Density)
         enddo
         hmix = scal_mp1_cc(i,RhoH) / scal_mp1_cc(i,Density)
         
         scal_mp1_cc(i,Temp) = scal_m_cc(i, Temp)
         ! get the new value for the temperature
         call FORT_TfromHYpt(scal_mp1_cc(i,Temp), hmix, Y, &
                             Nspec, errMax, NiterMax, res, Niter)
         if (Niter.lt.0) then
            print *,'strang_chem: H to T solve failed'
            print *,'Niter=',Niter
            stop
         endif

      enddo
      
      call fill_scal_cc_ghost_cells(scal_mp1_cc)
      
      end subroutine reaction_correction
