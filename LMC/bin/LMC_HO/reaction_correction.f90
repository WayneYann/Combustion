      subroutine reaction_correction(scal_mp1_cc, scal_m_cc, advection_kp1, &
                                     advection_k, diffusion_kp1, diffusion_k, &
                                     I_k, dtm, lo, hi)
      implicit none
      include 'spec.h'
      double precision, intent(inout) ::  scal_mp1_cc(-2:nx+1, nscal)
      double precision, intent(in   ) ::    scal_m_cc(-2:nx+1, nscal)
      double precision, intent(in   ) :: advection_kp1(0 :nx-1,nscal)
      double precision, intent(in   ) ::   advection_k(0 :nx-1,nscal)
      double precision, intent(in   ) :: diffusion_kp1(-1:nx,  nscal)
      double precision, intent(in   ) ::   diffusion_k(-1:nx,  nscal)
      double precision, intent(in   ) ::           I_k(0:nx-1, nscal)
      double precision, intent(in   ) ::           dtm
      integer,          intnet(in   ) ::        lo, hi
      
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
      
      do i=lo,hi
         ! need to convert advection and diffusion to cell-centered 
         ! quantities here
         ! set up the parameters used for the VODE and BE solve
         do n = 1,Nspec
            is = FirstSpec+n-1
            
            rhs(n) = scal_m_cc(i, is) &
               + dtm*(advection_kp1(i, is) - advection_k(i, is)  &
                    + diffusion_kp1(i, is) - diffusion_k(i, is)  &
                                           -       wdot_k(i, n)) &
               + I_k(i, is)
            
            c_0(n) = rhs(n)
            c_1(n) = 0.d0
            
            rhoYold(n) = scal_m_cc(i, is)
         end do
         
         rhs(Nspec+1) = scal_m_cc(i, RhoH) &
            + dtm*(advection_kp1(i, RhoH) - advection_k(i, RhoH)  &
                 + diffusion_kp1(i, RhoH) - diffusion_k(i, RhoH)) &
            + I_k(i, is)
         c_0(0) = rhs(Nspec+1)
         c_1(0) = 0.d0
         
         rhohold = scal_m_cc(i, RhoH)
         
         rhoh_INIT = scal_m_cc(i, RhoH)
         T_INIT    = scal_m_cc(i, Temp)
         
         ! call VODE to solve the ODE, to get a guess for the BE Newton solve
         call chemsolve(guess, rhohguess, rhoYold, rhohold, FuncCount, &
                        dtm, diag, do_diag, ifail, i)
         
         ! use the result from VODE as the intial guess for the Newton solve
         guess(Nspec+1) = rhguess
         
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
         
         ! get the new value for the temperature
         call FORT_TfromHYpt(scal_mp1_cc(i,Temp), hmix, Y, &
                             Nspec, errMax, NiterMax, res, Niter)
         if (Niter.lt.0) then
            print *,'strang_chem: H to T solve failed'
            print *,'Niter=',Niter
            stop
         endif

      enddo
      
      end subroutine reaction_correction
