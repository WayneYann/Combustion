      subroutine strang_chem(scal_old,scal_new,
     $                       provide_wdot,
     $                       const_src,
     $                       I_R,dt,lo,hi,bc)
     
      use wchem
     
      implicit none
      include 'spec.h'
      real*8     scal_old(-2:nfine+1,nscal)
      real*8     scal_new(-2:nfine+1,nscal)
      logical    provide_wdot
      real*8    const_src( 0:nfine-1,nscal)
      real*8          I_R(-1:nfine  ,0:Nspec)
      real*8           dt
      
      integer lo,hi,bc(2)
      integer i,is,n,j
      
      real*8 Tnew
      
      real*8 rho_new
      real*8 hmix
      real*8 Y(Nspec)
      
      real*8 YTold(Nspec+1)
      real*8 YTnew(Nspec+1)
      real*8 YTnext(Nspec+1)
      real*8 guess(Nspec+1)
      real*8 rhguess
      real*8 wnew(Nspec), wold(Nspec)
      
      real*8 cp, hnew(Nspec), hold(Nspec), cp_old
      
      real*8 Qry(Nspec)
      
      double precision RWRK
      integer IWRK
      
      integer NiterMAX, Niter, ifail
      parameter (NiterMAX = 30)
      double precision res(NiterMAX), errMAX
      
      integer FuncCount, do_diag
      real*8 diag(Nreac)
      
      do_diag = 0
      
      errMax = hmix_TYP*1.e-20
      
      
c     Evolve chem over grid
      do i=lo,hi
         rho_new = 0.d0
         do n = 1,Nspec
            is = FirstSpec+n-1
            
            rho_new = rho_new + scal_new(i,is)
            
            YTnew(n) = scal_new(i,is)
            YTold(n) = scal_old(i,is)
            
            c_0(n) = const_src(i, is)
            c_1(n) = 0.d0
         end do
         c_0(0) = const_src(i,RhoH)
         c_1(0) = 0.d0
                  
         rhoh_INIT = scal_old(i, RhoH)
         T_INIT    = scal_old(i, Temp)
         YTnew(Nspec + 1) = scal_new(i,RhoH)
         YTold(Nspec + 1) = scal_old(i,RhoH)
         
         
         do n = 1,Nspec
           Y(n) = YTold(n)/scal_old(i,Density)
         enddo
         ! compute the production rate from the previous time step
         call CKWYR(scal_old(i,Density), scal_old(i,Temp), Y, IWRK, RWRK, wold)
         
         if (.not. provide_wdot) then
            Tnew = scal_new(i,Temp)
            do n = 1,Nspec
              Y(n) = YTnew(n)/rho_new
            enddo
            ! compute the production rate from the current time step, previous k
            call CKWYR(rho_new, Tnew, Y, IWRK, RWRK, wnew)
         end if
         
         YTnext(Nspec+1) = YTold(Nspec+1) + dt*const_src(i,RhoH)
         do n = 1, Nspec
            is = FirstSpec+n-1
            
            YTnext(n) = YTold(n) + dt*const_src(i, is)! - dt*wnew(n)*mwt(n) + dt*0.5d0*mwt(n)*(wnew(n)+wold(n))
         end do
         
         ! call VODE to solve the ODE, to get a guess for the BE Newton solve
            call chemsolve(guess, rhguess, YTold, scal_old(i,RhoH), FuncCount,
     &               dt, diag, do_diag, ifail, i)
            
            !rho_new = 0.d0
            !do n = 1, Nspec
            !   rho_new = rho_new + guess(n)
            !enddo
            rho_new = scal_new(i,Density)
            
            guess(Nspec+1) = rhguess
         
         !guess = YTold
         
         ! call the nonlinear backward Euler solver
         call bechem(guess, rho_new, YTnext, dt)
         !YTnext = guess
         
         scal_new(i,Density) = 0.d0
         do n = 1,Nspec
            scal_new(i,Density) = scal_new(i,Density) + YTnext(n)
            scal_new(i,FirstSpec+n-1) = YTnext(n)
         enddo
         scal_new(i, RhoH) = YTnext(Nspec + 1)
         
         do n = 1,Nspec
            Y(n) = scal_new(i,FirstSpec+n-1)/scal_new(i,Density)
         enddo
         
         hmix = scal_new(i,RhoH) / scal_new(i,Density)
         
         
         call FORT_TfromHYpt(scal_new(i,Temp),hmix,Y,
     &                          Nspec,errMax,NiterMAX,res,Niter)
         if (Niter.lt.0) then
            print *,'strang_chem: H to T solve failed'
            print *,'Niter=',Niter
            stop
         endif
         
         call CKWYR(scal_new(i,Density), scal_new(i,Temp), Y, iwrk, rwrk, wnew)
         
         do n = 1,Nspec
            is = FirstSpec + n - 1
!            I_R(i,n) = (scal_new(i,is)-scal_old(i,is)) / dt
!     $           - const_src(i,is)
!     $           - 0.5d0*(lin_src_old(i,is)+lin_src_new(i,is))
  
            I_R(i,n) = 0.5*mwt(n)*(wold(n)+wnew(n))
         enddo
         
         call set_bc_s(scal_new,lo,hi,bc)

      enddo
      
      print *,'... done chemistry'

      end
