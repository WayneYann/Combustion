      subroutine strang_chem(scal_old,scal_new,
     $                       const_src,lin_src_old,lin_src_new,
     $                       I_R,dt,lo,hi,bc)
      implicit none
      include 'spec.h'
      real*8     scal_old(-2:nfine+1,nscal)
      real*8     scal_new(-2:nfine+1,nscal)
      real*8    const_src( 0:nfine-1,nscal)
      real*8  lin_src_old( 0:nfine-1,nscal)
      real*8  lin_src_new( 0:nfine-1,nscal)
      real*8          I_R(-1:nfine  ,0:Nspec)
      real*8  dt
      integer lo,hi,bc(2)
      
      integer i,is,n,ifail
      real*8 RYold(Nspec), RYnew(Nspec), Told, Tnew
      real*8 rho_old
      integer FuncCount, do_diag
      real*8 diag(Nreac),hmix
      real*8 Y(Nspec)
      
      integer NiterMAX, Niter
      parameter (NiterMAX = 30)
      double precision res(NiterMAX), errMAX

c     Shut off diagnostics
      do_diag = 0

      print *,'... chemistry'
c     Evolve chem over grid
      do i=lo,hi

         rho_old = 0.d0
         do n = 1,Nspec
            RYold(n) = scal_old(i,FirstSpec+n-1)
            rho_old = rho_old + RYold(n)
         enddo
         
         Told = scal_old(i,RhoH)

c     Set linear source terms in common for ode integrators access
         do n = 1,Nspec
            is = FirstSpec + n - 1
            c_0(n) = const_src(i,is) + lin_src_old(i,is)
            c_1(n) = (lin_src_new(i,is) - lin_src_old(i,is))/dt
         enddo
         c_0(0) = const_src(i,RhoH) + lin_src_old(i,RhoH)
         c_1(0) = (lin_src_new(i,RhoH) - lin_src_old(i,RhoH))/dt
         rhoh_INIT = scal_old(i,RhoH)
         T_INIT = scal_old(i,Temp)

         call chemsolve(RYnew, Tnew, RYold, Told, FuncCount, dt,
     &                  diag, do_diag, ifail, i)
         if (ifail.ne.0) then
            print *,'solve failed, i=',i
            stop
         endif

         scal_new(i,Density) = 0.d0
         do n = 1,Nspec
            scal_new(i,Density) = scal_new(i,Density) + RYnew(n)
         enddo
         do n = 1,Nspec
            scal_new(i,FirstSpec+n-1) = RYnew(n)
            Y(n) = RYnew(n)/scal_new(i,Density)
         enddo

         scal_new(i,RhoH) = Tnew
         hmix = scal_new(i,RhoH) / scal_new(i,Density)
         errMax = hmix_TYP*1.e-20

         call FORT_TfromHYpt(scal_new(i,Temp),hmix,Y,
     &        Nspec,errMax,NiterMAX,res,Niter)

         if (Niter.lt.0) then
            print *,'strang_chem: H to T solve failed'
            print *,'Niter=',Niter
            stop
         endif

         do n = 1,Nspec
            is = FirstSpec + n - 1
            I_R(i,n) =
     $           (scal_new(i,is)-scal_old(i,is)) / dt
     $           - const_src(i,is)
     $           - 0.5d0*(lin_src_old(i,is)+lin_src_new(i,is))
         enddo

         call set_bc_s(scal_new,lo,hi,bc)

      enddo

      end
