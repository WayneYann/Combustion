      subroutine strang_chem(nx,scal_old,scal_new,
     $                       const_src,lin_src_old,lin_src_new,
     $                       intra,dt)
      implicit none
      include 'spec.h'
      integer nx
      real*8     scal_old(-1:nx  ,nscal)
      real*8     scal_new(-1:nx  ,nscal)
      real*8    const_src( 0:nx-1,nscal)
      real*8  lin_src_old( 0:nx-1,nscal)
      real*8  lin_src_new( 0:nx-1,nscal)
      real*8        intra( 0:nx-1,nscal)
      real*8  dt
      
      integer i,n,ifail
      integer ispec
      real*8 RYold(maxspec), RYnew(maxspec), Told, Tnew
      real*8 linSrcOLD(nscal),linSrcNEW(nscal)
      integer FuncCount, do_diag, IWRK
      real*8 diag(maxreac),Y(maxspec),hmix,RWRK
      
      integer NiterMAX, Niter
      parameter (NiterMAX = 30)
      double precision res(NiterMAX), errMAX

c     Find a good typical value for hmix
      hmix_TYP = ABS(scal_old(-1,RhoH)/scal_old(-1,Density))
      do i = 0,nx-1
         hmix_TYP=MAX(hmix_TYP,
     &                ABS(scal_old(i,RhoH)/scal_old(i,Density)))
      enddo

c     Evolve chem over grid
      do i = 0,nx-1
         do n = 1,Nspec
            RYold(n) = scal_old(i,FirstSpec+n-1)
         enddo
         Told = scal_old(i,Temp)

c     Set linear source terms in common for ode integrators access
         do n = 1,Nspec
            ispec = FirstSpec + n - 1
            c_0(n+1) = const_src(i,ispec) + lin_src_old(i,ispec)
            c_1(n+1) = (lin_src_new(i,ispec) - lin_src_old(i,ispec))/dt
         enddo
         c_0(1) = const_src(i,RhoH) + lin_src_old(i,RhoH)
         c_1(1) = (lin_src_new(i,RhoH) - lin_src_old(i,RhoH))/dt
         rhoh_INIT = scal_old(i,RhoH)

c         print *,'sc:',i
         call chemsolve(RYnew, Tnew, RYold, Told, FuncCount, dt,
     &                  diag, do_diag, ifail)
         if (ifail.ne.0) then
            print *,'solve failed, i=',i
            stop
         endif
         
         scal_new(i,Density) = 0.d0
         do n = 1,Nspec
            ispec = FirstSpec+n-1
            scal_new(i,Density) = scal_new(i,Density)+scal_new(i,ispec)
         enddo
         do n = 1,Nspec
            scal_new(i,FirstSpec+n-1)=RYnew(n)
            Y(n) = RYnew(n)/scal_new(i,Density)
         enddo

         scal_new(i,RhoH) = scal_old(i,RhoH)+
     &        + dt*const_src(i,RhoH)
     &        + 0.5d0*dt*(lin_src_old(i,RhoH)+lin_src_new(i,rhoH))
         hmix = scal_new(i,RhoH) / scal_new(i,Density)

         errMax = hmix_TYP * 1.e-10
         call FORT_TfromHYpt(Tnew,hmix,Y,errMax,NiterMAX,res,Niter)
         if (Niter.lt.0) then
            print *,'SC: H to T solve failed in F, Niter=',Niter
            print *,'hmix_TYP:',hmix_TYP
            stop
         endif
         scal_new(i,Temp) = Tnew

      enddo

c     Define change in state due to chemistry.
      do i = 0,nx-1
         do n = 1,nscal
            intra(i,n) =
     $           (scal_new(i,n)-scal_old(i,n)) / dt
     $           - const_src(i,n)
     $           - 0.5d0*(lin_src_old(i,n)+lin_src_new(i,n))
         enddo
      enddo
      
      print *,' '
      print *,' '
        
      end
