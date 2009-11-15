      subroutine driver()
      implicit none
      include 'spec.h'
      double precision RWRK
      integer n, k, IWRK, kname(0:maxspnml*maxspec-1), len
      character*(maxspnml) name

      integer NiterMAX, Niter
      parameter (NiterMAX = 30)
      double precision RY(maxspec), Y(maxspec), T, Patm, rho
      double precision rhoD(maxspec), kappa, mu, hmix
      double precision res(NiterMAX), errMAX
      double precision RYnew(maxspec), Ynew(maxspec), Tnew
      integer FC, do_diag
      double precision diag(maxreac), dt, sum
      
      call initchem()

      do n=1,Nspec
         Y(n) = 0
      enddo

c     H2 at phi=0.37
      Y(iH2) = 0.0107420536657d0
      Y(iO2) = 0.23041550611d0
      Y(iN2) = 1.d0 - Y(iH2) - Y(iO2)

      Patm = 1.d0
      T = 298.d0
      T = 1500.d0

      Pcgs = Patm*P1ATM

      call CKRHOY(Pcgs,T,Y,IWRK,RWRK,rho)
      print *,'rho: ',rho

c      call calc_diffusivities(T, Y, Patm, rhoD, kappa, mu)
c      print *,'mu: ',mu
c      print *,'kappa: ',kappa
c      do n=1,Nspec
c         name = specNames(n)
c         print *,'rhoD(',name(1:specNameLen),'): ',rhoD(n)
c      enddo


      call CKHBMS(T,Y,IWRK,RWRK,hmix)
      T = T*1.3

      errMax = ABS(hmix*1.e-20)
      call FORT_TfromHYpt(T,hmix,Y,errMax,NiterMAX,res,Niter)
      if (Niter.ge.0) then
         print *,'H to T solve converged in ',Niter,' iterations'
      else
         print *,'H to T solve failed, Niter=',Niter
      endif
      print *,'T new is ',T


      print *,'pre-chem state'
      print *,'T:',T
      print *,'Y:',(Y(n),n=1,Nspec)
      dt = 1.d-5
      do n=1,Nspec
         RY(n) = Y(n) * rho
      enddo
      verbose_vode = 1
      do n=1,Nspec+1
         c_0(n) = 0.d0
         c_1(n) = 0.d0
      enddo
      hmix_INIT = hmix
      hmix_TYP = 1.d5

      do_diag=0
      call chemsolve(RYnew, Tnew, RY, T, FC, Patm, dt, diag, do_diag)

      sum = 0.d0
      do n=1,Nspec
         sum = sum + RYnew(n)
      enddo
      do n=1,Nspec
         Ynew(n) = RYnew(n)/sum
      enddo

      errMax = ABS(hmix*1.e-10)
      print *,'errMax:',errMax
      call FORT_TfromHYpt(Tnew,hmix,Ynew,errMax,NiterMAX,res,Niter)
      if (Niter.gt.0) then
         print *,'H to T solve converged in ',Niter,' iterations'
      else
         print *,'H to T solve failed in DRIVER, Niter=',Niter
      endif

      print *,'post-chem state'

      print *,'T new is ',Tnew
      print *,'Y:',(Ynew(n),n=1,Nspec)
      end
