      subroutine driver()
      implicit none
      include 'spec.h'
      double precision RWRK
      integer n, k, IWRK, kname(0:maxspnml*maxspec-1), len
      character*(maxspnml) name

      integer NiterMAX, Niter
      parameter (NiterMAX = 30)
      double precision Y(maxspec), T, Patm, rho, Perg
      double precision rhoD(maxspec), kappa, mu, hmix
      double precision res(NiterMAX), errMAX
      double precision Ynew(maxspec), Tnew
      integer FC, do_diag
      double precision diag(maxreac), dt, sum
      
      call initeg()

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

      Perg = Patm*P1ATM

      call CKRHOY(Perg,T,Y,IWRK,RWRK,rho)
      print *,'rho: ',rho

      call calcDiffusivity(T, Y, Patm, rhoD, kappa, mu)
      print *,'mu: ',mu
      print *,'kappa: ',kappa
      do n=1,Nspec
         name = specNames(n)
         print *,'rhoD(',name(1:specNameLen),'): ',rhoD(n)
      enddo


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
      dt = 1.d-4
      call chemsolve(Ynew, Tnew, Y, T, FC, Patm, dt, diag, do_diag)

      print *,'post-chem state'
      print *,'T:',Tnew
      print *,'Y:',(Ynew(n),n=1,Nspec)
      sum = 0.d0
      do n=1,Nspec
         sum = sum + Ynew(n)
      enddo
      print *,'new sum:',sum
      print *,'T new is ',Tnew

      call CKHBMS(Tnew,Ynew,IWRK,RWRK,hmix)
      Tnew = Tnew*1.05

      errMax = ABS(hmix*1.e-20)
      call FORT_TfromHYpt(Tnew,hmix,Ynew,errMax,NiterMAX,res,Niter)
      if (Niter.gt.0) then
         print *,'H to T solve converged in ',Niter,' iterations'
      else
         print *,'H to T solve failed, Niter=',Niter
      endif

      end
