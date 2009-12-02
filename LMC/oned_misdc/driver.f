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
      integer FC, do_diag, i, NEVAL
      parameter (NEVAL=100)
      double precision diag(maxreac), dt, sum, eps
      double precision cp(0:NEVAL-1), cpp(0:NEVAL-1), Teval(0:NEVAL-1)
      
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
      call FORT_TfromHYpt(T,hmix,Y,Nspec,errMax,NiterMAX,res,Niter)
      if (Niter.ge.0) then
         print *,'H to T solve converged in ',Niter,' iterations'
      else
         print *,'H to T solve failed, Niter=',Niter
      endif
      print *,'T new is ',T



      print *,'debugging test:'

      hmix = 0.868214818423d+07;
      T=0.999984396440d+03;
      T=1.00000001d+03;
      Y(1) = 0.889979543278d-2
      Y(2) = 0.176897481907d+0
      Y(3) = 0.682771928057d-1
      Y(4) = 0.103061796579d-3
      Y(5) = 0.320304190562d-3
      Y(6) = 0.668446420264d-3
      Y(7) = 0.758079131426d-3
      Y(8) = 0.135677506994d-3
      Y(9) = 0.743939960808d+0

      errMax = ABS(hmix*1.e-20)
      call FORT_TfromHYpt(T,hmix,Y,Nspec,errMax,NiterMAX,res,Niter)
      if (Niter.ge.0) then
         print *,'H to T solve converged in ',Niter,' iterations'
      else
         print *,'H to T solve failed, Niter=',Niter
      endif
      print *,'T new is ',T

      T = 1.d3
      eps = 5.d-3
      open(unit=34,file='cp.dat',status='unknown')
      do i=0,NEVAL-1
         Teval(i) = T + 2*eps*(DBLE(i)-(NEVAL-1)/2.d0)/NEVAL
         CALL CKCPBS(Teval(i),Y,IWRK,RWRK,cp(i))
      enddo
      do i=0,NEVAL-2
         write(34,*)
c     &        0.5*(Teval(i+1)+Teval(i)),
     &        Teval(i),
     &        cp(i)
c     &        0.5*(cp(i+1)+cp(i)),
c     &        (cp(i+1)-cp(i))/(Teval(i+1)-Teval(i))
      enddo
      close(34)
      stop


      print *,'pre-chem state'
      print *,'T:',T
      print *,'Y:',(Y(n),n=1,Nspec)
      dt = 1.d-5
      do n=1,Nspec
         RY(n) = Y(n) * rho
      enddo
      verbose_vode = 1
      do n=0,Nspec
         c_0(n) = 0.d0
         c_1(n) = 0.d0
      enddo
      rhoh_INIT = hmix * rho
      hmix_TYP = 1.d5

      T=437.609349445155
      do_diag=0
      dt =  3.896686739246933d-5
      Y(1) = 5.241468896616197D-003
      Y(2) = 1.046629392866622D-010
      Y(3) = 5.510364915755311D-008
      Y(4) = 0.222367303274566D0    
      Y(5) = 1.259485524742000D-008
      Y(6) = 1.426435811729874D-002
      Y(7) = 2.781504475068043D-005
      Y(8) = 2.958869335307645D-005
      Y(9) = 0.758069398170248D0    
      call CKHBMS(T,Y,IWRK,RWRK,hmix)
      call CKRHOY(Pcgs,T,Y,IWRK,RWRK,rho)
      sum = 0.d0
      do n=1,Nspec
         sum = sum + Y(n)
         RY(n) = Y(n) * rho
      enddo
      rhoh_INIT = hmix * rho
      print *,'HACK   pre-chem state'
      print *,'T,sum:',T,sum
      print *,'Y:',(Y(n),n=1,Nspec)

      
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
      call FORT_TfromHYpt(Tnew,hmix,Ynew,
     &     Nspec,errMax,NiterMAX,res,Niter)
      if (Niter.ge.0) then
         print *,'H to T solve converged in ',Niter,' iterations'
      else
         print *,'H to T solve failed in DRIVER, Niter=',Niter
         do n=1,Nspec
            print *,res(n)
         enddo
      endif

      print *,'post-chem state'

      print *,'T new is ',Tnew
      print *,'Y:',(Ynew(n),n=1,Nspec)
      print *,'hmix:',hmix


      end
