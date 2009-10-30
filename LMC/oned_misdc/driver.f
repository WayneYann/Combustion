      subroutine driver()
      implicit none
      include 'spec.h'
      double precision RWRK
      integer n, k, IWRK, kname(0:maxspnml*maxspec-1), len
      integer iH2, iO2, iN2
      character*(maxspnml) specNames(maxspec), nameTMP

      double precision Y(maxspec), T, Patm, rho, Perg
      double precision rhoD(maxspec), kappa, mu
      
      call initeg()

      call cksyms(kname,maxspnml)
      do n=1,Nspec
         call convStr(kname((n-1)*maxspnml),maxspnml,specNames(n),len)
         if (specNames(n).eq.'H2') iH2=n
         if (specNames(n).eq.'O2') iO2=n
         if (specNames(n).eq.'N2') iN2=n
      enddo

      do n=1,Nspec
         Y(n) = 0
      enddo

c     H2 at phi=0.37
      Y(iH2) = 0.0107420536657d0
      Y(iO2) = 0.23041550611d0
      Y(iN2) = 1.d0 - Y(iH2) - Y(iO2)

      Patm = 1.d0
      T = 298.d0

      Perg = Patm*P1ATM

      call CKRHOY(Perg,T,Y,IWRK,RWRK,rho)
      print *,'rho: ',rho

      call calcDiffusivity(T, Y, Patm, rhoD, kappa, mu)
      print *,'mu: ',mu
      print *,'kappa: ',kappa
      do n=1,Nspec
         print *,'rhoD(',specNames(n),'): ',rhoD(n)
      enddo

      end
