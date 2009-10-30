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
      print *,'rho [CGS]: ',rho

      call calcDiffusivityMKS(T, Y, Patm, rhoD, kappa, mu)
      print *,'mu [MKS]: ',mu
      print *,'kappa [MKS]: ',kappa
      do n=1,Nspec
         print *,'rhoD(',specNames(n),')[MKS]: ',rhoD(n)
      enddo

      end



      block data transtat
      include 'spec.h'
      data traninit / -1 /
      data tranfile / 'tran.asc.chem-H' /
      data TMIN_TRANS / 0.d0 /
      data Pr / 0.7d0 /
      data Sc / 0.7d0 /
      data LeEQ1 / 0 /
      data thickFacTR / 1.d0 /
      end



      subroutine calcDiffusivityMKS(T, Y, Patm, rhoD, kappa, mu)
      implicit none
      include 'spec.h'
      double precision T, Y(*), Patm
      double precision rhoD(*), kappa, mu

      double precision Ptmp, Dt(maxspec), CPMS(maxspec)
      double precision RHO, Tt, Wavg
      double precision X(maxspec), alpha, l1, l2, cpmix, RWRK
      integer n, IWRK

      if (LeEQ1 .eq. 0) then
         
c     Ensure chem/tran initialized
         if (traninit.lt.0) call initeg()
         
         Ptmp = Patm * P1ATM
         Tt = MAX(T,TMIN_TRANS) 
         
         CALL CKMMWY(Y,IWRK,RWRK,Wavg)
         CALL CKCPMS(Tt,IWRK,RWRK,CPMS)
         CALL CKYTX(Y,IWRK,RWRK,X)
         CALL EGSPAR(Tt,X,Y,CPMS,EGRWRK,EGIWRK)
         CALL EGSV1(Ptmp,Tt,Y,Wavg,EGRWRK,Dt)
         CALL CKRHOY(Ptmp,Tt,Y,IWRK,RWRK,RHO)

         do n=1,Nspec
            rhoD(n) = RHO * Wavg * invmwt(n) * Dt(n) * 0.1d0
         end do
         
         alpha = 1.0D0
         CALL EGSL1(alpha, Tt, X, EGRWRK, l1)
         alpha = -1.0D0
         CALL EGSL1(alpha, Tt, X, EGRWRK, l2)
         kappa = .5 * (l1 + l2) * 1.d-5
         
         CALL EGSE3(Tt, Y, EGRWRK, mu)
         mu = mu * 1.d-1

         
         if (thickFacTR.ne.1.d0) then
            
            do n=1,Nspec
               rhoD(n) = rhoD(n) * thickFacTR
            end do
            kappa = kappa * thickFacTR
            
         endif
         
      else

         mu = 1.85e-5*(T/298.0)**.7         
         do n=1,Nspec
            rhoD(n) = mu * thickFacTR / Sc
         enddo
         
         CALL CKCPBS(T,Y,IWRK,RWRK,CPMIX)
         kappa = rhoD(1) * CPMIX * 1.d-4 * thickFacTR / Pr

      endif
      end

      subroutine initeg
      implicit none
      include 'spec.h'
      double precision RU, RUC, RWRK
      integer lout, IWRK, n
      if (traninit.lt.0) then
         open(unit=51,status='old',form='formatted',file=tranfile)
         call CNVTTR(51)
         lout = 6
         call EGINI(eg_nodes, lout, eg_IFLAG, eg_ITLS,
     &        EGRWRK, egr, EGIWRK, egi)

c     Other useful things
         call CKRP(IWRK, RWRK, RU, RUC, P1ATM)
         call CKWT(IWRK, RWRK, invmwt)         
         do n=1,Nspec
            invmwt(n) = 1.d0 / invmwt(n)
         end do

         traninit = 1
      endif
      end

      subroutine CNVTTR (LINKMC)      
      implicit none
      include 'spec.h'
      character*(maxspnml) CDUMMY
      character*16 CFMT, IFMT, LFMT, RFMT
      PARAMETER
     1     (CFMT='(8A16)', IFMT='(10I12)', LFMT='(L8)',
     2     RFMT='(1P,5E24.16)')
      integer MaxOrder, NOrd, NS, LINKEG, LINKMC, IDUMMY, K, NONS,
     &     NONSNS, N
      parameter (MaxOrder = 4, LINKEG=30)
      double precision RDUMMY, PATM
      logical LDUMMY

      double precision WT(maxspec), POLA(maxspec), DIPO(maxspec),
     &     SIGM(maxspec),
     &     EPSI(maxspec), ZROT(maxspec), COFE(maxspec*MaxOrder),
     &     COFL(maxspec*MaxOrder), COFD(maxspec*maxspec*MaxOrder)
      integer LINA(maxspec)
C-----------------------------------------------------------------------
C     CHEMKIN-III TRANSPORT LINKING FILE (Warning: CK I and II different)
C-----------------------------------------------------------------------
      READ (LINKMC,CFMT) CDUMMY
      READ (LINKMC,CFMT) CDUMMY
      READ (LINKMC,CFMT) CDUMMY
      READ (LINKMC,LFMT) LDUMMY
      READ (LINKMC,IFMT) IDUMMY, IDUMMY, NOrd, NS, IDUMMY


      if (NOrd .GT. MaxOrder) then
c         call bl_abort('Polyfit too large for EGLib!')
      end if
c
c     Make sure the transport & chemistry files are consistent.
c
      CALL CKINDX(IDUMMY,RDUMMY,Nelt,Nspec,Nreac,Nfit)

      if (NS .NE. Nspec) then
         print*, 'transport database thinks Nspec = ', NS
         print*, 'chemistry database thinks Nspec = ', Nspec
c         call bl_abort("chemistry & transport are not consistent")
      endif

      NONS = NOrd * NS
      NONSNS = NONS * NS
      READ (LINKMC,RFMT) PATM
      READ (LINKMC,RFMT) (WT(K),   K=1, NS)
      READ (LINKMC,RFMT) (EPSI(K), K=1, NS) 
      READ (LINKMC,RFMT) (SIGM(K), K=1, NS)
      READ (LINKMC,RFMT) (DIPO(K), K=1, NS)
      READ (LINKMC,RFMT) (POLA(K), K=1, NS)
      READ (LINKMC,RFMT) (ZROT(K), K=1, NS) 
      READ (LINKMC,IFMT) (LINA(K), K=1, NS)
      READ (LINKMC,RFMT) (COFL(N), N=1, NONS)
      READ (LINKMC,RFMT) (COFE(N), N=1, NONS)
      READ (LINKMC,RFMT) (COFD(N), N=1, NONSNS)
      rewind(unit=LLINKMC)
C-----------------------------------------------------------------------
C     WRITE EGLIB TRANSPORT LINKING FILE
C-----------------------------------------------------------------------
      OPEN (UNIT=LINKEG,STATUS='unknown',FORM='UNFORMATTED',
     &     FILE='Linkeg')
      WRITE (LINKEG) NS, NOrd, (WT(K), K=1, NS), (EPSI(K), K=1, NS), 
     &                (SIGM(K), K=1, NS), (DIPO(K), K=1, NS), 
     &                (POLA(K), K=1, NS), (ZROT(K), K=1, NS), 
     &                (LINA(K), K=1, NS), (COFE(N), N=1, NONS), 
     &                (COFL(N), N=1, NONS), (COFD(N), N=1, NONSNS)
      CLOSE(UNIT=LINKEG)
      end

      subroutine convStr(codedString,maxLen,string,strLen)
      implicit none
      integer codedString(*), maxLen, k, strLen
      character*(*) string
      do k=1,maxLen
         string(k:k) = ' '
         if (codedString(k).ne. ICHAR(' ')) then
            string(k:k) = CHAR(codedString(k))
         endif
      enddo
      strLen = k
      end


