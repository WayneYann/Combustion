
      block data chemdat
      include 'spec.h'
      data traninit / -1 /
      data tranfile / 'tran.asc.chem-H' /
      data TMIN_TRANS / 0.d0 /
      data Pr / 0.7d0 /
      data Sc / 0.7d0 /
      data LeEQ1 / 1 /
      data thickFacTR / 1.d0 /
      data thickFacCH / 1.d0 /
      data max_vode_subcycles / 15000 /
      data min_vode_timestep / 1.e-14 /
      data Pcgs / -1 /
      data iH2  / -1 /
      data iO2  / -1 /
      data iN2  / -1 /
      data iCH4 / -1 /
      data dvd_debug / 0 /
      end



      subroutine calc_diffusivities(nx, scal, beta, mu, dx, time)
      implicit none
      include 'spec.h'
      integer nx
      double precision scal(-1:nx,*)
      double precision beta(-1:nx,*)
      double precision mu(-1:nx)
      double precision time, dx

      double precision Dt(maxspec), CPMS(maxspec), Y(maxspec)
      double precision Tt, Wavg, rho
      double precision X(maxspec), alpha, l1, l2, cpmix, RWRK
      integer n, i, IWRK

      double precision cpi(9)

c     Ensure chem/tran initialized
      if (traninit.lt.0) call initchem()

      call set_bc_s(nx,scal,dx,time)
      if (LeEQ1 .eq. 0) then
         
         do i=-1, nx         
            Tt = MAX(scal(i,Temp),TMIN_TRANS) 
            do n=1,Nspec
               Y(n) = scal(i,FirstSpec+n-1) / scal(i,Density)
            enddo
            
            CALL CKMMWY(Y,IWRK,RWRK,Wavg)
            CALL CKCPMS(Tt,IWRK,RWRK,CPMS)
            CALL CKYTX(Y,IWRK,RWRK,X)
            CALL EGSPAR(Tt,X,Y,CPMS,EGRWRK,EGIWRK)
            CALL EGSV1(Pcgs,Tt,Y,Wavg,EGRWRK,Dt)
c            CALL CKRHOY(Pcgs,Tt,Y,IWRK,RWRK,RHO)

            do n=1,Nspec
               beta(i,FirstSpec+n-1)
     &              = scal(i,Density) * Wavg * invmwt(n) * Dt(n)
            end do
         
            alpha = 1.0D0
            CALL EGSL1(alpha, Tt, X, EGRWRK, l1)
            alpha = -1.0D0
            CALL EGSL1(alpha, Tt, X, EGRWRK, l2)
            beta(i,Temp) = .5 * (l1 + l2)
            CALL CKCPBS(scal(i,Temp),Y,IWRK,RWRK,CPMIX)
            beta(i,RhoH) = beta(i,Temp) / CPMIX
            
            CALL EGSE3(Tt, Y, EGRWRK, mu(i))            
         enddo
         
      else

         do i=-1, nx
c     Kanuary, Combustion Phenomena (Wiley, New York) 1982:  mu [g/(cm.s)] = 10 mu[kg/(m.s)]
            mu(i) = 10.d0 * 1.85e-5*(MAX(scal(i,Temp),1.d0)/298.0)**.7
c     For Le=1, rho.D = lambda/cp = mu/Pr
            rho = 0.d0
            do n=1,Nspec
               beta(i,FirstSpec+n-1) = mu(i) / Sc
               rho = rho + scal(i,FirstSpec+n-1)
            enddo
            
            do n=1,Nspec
               Y(n) = scal(i,FirstSpec+n-1) / rho
            enddo
            CALL CKCPBS(scal(i,Temp),Y,IWRK,RWRK,CPMIX)
            beta(i,RhoH) = mu(i) / Pr
            beta(i,Temp) = beta(i,RhoH) * CPMIX
         enddo
      endif

      if (thickFacTR.ne.1.d0) then
         do i=-1, nx
            do n=1,Nspec
               beta(i,FirstSpec+n-1) = beta(i,FirstSpec+n-1)*thickFacTR
            end do
            beta(i,Temp) = beta(i,Temp) * thickFacTR
            beta(i,RhoH) = beta(i,RhoH) * thickFacTR
         enddo
      endif

      end

      subroutine initchem
      implicit none
      include 'spec.h'
      double precision RU, RUC, RWRK
      integer lout, IWRK, n, len, offset,i
      integer kname(maxspnml*maxspec)

      if (traninit.lt.0) then
         open(unit=51,status='old',form='formatted',file=tranfile)
         call CNVTTR(51)
         lout = 6
         call EGINI(eg_nodes, lout, eg_IFLAG, eg_ITLS,
     &        EGRWRK, egr, EGIWRK, egi)

c     Other useful things
         call CKRP(IWRK, RWRK, RU, RUC, P1ATM)
         call CKWT(IWRK, RWRK, mwt)         
         do n=1,Nspec
            invmwt(n) = 1.d0 / mwt(n)
         end do

         call cksyms(kname,maxspnml)
         specNameLen = 0
         do n=1,Nspec
            offset = (n-1)*maxspnml+1
            call convStr(kname(offset),maxspnml,specNames(n),len)
            specNamelen = MAX(specNameLen,len)
            if (specNames(n).eq.'H2')  iH2=n
            if (specNames(n).eq.'O2')  iO2=n
            if (specNames(n).eq.'N2')  iN2=n
            if (specNames(n).eq.'CH4') iCH4=n
         enddo
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
      if (Nelt.gt.maxelts) then
         print *,'Too many elements, increase maxelts'
      endif
      if (Nspec.gt.maxspec) then
         print *,'Too many species, increase maxspec'
      endif
      if (Nreac.gt.maxreac) then
         print *,'Too many reactions, increase maxreac'
      endif

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
         if (codedString(k) .ne. ICHAR(' ')) then
            string(k:k) = CHAR(codedString(k))
            strLen = k
         endif
      enddo
      end


      subroutine vodeF_T_RhoY(N, TIME, Z, ZP, RPAR, IPAR)
      implicit none
      include 'spec.h'
      double precision TIME, Z(maxspec+1), ZP(maxspec+1), RPAR(*)
      integer N, IPAR(*)
      
      double precision RHO, CPB, SUM, Y(maxspec)
      double precision HK(maxspec), WDOTK(maxspec), C(maxspec), RWRK
      integer K, IWRK

      integer NiterMAX, Niter
      parameter (NiterMAX = 30)
      double precision res(NiterMAX), errMAX, hmix, T, DEN

      if (Pcgs.lt.0.d0) then
         print *,'vodeF_T_RhoY: Must set Pcgs before calling vode'
         stop
      endif

      RHO = 0.d0
      do K=1,Nspec
         RHO = RHO + Z(1+K)
      enddo

      do K=1,Nspec
         C(K) = Z(1+K)*invmwt(K)
         Y(K) = Z(1+K)/RHO
      enddo

      T = Z(1)
      hmix = (rhoh_INIT + c_0(1)*TIME + c_1(1)*TIME*TIME)/RHO
      errMax = ABS(hmix_TYP*1.e-6)
      call FORT_TfromHYpt(T,hmix,Y,errMax,NiterMAX,res,Niter)
      if (Niter.lt.0) then
         print *,'F: H to T solve failed in F, Niter=',Niter
         stop
      endif

      call CKHMS(T,IWRK,RWRK,HK)
      call CKCPBS(T,Y,IWRK,RWRK,CPB)
      call CKWC(T,C,IWRK,RWRK,WDOTK)
      SUM = 0.d0
      DO K = 1, Nspec
         ZP(K+1) = WDOTK(K)*mwt(K)/thickFacCH
     &               + c_0(1+K) + c_1(1+K)*TIME
         SUM = SUM + HK(K)*ZP(K+1)
      END DO
      ZP(1) = (c_0(1) + c_1(1)*TIME - SUM) / (RHO*CPB)
      END

      subroutine vodeJ(NEQ, T, Y, ML, MU, PD, NRPD, RPAR, IPAR)
      implicit none
      integer NEQ, NRPD, ML, MU, IPAR(*)
      double precision T, Y(NEQ), PD(NRPD,NEQ), RPAR(*)
      print *,'Should not be in vodeJ'
      stop
      end




      subroutine chemsolve(RYnew, Tnew, RYold, Told, FuncCount, dt,
     &     diag, do_diag, ifail)
      implicit none
      include 'spec.h'

      double precision YJ_SAVE(80)
      logical FIRST
      common /VHACK/ YJ_SAVE, FIRST
      save   /VHACK/

      integer do_diag, ifail
      double precision RYold(*), RYnew(*), Told, Tnew, FuncCount
      double precision dt, diag(*)
   
      integer NEQ, ITOL, IOPT, ITASK, open_vode_failure_file
      parameter (ITOL=1, IOPT=1, ITASK=1)
      double precision RTOL, ATOL(maxspec+1), ATOLEPS, TT1, TT2
      parameter (RTOL=1.0E-8, ATOLEPS=1.0E-8)
      external vodeF_T_RhoY, vodeJ, open_vode_failure_file
      integer n, MF, ISTATE, lout
      character*(maxspnml) name

      integer nsubchem, nsub, node
      double precision dtloc, weight, TT1save
      double precision C(maxspec),Q(maxreac), scale

      double precision dY(maxspec), Ytemp(maxspec),Yres(maxspec),sum
      double precision zp(maxspec+1)
      logical bad_soln


c     DVODE workspace requirements      
      integer dvr, dvi
      parameter (dvr = 22 + 9*(maxspec+1) + 2*(maxspec+1)**2)
      parameter (dvi = 30 + maxspec + 1)
      
      double precision DVRWRK(dvr)
      integer DVIWRK(dvi)

      double precision Z(maxspec+1)
      double precision RPAR, RWRK
      integer IPAR, IWRK

c     IOPT=1 parameter settings for VODE
      DVRWRK(4) = 0
      DVRWRK(5) = 0
      DVRWRK(7) = min_vode_timestep
      DVIWRK(5) = 0
      DVIWRK(6) = max_vode_subcycles
      DVIWRK(7) = 0

      if (do_diag.eq.1) nsubchem = nchemdiag
      
      MF = 22
      ATOL(1) = ATOLEPS
      TT1 = 0.d0
      TT2 = dt
      if (do_diag.eq.1) then
         nsub = nsubchem
         dtloc = dt/nsubchem
      else
         nsub = 1
         dtloc = dt
      endif
      ISTATE = 1
      NEQ = Nspec + 1


      Z(1) = Told
      do n=1,Nspec
         Z(1+n) = RYold(n)
      end do

c     Always form Jacobian to start
      FIRST = .TRUE.

      if (do_diag.eq.1) then
         FuncCount = 0
         do n=1,Nspec
            C(n) = Z(n+1)*invmwt(n)
         enddo
         CALL CKQC(Z(1),C,IWRK,RWRK,Q)
         do n=1,Nreac
            diag(n) = diag(n) + 0.5*dtloc*Q(n)
         enddo
      endif

      ifail = 0
      do node = 1,nsub
         if (node.lt.nsub) then
            weight = 1.d0
         else
            weight = 0.5d0
         endif

         TT1save = TT1
         TT2 = TT1 + dtloc

         CALL DVODE
     &        (vodeF_T_RhoY, NEQ, Z, TT1, TT2, ITOL, RTOL, ATOL,
     &        ITASK, ISTATE, IOPT, DVRWRK, dvr, DVIWRK,
     &        dvi, vodeJ, MF, RPAR, IPAR)

         TT1 = TT2

         if (do_diag.eq.1) then
            do n=1,Nspec
               C(n) = Z(n+1)*invmwt(n)
            enddo
            CALL CKQC(Z(1),C,IWRK,RWRK,Q)
            do n=1,Nreac
               diag(n) = diag(n) + weight*dtloc*Q(n)
            enddo
            FuncCount = FuncCount + DVIWRK(11)
         else
            FuncCount = DVIWRK(11)
         endif

         if (ISTATE.LE.-1  .or.  verbose_vode .eq. 1) then
            write(6,*) '......dvode done:'
            write(6,*) ' last successful step size = ',DVRWRK(11)
            write(6,*) '          next step to try = ',DVRWRK(12)
            write(6,*) '   integrated time reached = ',DVRWRK(13)
            write(6,*) '      number of time steps = ',DVIWRK(11)
            write(6,*) '              number of fs = ',DVIWRK(12)
            write(6,*) '              number of Js = ',DVIWRK(13)
            write(6,*) '    method order last used = ',DVIWRK(14)
            write(6,*) '   method order to be used = ',DVIWRK(15)
            write(6,*) '            number of LUDs = ',DVIWRK(19)
            write(6,*) ' number of Newton iterations ',DVIWRK(20)
            write(6,*) ' number of Newton failures = ',DVIWRK(21)
            write(6,*) '     comp with largest err = ',DVIWRK(16)
         end if

         Tnew = Z(1)
         do n=1,Nspec
            RYnew(n) = Z(n+1)
         end do
         
         if (ISTATE.LE.-1) ifail = 1

      enddo
      end


      subroutine FORT_TfromHYpt(T,Hin,Y,errMax,NiterMAX,res,Niter)
      implicit none
      include 'spec.h'
      double precision T,Y(*),H,Hin
      double precision TMIN,TMAX,errMAX
      integer NiterMAX,Niter,n,NiterDAMP
      parameter (TMIN=250, TMAX=5000)
      double precision  T0,cp,cv,dH,Tt,RoverWbar,Wbar,RU,RUC
      double precision res(0:NiterMAX-1),dT, Htarg
      logical out_of_bounds, converged, soln_bad, stalled
      double precision h300,cp300,h6500,cp6500
      integer ihitlo,ihithi,j,IWRK
      double precision RWRK

      out_of_bounds(Tt) = (Tt.lt.TMIN) .or. (Tt.gt.TMAX)

      NiterDAMP = NiterMAX
      if ((T.GE.TMIN).and.(T.LE.TMAX)) then
         T0 = T
      else
         T0 = 0.5*(TMIN+TMAX)
         T = T0
      end if
      Niter = 0
      dH = 0.d0
      soln_bad = .FALSE.
c     Hin in MKS, convert to CGS
c      Htarg = Hin * 1.d4
      Htarg = Hin
      ihitlo = 0
      ihithi = 0

      CALL CKHBMS(T,Y,IWRK,RWRK,H)
      dH = 2*ABS(H - Htarg)/(1.d0 + ABS(H) + ABS(Htarg))

      res(Niter) = dH
      converged = dH.le.errMAX

      do while ((.not.converged) .and. (.not.soln_bad))

         CALL CKCPBS(T,Y,IWRK,RWRK,cp)
         dT = (Htarg - H)/cp
         if ((Niter.le.NiterDAMP).and.(T+dT.ge.TMAX)) then
            T = TMAX
            ihithi = 1
         else if ((Niter.le.NiterDAMP).and.(T+dT.le.TMIN)) then
            T = TMIN
            ihitlo = 1
         else
            T = T + dT
         end if
         soln_bad = out_of_bounds(T)
         if (soln_bad) then
            Niter = -1
            goto 100
         else
            CALL CKHBMS(T,Y,IWRK,RWRK,H)
            dH = 2*ABS(H - Htarg)/(1.d0 + ABS(H) + ABS(Htarg))
            res(Niter) = dH
            Niter = Niter + 1
         end if
         if (Niter .ge. NiterMAX) then
            Niter = -2
            goto 100
         endif
         converged = (ABS(dH).le.errMAX) .or. (ABS(dT).le.errMAX)

         if ((ihitlo.eq.1).and.(H.gt.Htarg)) then
            T = TMIN
            CALL CKHBMS(T,Y,IWRK,RWRK,h300)
            CALL CKCPBS(T,Y,IWRK,RWRK,cp300)
            T=TMIN+(Htarg-h300)/cp300
            converged = .true.
         endif
         if ((ihithi.eq.1).and.(H.lt.Htarg)) then
            T = TMAX
            CALL CKHBMS(T,Y,IWRK,RWRK,h6500)
            CALL CKCPBS(T,Y,IWRK,RWRK,cp6500)
            T=TMAX+(Htarg-h6500)/cp6500
            converged = .true.
         endif
      end do

      return
c
c     Error condition....dump state and bail out
c
 100  continue

      write(6,997) 'T from (H,Y): failed'
      write(6,997) 'iterations tried = ',Niter
      write(6,998) 'initial T = ',T0
      write(6,998) 'current T = ',T
      write(6,998) 'species mass fracs:'
      do n = 1,Nspec
         write(6,998) '  ',Y(n)
      end do
      write(6,998)
      write(6,998) 'residual:'
      do n = 0,NiterMAX-1
         write(6,998) '  ',res(n)
      end do

 997  format(a,3(i4,a))
 998  format(a,d21.12)
      end
  
      integer function open_vode_failure_file ()
      implicit none
      character*30 name, myproc
      integer lout,i,j,k,idx

c     Hardwire the unit number to 26 for the moment
      lout = 26 
c      call bl_pd_myproc(i)
      i = 0
      write(myproc, *) i
      idx = 1 
      do j = 1, 30
         if (myproc(j:j) .ne. ' ') then
            idx = j
            goto 1
         end if 
      end do
 1    continue
      do k = 30, j+1, -1
         if (myproc(k:k) .ne. ' ') then
            goto 2
         end if
      end do
 2    continue
      write(name, '(2a)') 'vode.failed.', myproc(idx:k)
c      write(name, '(2a)') 'vode.failed.', myproc(idx:30)
      open(unit=lout, file=name, form='formatted', status='replace')
      open_vode_failure_file = lout
      end


      subroutine get_hmix_given_T_RhoY(nx,scal,dx)
      implicit none
      include 'spec.h'
      integer nx
      real*8 scal(-1:nx,*)
      real*8 dx
      
      integer i,n,IWRK
      real*8 Y(maxspec), rho, RWRK, hmix

      do i = 0,nx-1
         rho = 0.d0
         do n=1,Nspec
            rho = rho + scal(i,FirstSpec+n-1)
         enddo
         do n=1,Nspec
            Y(n) = scal(i,FirstSpec+n-1)/rho
         enddo
         call CKHBMS(scal(i,Temp),Y,IWRK,RWRK,hmix)
         scal(i,RhoH) = hmix * rho
      enddo
      end
