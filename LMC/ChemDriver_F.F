#include "ChemDriver_F.H"
#include "CONSTANTS.H"

#if defined(BL_USE_FLOAT) || defined(BL_T3E) || defined(BL_CRAY)
#define three4th    0.75
#define onepoint27  1.27
#define point4      0.4
#define point67     0.67
#define point14     0.14
#define onetenthsnd 0.0001
#define ten2minus18 1.0e-18
#else
#define three4th    0.75d0
#define onepoint27  1.27d0
#define point4      0.4d0
#define point67     0.67d0
#define point14     0.14d0
#define onetenthsnd 0.0001d0
#define ten2minus18 1.0d-18
#endif
      
      subroutine FORT_SETUPCKIO (nametr,nametrlen)
      implicit none
      integer nametrlen
      integer nametr(nametrlen)
#include "cdwrk.H"
#include "conp.H"
      integer maxlen, i
      parameter (maxlen=256)
      character filet*(maxlen)
      integer isioproc
c
c     Convert encoded name to a Fortran string.
c
      if (nametrlen.GT.maxlen) then
         call bl_abort('FORT_CKINTERP: tran file name too long')
      else if (nametrlen.NE.0) then
         do i = 1, nametrlen
            filet(i:i) = char(nametr(i))
         end do
c
c        Open problem-specific transport file for EGlib.
c
         LLINKMC = 51
         open(unit=LLINKMC,status='OLD',form='FORMATTED',
     &        file=filet(1:nametrlen), err=300)
      else
         call bl_abort('FORT_CKINTERP: transport file must be provided')
      endif
      return
 300  call bl_pd_is_ioproc(isioproc)
      if (isioproc .eq. 1) then
         print *,'Transfile = ',filet(1:nametrlen)
         call bl_abort('File not found')
      endif
      end

      subroutine FORT_LINKCKDB()
      implicit none
#include "cdwrk.H"
      end

      subroutine FORT_SETTMINTRANS(TminTRANS)
      implicit none
      REAL_T TminTRANS
#include "cdwrk.H"
      TMIN_TRANS = TminTRANS
      end

      subroutine FORT_SETVERBOSEVODE()
      implicit none
#include "cdwrk.H"
      verbose_vode = 1
      end

      subroutine FORT_SETVODESUBCYC(maxcyc)
      implicit none
      integer maxcyc
#include "cdwrk.H"
      max_vode_subcycles = maxcyc
      end

      subroutine FORT_SETSPECSCALY(name, nlength)
      implicit none
#include "cdwrk.H"
      integer nlength, name(nlength), i, j, maxlen
      REAL_T val
      parameter (maxlen=256)
      character filet*(maxlen)
      character*(maxspnml) spname, spinname
c      
c     Convert encoded names to strings, and open file
c
      if (nlength.GT.maxlen) then
         call bl_abort('FORT_SETSPECSCAL: scale file name too long')
      end if
      
      do i = 1, nlength
         filet(i:i) = char(name(i))
      end do
      open(unit=51,status='OLD',form='FORMATTED',
     &     file=filet(1:nlength),err=30)
      
 10   continue 
      read(51,*,end=20) spinname, val
      do j = 1,Nspec
         call get_spec_name(spname,j)
         if (spname .eq. spinname) then
            spec_scalY(j) = ABS(val)
         end if
      end do
      goto 10
 20   close(51)
      goto 40
 30   write(6,*) 'Trouble opening file = ',filet(1:nlength)
      call bl_abort(" ")
 40   continue 
      end

      subroutine FORT_CLOSECKIO()
      implicit none
#include "cdwrk.H"
      close(LLINKMC)
      end

      subroutine FORT_INTERPCHEM()
      implicit none
#include "cdwrk.H"
#include "conp.H"
      logical error
      integer ioproc, myid, ierr, n, RTOT, lout, namlen, i
      integer FORT_GETCKSPECNAME
      parameter (ioproc = 0)
      character*(maxspnml) name
      integer coded(maxspnml)

      call BL_PD_MYPROC(myid)
      
      if (myid .eq. ioproc) then
c
c        Suck in mc linking file, spit out EG linking file.
c
         call CNVTTR (LLINKMC)
      endif
      
      call BL_PD_BARRIER()
      
      lout = 6
      call EGINI(eg_nodes, lout, eg_IFLAG, eg_ITLS,
     &     RWRK(egbr), egr, IWRK(egbi), egi)
c      
c     Set pointers in conp common blocks (used in conpF evaluations)
c
      CALL CKINDX(IWRK(ckbi),RWRK(ckbr),Nelt,Nspec,Nreac,Nfit)
      NEQ   = Nspec + 1
      NP    = dvdbr
      NWT   = NP  + 1
      NZ    = NWT + Nspec
      RTOT  = NZ  + NEQ - 1

      if (RTOT .GT. dvder) then
         write(6,*) 'Memory layout bust, dvdr not big enough'
         write(6,*) RTOT, dvder
         call bl_abort(" ")
      end if
c
c     Set molecular weights
c
      CALL CKWT(IWRK(ckbi), RWRK(ckbr), RWRK(NWT))
c
c     Find N2 in the list.
c
      iN2 = -1
      do n = 1,Nspec
         call get_spec_name(name,n)
         if (name .eq. 'N2' ) iN2 = n
      end do
      if (iN2.eq.-1)
     &     write(6,*) '.....warning: no N2 in chemistry species list'
      end

      subroutine CNVTTR (LINKMC)
      implicit none
#include "cdwrk.H"
      character*(maxspnml) CDUMMY
      character*16 CFMT, IFMT, LFMT, RFMT
      PARAMETER
     1     (CFMT='(8A16)', IFMT='(10I12)', LFMT='(L8)', RFMT='(1P,5E24.16)')
      integer MaxOrder, NOrd, NS, LINKEG, LINKMC, IDUMMY, K, NONS, NONSNS, N
      parameter (MaxOrder = 4, LINKEG=30)
      REAL_T RDUMMY, PATM
      logical LDUMMY

      REAL_T WT(maxspec), POLA(maxspec), DIPO(maxspec), SIGM(maxspec),
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
         call bl_abort('Polyfit too large for EGLib!')
      end if
c
c     Make sure the transport & chemistry files are consistent.
c
      CALL CKINDX(IWRK(ckbi),RWRK(ckbr),Nelt,Nspec,Nreac,Nfit)

      if (NS .NE. Nspec) then
         print*, 'transport database thinks Nspec = ', NS
         print*, 'chemistry database thinks Nspec = ', Nspec
         call bl_abort("chemistry & transport are not consistent")
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

      integer function FORT_GETCKMAXNAMELEN()
      implicit none
#include "cdwrk.H"
      FORT_GETCKMAXNAMELEN = maxspnml
      end

      subroutine FORT_GETCKDIMPARAMS(imaxreac, imaxspec, imaxelts,
     &                               imaxord, imaxthrdb, imaxtp, imaxsp,
     &                               imaxspnml)
      implicit none
      integer imaxreac, imaxspec, imaxelts, imaxord
      integer imaxthrdb, imaxtp, imaxsp, imaxspnml
#include "cdwrk.H"
      imaxreac = maxreac
      imaxspec = maxspec
      imaxelts = maxelts
      imaxord = 10
      imaxthrdb = maxthrdb
      imaxtp = maxtp
      imaxsp = maxsp
      imaxspnml = maxspnml
      end

      subroutine FORT_FINDLHS(reactions, Nreacs, id)
      implicit none
#include "cdwrk.H"
      integer reactions(*), Nreacs, id
      integer j, n, Ndim, Nids, KI(maxsp), NU(maxsp)
#ifdef MIKE
      Ndim = maxsp
      if ((id.le.0).or.(id.gt.Nspec)) then
         write(6,*) 'FORT_FINDLHS:  species id out of range: ',id
         call bl_abort(" ")
      end if
      Nreacs = 0
      do j=1,Nreac
         CALL CKINU(j, Ndim, IWRK(ckbi), RWRK(ckbr), Nids, KI, NU)
         do n=1,Nids
            if ((KI(n).eq.id).and.(NU(n).lt.0)) then
               Nreacs = Nreacs + 1
               reactions(Nreacs) = j
            endif
         end do
      end do
#else
      call bl_abort("FORT_FINDLHS not implemented")
#endif
      end

      subroutine FORT_FINDRHS(reactions, Nreacs, id)
      implicit none
#include "cdwrk.H"
      integer reactions(*), Nreacs, id
      integer j, n, Ndim, Nids, KI(maxsp), NU(maxsp)
#ifdef MIKE
      Ndim = maxsp
      if ((id.le.0).or.(id.gt.Nspec)) then
         write(6,*) 'FORT_FINDRHS:  species id out of range: ',id
         call bl_abort(" ")
      end if
      Nreacs = 0
      do j=1,Nreac
         CALL CKINU(j, Ndim, IWRK(ckbi), RWRK(ckbr), Nids, KI, NU)
         do n=1,Nids
            if ((KI(n).eq.id).and.(NU(n).gt.0)) then
               Nreacs = Nreacs + 1
               reactions(Nreacs) = j
            endif
         end do
      end do
#else
      call bl_abort("FORT_FINDRHS not implemented")
#endif
      end

      subroutine FORT_SETNU(nu,lenNU)
      implicit none
#include "cdwrk.H"
      integer lenNU
      integer nu(maxreac,maxspec)
      integer i

      if (lenNU .lt. maxreac*maxspec) then
         write(6,*) 'FORT_CKNU:  nu work array too small: '
         call bl_abort(" ")
      endif
      call CKNU(maxreac, IWRK(ckbi), RWRK(ckbr), nu)
      end

      subroutine FORT_CKINU(Nids,KI,lenKI,NU,lenNU,rxnID,nuAll)
      implicit none
#include "cdwrk.H"
      integer lenKI,lenNU,NDIM1
      integer rxnID, Nids, KI(lenKI), NU(lenNU), nuAll(maxreac,maxspec)
      integer Ndim, k
      Ndim = MIN(lenKI,lenNU)
      if ((rxnID.le.0).or.(rxnID.gt.Nreac)) then
         write(6,*) 'FORT_CKINU:  reaction id out of range: ',rxnID
         call bl_abort(" ")
      end if
      if (Ndim.lt.maxsp) then
         call bl_abort('FORT_CKINU:  KI or NU not long enough')
      end if
      Nids = 0
      do k=1,Nspec
         if (nuAll(rxnID,k).ne.0) then
            Nids = Nids + 1
            KI(Nids) = k
            NU(Nids) = nuAll(rxnID,k)
         endif
      enddo
      end

      integer function FORT_CKELTXINSPY(eltID, spID)
      implicit none
#include "cdwrk.H"
      integer eltID, spID, i, j
      integer NCF(maxelts,maxspec)
      CALL CKNCF(maxelts, IWRK(ckbi), RWRK(ckbr), NCF)
      FORT_CKELTXINSPY = NCF(eltID+1,spID+1)
      end

      integer function FORT_GETCKNUMSPEC()
      implicit none
#include "cdwrk.H"
      FORT_GETCKNUMSPEC = Nspec
      end

      integer function FORT_GETCKNUMELT()
      implicit none
#include "cdwrk.H"
      FORT_GETCKNUMELT = Nelt
      end

      integer function FORT_GETCKNUMREAC()
      implicit none
#include "cdwrk.H"
      FORT_GETCKNUMREAC = Nreac
      end

      double precision function FORT_RUNIV()
      implicit none
      double precision Ruc, Pa
#include "cdwrk.H"
      call CKRP(IWRK(ckbi),RWRK(ckbr),FORT_RUNIV,Ruc,Pa)
c     1 erg/(mole.K) = 1.e-4 J/(kmole.K)
      FORT_RUNIV = FORT_RUNIV*1.d-4
      end

      double precision function FORT_P1ATMMKS()
      implicit none
      double precision Ru, Ruc, Pa
#include "cdwrk.H"
      call CKRP(IWRK(ckbi),RWRK(ckbr),Ru,Ruc,Pa)
c     1 N/(m.m) = 0.1 dyne/(cm.cm)
      FORT_P1ATMMKS = Pa*1.d-1
      end

      integer function FORT_GETCKELTNAME(i, coded)
      implicit none
#include "cdwrk.H"
      integer i
      integer coded(*)
      integer names(0:maxelts*2)
      integer ctr, nlen
      nlen = 2
      call CKSYME(names,nlen)
      coded(1) = names(2*(i-1)  )
      coded(2) = names(2*(i-1)+1)
      if (coded(2).eq.ICHAR(' ')) then
         FORT_GETCKELTNAME = 1
      else
         FORT_GETCKELTNAME = 2
      endif
      end

      integer function FORT_GETCKSPECNAME(i, coded)
      implicit none
#include "cdwrk.H"
      integer i
      integer coded(*)
      integer names(maxspec*maxspnml)
      integer j, str_len
      call CKSYMS(names, maxspnml)
      do j = 1, maxspnml
         coded(j) = names(maxspnml*(i-1)+j)
      end do
      do j = 1, maxspnml
         if (coded(j).eq.ICHAR(' ')) then
            str_len = j
            exit
         endif 
      end do
      FORT_GETCKSPECNAME = str_len - 1
      end

      integer function FORT_CKSYMR(fortReacIdx, coded)
      implicit none
#include "cdwrk.H"
      integer fortReacIdx
      integer coded(*)
      character*(72) line 
      integer j, str_len, istr, iend, lout
      logical error
      data error /.false./
#ifdef MIKE
      lout = 6
      call CKSYMR(fortReacIdx,lout,IWRK(ckbi),RWRK(ckbr),
     &     CWRK(ckbc),str_len,line,error)
      if (error) then
         write(lout,*) 'Could not get reaction name for ',fortReacIdx
         call bl_abort(" ")
      end if
c      
c     Encode the name for transfer to C++
c
      istr = 1
      do while (line(istr:istr) .EQ. ' ')
         istr = istr + 1
      end do
      do j = 0, str_len-1
         coded(j+1) = ICHAR(line(istr+j:istr+j))
      end do
      FORT_CKSYMR = str_len
#else
      FORT_CKSYMR = 0
      call bl_abort("FORT_CKSYMR not implemented")
#endif
      end

      subroutine get_spec_name(name, j)
      implicit none
#include "cdwrk.H"
      integer i, j, FORT_GETCKSPECNAME
      integer coded(maxspnml), len
      character*(maxspnml) name
      len = FORT_GETCKSPECNAME(j, coded)
      do i = 1, maxspnml
         name(i:i) = ' '
      end do
      do i = 1, len
         name(i:i) = char(coded(i))
      end do
      end

      subroutine get_spec_number(name, j)
      implicit none
#include "cdwrk.H"
      integer j, n
      character*(*) name
      character*(maxspnml) locName
      
      j = -1
      do n = 1, Nspec
         call get_spec_name(locName, n)
         if (locName .EQ. name) j = n
      end do
      end

      subroutine FORT_GETCKMWT(mwt)
      implicit none
#include "cdwrk.H"
      REAL_T mwt(*)
      integer n
c     Result in kg/kmole
      call CKWT(IWRK(ckbi),RWRK(ckbr),mwt)
      end

      subroutine FORT_GETCKAWT(awt)
      implicit none
#include "cdwrk.H"
      REAL_T awt(*)
#ifdef MIKE
      integer n
c     Result in kg/kmole
      call CKAWT(IWRK(ckbi),RWRK(ckbr),awt)
#else
      call bl_abort("FORT_GETCKAWT not implemented")
#endif
      end

      subroutine conpFY(N, TIME, Z, ZP, RPAR, IPAR)
      implicit none
#include "cdwrk.H"
#include "conp.H"
#include "visc.H"
      REAL_T TIME, Z(NEQ), ZP(NEQ), RPAR(*)
      integer N, IPAR(*)
      
      REAL_T RHO, CPB, SUM, H, WDOT, WT, THFAC
      integer K

      REAL_T CONC(maxspec), WDOTS(maxspec), ENTHALPY(maxspec)
C
C     Variables in Z are:  Z(1)   = T
C                          Z(K+1) = Y(K)
      
      CALL CKRHOY(RPAR(NP),Z(1),Z(2),IPAR(ckbi),RPAR(ckbr),RHO)
      CALL CKCPBS(Z(1),Z(2),IPAR(ckbi),RPAR(ckbr),CPB)
      CALL CKYTCP(RPAR(NP),Z(1),Z(2),IPAR(ckbi),RPAR(ckbr),CONC)
c
c     Get net production rates.  Compute such that production from -ve
c     reactants gives zero contrib.
c      
c      do k=1,Nspec
c         if (CONC(k) .lt. zero) write(6,*) '.....negative C',k
c      end do
      
#define MAKE_C_POS
#undef MAKE_C_POS
#ifdef MAKE_C_POS
      do k=1,Nspec
         CONC(k) = MAX(CONC(k),zero)
      end do 
#endif
      
      CALL CKWC(Z(1), CONC, IPAR(ckbi), RPAR(ckbr), WDOTS)
      CALL CKHMS(Z(1), IPAR(ckbi), RPAR(ckbr), ENTHALPY)
C
C     Form governing equation
C
      THFAC = one / thickFacCH
      SUM = zero
      DO K = 1, Nspec
         H    = ENTHALPY(K)
         WDOT = WDOTS(K) * THFAC
         WT   = RPAR(NWT+K-1)
         ZP(K+1) = WDOT * WT / RHO
         SUM = SUM + H * WDOT * WT
      END DO
      ZP(1) = -SUM / (RHO*CPB)

#if 0
      print*, 'Z:'
      do k = 1, Nspec+1
         write(6,996) Z(K)
      end do
      print*, 'ZP:'
      do k = 1, Nspec+1
         write(6,996) ZP(K)
      end do
      print*, 'WDOT:'
      do k = 1, Nspec
         write(6,996) WDOTS(K)
      end do

996   format(e30.22)
#endif

      END

      subroutine conpJY(N, TN, Y, SAVF, NFE, FTEM, ML,
     &                  MU, PD, NRPD, RPAR, IPAR)
      implicit none
      REAL_T SAVF
      REAL_T PD, RPAR(*), TN, Y
      dimension SAVF(*)
      integer N, NRPD, ML, MU, IPAR(*)
      dimension Y(N), PD(NRPD,N)
      REAL_T FTEM
      dimension FTEM(*)
      integer NFE
#include "cdwrk.H"
#include "conp.H"
      call bl_abort("conpJY: SHOULD NOT BE HERE!")
      END

      integer function TfromeYpt(T,ein,Y,errMax,NiterMAX,res)
      implicit none
#include "cdwrk.H"
      integer NiterMAX,Niter,n,NiterDAMP
      REAL_T T,ein,Y(*),errMAX,res(0:NiterMAX-1)
      REAL_T TMIN,TMAX,e

      parameter (TMIN=200, TMAX=5000)
      REAL_T  T0,h,cp,cv,de,temp,RoverWbar,Wbar,RU,RUC,P1ATM
      REAL_T dT, etarg
      logical out_of_bounds, converged, soln_bad, stalled
      REAL_T e300,cv300,e6500,cv6500
      integer ihitlo,ihithi

      out_of_bounds(temp) = (temp.lt.TMIN) .or. (temp.gt.TMAX)

      NiterDAMP = NiterMAX
      if ((T.GE.TMIN).and.(T.LE.TMAX)) then
         T0 = T
      else
         T0 = half*(TMIN+TMAX)
         T = T0
      end if
      Niter = 0
      de = zero
      soln_bad = .FALSE.
      etarg = ein * BL_REAL_E(1.0,4)
      ihitlo = 0
      ihithi = 0

      CALL CKUBMS(T,Y,IWRK(ckbi),RWRK(ckbr),e)

      de = two*ABS(e - etarg)/(one + ABS(e) + ABS(etarg))
      res(Niter) = de
      converged = de.le.errMAX

      do while ((.not.converged) .and. (.not.soln_bad))
         CALL CKCVBS(T,Y,IWRK(ckbi),RWRK(ckbr),cv)
         dT = (etarg - e)/cv
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
            TfromeYpt = -1
            goto 100
         else
            CALL CKUBMS(T,Y,IWRK(ckbi),RWRK(ckbr),e)
            de = two*ABS(e - etarg)/(one + ABS(e) + ABS(etarg))
            res(Niter) = de
            Niter = Niter + 1
         end if
         if (Niter .ge. NiterMAX) then
            TfromeYpt = -2
            goto 100
         endif
         converged = (de.le.errMAX) .or. (ABS(dT).le.errMAX)

         if((ihitlo.eq.1).and.(e.gt.etarg))then
            T = 300.d0
            CALL CKUBMS(T,Y,IWRK(ckbi),RWRK(ckbr),e300)
            CALL CKCVBS(T,Y,IWRK(ckbi),RWRK(ckbr),cv300)
            T=300.d0+(etarg-e300)/cv300
            converged = .true.
         endif
         if((ihithi.eq.1).and.(e.lt.etarg))then
            T = 6500.d0
            CALL CKUBMS(T,Y,IWRK(ckbi),RWRK(ckbr),e6500)
            CALL CKCVBS(T,Y,IWRK(ckbi),RWRK(ckbr),cv6500)
            T=6500.d0+(etarg-e6500)/cv6500
            converged = .true.
         endif

      end do

c     Set max iters taken during this solve, and exit
      TfromeYpt = Niter
      return

c     Error condition....dump state and bail out
 100  continue

      write(6,997) 'T from (e,Y): failed'
      write(6,997) 'iterations tried = ',Niter
      write(6,998) 'initial T = ',T0
      write(6,998) 'current T = ',T
      write(6,998) 'species mass fracs:'
      do n = 1,Nspec
         write(6,998) '  ',Y(n)
      end do
      write(6,998)
      write(6,998) 'residual = e - h + RT/Wbar [cgs]'
      do n = 0,Niter-1
         write(6,998) '  ',res(n)
      end do

 997  format(a,3(i4,a))
 998  format(a,d21.12)

      end

      subroutine FORT_TfromHYpt(T,Hin,Y,errMax,NiterMAX,res,Niter)
      implicit none
#include "cdwrk.H"
      REAL_T T,Y(*),H,Hin
      REAL_T TMIN,TMAX,errMAX
      integer NiterMAX,Niter,n,NiterDAMP
      parameter (TMIN=250, TMAX=5000)
      REAL_T  T0,cp,cv,dH,temp,RoverWbar,Wbar,RU,RUC,P1ATM
      REAL_T res(0:NiterMAX-1),dT, Htarg
      logical out_of_bounds, converged, soln_bad, stalled
      REAL_T h300,cp300,h6500,cp6500
      integer ihitlo,ihithi,j

      out_of_bounds(temp) = (temp.lt.TMIN-one) .or. (temp.gt.TMAX)

      NiterDAMP = NiterMAX
      if ((T.GE.TMIN).and.(T.LE.TMAX)) then
         T0 = T
      else
         T0 = half*(TMIN+TMAX)
         T = T0
      end if
      Niter = 0
      dH = zero
      soln_bad = .FALSE.
      Htarg = Hin * 1.d4
      ihitlo = 0
      ihithi = 0

      CALL CKHBMS(T,Y,IWRK(ckbi),RWRK(ckbr),H)
      dH = two*ABS(H - Htarg)/(one + ABS(H) + ABS(Htarg))
      res(Niter) = dH
      converged = dH.le.errMAX

      do while ((.not.converged) .and. (.not.soln_bad))

         CALL CKCPBS(T,Y,IWRK(ckbi),RWRK(ckbr),cp)
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
            CALL CKHBMS(T,Y,IWRK(ckbi),RWRK(ckbr),H)
            dH = two*ABS(H - Htarg)/(one + ABS(H) + ABS(Htarg))
            res(Niter) = dH
            Niter = Niter + 1
         end if
         if (Niter .ge. NiterMAX) then
            Niter = -2
            goto 100
         endif
         converged = (dH.le.errMAX) .or. (ABS(dT).le.errMAX)

         if ((ihitlo.eq.1).and.(H.gt.Htarg)) then
            T = TMIN
            CALL CKHBMS(T,Y,IWRK(ckbi),RWRK(ckbr),h300)
            CALL CKCPBS(T,Y,IWRK(ckbi),RWRK(ckbr),cp300)
            T=TMIN+(Htarg-h300)/cp300
            converged = .true.
         endif
         if ((ihithi.eq.1).and.(H.lt.Htarg)) then
            T = TMAX
            CALL CKHBMS(T,Y,IWRK(ckbi),RWRK(ckbr),h6500)
            CALL CKCPBS(T,Y,IWRK(ckbi),RWRK(ckbr),cp6500)
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
      call bl_pd_myproc(i)
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
    
      block data tranjunk
#include "cdwrk.H"
#include "conp.H"
      
      data verbose_vode / 0 /
      data max_vode_subcycles / 15000 /
      data spec_scalY / maxspec*one /
      data thickFacCH / 1.d0 /
      data nstiff / 1 /

      end
