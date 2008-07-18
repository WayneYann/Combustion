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
         open(unit=LLINKMC,status='UNKNOWN',form='FORMATTED',
     &        file=filet(1:nametrlen))
      else
         call bl_abort('FORT_CKINTERP: transport file must be provided')
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
      integer ioproc, myid, ierr, n, RTOT, lout
      parameter (ioproc = 0)
      character*(maxspnml) name

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
      NRHO  = NP    + 1
      NTEMP = NRHO  + 1
      NHMIX = NTEMP + 1
      NWT   = NHMIX + 1
      NC    = NWT   + Nspec
      NH    = NC    + Nspec
      NCPK  = NH    + Nspec
      NWDOT = NCPK  + Nspec
      NZ    = NWDOT + NEQ
      RTOT  = NZ    + NEQ - 1
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
      if (iN2.eq.-1) write(6,*) '.....warning: no N2 in CKD_2D'
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
      imaxord = maxord
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
      integer names(maxelts*2)
      integer ctr, nlen
      nlen = 2
      call CKSYME(names,nlen)
      coded(1) = ICHAR(names(2*i  ))
      coded(2) = ICHAR(names(2*i+1))
      FORT_GETCKELTNAME = 2
      end

      integer function FORT_GETCKSPECNAME(i, coded)
      implicit none
#include "cdwrk.H"
      integer i
      integer coded(*)
      int names(maxspec*maxspnml)
      integer j, str_len
      call CKSYMS(names, maxspnml)
      str_len = INDEX(names(maxspnml*i:maxspnml*(i+1)-1),' ')
      do j = 1, maxspnml
         coded(j) = ICHAR(names(maxspnml*i+j))
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
      integer j, lout
      character*(maxspnml) name
      character*(maxspnml) dummyCh(maxspec)
      
      logical error
      data error /.false./

      lout = 6
      call CKSYMS(CWRK(ckbc),lout,dummyCh,error)
      if (error) then
         call bl_abort('Could not get species names')
      end if
      name = dummyCh(j)
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
C
C     Variables in Z are:  Z(1)   = T
C                          Z(K+1) = Y(K)
      
      CALL CKRHOY(RPAR(NP),Z(1),Z(2),IPAR(ckbi),RPAR(ckbr),RHO)
      CALL CKCPBS(Z(1),Z(2),IPAR(ckbi),RPAR(ckbr),CPB)
      CALL CKYTCP(RPAR(NP), Z(1), Z(2), IPAR(ckbi), RPAR(ckbr), RPAR(NC))
c
c     Get net production rates.  Compute such that production from -ve
c     reactants gives zero contrib.
c      
c      do k=1,Nspec
c         if (RPAR(NC-1+k) .lt. zero) write(6,*) '.....negative C',k
c      end do
      
#define MAKE_C_POS
#undef MAKE_C_POS
#ifdef MAKE_C_POS
      do k=1,Nspec
         RPAR(NC-1+k) = MAX(RPAR(NC-1+k),zero)
      end do 
#endif
      
      CALL CKWC(Z(1), RPAR(NC), IPAR(ckbi), RPAR(ckbr), RPAR(NWDOT))
      CALL CKHMS(Z(1), IPAR(ckbi), RPAR(ckbr), RPAR(NH))
C
C     Form governing equation
C
      THFAC = one / thickFacCH
      SUM = zero
      DO K = 1, Nspec
         H    = RPAR(NH    + K - 1)
         WDOT = RPAR(NWDOT + K - 1) * THFAC
         WT   = RPAR(NWT   + K - 1)
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
         write(6,996) RPAR(NWDOT + K - 1)
      end do

996   format(e30.22)
#endif

      END

      subroutine conpFC(N, TIME, Z, ZP, RPAR, IPAR)
      implicit none
#include "cdwrk.H"
#include "conp.H"
      REAL_T TIME, Z(NEQ), ZP(NEQ), RPAR(*)
      integer N, IPAR(*)
      
      REAL_T RHO, CPB, CP,SUM, H, WDOT, WT
      integer K,NY
      
      NY = NC
C
C     Variables in Z are:  Z(1)   = T, Z(K+1) = C(K)
C     Variables in RPAR:  RPAR(NP) = P
C                         RPAR(NWDOT) = WDOT 
c                         RPAR(NH)    = H
c                         RPAR(NY)    = Y
c                         RPAR(NRHO)  = RHO
c
      CALL CKCTY(Z(2), IPAR(ckbi), RPAR(ckbr), RPAR(NY))
      CALL CKCPBS(Z(1),RPAR(NY),IPAR(ckbi),RPAR(ckbr),CPB)
c
c     Get net production rates.  Compute such that production from -ve
c     reactants gives zero contrib.
c      
      CALL CKWC(Z(1), Z(2), IPAR(ckbi), RPAR(ckbr), RPAR(NWDOT))
      CALL CKHMS(Z(1), IPAR(ckbi), RPAR(ckbr), RPAR(NH))
      CALL CKCPMS(Z(1),IPAR(ckbi), RPAR(ckbr), RPAR(NCPK))
C
C     Form governing equation
C
      SUM = zero
      DO K = 1, Nspec
         WDOT = RPAR(NWDOT + K - 1)
         H    = RPAR(NH    + K - 1)
         WT   = RPAR(NWT   + K - 1)
         CP   = RPAR(NCPK  + K - 1)
         ZP(K+1) = WDOT 
         SUM = SUM + H * WDOT * WT
      ENDDO
      ZP(1) = -SUM / (RPAR(NRHO)*CPB)
      END

      subroutine conpFCnoT(N, TIME, Z, ZP, RPAR, IPAR)
      implicit none
#include "cdwrk.H"
#include "conp.H"
      REAL_T TIME, Z(NEQ), ZP(NEQ), RPAR(*)
      integer N, IPAR(*)
      
      REAL_T Hmix, Temp0
      integer NY, k
      
      REAL_T HtoTerrMAX
      integer HtoTiterMAX, Niter, FORT_TfromHYpt
      parameter (HtoTerrMAX = BL_REAL_E(1.0,-8))
      parameter (HtoTiterMAX = 20)
      REAL_T res(0:HtoTiterMAX-1)

      NY = NC
C
C     Variables in Z are:  Z(K) = C(K) in moles/cm**3
C     Variables in RPAR:  RPAR(NHMIX) = HMIX (erg/g = 1.e-4 J/kg)
c                         RPAR(NY)    = Y
c                         RPAR(NTEMP) = Temp
c
      CALL CKCTY(Z,IPAR(ckbi),RPAR(ckbr),RPAR(NY))

      Temp0 = RPAR(NTEMP)
      Hmix = RPAR(NHMIX)*BL_REAL_E(1.0,-4)
      Niter = FORT_TfromHYpt(RPAR(NTEMP),Hmix,RPAR(NY),HtoTerrMax,HtoTiterMAX,res)
      if (Niter .lt. 0) goto 100
      CALL CKWC(RPAR(NTEMP),Z,IPAR(ckbi),RPAR(ckbr),ZP)
      RETURN

c     Error condition....dump state and bail out
 100  continue
      
      write(6,997) 'T from (Hmix,Y): failed in conpFCnoT()'
      write(6,997) 'iterations tried = ',Niter
      write(6,998) 'initial T = ',Temp0
      write(6,998) 'current T = ',RPAR(NTEMP)
      write(6,998) 'species mass fracs:'
      do k = 1,Nspec
         write(6,998) '  ',RPAR(NY+k-1)
      end do
      write(6,998)
      write(6,998) 'residual = HMIX - sum(H(T)_i.Y_i):'
      do k = 0,HtoTiterMAX-1
         write(6,998) '  ',res(k)
      end do
      
 997  format(a,2(i4,a))
 998  format(a,d21.12)
      end

      subroutine conpJC(N, TN, C, SAVF, NFE, FTEM, ML,
     &                  MU, PD, NRPD, RPAR, IPAR)
      implicit none
      REAL_T SAVF
      REAL_T PD, RPAR(*), TN, C
      dimension SAVF(*)
      integer N, NRPD, ML, MU, IPAR(*)
      dimension C(N), PD(NRPD,N)
      REAL_T FTEM
      dimension FTEM(*)
      integer NFE
#include "cdwrk.H"
#include "conp.H"
      call bl_abort('conpJC: what are you doing calling me?')
      end

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

      integer function FORT_TfromHYpt(T,Hin,Y,errMax,NiterMAX,res)
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
            FORT_TfromHYpt = -1
            goto 100
         else
            CALL CKHBMS(T,Y,IWRK(ckbi),RWRK(ckbr),H)
            dH = two*ABS(H - Htarg)/(one + ABS(H) + ABS(Htarg))
            res(Niter) = dH
            Niter = Niter + 1
         end if
         if (Niter .ge. NiterMAX) then
            FORT_TfromHYpt = -2
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
c
c     Set max iters taken during this solve and exit
c
      FORT_TfromHYpt = Niter
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
      do n = 0,Niter-1
         write(6,998) '  ',res(n)
      end do

 997  format(a,3(i4,a))
 998  format(a,d21.12)
      end

      subroutine chemeqpt(Perg,Tin,Tout,Yin,Yout,dt)
      implicit none
#include "cdwrk.H"
#include "conp.H"
      REAL_T Perg,Tin,Tout,Yin(maxspec),Yout(maxspec),dt
      REAL_T hmix,e,c(maxspec)
      integer n
      
      REAL_T epsmn,epsmx,dtmn,tnot,ymn(maxspec),prt
      integer itermx,ns,TfromeYpt,FORT_TfromHYpt,Niter

      REAL_T HtoTerrMAX
      integer HtoTiterMAX
      parameter (HtoTerrMAX = BL_REAL_E(1.0,-8))
      parameter (HtoTiterMAX = 20)
      REAL_T res(0:HtoTiterMAX-1)

      logical const_h
      external ckdfe_consth,ckdfe_conste

      const_h = .false.
      const_h = .true.

      CALL CKYTCP(Perg,Tin,Yin,IWRK(ckbi),RWRK(ckbr),c)
      RWRK(NTEMP) = Tin
      if (const_h) then
         CALL CKHBMS(Tin,Yin,IWRK(ckbi),RWRK(ckbr),hmix)
         RWRK(NHMIX) = hmix*BL_REAL_E(1.0,-4)
      else
         CALL CKUBMS(Tin,Yin,IWRK(ckbi),RWRK(ckbr),e)
         RWRK(NHMIX) = e*BL_REAL_E(1.0,-4)
      end if

c     Reset chemeq2 defaults
      epsmn = zero
      epsmx = zero
      dtmn = zero
      tnot = zero
      itermx = 1
      ns = Nspec
      do n=1,Nspec
         ymn(n) = BL_REAL_E(0.,0)
      end do
      prt = 1
      call chemsp(epsmn, epsmx, dtmn, tnot, itermx, ns, ymn, prt)

c     Call solver
      if (const_h) then
         call chemeq2(dt,ckdfe_consth,Nspec,c)
      else
         call chemeq2(dt,ckdfe_conste,Nspec,c)
      end if

c     Pull output data from interface
      call CKCTY(c, IWRK(ckbi), RWRK(ckbr), Yout)
      Tout = RWRK(NTEMP)

c     Update temperature field
      if (const_h) then
         Niter = FORT_TfromHYpt(Tout,RWRK(NHMIX),Yout,HtoTerrMax,HtoTiterMAX,res)
         if (Niter .lt. 0) then
            write(6,*) 'T from h,y solve in chemeqpt() failed',Niter
            call bl_abort(" ")
         end if
      else
         Niter = TfromeYpt(Tout,RWRK(NHMIX),Yout,HtoTerrMax,HtoTiterMAX,res)
         if (Niter .lt. 0) then
            write(6,*) 'T from e,y solve in chemeqpt() failed',Niter
            call bl_abort(" ")
         end if
      end if

      RWRK(NTEMP) = Tout
      end

      subroutine ckdfe_consth(c, q, d, t)
c
c     Derivative function evaluator based on ChemKin for use in ChemEQ2
c
      implicit none
#include "cdwrk.H"
#include "conp.H"
      REAL_T c(*), q(*), d(*), t
      integer FORT_TfromHYpt, Niter
      REAL_T HtoTerrMax
      integer HtoTiterMAX
      parameter (HtoTerrMax = BL_REAL_E(1.0,-8))
      parameter (HtoTiterMAX = 20)
      REAL_T temp, hmix, res(0:HtoTiterMAX-1), y(maxspec)
      integer n

      temp = RWRK(NTEMP)
      hmix = RWRK(NHMIX)

c     Update temperature to agree with current species/enthalpy values
      call CKCTY(c, IWRK(ckbi), RWRK(ckbr), y)
      Niter = FORT_TfromHYpt(temp,hmix,y,HtoTerrMax,HtoTiterMAX,res)
      if (Niter .lt. 0) then
         write(6,*) 'T from h,y solve in ckdfe_consth() failed',Niter
         call bl_abort(" ")
      end if
      RWRK(NTEMP) = temp
#ifdef MIKE
      call CKCDC(temp, c, IWRK(ckbi), RWRK(ckbr), q, d)
#else
      call bl_abort("ckdfe_consth not implemented")
#endif
      end

      subroutine ckdfe_conste(c, q, d, t)
c
c     Derivative function evaluator based on ChemKin for use in ChemEQ2
c
      implicit none
#include "cdwrk.H"
#include "conp.H"
      REAL_T c(*), q(*), d(*), t
      integer TfromeYpt, Niter
      REAL_T HtoTerrMax
      integer HtoTiterMAX
      parameter (HtoTerrMax = BL_REAL_E(1.0,-8))
      parameter (HtoTiterMAX = 20)
      REAL_T temp, e, res(0:HtoTiterMAX-1), y(maxspec)
      integer n

      temp = RWRK(NTEMP)
      e = RWRK(NHMIX)

c     Update temperature to agree with current species/enthalpy values
      call CKCTY(c, IWRK(ckbi), RWRK(ckbr), y)
      Niter = TfromeYpt(temp,e,y,HtoTerrMax,HtoTiterMAX,res)
      if (Niter .lt. 0) then
         call bl_abort('T from e,y solve in ckdfe_conste() failed')
      end if

      RWRK(NTEMP) = temp
#ifdef MIKE
      call CKCDC(temp, c, IWRK(ckbi), RWRK(ckbr), q, d)
#else
      call bl_abort("ckdfe_consth not implemented")
#endif
      end

      subroutine chemeq2(dtg, gsub, n, y)
c
cd* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
cd
cd    chemeq2(dtg, gsub, n, y)
cd
cd    original chemeq development:
cd    originators:  t.r. young                                  nrl 1982
cd    vax version:  t.r. young            nrl code 4040         may 1983
cd    workstation:  g. patnaik            berkeley research     jun 1995
cd
cd    chemeq2 development: d.r. mott      nrl code 6404         may 1999
cd
cd
cd    Description:  Subroutine chemeq2 solves a class of "stiff" ODEs
cd    associated with reactive flow problems that cannot be readily
cd    solved by the standard classical methods. In contrast to the
cd    original chemeq subroutine, this version uses the same
cd    quasi-steady-state update for every species regardless of the
cd    timescale for  that species. An adaptive stepsize is  chosen to
cd    give accurate results for the fastest changing quantity, and a
cd    stability check on the timestep is also available when the
cd    corrector is iterated.
cd
cd    NOTE:  The accuracy-based timestep calculation can be augmented
cd    with a stability-based check when at least three corrector
cd    iterations are performed.  To include this check, "uncomment"
cd    the lines that start with "D", or use the compiler flag "-d_lines"
cd    if available to compile the code including these lines.  If the
cd    lines are manually uncommented,  the continuation characters
cd    must be placed in the correct column.  For most problems, the
cd    stability check is not needed, and eliminating the calculations
cd    and logic associated with the check enhances performance.
cd
cd    The routine assumes that all of the integrated quantites and the
cd    time step are positive.
cd
cd    argument list definition (name, type, description, input vs. output):
cd    dtg            real         the interval of integration or the   i
cd                                range of the independent variable.
cd                                0.0 <= t <= dtg. (global timestep)
cd    gsub           real         the name of the derivitive function  i
cd                                evaluator subroutine.
cd    n              integer      the number of equations to be        i
cd                                integrated. an error exisis if n is
cd                                greater than nd set by the parameter
cd                                statement.
cd    y(n)           real         the initial values at call time    i/o
cd                                and the final values at return time.
cd
cd    Language and limitations:  This subroutine is written in standard
cd    FORTRAN 77.  For high accuracy, this routine should be compiled
cd    using whatever "double precision" flag is appropriate for the
cd    platform being used (such as "f77 -r8 . . . .")
cd
cd    Entry points:  Four entry points are provided for flexibility and
cd    optimum control.  This structure was maintained from the original
cd    chemeq subroutine to ensure compatiability with previous
cd    applications that use chemeq.
cd
cd    chemeq2:  advances the equations the given increment 'dtg'.
cd
cd    chemct:  informative, prints the values of the indicative
cd             counters listed below;
cd             1.  the number of derivative function evaluations.
cd             2.  the number times the integration step was restarted
cd                 due to nonconvergence of the predictor-corrector
cd                 scheme.
cd
cd    chemsp:  provides the user with the option to reset the most
cd             important control parameters.
cd
cd    chempr:  informative, prints out internal variables for diagnostic
cd             purposes.
cd
cd    subroutines referenced:
cd
cd    gsub;    whose actual name and definition are supplied by the user
cd             is called to obtain the derivitive functions.
cd
cd    call gsub(y, q, d, t)
cd    argument list to gsub;
cd    y(n)              real      current values of the dependent      i
cd                                variable.
cd    q(n)              real      calculated formation rates.          o
cd    d(n)              real      calculated loss rates.               o
cd    t                 real      current value of the independent     i
cd                                variable.
cd
cd    chemer: Called whenever an error is detected. Currently the
cd            only error recognized is a time step that is too small.
cd
cd    call chemer(y, n)
cd    argument list to chemer;   (same definition as "chemeq2").
cd
cd* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      implicit none
      integer nd
#include "cdwrk.H"
cHACK      parameter    (nd = 10)
      parameter    (nd = maxspec)
c
      external gsub
      integer      n, ns, lo, i
      integer      itermax, iter, itermx
c
c     the following are counters (this call & total) for gsub calls
c      and timestep repeats
c
      integer      gcount,  rcount, tgcnt, trcnt
      integer      gc,rc,tgc,trc
c
      REAL_T       ts, tn, tfd, tmk
      REAL_T       y(n)
      REAL_T       ymin(nd), ymn(nd)
      REAL_T       q(nd), d(nd), rtaus(nd), y1(nd)
      REAL_T       ys(nd), y0(nd), rtau(nd)
      REAL_T       alpha, qs(nd)
      REAL_T       scr1, scr2, scrarray(nd)
      REAL_T       epscl, dtg, dtmin, sqreps, tstart, dt, dto
      REAL_T       epsmax, epsmin, rswitch
      REAL_T       epsmx, epsmn, dtmn, tnot, prt
      REAL_T       scrtch, ascr, eps
      REAL_T       rtaui, rtaub, qt, pb, dtc, rteps
c
c   ym1, ym2, and stab are used only for the stability check on dt
#undef CHEMEQ2_STABILITY_CHECK
#define CHEMEQ2_STABILITY_CHECK
#ifdef CHEMEQ2_STABILITY_CHECK
      REAL_T       ym1(nd), ym2(nd), stab
#endif
c
      save

      data         gcount,  rcount, tgcnt, trcnt/4*0/
      data          itermax/1/, epscl/100.0/
      data         tfd/1.000008/, dtmin/1.0e-15/, sqreps/0.50/
      data         tstart, dt/2*0.0/, tn/0.0e+00/, q/nd*0.0/
      data         epsmax/10.0/, lo/16/, epsmin/1.0e-02/, d/nd*0.0/
      data         rswitch/  5.965900 /
c
c rswitch for 4-4 pade: 5.9659
c
cd     check input parameters.
          if(n .gt. nd) then
             write(lo, 1002) n, nd
 1002     format(5(/),'from -chemeq2-  :  no. of eq.s requested is too',
     .    ' large'/' requested (',i5,'), max. allowed (',i5,')')
             call bl_abort(" ")
          end if
c
c     initialize the control parameters.
  110     tn = 0.0e+00
c
c     store and limit to 'ymin' the initial values.
          do i = 1, n
             q(i) = 0.0
             d(i) = 0.0
             y0(i) = y(i)
             y(i) = max(y(i), ymin(i))
          end do

c     evaluate the derivitives of the initial values.

          call gsub(y, q, d, tn + tstart)
          gcount = gcount + 1
c
c     estimate the initial stepsize.
c
c     strongly increasing functions(q >>> d assumed here) use a step-
c     size estimate proportional to the step needed for the function to
c     reach equilibrium where as functions decreasing or in equilibrium
c     use a stepsize estimate directly proportional to the character-
c     istic stepsize of the function. convergence of the integration
c     scheme is likely since the smallest estimate is chosen for the
c     initial stepsize.
          scrtch = 1.0e-25
          do  i = 1, n
             ascr = abs(q(i))
c             scr2 = sign(1./y(i),.1*epsmin*ascr - d(i) )
c  HACK
             scr2 = sign(1./(y(i)+1.0e-30),.1*epsmin*ascr - d(i) )
             scr1 = scr2 * d(i)
             scrtch = max(scr1,-abs(ascr-d(i))*scr2,scrtch)
          end do
          dt = min(sqreps/scrtch,dtg)
c
c     the starting values are stored.
 100      ts = tn
c
          do i=1,n
c             rtau(i) = dt*d(i)/y(i)
c  HACK
             rtau(i) = dt*d(i)/(y(i)+1.0e-30)
             ys(i) = y(i)
             qs(i) = q(i)
             rtaus(i) = rtau(i)
          end do
c
c
c     find the predictor terms.
  101 continue
c
      do i = 1,n
c
c     prediction
c
         rtaui = rtau(i)
c
c        note that one of two approximations for alpha is chosen:
c        1) Pade b for all rtaui (see supporting memo report)
c             or
c        2) Pade a for rtaui<=rswitch,
c           linear approximation for rtaui > rswitch
c           (again, see supporting NRL memo report (Mott et al., 2000) )
c
c   Option 1): Pade b
c
         alpha = (180.+rtaui*(60.+rtaui*(11.+rtaui)))
     &        / (360. + rtaui*(60. + rtaui*(12. + rtaui)))
c
c   Option 2): Pade a or linear
c
c         if(rtaui.le.rswitch) then
c            alpha = (840.+rtaui*(140.+rtaui*(20.+rtaui)))
c     &           / (1680. + 40. * rtaui*rtaui)
c         else
c            alpha = 1.-1./rtaui
c         end if
c
         scrarray(i) = (q(i)-d(i))/(1.0 + alpha*rtaui)
      end do
c
      iter = 1
      do while(iter.le.itermax)

c     limit decreasing functions to their minimum values.
         do i = 1,n
#ifdef CHEMEQ2_STABILITY_CHECK
            ym2(i) = ym1(i)
            ym1(i) = y(i)
#endif
            y(i) = max(ys(i) + dt*scrarray(i), ymin(i))
         end do

         if(iter.eq.1) then
c
c  the first corrector step advances the time (tentatively) and
c   saves the initial predictor value as y1 for the timestep check later
            tn = ts + dt
            do i=1,n
               y1(i) = y(i)
            end do
         end if
c
c     evaluate the derivitives for the corrector.
c
          call gsub(y, q, d, tn + tstart)
          gcount = gcount + 1
          eps = 1.0e-10

          do i = 1,n
c
c             rtaub = .5*(rtaus(i)+dt*d(i)/y(i))
c  HACK
             rtaub = .5*(rtaus(i)+dt*d(i)/(y(i)+1.e-30))
c
c  Same options for calculating alpha as in predictor:
c
c  Option 1): Pade b
c
                alpha = (180.+rtaub*(60.+rtaub*(11.+rtaub)))
     &               / (360. + rtaub*(60. + rtaub*(12. + rtaub)))

c  Option 2):  Pade a or linear
c
c             if(rtaub.le.rswitch) then
c                alpha = (840.+rtaub*(140.+rtaub*(20.+rtaub)))
c     &               / (1680. + 40.*rtaub*rtaub)
c             else
c                alpha = 1.-1./rtaub
c             end if

             qt = qs(i)*(1. - alpha) + q(i)*alpha
             pb = rtaub/dt
             scrarray(i) = (qt - ys(i)*pb) / (1.0 + alpha*rtaub)
c
          end do
c
          iter = iter + 1
c
       end do
c
c     calculate new f, check for convergence, and limit decreasing
c     functions. the order of the operations in this loop is important.
       do i = 1,n
          scr2 = max(ys(i) + dt*scrarray(i), zero)
          scr1 = abs(scr2 - y1(i))
          y(i) = max(scr2, ymin(i))
#ifdef CHEMEQ2_STABILITY_CHECK
         ym2(i) = ym1(i)
         ym1(i) = y(i)
#endif
c
          if(.25*(ys(i) + y(i)).gt.ymin(i)) then
c             scr1 = scr1/y(i)
c  HACK
             scr1 = scr1/(y(i)+1.0e-30)
             eps = max(.5*(scr1+
     &            min(abs(q(i)-d(i))/(q(i)+d(i)+1.0e-30),scr1)),eps)
c
          end if
       end do
       eps = eps*epscl
c
c     print out dianostics if stepsize becomes too small.

c       if(dt .le. dtmin + 1.0e-16*tn) then
       if(dt .le. 1.0e-16*tn) then
          write(lo, 1003) dt, tn, dtmin
          do i = 1,n
             dtc = epsmin*y(i)/(abs(q(i)-d(i)) + 1.0e-30)
             write(lo, 1004) q(i), d(i), y(i), rtau(i), dtc,
     &            q(i)-d(i),ys(i), y0(i), ymin(i)
          end do

 1003     format('1    chemeq error;   stepsize too small ! ! !', /,
     1    '     dt = ', 1pe10.3, ' tn = ', d25.15,
     2    ' dtmin = ',e10.3, //, 11x, 'q', 10x, 'd', 10x, 'y',
     3         8x, 'rtau', 8x, 'dtc', 7x, 'q - d',7x, 'ys',
     4         9x, 'y0', 8x, 'ymin')
c 1004     format(5x, 1p12e11.3)
 1004     format(5x, 12(1pe11.3))
             dt = dtg - ts
             dt = min(dtmin, abs(dt))
c
c         call error diagnostic routine
             call chemer
c
          end if
c
c         check for convergence.
c
c  The following section is used for the stability check
#ifdef CHEMEQ2_STABILITY_CHECK
          stab = 0.01
          if(itermax.ge.3) then
             do i=1,n
                stab = max(stab, abs(y(i)-ym1(i))/
     &               (abs(ym1(i)-ym2(i))+1.e-20*y(i)))
             end do
          endif
#endif

#ifdef CHEMEQ2_STABILITY_CHECK
       if(eps .le. epsmax .and. stab .le. 1.) then
#else
       if(eps .le. epsmax) then
#endif
c
c     Valid step.  Return if dtg has been reached.
c
          if(dtg .le. tn*tfd) return
       else
c
c     Invalid step; reset tn to ts
c
          tn = ts
       end if
C
c      perform stepsize modifications.
c      estimate sqrt(eps) by newton iteration.
c
       rteps = 0.5*(eps + 1.0)
       rteps = 0.5*(rteps + eps/rteps)
       rteps = 0.5*(rteps + eps/rteps)
c
       dto = dt
#ifdef CHEMEQ2_STABILITY_CHECK
       dt = min(dt*(1.0/rteps + .005), tfd*(dtg - tn), dto/(stab+.001))
#else
       dt = min(dt*(1.0/rteps + .005), tfd*(dtg - tn))
#endif
c
c     begin new step if previous step converged.
c
#ifdef CHEMEQ2_STABILITY_CHECK
       if(eps .gt. epsmax .or. stab. gt. 1) then
#else
       if(eps .gt. epsmax) then
#endif
          rcount = rcount + 1
c
c         After an unsuccessful step the initial timescales don't
c         change, but dt does, requiring rtaus to be scaled by the
c         ratio of the new and old timesteps.
c
          dto = dt/dto
          i = 1
          do while(i.le.n)
             rtaus(i) = rtaus(i)*dto
             i = i+1
          end do
c
c         Unsuccessful steps return to line 101 so that the initial
c         source terms do not get recalculated.
c
          goto 101
       end if
c
c      Successful step; get the source terms for the next step
c      and continue back at line 100
c
       call gsub(y, q, d,tn + tstart)
       gcount = gcount + 1
       go to 100
c
c
c
       entry chemct (tmk,gc,rc,tgc,trc)
c     ----- ------
c
cd* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
cd
cd    chemct (tmk)
cd    write out the values of the various indicative counters that the
cd    program keeps.
cd
cd    argument list definition:
cd    tmk               REAL_T    a floating point number printed      i
cd                                to identify the call.
cd
cd    output variable definition:
cd    tmk               REAL_T    floating point identifier.
cd    gcount            integer   number of derivative subroutine calls
cd                                since the last call.
cd    rcount            integer   number of times stepsize was reduced
cd                                since last call.
cd    tgcnt             integer   total of gcount to this call.
cd    trcnt             integer   total of rcount to this call.
cd
cd* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
          tgcnt = tgcnt + gcount
          trcnt = trcnt + rcount
c
c     print out indicative counters.
c          write(lo, 1000) tmk, gcount, rcount, tgcnt,
c     .    trcnt
c 1000     format(' chemeq indices; tmk = ', 1pe10.3,
c     .    ' gcount, rcount = ', 2i7, ' totals: ', 2i7)

          gc = gcount
          rc = rcount
          tgc = tgcnt
          trc = trcnt
c
c     reset counters.
          gcount = 0
          rcount = 0

          return
c
c
c
      entry chemsp(epsmn, epsmx, dtmn, tnot, itermx, ns, ymn, prt)
c     ----- ------
c
cd* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
cd
cd    chemsp(epsmn, epsmx, dtmn, tnot, itermx, prt)
cd
cd    reset any local control parameters if their respective input
cd    values are greater than zero. default values are used if the
cd    input values are zero or less repectively.
cd
cd    argument list definition:
cd    epsmn             REAL_T    the maximum relative error allowed   i
cd                                for convergence of the corrector step.
cd                                default value:  1.0e-02
cd    epsmx             REAL_T    this number provides the basis for   i
cd                                deciding weather convergence can be
cd                                achieved with out added stepsize
cd                                reduction. if eps/epsmin is greater
cd                                than epsmx further reduction is
cd                                applied.
cd                                default value : 10.0
cd    dtmn              REAL_T    the smallest stepsize allowed.       i
cd                                default value:  1.0e-15
cd    tnot              REAL_T    the initial value of the independent i
cd                                variable t.
cd                                default value:  0.0
cd    itermx            i         number of times the corrector is applied
cd                                default value:  1
cd    ns                integer   number of entries in ymin to reset   i
cd
cd    ymn(nd)           REAL_T    minimum values allowed for y         i
cd                                default value: 1.0e-20
cd    prt               REAL_T    controls the output of chemsp.  any  i
cd                                non zero value suppresses all print
cd                                output from this entry.
cd
cd* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
          epsmin = 1.0e-02
          if(epsmn .gt. 0.0)epsmin = epsmn
          if(epsmn .gt. 0.0)sqreps = 5.0*sqrt(epsmin)
          epscl = 1.0/epsmin
          epsmax = 10.0
          if(epsmx .gt. 0.0)epsmax = epsmx
          dtmin = 1.0e-15
          if(dtmn .gt. 0.0)dtmin = dtmn
          tstart = tnot
          itermax = 1
          if(itermx.gt. 0) itermax = itermx
          do i=1,ns
             ymin(i) = 1.e-20
c             if(ymn(i).gt.0.) ymin(i) = ymn(i)
             if(ymn(i).ge.0.) ymin(i) = ymn(i)
          end do
c
c     print new values of control parameters.
          if(prt .eq. 0.0) then
             write(lo, 1001) epsmn, epsmx, dtmn, tnot, itermx
             write(lo, 1005) ns
             if (ns.gt.0) write(lo,1006) (ymin(i), i=1,ns)
          end if

 1001     format(' initalize "chemeq2" via "chemsp"', /,
     .    ' epsmn, epsmx, dtmn, tnot, itermx = ', 5(1pg10.3))
 1005     format(' ns = ',I5)
 1006     format(' ymin: ',50e12.3)

          return
c
c
c
      entry chempr (y, n)
c     ----- ------
c
cd* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
cd
cd    chempr (y, n)
cd
cd    chempr may be called whenever an error occurs that can be
cd    attributed to the results of chemeq. a partial set of the internal
cd    variables is printed as a diagnostic.
cd
cd    argument list definition:
cd    y(n)         r  current values of the dependent variable.     i
cd    n            i  the number of entries in y and ymin.          i
cd
cd* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
          write(lo, 1003) dt, tn, dtmin
          do 45 i = 1,n
          dtc = epsmin*y(i)/(abs(q(i) - d(i)) + 1.0e-30)
   45     write(lo, 1004) q(i), d(i), y(i), rtau(i),
     &         dtc, q(i)-d(i),  ys(i), y0(i), ymin(i)
c
      end

      subroutine chemer
c     ---------- ------
c
c     diagnostic routine for stiff o.d.e. solver  -chemeq-
c
      print 1001
      print *,' chemer: time step smaller than small'
      print *,' fraction of total integration interval...bailing'
      print *,' (stop 69) executed from -chemer-'
      print 1001
      call bl_abort("69")
 1001 format(2(/))
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
      REAL_T TSAVE, SMH(0:maxspec-1), EQK(0:maxreac-1)
      REAL_T RKFT(0:maxreac-1), RKRT(0:maxreac-1)
      common / junk / TSAVE, SMH, EQK, RKFT, RKRT
      
      data verbose_vode / 0 /
      data max_vode_subcycles / 15000 /
      data spec_scalY / maxspec*one /
      data thickFacCH / 1.d0 /
      data nstiff / 1 /

      end
