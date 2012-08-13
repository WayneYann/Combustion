      !
      ! Hardwired sizes.
      !
      integer maxreac, maxspec, maxelts, maxthrdb, maxtp, maxsp
      integer maxspnml
      parameter (maxreac = 472, maxspec=72, maxelts=6)
      parameter (maxthrdb=20, maxtp=3, maxsp=12, maxspnml = 16)
      !
      ! The only ChemKin stuff we need.
      !
      integer, parameter :: ckbi=1, ckbr=1
      !
      ! EGLib package for transport properties.
      !
      ! For normal LMC usage:
      !
      integer, parameter :: eg_nodes=1, eg_IFLAG=3, eg_ITLS=1
      !
      ! If you're using the DD* routines:
      !
      ! integer, parameter :: eg_nodes=1, eg_IFLAG=7, eg_ITLS=3

      !
      ! DVODE & DVODE driver workspace requirements.
      !
      ! dvbr & dvbi index to start of DVODE real & integer workspace.
      !
      integer, parameter :: dvbr=1, dvbi=1
      !
      ! These are the sizes of the DVODE workspaces. The two double
      ! precision workspaces are grouped together with the driver
      ! space coming after the DVODE space proper.
      !
      integer dvr, dvi, dvdr
      !
      ! dver  indexes to end   of DVODE        workspace.
      ! dvdbr indexes to start of DVODE driver workspace.
      ! dvder indexes to end   of DVODE driver workspace.
      !
      integer dver, dvdbr, dvder

      parameter (dvr = 22 + 9*(maxspec+1) + 2*(maxspec+1)**2)
      parameter (dvi = 30 + maxspec + 1)
      parameter (dvdr = 2*maxspec+3)

      parameter (dver =dvbr+dvr-1)
      parameter (dvdbr=dver+1)
      parameter (dvder=dvdbr+dvdr-1)
      !
      ! Total workspace requirements - add 1 to guard against 0 size.
      !
      integer lenr_fixed, leni_fixed
      parameter (lenr_fixed = dvr + dvdr + 1)
      parameter (leni_fixed = dvi + 1)

      double precision RWRK(lenr_fixed)
      integer          IWRK(leni_fixed)

      common / ckdrwork / RWRK
      common / ckdiwork / IWRK
      save /ckdrwork/, /ckdiwork/
      !
      ! EGlib workspace -- allocated once in FORT_INITCHEM()
      !
      integer,          pointer :: EGIWRK(:)
      double precision, pointer :: EGRWRK(:)
      common / EGWORK / EGIWRK, EGRWRK
      save   / EGWORK /

!$omp threadprivate(/EGWORK/)
      !
      ! Actual sizes
      !
      integer            Nelt, Nspec, Nreac, Nfit
      common / ckdsize / Nelt, Nspec, Nreac, Nfit
      save /ckdsize/
      !
      ! I/O
      !
      integer          LLINKMC
      common / ckdio / LLINKMC
      save /ckdio/

      integer verbose_vode, max_vode_subcycles
      common /ckdio1/ verbose_vode, max_vode_subcycles
      save /ckdio1/
      !
      ! Integrator
      !
      double precision  spec_scalY(maxspec)
      common /vode1/ spec_scalY
      save /vode1/

      double precision  thickFacCH
      common /vode2/ thickFacCH
      save /vode2/
      !            
      ! Transport library
      !
      double precision  TMIN_TRANS
      common /tran2/ TMIN_TRANS
      save /tran2/

      integer Ncoefs,Doffset,PTCoffset,DMIXoffset,TDoffset
      common /tran3/ Ncoefs,Doffset,PTCoffset,DMIXoffset,TDoffset
      save /tran3/
