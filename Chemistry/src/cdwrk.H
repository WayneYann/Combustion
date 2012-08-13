
!     Hardwired sizes
      integer maxreac, maxspec, maxelts, maxthrdb, maxtp, maxsp
      integer   maxspnml
      parameter (maxreac = 472, maxspec=72, maxelts=6)
      parameter (maxthrdb=20, maxtp=3, maxsp=12, maxspnml = 16)

!     ChemKin workspace requirements
      integer ckr, cki, ckl, ckc
      parameter (ckr = 0)
      parameter (cki = 0)
      parameter (ckl = 0)
      parameter (ckc = 0)

!     EGLib package for transport properties (using scalar "S" mode):
!     1) the values IFLAG=7 and ITLS=3 enable all the "S" library routines
!     2) the original coding of egr (which was for the ITLS=3) was in error
!      integer egr, egi, egl, egc, eg_nodes, eg_IFLAG, eg_ITLS
      integer eg_nodes, eg_IFLAG, eg_ITLS
      parameter (eg_nodes = 1)
!      parameter (eg_IFLAG = 7) ! Use this if using DD* routines.
!      parameter (eg_ITLS  = 3) ! Use this if using DD* routines.
      parameter (eg_IFLAG = 3) ! Use this for normal LMC usage.
      parameter (eg_ITLS  = 1) ! Use this for normal LMC usage.
!      parameter (egr = 23 + 14*maxspec + 32*maxspec**2 + 13*eg_nodes 
!     &         + 30*eg_nodes*maxspec + 5*eg_nodes*maxspec**2)
!      parameter (egi = maxspec)
!      parameter (egl = 0)
!      parameter (egc = 0)

!     DVODE workspace requirements      
      integer dvr, dvi, dvl, dvc
      parameter (dvr = 22 + 9*(maxspec+1) + 2*(maxspec+1)**2)
      parameter (dvi = 30 + maxspec + 1)
      parameter (dvl = 0)
      parameter (dvc = 0)
      
!     DVODE driver (for transient ODE integration) workspace requirements
      integer dvdr, dvdi, dvdl, dvdc
      parameter (dvdr = 2*maxspec+3)
      parameter (dvdi = 0)
      parameter (dvdl = 0)
      parameter (dvdc = 0)

!     Workspace memory layout
      integer ckbr, dvbr, dvdbr, cker, dver, dvder
      integer ckbi, dvbi, dvdbi, ckei, dvei, dvdei
      integer ckbl, dvbl, dvdbl, ckel, dvel, dvdel
      integer ckbc, dvbc, dvdbc, ckec, dvec, dvdec
      parameter (ckbr=1,ckbi=1,ckbl=1,ckbc=1)
      parameter (cker=ckbr+ckr-1,ckei=ckbi+cki-1)
      parameter (ckel=ckbl+ckl-1,ckec=ckbc+ckc-1)
      parameter (dvbr=1, dvbi=1, dvbl=1, dvbc=1)
      parameter (dver =dvbr+dvr-1,  dvei=dvbi+dvi-1)
      parameter ( dvel =dvbl+dvl-1,  dvec=dvbc+dvc-1)
      parameter (dvdbr=dver+1,dvdbi=dvei+1,dvdbl=dvel+1,dvdbc=dvec+1)
      parameter (dvder=dvdbr+dvdr-1,dvdei=dvdbi+dvdi-1)
      parameter ( dvdel=dvdbl+dvdl-1,dvdec=dvdbc+dvdc-1)
      
!     Total workspace requirements -- add 1 to guard against null arrays.
      integer lenr_fixed, leni_fixed, lenl_fixed, lenc_fixed
      parameter (lenr_fixed = ckr + dvr + dvdr + 1)
      parameter (leni_fixed = cki + dvi + dvdi + 1)
      parameter (lenl_fixed = ckl + dvl + dvdl + 1)
      parameter (lenc_fixed = ckc + dvc + dvdc + 1)
      
      double precision     RWRK(lenr_fixed)
      integer              IWRK(leni_fixed)
      logical              LWRK(lenl_fixed)
      character*(maxspnml) CWRK(lenc_fixed)
      
      common / ckdrwork / RWRK
      common / ckdiwork / IWRK
      common / ckdlwork / LWRK
      common / ckdcwork / CWRK
      save /ckdrwork/, /ckdiwork/, /ckdlwork/, /ckdcwork/
      !
      ! EGlib workspace -- allocated once in FORT_INITCHEM()
      !
      integer,          pointer :: EGIWRK(:)
      double precision, pointer :: EGRWRK(:)
      common / EGWORK / EGIWRK, EGRWRK
      save   / EGWORK /

!$omp threadprivate(/EGWORK/)

!     Actual sizes
      integer            Nelt, Nspec, Nreac, Nfit
      common / ckdsize / Nelt, Nspec, Nreac, Nfit
      save /ckdsize/

!     I/O
      integer          LLINKMC
      common / ckdio / LLINKMC
      save /ckdio/

      integer verbose_vode, max_vode_subcycles
      common /ckdio1/ verbose_vode, max_vode_subcycles
      save /ckdio1/

!     Integrator
      double precision  spec_scalY(maxspec)
      common /vode1/ spec_scalY
      save /vode1/

      double precision  thickFacCH
      common /vode2/ thickFacCH
      save /vode2/
            
!     Transport library
      double precision  TMIN_TRANS
      common /tran2/ TMIN_TRANS
      save /tran2/

      integer Ncoefs,Doffset,PTCoffset,DMIXoffset,TDoffset
      common /tran3/ Ncoefs,Doffset,PTCoffset,DMIXoffset,TDoffset
      save /tran3/
