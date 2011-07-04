      module cdwrk_module

!     Hardwired sizes
!     integer, parameter ::  maxreac  = 1000, maxspec = 200, maxelts = 20
!     integer, parameter ::  maxthrdb =   10, maxtp   =   3, maxsp = 12
!     integer, parameter ::  maxspnml =   16

      integer, parameter ::  maxreac  =  472, maxspec =  72, maxelts = 6
      integer, parameter ::  maxthrdb =   20, maxtp   =   3, maxsp = 12
      integer, parameter ::  maxspnml =   16

!     integer, parameter ::  maxreac  =    3, maxspec =   6, maxelts = 4 
!     integer, parameter ::  maxthrdb =   10, maxtp   =   3, maxsp = 12
!     integer, parameter ::  maxspnml =   16

!     ChemKin workspace requirements
      integer, parameter ::  ckr = 0
      integer, parameter ::  cki = 0
      integer, parameter ::  ckl = 0
      integer, parameter ::  ckc = 0

!     Multicomponent transport properties package workspace requirements
      integer, parameter ::  MAXFIT = 7, NO = 4, NFDIM = 165, NT = 50
      integer, parameter ::  NRANGE = MAXTP-1, NLITEMAX = 2

      integer, parameter :: mcrTranfit = 0, mcrTranlib = 0
      integer, parameter ::  mcr = mcrTranfit+mcrTranlib

      integer, parameter :: mciTranfit = 0, mciTranlib = 0
      integer, parameter ::  mci = mciTranfit+mciTranlib

      integer, parameter :: mccTranfit = 0, mccTranlib = 0
      integer, parameter ::  mcc = mccTranfit+mccTranlib

      integer, parameter :: mclTranfit = 0, mclTranlib = 0
      integer, parameter ::  mcl = mclTranfit+mclTranlib

!     EGLib package for transport properties (using scalar "S" mode): 
!     1) the values IFLAG=7 and ITLS=3 enable all the "S" library routines 
!     2) the original coding of egr (which was for the ITLS=3) was in error

      integer, parameter ::  egi = maxspec, egl = 0, egc = 0
      integer, parameter ::  eg_nodes = 1, eg_IFLAG = 7, eg_ITLS = 3

      integer, parameter ::  egr = 23 + 14*maxspec + 32*maxspec**2 + 13*eg_nodes &
                                 + 30*eg_nodes*maxspec + 5*eg_nodes*maxspec**2

!     DVODE workspace requirements      
      integer, parameter ::  dvl = 0, dvc = 0
      integer, parameter ::  dvr = 22 + 9*(maxspec+1) + 2*(maxspec+1)**2
      integer, parameter ::  dvi = 30 + maxspec + 1
      
!     DVODE driver (for transient ODE integration) workspace requirements
      integer, parameter ::dvdi = 0, dvdl = 0, dvdc = 0
      integer, parameter :: dvdr = 2*(maxspec+1)

!     Workspace memory layout
      integer, parameter :: ckbr=1,ckbi=1,ckbl=1,ckbc=1 
      integer, parameter :: cker=ckbr+ckr-1,ckei=ckbi+cki-1, &
                 ckel=ckbl+ckl-1,ckec=ckbc+ckc-1
      integer, parameter :: mcbr=cker+1, mcbi=ckei+1, mcbl=ckel+1, mcbc=ckec+1
      integer, parameter :: mcer =mcbr+mcr-1,  mcei=mcbi+mci-1, &
                 mcel =mcbl+mcl-1,  mcec=mcbc+mcc-1
      integer, parameter :: egbr=mcer+1, egbi=mcei+1, egbl=mcel+1, egbc=mcec+1 
      integer, parameter :: eger =egbr+egr-1,  egei=egbi+egi-1, &
                 egel =egbl+egl-1,  egec=egbc+egc-1
      integer, parameter :: dvbr=eger+1, dvbi=egei+1, dvbl=egel+1, dvbc=egec+1 
      integer, parameter :: dver =dvbr+dvr-1,  dvei=dvbi+dvi-1, &
                 dvel =dvbl+dvl-1,  dvec=dvbc+dvc-1
      integer, parameter :: dvdbr=dver+1,dvdbi=dvei+1,dvdbl=dvel+1,dvdbc=dvec+1 
      integer, parameter :: dvder=dvdbr+dvdr-1,dvdei=dvdbi+dvdi-1, &
                 dvdel=dvdbl+dvdl-1,dvdec=dvdbc+dvdc-1
      
      ! Total workspace requirements -- add 1 to guard against null arrays.
      integer, parameter :: lenr_fixed = ckr + mcr + egr + dvr + dvdr + 1
      integer, parameter :: leni_fixed = cki + mci + egi + dvi + dvdi + 1
      integer, parameter :: lenl_fixed = ckl + mcl + egl + dvl + dvdl + 1
      integer, parameter :: lenc_fixed = ckc + mcc + egc + dvc + dvdc + 1
      
      double precision, save :: RWRK(lenr_fixed)
      integer         , save :: IWRK(leni_fixed)
      logical         , save :: LWRK(lenl_fixed)
      character*(maxspnml), save :: CWRK(lenc_fixed)
      
      ! Actual sizes
      integer, save :: Nelt, Nspec, Nreac, Nfit

      ! I/O
      integer          LLINKMC
      common / ckdio / LLINKMC
      save /ckdio/

      integer, save ::  verbose_vode, max_vode_subcycles

      ! Integrator
      double precision, save ::  spec_scalY(maxspec)

      double precision, save :: thickFacCH
            
      ! Transport library
      double precision, save ::  TMIN_TRANS

      integer, save ::  Ncoefs,Doffset,PTCoffset,DMIXoffset,TDoffset

      end module cdwrk_module
