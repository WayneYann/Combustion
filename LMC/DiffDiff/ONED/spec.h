
c     Chem species, etc
      integer maxreac, maxspec, maxelts, maxthrdb, maxtp, maxsp,
     &     maxspnml
      parameter (maxreac = 100, maxspec=30, maxelts=6,
     &     maxthrdb=10, maxtp=3, maxsp=12, maxspnml = 16)

c     nscal: room for rho, rhoH, Temp + species (rho.Y)
      integer maxscal, nx
c      parameter (maxscal = maxspec + 3, nx = 256)
      parameter (maxscal = maxspec + 3, nx = 64)

      integer Nelt, Nspec, Nreac, Nfit, iH2, iO2, iCH4,
     &     iN2, specNameLen, Density, Temp, RhoH, 
     &     RhoRT, FirstSpec, LastSpec, nscal
      logical nochem_hack
      common / speci / Nelt, Nspec, Nreac, Nfit, iH2, iO2, iCH4,
     &     iN2, specNameLen, Density, Temp, RhoH, 
     &     RhoRT, FirstSpec, LastSpec, nscal
      common / specL / nochem_hack
      save /speci/, /specL/

      character*(maxspnml) specNames(maxspec)
      common / specc / specNames
      save /specc/

c     EGLib stuff
      integer eg_nodes, eg_IFLAG, eg_ITLS, egr, egi, LLINKMC
      parameter (eg_nodes = 1)
      parameter (eg_IFLAG = 7)
      parameter (eg_ITLS  = 3)
      parameter ( egr = 23 + 14*maxspec + 32*maxspec**2 + 13*eg_nodes
     &     + 30*eg_nodes*maxspec + 5*eg_nodes*maxspec**2)
      parameter (egi = maxspec)
      parameter (LLINKMC = 56)
      
      integer traninit, EGIWRK(egi), LeEQ1,
     &     max_vode_subcycles, verbose_vode
      common / trani / traninit, EGIWRK, LeEQ1,
     &     max_vode_subcycles, verbose_vode
      save /trani/

      real*8 EGRWRK(egr), invmwt(maxspec), mwt(maxspec),
     &     RU, RUC, P1ATM, TMIN_TRANS, min_vode_timestep
      common / tranr / EGRWRK, invmwt, mwt, RU, RUC, P1ATM, TMIN_TRANS,
     &     min_vode_timestep
      save /tranr/

      character*256 tranfile
      common / tranc / tranfile
      save /tranc/

c     DVODE driver stuff
      integer nchemdiag, dvd_debug
      common / dvdi / nchemdiag, dvd_debug
      save /dvdi/

      real*8 c_0(0:maxspec), c_1(0:maxspec), rhoh_INIT,
     &     hmix_TYP
      common / dvdr / c_0, c_1, rhoh_INIT, hmix_TYP
      save /dvdr/


c     Driver stuff
      real*8 Pcgs,errMax,dtRedFac,big,small,smallDt,typVal(maxscal)
      integer NiterMAX,setTfromH,rhoInTrans,advance_RhoH,alt_spec_update
      integer probtype,Ncorrect
      real*8 problo,probhi,flame_offset
      parameter (errMAX=1.d-8)
      parameter (NiterMAX=20)
      parameter (big=1.d30)
      parameter (small=1.d-30)
      parameter (smallDt=1.d-30)
      integer N1dMAX
      parameter (N1dMAX=(maxspec+1)*nx)
      common / drvcomr / Pcgs,dtRedFac,typVal,problo,probhi,flame_offset
      common / drvcomi / setTfromH,rhoInTrans,advance_RhoH,alt_spec_update,
     &     probtype,Ncorrect
      save /drvcomr/, /drvcomi/

c     Dummy arrays for CK calls
      real*8 RWRK
      integer IWRK
      common /ckcomr/ RWRK
      common /ckcomi/ IWRK


