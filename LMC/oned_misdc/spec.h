
c     Chem species, etc
      integer maxreac, maxspec, maxelts, maxthrdb, maxtp, maxsp,
     &     maxspnml
      parameter (maxreac = 100, maxspec=30, maxelts=3,
     &     maxthrdb=10, maxtp=3, maxsp=12, maxspnml = 16)

      integer iH2, iO2, iCH4, iN2, specNameLen
      common / speci / iH2, iO2, iCH4, iN2, specNameLen

      character*(maxspnml) specNames(maxspec)
      common / specc / specNames

c     EGLib stuff
      integer eg_nodes, eg_IFLAG, eg_ITLS, egr, egi, LLINKMC
      parameter (eg_nodes = 1)
      parameter (eg_IFLAG = 7)
      parameter (eg_ITLS  = 3)
      parameter ( egr = 23 + 14*maxspec + 32*maxspec**2 + 13*eg_nodes
     &     + 30*eg_nodes*maxspec + 5*eg_nodes*maxspec**2)
      parameter (egi = maxspec)
      parameter (LLINKMC = 56)
      
      integer traninit, EGIWRK(egi), Nelt, Nspec, Nreac, Nfit,
     &     LeEQ1, max_vode_subcycles, verbose_vode
      common / trani / EGIWRK, traninit, Nelt, Nspec, Nreac, Nfit,
     &     LeEQ1, max_vode_subcycles, verbose_vode

      double precision EGRWRK(egr), invmwt(maxspec), mwt(maxspec),
     &     P1ATM, TMIN_TRANS, Pr, Sc, thickFacTR, thickFacCH,
     &     min_vode_timestep
      common / tranr / EGRWRK, invmwt, mwt, P1ATM, TMIN_TRANS, Pr, Sc,
     &     thickFacTR, thickFacCH, min_vode_timestep

      character*256 tranfile
      common / tranc / tranfile


c     DVODE driver stuff
      integer nchemdiag
      common / dvdi / nchemdiag

      double precision Pcgs_dvd
      common / dvdr / Pcgs_dvd

