
c     Chem species, etc
      integer maxreac, maxspec, maxelts, maxthrdb, maxtp, maxsp,
     &     maxspnml
      parameter (maxreac = 100, maxspec=30, maxelts=3,
     &     maxthrdb=10, maxtp=3, maxsp=12, maxspnml = 16)

      integer Nelt, Nspec, Nreac, Nfit, iH2, iO2, iCH4,
     &     iN2, specNameLen, nscal, Density, Temp, RhoH, 
     &     RhoRT, FirstSpec, LastSpec
      common / speci / Nelt, Nspec, Nreac, Nfit, iH2, iO2, iCH4,
     &     iN2, specNameLen, nscal, Density, Temp, RhoH, 
     &     RhoRT, FirstSpec, LastSpec

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
      
      integer traninit, EGIWRK(egi), LeEQ1,
     &     max_vode_subcycles, verbose_vode
      common / trani / traninit, EGIWRK, LeEQ1,
     &     max_vode_subcycles, verbose_vode

      double precision EGRWRK(egr), invmwt(maxspec), mwt(maxspec),
     &     P1ATM, TMIN_TRANS, min_vode_timestep
      common / tranr / EGRWRK, invmwt, mwt, P1ATM, TMIN_TRANS,
     &     min_vode_timestep

      character*256 tranfile
      common / tranc / tranfile


c     DVODE driver stuff
      integer nchemdiag
      common / dvdi / nchemdiag

      double precision c_0(maxspec+1), c_1(maxspec+1)
      common / dvdr / c_0, c_1

c     LMC alg stuff
      integer probtype, misdc_iterMAX
      common / lmci / probtype, misdc_iterMAX

      double precision dpdt_factor, Pcgs, T_bc(0:1), rho_bc(0:1),
     &     Y_bc(maxspec,0:1), h_bc(0:1), u_bc(0:1), flame_offset, 
     &     coef_avg_harm, Pr, Sc, thickFacTR, thickFacCH
      common / lmcr / dpdt_factor, Pcgs, T_bc, rho_bc,
     &     Y_bc, h_bc, u_bc, flame_offset,
     &     coef_avg_harm, Pr, Sc, thickFacTR, thickFacCH
