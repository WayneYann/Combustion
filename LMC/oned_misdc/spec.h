
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
      save /speci/

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

      double precision EGRWRK(egr), invmwt(maxspec), mwt(maxspec),
     &     P1ATM, TMIN_TRANS, min_vode_timestep
      common / tranr / EGRWRK, invmwt, mwt, P1ATM, TMIN_TRANS,
     &     min_vode_timestep
      save /tranr/

      character*256 tranfile
      common / tranc / tranfile
      save /tranc/

c     DVODE driver stuff
      integer nchemdiag, dvd_debug
      common / dvdi / nchemdiag, dvd_debug
      save /dvdi/

      double precision c_0(0:maxspec), c_1(0:maxspec), rhoh_INIT,
     &     hmix_TYP
      common / dvdr / c_0, c_1, rhoh_INIT, hmix_TYP
      save /dvdr/

c     LMC alg stuff
      integer probtype, misdc_iterMAX,on_lo,on_hi,max_order,
     &     divu_ceiling_flag, predict_temp_for_coeffs
      parameter (on_lo = 0, on_hi = 1, max_order = 3)
      common / lmci / probtype, misdc_iterMAX, divu_ceiling_flag,
     &     predict_temp_for_coeffs
      save /lmci/

      double precision dpdt_factor, Pcgs, T_bc(0:1), rho_bc(0:1),
     &     Y_bc(maxspec,0:1), h_bc(0:1), u_bc(0:1), flame_offset, 
     &     coef_avg_harm, Pr, Sc, thickFacTR, thickFacCH,
     &     rho_divu_ceiling, divu_dt_factor
      common / lmcr / dpdt_factor, Pcgs, T_bc, rho_bc,
     &     Y_bc, h_bc, u_bc, flame_offset,
     &     coef_avg_harm, Pr, Sc, thickFacTR, thickFacCH,
     &     rho_divu_ceiling, divu_dt_factor
      save /lmcr/
