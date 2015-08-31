!     Chem species, etc
      integer maxreac, maxspec, maxelts, maxthrdb, maxspnml, maxlev
      parameter (maxreac = 325, maxspec=53, maxelts=6, maxthrdb=10, maxspnml=16, maxlev=3)

!     nscal: room for rho, rhoH, Temp, RhoRT + species (rho.Y)
      integer maxscal
      parameter (maxscal = maxspec + 4)

      integer Nelt, Nspec, Nreac, Nfit, iH2, iO2, iCH4, &
            iCH3OCH3, iCO2, iH2O, iN2, Density, Temp, RhoH, &
            RhoRT, FirstSpec, LastSpec, nscal, nnodes
!     todo: get rid of dodiff
      logical dodiff
      common / speci / Nelt, Nspec, Nreac, Nfit, iH2, iO2, iCH4, &
            iCH3OCH3, iCO2, iH2O, iN2, Density, Temp, RhoH, &
           RhoRT, FirstSpec, LastSpec, nscal, nnodes, dodiff
      save /speci/

      character*(maxspnml) specNames(maxspec)
      common / specc / specNames
      save /specc/

!     EGLib stuff
      integer eg_nodes, eg_IFLAG, eg_ITLS, egr, egi, LLINKMC
      parameter (eg_nodes = 1)
      parameter (eg_IFLAG = 3)
      parameter (eg_ITLS  = 1)
      parameter ( egr = 23 + 14*maxspec + 32*maxspec**2 + 13*eg_nodes &
          + 30*eg_nodes*maxspec + 5*eg_nodes*maxspec**2)
      parameter (egi = maxspec)
      parameter (LLINKMC = 56)
      
      integer traninit, EGIWRK(egi), LeEQ1, max_vode_subcycles, verbose_vode
      common / trani / traninit, EGIWRK, LeEQ1, max_vode_subcycles, verbose_vode
      save /trani/

      double precision EGRWRK(egr), invmwt(maxspec), mwt(maxspec), &
           P1ATM, TMIN_TRANS, min_vode_timestep
      common / tranr / EGRWRK, invmwt, mwt, P1ATM, TMIN_TRANS, &
           min_vode_timestep
      save /tranr/

      character*256 tranfile
      common / tranc / tranfile
      save /tranc/

!     DVODE driver stuff
      integer nchemdiag
      common / dvdi / nchemdiag
      save /dvdi/

      double precision c_0(0:maxspec), c_1(0:maxspec), rhoh_INIT, T_INIT, hmix_TYP
      common / dvdr / c_0, c_1, rhoh_INIT, T_INIT, hmix_TYP
      save /dvdr/

!     LMC alg stuff
      integer misdc_iterMAX,on_lo,on_hi,max_order, &
           divu_ceiling_flag, is_first_initial_iter, doing_init_iters, &
           unlim, lim_rxns, coef_avg_harm, fancy_predictor, &
           initial_S_type, recompute_S, probtype
      parameter (on_lo = 0, on_hi = 1, max_order = 3)
      common / lmci / misdc_iterMAX, divu_ceiling_flag, &
           is_first_initial_iter, doing_init_iters, unlim, lim_rxns, &
           coef_avg_harm, fancy_predictor, &
           initial_S_type, recompute_S, probtype
      save /lmci/

      double precision dpdt_factor, Pcgs, T_bc(0:1), rho_bc(0:1), &
           Y_bc(maxspec,0:1), h_bc(0:1), u_bc(0:1), flame_offset, &
           Pr, Sc, rho_divu_ceiling, divu_dt_factor, V_in, vel_TYP
      common / lmcr / dpdt_factor, Pcgs, T_bc, rho_bc, &
           Y_bc, h_bc, u_bc, flame_offset, &
           Pr, Sc, rho_divu_ceiling, divu_dt_factor, V_in, vel_TYP
      save /lmcr/

      integer nx
      common / amri / nx
      save /amri/

      logical subcycling
      common /amrl / subcycling
      save /amrl/
