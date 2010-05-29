      
      REAL_T vorterr, temperr, adverr, tempgrad, flametracval
      REAL_T domnlo(SDIM), domnhi(SDIM)
      REAL_T nozzlerefwidth, nozzlerefheight, nozzlerefloc
      
      integer probtype, max_vort_lev, max_temp_lev
      integer maximum_nozzle_level
      integer fuelID, oxidID, prodID
      integer refine_nozzle, noInletSpecDiff

      common /prob1/ vorterr, temperr, adverr, tempgrad, flametracval
      common /prob2/ domnlo, domnhi
      common /prob3/ probtype, max_vort_lev, max_temp_lev, maximum_nozzle_level
      common /prob4/ nozzlerefwidth,nozzlerefheight,nozzlerefloc
      common /prob5/ refine_nozzle, noInletSpecDiff

      REAL_T flowrate_fu_AR
      REAL_T flowrate_fu_CH4
      REAL_T flowrate_fu_H2
      REAL_T flowrate_fu_H2O
      REAL_T flowrate_fu_N2
      REAL_T flowrate_fu_NH3
      REAL_T flowrate_fu_NO
      REAL_T flowrate_fu_O2
      REAL_T flowrate_ox_AR
      REAL_T flowrate_ox_H2O
      REAL_T flowrate_ox_N2
      REAL_T flowrate_ox_O2
      REAL_T T_ox, R_ox_1, R_ox_2, T_fu, R_fu
      REAL_T x_hot, y_hot_lo, y_hot_hi, T_hot
      REAL_T a_ox, b_ox, c_ox, Trad_bg, time_rad, dt_rad
      REAL_T y_wall2, y_wall3, T_wall1, T_wall2, T_wall3

      common / myprob1 / 
     &   flowrate_fu_AR, flowrate_fu_CH4, flowrate_fu_H2, flowrate_fu_H2O,
     &   flowrate_fu_N2, flowrate_fu_NH3, flowrate_fu_NO, flowrate_fu_O2
      common / myprob2 / 
     &   flowrate_ox_AR, flowrate_ox_H2O, flowrate_ox_N2, flowrate_ox_O2
      common / myprob3 / T_ox, R_ox_1, R_ox_2, T_fu, R_fu
      common / myprob4 / x_hot, y_hot_lo, y_hot_hi, T_hot
      common / myprob5 / a_ox, b_ox, c_ox, Trad_bg, time_rad, dt_rad
      common / myprob6 / y_wall2, y_wall3, T_wall1, T_wall2, T_wall3
      
      save /prob1/, /prob2/, /prob3/, /prob4/, /prob5/, /myprob1/, /myprob2/,
     &      /myprob3/,/myprob4/, /myprob5/, /myprob6/
