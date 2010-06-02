      
      REAL_T V_fu, V_ox, T_fu, T_ox, R_fu, R_hot, T_hot, T_wall
      REAL_T y_hot_lo, y_hot_hi, temperr, vorterr, tempgrad
      REAL_T flametracval, nozzle_height, nozzle_radius
      REAL_T domnlo(SDIM), domnhi(SDIM), V_in_rampdt

      integer max_temp_lev, max_vort_lev, max_trac_lev, max_nozzle_lev
      integer probtype, fuelID, oxidID, prodID, noInletSpecDiff

      common /prob1/ V_fu, V_ox, T_fu, T_ox, R_fu, R_hot, T_hot, T_wall
      common /prob2/ y_hot_lo, y_hot_hi, temperr, vorterr, tempgrad
      common /prob3/ flametracval, nozzle_height, nozzle_radius
      common /prob4/ domnlo, domnhi, V_in_rampdt
      common /prob5/ max_temp_lev, max_vort_lev, max_trac_lev, max_nozzle_lev
      common /prob6/ probtype, fuelID, oxidID, prodID, noInletSpecDiff

      save /prob1/, /prob2/, /prob3/, /prob4/, /prob5/, /prob6/
