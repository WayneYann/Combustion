      
      REAL_T vorterr, temperr, adverr, tempgrad, flametracval, twall
      REAL_T splitx, splity, traceSpecVal, xfrontw, yfrontw
      REAL_T domnlo(SDIM), domnhi(SDIM)
      
      integer probtype, max_vort_lev, max_temp_lev, max_trac_lev, max_nozzle_lev
      integer fuelID, oxidID, prodID, nspecies

      common /prob1/ vorterr, temperr, adverr, tempgrad, flametracval, twall,
     &               splitx, splity, traceSpecVal, xfrontw, yfrontw
      common /prob2/ domnlo, domnhi
      common /prob3/ probtype, max_vort_lev, max_temp_lev, max_trac_lev,
     &               max_nozzle_lev
      common /prob4/ fuelID, oxidID, prodID, nspecies
      save /prob1/, /prob2/, /prob3/, /prob4/

      integer refine_nozzle
      REAL_T refine_nozzle_x, refine_nozzle_y, blobx, bloby, blobr, xcen, y_offset
      common /myprob1/ refine_nozzle
      common /myprob2/ refine_nozzle_x, refine_nozzle_y, blobx, bloby,
     &     blobr, xcen, y_offset

      REAL_T v_strength,v_width,v_xcen,v_ycen, v_cl_x
      common /myprob3/ v_strength,v_width,v_xcen,v_ycen,v_cl_x

      REAL_T v_blob_r, v_blob_T, v_blob_airfrac, turb_scale,phi_in
      common /myprob4/ v_blob_r, v_blob_T, v_blob_airfrac, turb_scale,phi_in

      REAL_T stTh, Rf, V_in, V_co, T_in, R_hot, R_hotBL, stBL, pipeTh, pipeBL
      REAL_T tV_in_l, tV_in_r, V_in_new, tV_co_l, tV_co_r, V_co_new
      REAL_T T_stick, rho_stick, h_stick
      common /myprob5/ stTh, Rf, V_in, V_co, T_in, R_hot, R_hotBL, stBL,
     &     pipeTh, pipeBL, tV_in_l, tV_in_r, V_in_new, tV_co_l, tV_co_r,
     &     V_co_new, T_stick, rho_stick, h_stick
      save /myprob1/, /myprob2/, /myprob3/, /myprob4/, /myprob5/ 

      logical forceInflow
      REAL_T probSizeFile(3), dxFile(3)
      integer dimFile
      common /inflow1/forceInflow
      common /inflow2/probSizeFile, dxFile
      common /inflow3/dimFile
      save /inflow1/, /inflow2/, /inflow3/

      integer waves
      REAL_T amplitude, coefficient(2, 100)
      common / wavesI / waves
      common / wavesR / amplitude, coefficient
      save /wavesI/, /wavesR/
