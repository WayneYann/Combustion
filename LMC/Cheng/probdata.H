      
      REAL_T vorterr, temperr, adverr, tempgrad, flametracval,
     &               traceSpecVal
      REAL_T domnlo(SDIM), domnhi(SDIM)

      integer probtype, max_vort_lev, max_temp_lev
      integer fuelID, oxidID, prodID, nspecies
      character*1 direction

      common /prob1/ vorterr, temperr, adverr, tempgrad, flametracval,
     &               traceSpecVal
      common /prob2/ domnlo, domnhi
      common /prob3/ probtype, max_vort_lev, max_temp_lev
      common /prob4/ direction
      save /prob1/,/prob2/,/prob3/,/prob4/

      integer refine_nozzle, max_nozzle_lev,
     &         refine_stick, max_stick_lev
      REAL_T refine_nozzle_x, refine_nozzle_y, refine_nozzle_z, 
     &         nozzle_width, refine_stick_x, refine_stick_z
      common /myprob1/ refine_nozzle, max_nozzle_lev,
     &         refine_stick, max_stick_lev
      common /myprob2/ refine_nozzle_x, refine_nozzle_y, 
     &                 refine_nozzle_z, nozzle_width,
     &                 refine_stick_x, refine_stick_z

      save /myprob1/,/myprob2/

      logical add_turb
      REAL_T zturb_top, zturb_bot, turb_scale, anisotsc
      common /myprob4/ zturb_top, zturb_bot, turb_scale, anisotsc, add_turb
      save /myprob4/

      logical forceInflow
      REAL_T probSizeFile(3)
      common /inflow1/forceInflow
      common /inflow2/probSizeFile
      save /inflow1/, /inflow2/

      REAL_T Rf, wallTh, Ro, stTh, Vin, zBL,
     &            dBL, stBL, Thot, rhot, zhot, zstandoff, swK, swW,
     &            hotspot_r, hotspot_z, hotspot_dz, hotspot_T, hotspot_htc,
     &            hotspot_w, vheight, Vco_l, Vco_r, tVco_l, tVco_r
      common /myprob3/ Rf, wallTh, Ro, stTh, Vin, zBL,
     &            dBL, stBL, Thot, rhot, zhot, zstandoff, swK, swW,
     &            hotspot_r, hotspot_z, hotspot_dz, hotspot_T, hotspot_htc,
     &            hotspot_w, vheight, Vco_l, Vco_r, tVco_l, tVco_r
      save /myprob3/
