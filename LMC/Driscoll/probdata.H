      
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
     &         refine_stick, max_stick_lev, max_trac_lev, min_nogrid_lev
      REAL_T refine_nozzle_x, refine_nozzle_y, refine_nozzle_z, 
     &         nozzle_width, refine_stick_x, refine_stick_z, nogrid_radius
      common /myprob1/ refine_nozzle, max_nozzle_lev,
     &         refine_stick, max_stick_lev, max_trac_lev, min_nogrid_lev
      common /myprob2/ refine_nozzle_x, refine_nozzle_y, 
     &                 refine_nozzle_z, nozzle_width,
     &                 refine_stick_x, refine_stick_z, nogrid_radius

      save /myprob1/,/myprob2/

      REAL_T turb_scale, anisotsc
      common /myprob4/ turb_scale, anisotsc
      save /myprob4/

      logical forceInflow
      REAL_T probSizeFile(3),vheight
      common /inflow1/forceInflow
      common /inflow2/probSizeFile,vheight
      save /inflow1/, /inflow2/

      REAL_T Vin, zBL, dBL, stBL, Thot, xhot, yhot, zhot, dxhot, dyhot,
     &            dzhot, zstandoff, Vco_l, Vco_r, tVco_l, tVco_r, slotL,
     &            slotR, sampleLocFUEL, sampleLocPROD, Tinit, velBL
      common /myprob3/ Vin, zBL, dBL, stBL, Thot, xhot, yhot, zhot,
     &            dxhot, dyhot, dzhot, zstandoff, Vco_l, Vco_r, tVco_l,
     &            tVco_r, slotL, slotR, sampleLocFUEL, sampleLocPROD,
     &            Tinit, velBL
      save /myprob3/

      integer iN2, iCH4, iO2
      common /speccom/ iN2, iCH4, iO2
