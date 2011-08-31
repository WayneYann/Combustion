      
      REAL_T vorterr, temperr, adverr, tempgrad, flametracval, twall
      REAL_T splitx, splity, splitz, traceSpecVal
      REAL_T xfrontw, yfrontw, zfrontw
      REAL_T domnlo(SDIM), domnhi(SDIM)

      integer probtype, max_vort_lev, max_temp_lev
      integer fuelID, oxidID, prodID, nspecies
      character*1 direction

      common /prob1/ vorterr, temperr, adverr, tempgrad, flametracval, twall,
     &               splitx, splity, splitz, traceSpecVal, 
     &               xfrontw, yfrontw, zfrontw
      common /prob2/ domnlo, domnhi
      common /prob3/ probtype, max_vort_lev, max_temp_lev
      common /prob4/ direction
      save /prob1/,/prob2/,/prob3/,/prob4/

      integer refine_nozzle
      REAL_T refine_nozzle_x, refine_nozzle_y, refine_nozzle_z, 
     &         nozzle_width, blobx, bloby, blobr, blobz, xcen
      common /myprob1/ refine_nozzle
      common /myprob2/ refine_nozzle_x, refine_nozzle_y, 
     &                 refine_nozzle_z, nozzle_width,
     &                 blobx, bloby, blobr, blobz, xcen      

      REAL_T v_strength,v_width,v_xcen,v_ycen,v_zcen,v_cl_x, vel_frame,
     &  pmf_offset
      common /myprob3/ v_strength,v_width,v_xcen,v_ycen,v_zcen,v_cl_x,
     &  vel_frame,pmf_offset

      save /myprob1/,/myprob2/,/myprob3/

      logical add_turb
      REAL_T zturb_top, zturb_bot, turb_scale
      common /myprob4/ zturb_top, zturb_bot, turb_scale, add_turb
      save /myprob4/

      REAL_T T_in, XH2_in, XO2_in, T_out, XH2_out, XO2_out,
     &     blobrc, blobs, turbr, turbw
      common /myprob5/ T_in, XH2_in, XO2_in, T_out, XH2_out, XO2_out,
     &     blobrc, blobs, turbr, turbw
      save /myprob5/

      logical forceInflow
      REAL_T probSizeFile(3)
      common /inflow1/forceInflow
      common /inflow2/probSizeFile
      save /inflow1/, /inflow2/

      integer lmodemin, lmodemax, lpeak
      parameter (lmodemin=14, lmodemax=20, lpeak=17)
      REAL_T alphalm(0:lmodemax,lmodemin:lmodemax)
      REAL_T betalm(0:lmodemax,lmodemin:lmodemax)
      REAL_T gammalm(0:lmodemax,lmodemin:lmodemax)
      common /YlmRand/ alphalm, betalm, gammalm
      
