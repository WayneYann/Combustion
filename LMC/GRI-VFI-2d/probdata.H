      
      REAL_T vorterr, temperr, adverr, tempgrad, flametracval, twall
      REAL_T splitx, splity, traceSpecVal, xfrontw, yfrontw
      REAL_T domnlo(SDIM), domnhi(SDIM)
      
      integer probtype, max_vort_lev, max_temp_lev
      integer fuelID, oxidID, prodID, nspecies

      common /prob1/ vorterr, temperr, adverr, tempgrad, flametracval, twall,
     &               splitx, splity, traceSpecVal, xfrontw, yfrontw
      common /prob2/ domnlo, domnhi
      common /prob3/ probtype, max_vort_lev, max_temp_lev

      integer refine_nozzle
      REAL_T refine_nozzle_x, refine_nozzle_y, nozzle_width, blobx, bloby,
     &     blobr, xcen, y_offset
      common /myprob1/ refine_nozzle
      common /myprob2/ refine_nozzle_x, refine_nozzle_y, nozzle_width,
     &     blobx, bloby, blobr, xcen, y_offset

      REAL_T v_strength,v_width,v_xcen,v_ycen, v_cl_x
      common /myprob3/ v_strength,v_width,v_xcen,v_ycen,v_cl_x

      REAL_T v_blob_r, v_blob_T, v_blob_airfrac
      common /myprob4/ v_blob_r, v_blob_T, v_blob_airfrac
