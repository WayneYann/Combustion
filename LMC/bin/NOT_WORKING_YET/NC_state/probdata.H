#include "cdwrk.H"
      
      REAL_T vorterr, temperr, adverr, tempgrad, flametracval
      common /prob1/ vorterr, temperr, adverr, tempgrad, flametracval

      REAL_T domnlo(SDIM), domnhi(SDIM)
      common /prob2/ domnlo, domnhi
      
      integer probtype, max_vort_lev, max_temp_lev, max_trac_lev
      common /prob3/ probtype, max_vort_lev, max_temp_lev, max_trac_lev

      integer fuelID, oxidID, prodID, nspecies
      common /prob4/ fuelID, oxidID, prodID, nspecies

      save /prob1/, /prob2/, /prob3/, /prob4/


      REAL_T X_comp(maxspec), Y_comp(maxspec), phi_in, X_nitrogen,rshift
      common /myconc1/ X_comp, Y_comp, phi_in, X_nitrogen,rshift
      integer iH2, iO2
      common /myconc2/ iH2, iO2

      save /myconc1/, /myconc2/

      REAL_T v_strength,v_width,v_xcen,v_ycen, v_cl_x, x_hot,y_hot,r_hot,hot_wid,T_hot,T_amb
      common /myprob3/ v_strength,v_width,v_xcen,v_ycen,v_cl_x, x_hot,y_hot,r_hot,hot_wid,T_hot,T_amb

      save /myprob3/


      logical bcinit
      common / stats/ bcinit

      integer Nzones
      parameter (Nzones=5)      

      REAL_T u_bc(Nzones), v_bc(Nzones), rho_bc(Nzones)
      REAL_T Y_bc(0:maxspec-1, Nzones), T_bc(Nzones)
      REAL_T h_bc(Nzones)
      
      common / bcs / u_bc, v_bc, rho_bc, Y_bc, T_bc, h_bc
