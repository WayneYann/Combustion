c     Note: Requires ckdwrk.H to obtain the value of maxspec
      
      logical bcinit
      common / stats/ bcinit

      integer Nzones
      parameter (Nzones=4)      

      REAL_T u_bc(Nzones), v_bc(Nzones), w_bc(Nzones), rho_bc(Nzones)
      REAL_T Y_bc(0:maxspec-1, Nzones), T_bc(Nzones)
      REAL_T h_bc(Nzones)
      
      common / bcs / u_bc, v_bc, w_bc, rho_bc, Y_bc, T_bc, h_bc

      character(50) ac_hist_file
      REAL_T tau_control, cfix, coft_old, sest, V_in_old, corr,
     &     changeMax_control, tbase_control, dV_control, scale_control,
     &     zbase_control
      common / control1 / tau_control, cfix, coft_old, sest, V_in_old,
     &     corr, changeMax_control, tbase_control, dV_control, scale_control,
     &     zbase_control, ac_hist_file

      save /stats/, /bcs/, /control1/
