c     Note: Requires ckdwrk.H to obtain the value of maxspec
      
      logical bcinit
      common /stats/ bcinit

      integer Nzones
      parameter (Nzones=5)      

      REAL_T u_bc(Nzones), v_bc(Nzones), w_bc(Nzones), rho_bc(Nzones)
      REAL_T Y_bc(0:maxspec-1, Nzones), T_bc(Nzones)
      REAL_T h_bc(Nzones)
      
      common /bcs/ u_bc, v_bc, w_bc, rho_bc, Y_bc, T_bc, h_bc
      
      save /stats/,/bcs/      
