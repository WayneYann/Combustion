 &fortin
  probtype = 3
  flct_file = "swr_27830-27990.bin"
  forceInflow = .TRUE.
  nCompInflow = 3
  zstandoff = -0.0191
  zBL = .008
  zBL = .004
  vheight = 000.0

  Vin = 1.0
  tVco_l = 0.0
  tVco_r = 1.0
  Vco_l  = 5.0
  Vco_r  = 5.0

  turb_scale=1.

  flametracval = 1.e-9
  max_temp_lev = 0
  temperr = 100
  tempgrad = 300.0
  max_vort_lev = 0
  refine_stick = 0
  refine_nozzle = 1
  max_nozzle_lev = 1
  refine_nozzle_x = 0.03
  refine_nozzle_z = 0.12
  refine_nozzle_z = 0.04

  swK = 0.0
  dBL = .0001
  anisotsc = 1.
  wallTh = 0.00375
  Ro = .025
  Rf = .025
  stTh = -0.010
  rhot = .01
 /
 &heattransin
  pamb = 101325.
  dpdt_factor = .3
 /
