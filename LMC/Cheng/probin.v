 &fortin
  probtype = 3
  flct_file = ""
  flct_file = "../Turb.c0100.bin"
  forceInflow = .FALSE.
  forceInflow = .TRUE.
  zstandoff = -0.03
  zBL = .008
  vheight = 0.0

  Vco_l = 0.03
  Vco_r = 0.03
  tVco_l = 0.0
  tVco_r = 0.0

  rhot = .025
  Vin = 3.0
  Vin = 1.5
  Vin = 5.0

  turb_scale=1.667

  flametracval = 1.e-9
  max_temp_lev = 0
  temperr = 100
  max_vort_lev = 0
  refine_stick = 1
  max_stick_lev = 2
  refine_stick_x = 0.003
  refine_stick_z = 0.003
  refine_nozzle = 1
  max_nozzle_lev = 1
  refine_nozzle_x = 0.06
  refine_nozzle_z = 0.12
  refine_nozzle_z = 0.07
  tempgrad = 300.0

  swK = 0.0
  stBL = .001
  wallTh = 0.0
  Ro = .05
  stTh = 0.001
  stTh = 0.003
 /
 &heattransin
  pamb = 101325.
  dpdt_factor = .3
 /
