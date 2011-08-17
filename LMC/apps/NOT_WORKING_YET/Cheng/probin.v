 &fortin
  probtype = 3
  flct_file = ""
  flct_file = "Turb_plt0620"
  forceInflow = .FALSE.
  forceInflow = .TRUE.
  zstandoff = -0.03
  zBL = .008
  vheight = 0.0
  numInflPlanesStore = 16

  phiin = 0.7

  Vco_l  = 1.5
  Vco_r  = 1.5
  tVco_l = 0.0
  tVco_r = 1.0

  rhot = .025
  Vin = 3.9

  turb_scale=1.5

  flametracval = 1.e-9
  max_trac_lev = 3
  max_temp_lev = 0
  temperr = 100
  max_vort_lev = 0
  refine_stick = 1
  max_stick_lev = 3
  refine_stick_x = 0.003
  refine_stick_z = 0.003
  refine_nozzle = 1
  max_nozzle_lev = 2
  refine_nozzle_x = 0.03
  refine_nozzle_z = 0.07
  tempgrad = 300.0

  Tstick = 1000.0
  stL = .04
  stLw = .003

  swK = 0.0
  dBL = .17
  anisotsc = 1.0
  stBL = .002
  wallTh = 0.00375
  Ro = .025
  Rf = .025
  stTh = 0.002
  rhot = .006
 /
 &heattransin
  pamb = 101325.
  dpdt_factor = .3
 /
