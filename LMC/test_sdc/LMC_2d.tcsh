# CEG 6/15/2009
#
# This script run LMC with sdc 
#  varying # sdc iterations and dx for convergence studies


# calc does integer calculations, calcfx does double precision
alias calc 'awk "BEGIN{ print \!* }" '
alias calcfx ' awk -v CONVFMT="%12.2f" -v OFMT="%.9g" "BEGIN{  print \!* }" '

foreach DX ( 256 128 64 32 )
#foreach DX ( 64 )
foreach SDCITER ( 0 1 2 )

set FIXEDDT=`calcfx 0.0000032*(256/$DX)`
set NSTEP=`calcfx 0.000256/$FIXEDDT`

echo "dx is " $DX
echo "sdciter is " $SDCITER
echo "fixed_dt is " $FIXEDDT
echo "max_step is " $NSTEP

cat > input << EOF
ht.use_sdc = 1
ht.sdc_iters = $SDCITER

max_step  = $NSTEP
#stop_time = 0.000512 

ns.fixed_dt = $FIXEDDT 

amr.max_grid_size   = 64
amr.n_cell = $DX $DX

#ns.hack_nochem = 1

amr.plot_int = 10

ht.pltfile = plt00050
ht.plot_auxDiags = 0
ht.plot_rhoY=0
ht.plot_molefrac=1
ht.plot_ydot=0
ht.v = 0

mg.usecg = 2
mg.maxiter = 1000
mg.nu_0 = 1
mg.nu_1 = 4
mg.nu_2 = 4
mg.nu_f = 40
mg.v = 0
cg.maxiter = 1000
cg.v = 0

amr.max_level = 0
amr.ref_ratio       = 2 2 2 2   # refinement ratio
amr.regrid_int      = 2       # how often to regrid
amr.n_error_buf     = 2 2 2 2 # number of buffer cells in error est
amr.grid_eff        = 0.7     # what constitutes an efficient grid
amr.blocking_factor = 16      # block factor in grid generation
amr.check_file      = chk
amr.check_int       = 100     # number of timesteps between checkpoints
amr.plot_file       = plt
amr.grid_log        = grdlog  # name of grid logging file
amr.derive_plot_vars=rhoRT mag_vort
amr.probin_file = probin-control # This will default to file "probin" if not set

# ------------------  INPUTS TO PHYSICS CLASS -------------------
#ns.cfl            = 0.75      # cfl number for hyperbolic system
ns.cfl            = 0.5      # cfl number for hyperbolic system
ns.init_shrink    = 1.0        # scale down initial timestep
ns.change_max     = 1.0       # scale back initial timestep

ns.v   = 0

ns.init_iter      = 0        # number of init iters to def pressure
ns.num_divu_iters = 0

ns.dt_cutoff       = 5.e-10   # level 0 timestep below which we halt
ns.visc_tol        = 1.0e-14  # tolerence for viscous solves
ns.visc_abs_tol    = 1.0e-14  # tolerence for viscous solves
ns.vel_visc_coef   = 1.983e-5
ns.temp_cond_coef  = 2.6091e-5
ns.scal_diff_coefs = -0.01
ns.variable_vel_visc  = 1
ns.variable_scal_diff = 1
ns.gravity        = 0        # body force  (gravity in MKS units)
ns.sum_interval   = 1        # timesteps between computing mass
ns.do_reflux      = 1        # 1 => do refluxing
ns.do_mac_proj    = 1        # 1 => do MAC projection
ns.do_sync_proj   = 1        # 1 => do Sync Project
ns.do_MLsync_proj = 1
ns.divu_relax_factor   = 0.0
ns.be_cn_theta = 0.5
ns.S_in_vel_diffusion = 1
ns.do_temp = 1
ns.do_diffuse_sync = 1
ns.do_reflux_visc  = 1
ns.divu_ceiling = 1
ns.divu_dt_factor = .4
ns.min_rho_divu_ceiling = .01
# FIXME
#ns.unity_Le = 1
ns.unity_Le = 0
#ns.tranfile        = ./tran.asc.chem-H
#ns.fuelName        = H2
#ns.flameTracName   = H
#ns.consumptionName = H2
ns.oxidizerName    = O2

ns.tranfile        = ./tran.asc.CH4-2step
ns.fuelName        = CH4
ns.consumptionName = CH4
ns.flameTracName   = CO

ns.plot_consumption = 1
ns.plot_heat_release = 1
ns.do_active_control = false
ns.do_init_vort_proj = 0

# ----------------  PROBLEM DEPENDENT INPUTS
ns.lo_bc          = 0 1
ns.hi_bc          = 0 2

geometry.is_periodic = 1 0
geometry.coord_sys = 0  # 0 => cart, 1 => RZ
geometry.prob_lo   =   0.  0. # m
#geometry.prob_hi   =  .004 .004
geometry.prob_hi   =  .008 .008

# >>>>>>>>>>>>>  BC FLAGS <<<<<<<<<<<<<<<<
# 0 = Interior           3 = Symmetry
# 1 = Inflow             4 = SlipWall
# 2 = Outflow            5 = NoSlipWall


# ------------------  INPUTS TO GODUNOV CLASS ----------------------
godunov.slope_order = 4
godunov.use_unlimited_slopes = 1

# ------------------  INPUTS TO DIFFUSION CLASS --------------------
diffuse.use_cg_solve = 0
diffuse.max_order = 4
diffuse.tensor_max_order = 4
diffuse.use_tensor_cg_solve = 0
diffuse.v = 0
diffuse.Rhs_in_abs_tol = 1

# ------------------  INPUTS TO PROJECTION CLASS -------------------
proj.proj_tol       = 1.0e-10  # tolerence for projections
proj.sync_tol       = 1.0e-11 # tolerence for projections
proj.rho_wgt_vel_proj = 0      # 0 => const den proj, 1 => rho weighted
proj.Pcode          = 1
proj.filter_factor  = 0.0
proj.do_outflow_bcs = 1
proj.divu_minus_s_factor = 0.
proj.proj_2 = 1
proj.add_vort_proj = 0
proj.v = 0

# ------------------  INPUTS TO MACPROJ CLASS -------------------
mac.mac_tol        = 1.0e-12  # tolerence for mac projections
mac.mac_sync_tol   = 1.0e-9   # tolerence for mac SYNC projection
mac.mac_abs_tol    = 1.0e-14
mac.use_cg_solve   = 0
mac.do_outflow_bcs = 1

#
# Select form of FAB output: default is NATIVE
#
#   ASCII  (this is very slow)
#   NATIVE (native binary form on machine -- the default)
#   IEEE32 (useful if you want 32bit files when running in double precision)
#   8BIT   (eight-bit run-length-encoded)
#
fab.format = NATIVE

#
# Initializes distribution strategy from ParmParse.
#
# ParmParse options are:
#
#   DistributionMapping.strategy = ROUNDROBIN
#   DistributionMapping.strategy = KNAPSACK
#
# The default strategy is ROUNDROBIN.
#
DistributionMapping.strategy = KNAPSACK
DistributionMapping.do_not_minimize_comm_costs=0
DistributionMapping.verbose=0

EOF

mpiexec -n 2 main2d.Linux.Intel.Intel.MPI.ex input > output_${DX}_${SDCITER}

if ( $DX == 32 ) then
 mv plt00010 convergence/plt${SDCITER}032
 mv chk00010 convergence/chk${SDCITER}032
endif
if ( $DX == 64 ) then
 mv plt00020 convergence/plt${SDCITER}064
 mv chk00020 convergence/chk${SDCITER}064
endif
if ( $DX == 128 ) then
 mv plt00040 convergence/plt${SDCITER}128
 mv chk00040 convergence/chk${SDCITER}128
endif
if ( $DX == 256 ) then
 mv plt00080 convergence/plt${SDCITER}256
 mv chk00080 convergence/chk${SDCITER}256
endif

end
end
