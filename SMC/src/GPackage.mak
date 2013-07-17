f90sources += main.f90

f90sources += advance.f90
f90sources += build_info.f90
f90sources += checkpoint.f90
f90sources += cputime.f90
f90sources += derivative_stencil.f90
f90sources += initialize.f90

ifdef EXPAND
  f90sources += kernels_exp.f90
else
  f90sources += kernels.f90
endif

f90sources += kernels_2d.f90

f90sources += make_plot_variables.f90
f90sources += make_plotfile.f90
f90sources += nscbc.f90
f90sources += physbndry_reg.f90
f90sources += probin.f90
f90sources += sdcquad.f90
f90sources += smc.f90
f90sources += smc_bc.f90
f90sources += smcdata.f90
f90sources += threadbox.f90
f90sources += time.f90

ifdef CONVERGENCE
  f90sources += trans_prop_conv.f90
else
  f90sources += transport_properties.f90
endif

f90sources += variables.f90
f90sources += write_job_info.f90
