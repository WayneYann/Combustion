f90sources += chemistry_module.f90 vode_module.f90 tranlib_module.f90

ifndef USE_EGZ
   f90sources += eglib_module.f90
endif

