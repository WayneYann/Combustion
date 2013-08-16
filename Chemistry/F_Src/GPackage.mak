f90sources += chemistry_module.f90

ifdef USE_EGZ
   f90sources += egz_module.f90
else
   f90sources += eglib_module.f90

   fsources += EGini.f

   #fsources += EGSlib.f
   fsources += EGMlib.f
   #fsources += EGFlib.f
endif

f90sources += vode_module.f90
fsources += vode.f LinAlg.f

f90sources += tranlib_module.f90
fsources += tranlib_d.f math_d.f
