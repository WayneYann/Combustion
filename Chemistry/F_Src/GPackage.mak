f90sources += chemistry_module.f90

ifdef USE_EGZ
   f90sources += egz_module.f90
else
   f90sources += eglib_module.f90

   fsources += EGaux.f
   fsources += EGini.f

   #fsources += EGSlib.f
   fsources += EGMlib.f
   #fsources += EGFlib.f
endif
