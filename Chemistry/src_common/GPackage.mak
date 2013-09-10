
fsources += vode.f LinAlg.f math_d.f tranlib_d.f

ifdef USE_EGZ
  f90sources += egz_module.f90
else
  f90sources += eglib_module.f90
  fsources += EGini.f
  ifdef USE_EGM
     fsources += EGMlib.f
  else
     ifdef USE_EGF
        fsources += EGFlib.f
     else
        fsources += EGSlib.f
     endif
  endif  
endif

