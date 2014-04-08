f90sources += init_data.f90

ifeq ($(CHEMISTRY_MODEL),LIDRYER)
  fsources += LiDryer_040_01ATM.f
else ifeq ($(CHEMISTRY_MODEL),DRM19)
  fsources += drm19Soln_seed_0.50.f
else ifeq ($(CHEMISTRY_MODEL),GRI30)
  fsources += gri30_070.f
else ifeq ($(CHEMISTRY_MODEL), LUDME)
  fsources += LuDME_0700.f
endif
