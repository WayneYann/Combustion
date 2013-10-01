PRECISION       = DOUBLE
DEBUG           = TRUE
DEBUG           = FALSE
DIM             = 2
COMP            = g++
FCOMP           = gfortran
USE_MPI         = FALSE
USE_OMP         = FALSE

CHEMISTRY_MODEL = LUDME

EBASE = DMEjet

USE_EGZ := TRUE
USE_WORK_SPACE_MODULE := TRUE

# Set paths to BoxLib and Chemistry
TOP = ../../../../..
BOXLIB_DIR = ${TOP}/BoxLib
CHEMISTRY_DIR = Combustion/Chemistry

include $(BOXLIB_HOME)/Tools/C_mk/Make.defs

Bdirs := $(CHEMISTRY_DIR)/tools/bdf/DMEjetTest

Bdirs += $(CHEMISTRY_DIR)/src_common
Bdirs += $(CHEMISTRY_DIR)/src_f90

Bpack	+= $(foreach dir, $(Bdirs), $(TOP)/$(dir)/Make.package)
Blocs	+= $(foreach dir, $(Bdirs), $(TOP)/$(dir))

Pdirs 	:= C_BaseLib

Bpack	+= $(foreach dir, $(Pdirs), $(BOXLIB_HOME)/Src/$(dir)/Make.package)
Blocs	+= $(foreach dir, $(Pdirs), $(BOXLIB_HOME)/Src/$(dir))

MECHDIR = $(TOP)/$(CHEMISTRY_DIR)/data
ifeq (${CHEMISTRY_MODEL}, CHEMH)
  cEXE_sources += chem-H.c
  vpath %.f $(VPATH_LOCATIONS) ${MECHDIR}/chem-H/PMFs
  vpath %.c $(VPATH_LOCATIONS) ${MECHDIR}/chem-H
else ifeq (${CHEMISTRY_MODEL}, LIDRYER)
  cEXE_sources += LiDryer.c
  vpath %.f $(VPATH_LOCATIONS) ${MECHDIR}/LiDryer/PMFs
  vpath %.c $(VPATH_LOCATIONS) ${MECHDIR}/LiDryer
else ifeq ($(CHEMISTRY_MODEL), DRM19)
  cEXE_sources += drm19.c
  vpath %.f $(VPATH_LOCATIONS) ${MECHDIR}/gri/PMFs
  vpath %.c $(VPATH_LOCATIONS) ${MECHDIR}/gri
else ifeq ($(CHEMISTRY_MODEL), GRI30)
  cEXE_sources += grimech30.c
  vpath %.f $(VPATH_LOCATIONS) ${MECHDIR}/gri/PMFs
  vpath %.c $(VPATH_LOCATIONS) ${MECHDIR}/gri
else ifeq (${CHEMISTRY_MODEL}, LUDME)
  cEXE_sources += LuDME.c
  vpath %.f $(VPATH_LOCATIONS) ${MECHDIR}/Lu/PMFs
  vpath %.c $(VPATH_LOCATIONS) ${MECHDIR}/Lu
else ifneq (${CHEMISTRY_MODEL}, NULL)
  $(error "Unknown CHEMISTRY_MODEL: $(CHEMISTRY_MODEL)")
endif

include $(Bpack)

INCLUDE_LOCATIONS += $(Blocs)
VPATH_LOCATIONS   += $(Blocs)


all: $(executable) 
	@echo SUCCESS

vpath %.c   . $(VPATH_LOCATIONS)
vpath %.cpp . $(VPATH_LOCATIONS)
vpath %.h   . $(VPATH_LOCATIONS)
vpath %.H   . $(VPATH_LOCATIONS)
vpath %.F   . $(VPATH_LOCATIONS)
vpath %.f90 . $(VPATH_LOCATIONS)
vpath %.f   . $(VPATH_LOCATIONS)
vpath %.fi  . $(VPATH_LOCATIONS)


include $(BOXLIB_HOME)/Tools/C_mk/Make.rules


#-----------------------------------------------------------------------------
# for debugging.  To see the value of a Makefile variable,
# e.g. Fmlocs, simply do "make print-Fmlocs".  This will
# print out the value.
print-%: ; @echo $* is $($*)
