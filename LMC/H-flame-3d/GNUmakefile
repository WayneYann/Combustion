#
# $Id: GNUmakefile,v 1.5 2004-07-13 21:16:37 lijewski Exp $
#
PBOXLIB_HOME  = ../..
TOP           = $(PBOXLIB_HOME)
#
# User-settable variables.
#
PRECISION      = DOUBLE
DEBUG	       = FALSE
DEBUG	       = TRUE
PROFILE        = FALSE
DIM	       = 3
PRVERSION      = v7
COMP           = g++
FCOMP          = Intel
USE_MPI        = TRUE
USE_MPI        = FALSE
BUILD_IN_PLACE = TRUE
RUNDIR         = run
EBASE          = amr
USE_HGPROJ_SERIAL := FALSE
#
# Choose your reaction mechanism & thermal properties.
#
REACTION_MECHANISM=CH4-2STEP
REACTION_MECHANISM=INERT30
REACTION_MECHANISM=GRI30
REACTION_MECHANISM=GRI12
REACTION_MECHANISM=DRM19
REACTION_MECHANISM=CHEMH
#
# Packages below are of the form dir/lib. dir=location of the lib sources,
#  lib=library root name.  This info is separated into to corresponding
#  lists.
#
ProjPackage = $(TOP)/hgproj/proj
ifeq ($(USE_HGPROJ_SERIAL),TRUE)
ProjPackage = $(TOP)/hgproj-serial/proj
INCLUDE_LOCATIONS += $(dir $(ProjPackage))include/$(DIM)d.$(PRVERSION)
endif

Packages :=	$(TOP)/LMC/ht \
		$(TOP)/iamrlib/iamr  \
		$(TOP)/tensorMG/mcmg \
		$(ProjPackage) \
		$(TOP)/mglib/mg \
		$(TOP)/amrlib/amr \
		$(TOP)/bndrylib/bndry \
		$(TOP)/BoxLib/box 			
#
# Set standard defs 
#
include $(TOP)/mk/Make.defs
ifeq ($(USE_HGPROJ_SERIAL),TRUE)
DEFINES += -DBL_USE_HGPROJ_SERIAL
endif

ifeq ($(DIM),3)
  ifeq ($(findstring 7, $(PRVERSION)), 7)
      DEFINES += -DBL_PRVERSION=7
  endif
endif
#
# Setup build locations/libraries
#
MyLibRoots := $(notdir $(Packages))
MyPackageDirs := $(dir $(Packages))
MySrcDirs = . 
INCLUDE_LOCATIONS += ..
INCLUDE_LOCATIONS += $(TOP)/iamrlib/InflowForce
ifeq ($(BUILD_IN_PLACE),TRUE)
  MySrcDirs += $(MyPackageDirs)
  INCLUDE_LOCATIONS += $(MySrcDirs)
else
  LIBRARIES += $(addsuffix $(DIM)d, $(addprefix -l, $(MyLibRoots)))
  INCLUDE_LOCATIONS += $(TOP)/include
  ifeq ($(USE_HGPROJ_SERIAL),TRUE)
    LIBRARY_LOCATIONS += $(dir $(ProjPackage))lib/$(machineSuffix)
    LIBRARIES := $(patsubst -lproj$(DIM)d, -lproj$(DIM)d.$(PRVERSION), $(LIBRARIES))
  endif
  LIBRARY_LOCATIONS += $(TOP)/lib/$(machineSuffix)
endif
include $(addsuffix /Make.package, $(MySrcDirs))

vpath %.cpp $(MySrcDirs) ..
vpath %.F   $(MySrcDirs) .. 
vpath %.H   $(MySrcDirs) ..
vpath %.h   $(MySrcDirs) ..
vpath %.f   $(MySrcDirs) ..
vpath %.a   $(LIBRARY_LOCATIONS)
vpath %.H   $(TOP)/iamrlib/InflowForce
vpath %.F   $(TOP)/iamrlib/InflowForce


MyMakeLine = $(MAKE) PRECISION=$(PRECISION) PROFILE=$(PROFILE) COMP=$(COMP) \
	DEBUG=$(DEBUG) DIM=$(DIM) USE_MPI=$(USE_MPI) PRVERSION=$(PRVERSION) \
	USE_HGPROJ_SERIAL=$(USE_HGPROJ_SERIAL) \
	DEFINES="$(DEFINES)" EBASE= LBASE=$@

MyLibDeps := $(shell perl -e '$$,=" ";print reverse @ARGV ;' $(MyLibRoots))

all: $(executable)

$(executable): $(LIBRARIES)

libs: $(MyLibDeps)

$(MyLibDeps):
	cd $(dir $(filter %/$@, $(Packages))); $(MyMakeLine) install

cleanlibs:
	$(foreach lib,$(MyPackageDirs),cd $(lib); $(MyMakeLine) clean)

nodata:
	\rm -rf chk* plt* grdlog inspect.fab *.mfab *.fab *.asc *.out
#
# Create MS-IDE dsp file to build this project
#
DspOlevel = 1
DspFileRoot = $(EBASE)$(DIM)d
ifeq ($(DEBUG), TRUE)
  DspFileRoot := $(join $(DspFileRoot),.DEBUG)
  DspOlevel := 0
endif
dsp: $(DspFileRoot).dsp
.dspDepends: $(FEXE_sources) $(FEXE_headers) $(fEXE_sources) $(CEXE_headers) $(CEXE_sources) $(INCLUDE_LOCATIONS)
	@echo $^ > .dspDepends

$(DspFileRoot).dsp: .dspDepends
	@(echo Building DSP file = $(DspFileRoot).dsp;\
	$(TOP)/scripts/dsp.mak -p $(DspFileRoot) -t $(TOP) -d $(DIM) -o $(DspFileRoot).dsp -O $(DspOlevel) -f $^)

.PHONY: dsp .dspDepends libs cleanlibs nodata touch_hg

include $(TOP)/mk/Make.rules
