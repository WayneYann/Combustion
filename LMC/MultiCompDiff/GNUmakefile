#
# $Id: GNUmakefile,v 1.22 2007-09-10 19:47:56 sepp Exp $
#
PBOXLIB_HOME  = ../..
TOP           = $(PBOXLIB_HOME)
#
# User-settable variables.
#
PRECISION         = DOUBLE
DEBUG	          = TRUE
DEBUG	          = FALSE
PROFILE           = TRUE
PROFILE           = FALSE
DIM    	          = 2
PRVERSION         = v5
COMP              = mpCC
FCOMP             = xlf
COMP              = Intel
FCOMP             = Intel
RUNDIR            = run
BUILD_IN_PLACE    = TRUE
EBASE             = amr
USE_HGPROJ_SERIAL = FALSE
HGDEBUG           = FALSE
USE_FLCTS         = TRUE
USE_FLCTS         = FALSE
USE_MPI           = FALSE
USE_MPI           = TRUE
#
# Uncomment the following if you want RhoRT stored in Tracer slot.
#
#DEFINES += -DBL_RHORT_IN_TRACER
#
# If you want to plot reaction rates uncomment the following
#
#DEFINES += -DBL_PLOT_REACS
#
# If you want to plot fuel consumption rates uncomment the following
#
DEFINES += -DBL_PLOT_CONSUMPTION
#
# Choose your reaction mechanism & thermal properties.
#
REACTION_MECHANISM=CH4-2STEP
REACTION_MECHANISM=INERT30
REACTION_MECHANISM=GRI12
REACTION_MECHANISM=GRI30NON
REACTION_MECHANISM=PROPANE
REACTION_MECHANISM=GRI30HNO
REACTION_MECHANISM=GRI211HNO
REACTION_MECHANISM=GRI30
REACTION_MECHANISM=DRM19
REACTION_MECHANISM=CHEMH
#
# Packages below are of the form dir/lib.  dir=location of the lib sources,
# lib=library root name.  This info is separated into to corresponding lists.  
#
ProjPackage = $(TOP)/hgproj/proj
ifeq ($(USE_HGPROJ_SERIAL),TRUE)
  ProjPackage = $(TOP)/hgproj-serial/proj
  INCLUDE_LOCATIONS += $(dir $(ProjPackage))include/$(DIM)d.$(PRVERSION)
endif

Packages := $(TOP)/LMC/ht \
	    $(TOP)/iamrlib/iamr \
	    $(TOP)/tensorMG/mcmg \
	    $(ProjPackage) \
	    $(TOP)/mglib/mg \
	    $(TOP)/amrlib/amr \
	    $(TOP)/bndrylib/bndry \
	    $(TOP)/BoxLib/box
#
# Set standard definitions.
#
include $(TOP)/mk/Make.defs
ifeq ($(USE_HGPROJ_SERIAL),TRUE)
  DEFINES += -DBL_USE_HGPROJ_SERIAL
endif

ifeq ($(DIM),2)
  ifeq ($(findstring 5, $(PRVERSION)), 5)
    DEFINES += -DBL_PRVERSION=5
  endif
  ifeq ($(findstring 9, $(PRVERSION)), 9)
    DEFINES += -DBL_PRVERSION=9
  endif
endif
#
# Setup build locations/libraries.
#
MyLibRoots := $(notdir $(Packages))
MyPackageDirs := $(dir $(Packages))
MySrcDirs = .
INCLUDE_LOCATIONS += ..
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

ifeq ($(HGDEBUG),TRUE)
  DEFINES+=-DHG_DEBUG
endif
ifeq (${USE_FLCTS}, TRUE)
  DEFINES += -DBL_DO_FLCT
  INCLUDE_LOCATIONS += $(TOP)/iamrlib/InflowForce
  MySrcDirs += $(TOP)/iamrlib/InflowForce
endif

vpath %.cpp $(MySrcDirs) ..
vpath %.F   $(MySrcDirs) ..
vpath %.f90 $(MySrcDirs) ..
vpath %.H   $(MySrcDirs) ..
vpath %.h   $(MySrcDirs) ..
vpath %.f   $(MySrcDirs) ..
vpath %.a   $(LIBRARY_LOCATIONS)

MyMakeLine =  $(MAKE) PRECISION=$(PRECISION) PROFILE=$(PROFILE) COMP=$(COMP) \
	      DEBUG=$(DEBUG) DIM=$(DIM) USE_MPI=$(USE_MPI) PRVERSION=$(PRVERSION) \
	      USE_HGPROJ_SERIAL=$(USE_HGPROJ_SERIAL) \
	      DEFINES="$(DEFINES)" EBASE= LBASE=$@

MyLibDeps := $(shell perl -e '$$,=" ";print reverse @ARGV ;' $(MyLibRoots))

find_files = $(wildcard $(dir)*cpp) $(wildcard $(dir)*H) $(wildcard $(dir)*.F)
projFiles := $(foreach dir,$(dir $(ProjPackage)),$(find_files))
iamrFiles := $(foreach dir,$(dir $(filter %/iamr,$(Packages))),$(find_files))
StencilSpecificProjFiles = $(shell grep -l USE_HGPROJ_SERIAL $(iamrFiles)) \
                           $(shell grep -l BL_PRVERSION $(iamrFiles)) \
                           $(wildcard $(dir $(ProjPackage))*F) \
                           $(wildcard $(dir $(ProjPackage))*cpp) \
                           $(wildcard $(dir $(ProjPackage))*H)
ifeq ($(USE_HGPROJ_SERIAL),TRUE)
  StencilSpecificProjFiles += $(shell grep -l hg_version.H $(projFiles))
endif
#
# Targets
#
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
# Touch all stencil-specific files in serial projection
#
touch_hg:
	touch $(StencilSpecificProjFiles)

.PHONY: libs cleanlibs nodata touch_hg

include $(TOP)/mk/Make.rules
