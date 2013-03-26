
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
COMP           = Intel
FCOMP          = Intel
USE_MPI        = TRUE
USE_MPI        = FALSE
BUILD_IN_PLACE = TRUE
PMFDIR         = ../PMFfiles
USE_FLCTS      = FALSE
#
# Uncomment the following if you want RhoRT stored in Tracer slot.
#
#DEFINES += -DBL_RHORT_IN_TRACER
#
# Choose your reaction mechanism & thermal properties.
#
REACTION_MECHANISM=CH4-2STEP
REACTION_MECHANISM=INERT30
REACTION_MECHANISM=GRI30NON
REACTION_MECHANISM=GRI30
REACTION_MECHANISM=GRI12
REACTION_MECHANISM=DRM19
REACTION_MECHANISM=PROPANE
REACTION_MECHANISM=CHEMH
#
# Packages below are of the form dir/lib. dir=location of the lib sources,
#  lib=library root name.  This info is separated into to corresponding
#  lists.
#
Packages :=	$(TOP)/LMC/ht \
		$(TOP)/iamrlib/iamr  \
		$(TOP)/tensorMG/mcmg \
		$(TOP)/hgproj/proj \
		$(TOP)/mglib/mg \
		$(TOP)/amrlib/amr \
		$(TOP)/bndrylib/bndry \
		$(TOP)/BoxLib/box 			
#
# Set standard defs 
#
include $(TOP)/mk/Make.defs
ifeq ($(DIM),3)
  ifeq ($(findstring 7, $(PRVERSION)), 7)
      DEFINES += -DBL_PRVERSION=7
  endif
endif
ifeq (${USE_FLCTS}, TRUE)
  DEFINES += -DBL_DO_FLCT
  INCLUDE_LOCATIONS += $(TOP)/iamrlib/InflowForce
  MySrcDirs += $(TOP)/iamrlib/InflowForce
endif
#
# Setup build locations/libraries
#
MyLibRoots := $(notdir $(Packages))
MyPackageDirs := $(dir $(Packages))
MySrcDirs = . 
MySrcDirs += $(MyPackageDirs)
INCLUDE_LOCATIONS += $(MySrcDirs)
include $(addsuffix /Make.package, $(MySrcDirs))

vpath %.cpp $(MySrcDirs) ..
vpath %.F   $(MySrcDirs) .. 
vpath %.H   $(MySrcDirs) ..
vpath %.h   $(MySrcDirs) ..
vpath %.f   $(MySrcDirs) ${PMFDIR} ..
vpath %.a   $(LIBRARY_LOCATIONS)


all: $(executable)

$(executable): $(LIBRARIES)

include $(TOP)/mk/Make.rules
