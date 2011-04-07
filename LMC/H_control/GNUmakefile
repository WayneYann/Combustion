
PBOXLIB_HOME  = ../..
TOP           = $(PBOXLIB_HOME)
PMFdir        = ../PMFfiles
#
# User-settable variables.
#
PRECISION         = DOUBLE
DEBUG	          = FALSE
DEBUG	          = TRUE
PROFILE           = TRUE
PROFILE           = FALSE
DIM    	          = 3
DIM    	          = 2
FCOMP             = gfortran
COMP              = g++
BUILD_IN_PLACE    = TRUE
EBASE             = amr
HGDEBUG           = TRUE
HGDEBUG           = FALSE
USE_FLCTS         = FALSE
USE_VELOCITY      = FALSE
USE_MPI           = FALSE
USE_MPI           = TRUE
USE_PARTICLES     = FALSE
USE_PARTICLES     = TRUE

# This allows GetForce to call a Fortran routine
#DEFINES += -DGENGETFORCE

# This is for the time-dependent turbulence forcing term
#DEFINES += -DDO_LMC_FORCE

#NEEDS_FLUSH_F := TRUE
#
# Turn this on to restart from a plotfile.
#
#USE_NEWMECH = TRUE
#
# Uncomment the following if you want RhoRT stored in Tracer slot.
#
#DEFINES += -DBL_RHORT_IN_TRACER

DEFINES += -DCG_USE_OLD_CONVERGENCE_CRITERIA
#
# Choose your reaction mechanism & thermal properties.
#
#REACTION_MECHANISM=INERT30
#REACTION_MECHANISM=GRI12
#REACTION_MECHANISM=GRI30NON
#REACTION_MECHANISM=GRI30
#REACTION_MECHANISM=CH4-2STEP
#REACTION_MECHANISM=GRI30HNO
#REACTION_MECHANISM=GLARBORG
#REACTION_MECHANISM=LIDRYERFN2
#REACTION_MECHANISM=GRIMECH30-NON
REACTION_MECHANISM=DRM19

ProjPackage = $(TOP)/hgproj

ifeq ($(USE_HGPROJ_SERIAL),TRUE)
  ProjPackage = $(TOP)/hgproj-serial
  INCLUDE_LOCATIONS += $(dir $(ProjPackage))include/$(DIM)d.$(PRVERSION)
endif

ifeq ($(USE_VELOCITY), TRUE)
  DEFINES += -DBL_NOLINEVALUES -DBL_USE_VELOCITY -DBL_PARALLEL_IO
  INCLUDE_LOCATIONS += $(TOP)/pAmrvis
endif

ifeq ($(USE_NEWMECH), TRUE)
  DEFINES += -DBL_NOLINEVALUES -DBL_USE_NEWMECH -DBL_PARALLEL_IO
  INCLUDE_LOCATIONS += $(TOP)/pAmrvis
endif

Packages := $(TOP)/LMC \
	    $(TOP)/iamrlib \
	    $(TOP)/tensorMG \
	    $(ProjPackage) \
	    $(TOP)/mglib \
	    $(TOP)/amrlib \
	    $(TOP)/bndrylib \
	    $(TOP)/BoxLib

ifeq ($(USE_PARTICLES), TRUE)
  DEFINES += -DPARTICLES
  Packages += $(TOP)/particles
endif
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
MySrcDirs = .
INCLUDE_LOCATIONS += ..
MySrcDirs += $(Packages)
INCLUDE_LOCATIONS += $(MySrcDirs)

include $(addsuffix /Make.package, $(MySrcDirs))

ifeq ($(HGDEBUG),TRUE)
  DEFINES+=-DHG_DEBUG
endif
ifeq (${USE_FLCTS}, TRUE)
  DEFINES += -DBL_DO_FLCT
  INCLUDE_LOCATIONS += $(TOP)/iamrlib/InflowForce
  MySrcDirs += $(TOP)/iamrlib/InflowForce
endif

ifeq ($(USE_VELOCITY), TRUE)
  MySrcDirs += $(TOP)/pAmrvis
else
  MySrcDirs += $(TOP)/pAmrvis
endif

vpath %.cpp $(MySrcDirs) ${PMFdir}
vpath %.F   $(MySrcDirs) ${PMFdir}
vpath %.f90 $(MySrcDirs) ${PMFdir}
vpath %.H   $(MySrcDirs) ${PMFdir}
vpath %.h   $(MySrcDirs) ${PMFdir}
vpath %.f   $(MySrcDirs) ${PMFdir}
vpath %.a   $(LIBRARY_LOCATIONS)

all: $(executable)

include $(TOP)/mk/Make.rules
