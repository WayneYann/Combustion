
# $Id: GNUmakefile,v 1.19 2010-12-15 00:21:45 marc Exp $
#
PBOXLIB_HOME  = ../..
TOP           = $(PBOXLIB_HOME)
PMFdir        = ../PMFfiles
#
# User-settable variables.
#
PRECISION         = DOUBLE
DEBUG	          = FALSE
PROFILE           = TRUE
PROFILE           = FALSE
DIM    	          = 2
DIM    	          = 3
FCOMP             = Intel
COMP              = Intel
BUILD_IN_PLACE    = TRUE
EBASE             = amr
HGDEBUG           = TRUE
HGDEBUG           = FALSE
USE_FLCTS         = FALSE
USE_VELOCITY      = FALSE
USE_MPI           = TRUE
USE_MPI           = FALSE


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
# Turn this on to get the ability to set the velocities in initData() from a plotfile.
#

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
#REACTION_MECHANISM=DRM19
#REACTION_MECHANISM=GLARBORG
#REACTION_MECHANISM=LIDRYERFN2
REACTION_MECHANISM=GRIMECH30-NON
#
# Packages below are of the form dir/lib.  dir=location of the lib sources,
# lib=library root name.  This info is separated into to corresponding lists.  
#
ProjPackage = $(TOP)/hgproj/proj
ifeq ($(USE_HGPROJ_SERIAL),TRUE)
  ProjPackage = $(TOP)/hgproj-serial/proj
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

ifeq ($(USE_VELOCITY), TRUE)
  MySrcDirs += $(TOP)/pAmrvis
else
  MySrcDirs += $(TOP)/pAmrvis
endif

#CPPFLAGS += -DMPICH_SKIP_MPICXX

vpath %.cpp $(MySrcDirs) ${PMFdir}
vpath %.F   $(MySrcDirs) ${PMFdir}
vpath %.f90 $(MySrcDirs) ${PMFdir}
vpath %.H   $(MySrcDirs) ${PMFdir}
vpath %.h   $(MySrcDirs) ${PMFdir}
vpath %.f   $(MySrcDirs) ${PMFdir}
vpath %.a   $(LIBRARY_LOCATIONS)

all: $(executable)

o/3d.Linux.Intel.Intel.PROF.MPI.EXE/amr_real3d.2.o:
	ifort -module o/3d.Linux.Intel.Intel.PROF.MPI.EXE -O3 -hlo0 -ip -c ../../hgproj/amr_real3d.2.f -o o/3d.Linux.Intel.Intel.PROF.MPI.EXE/amr_real3d.2.o

o/3d.Linux.Intel.Intel.MPI.EXE/amr_real3d.2.o:
	ifort -module o/3d.Linux.Intel.Intel.MPI.EXE -O3 -hlo0 -ip -c ../../hgproj/amr_real3d.2.f -o o/3d.Linux.Intel.Intel.MPI.EXE/amr_real3d.2.o

o/3d.Linux.Intel.Intel.EXE/amr_real3d.2.o:
	ifort -module o/3d.Linux.Intel.Intel.EXE -O3 -hlo0 -ip -c ../../hgproj/amr_real3d.2.f -o o/3d.Linux.Intel.Intel.EXE/amr_real3d.2.o


include $(TOP)/mk/Make.rules
