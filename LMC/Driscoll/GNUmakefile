#
# $Id: GNUmakefile,v 1.6 2009-11-19 19:08:05 lijewski Exp $
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
COMP           = Intel
FCOMP          = Intel
USE_FLCTS      = FALSE
USE_FLCTS      = TRUE
USE_MPI        = TRUE
USE_MPI        = FALSE
EBASE          = amr
USE_HGPROJ_SERIAL = FALSE
USE_NEWMECH    = TRUE
USE_NEWMECH    = FALSE
#
#DEFINES      += -static
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
#DEFINES += -DBL_PLOT_CONSUMPTION
#
# Choose your reaction mechanism & thermal properties.
#
REACTION_MECHANISM=CHEMH
REACTION_MECHANISM=GRI30
REACTION_MECHANISM=GRI12
REACTION_MECHANISM=INERT30
REACTION_MECHANISM=PROPANE
REACTION_MECHANISM=CH4-2STEP
REACTION_MECHANISM=DRM19


DEFINES += -DCG_USE_OLD_CONVERGENCE_CRITERIA

include $(TOP)/mk/Make.defs

Bdirs := LMC LMC/Cheng BoxLib iamrlib amrlib tensorMG mglib bndrylib hgproj

ifeq ($(DIM),3)
  ifeq ($(findstring 7, $(PRVERSION)), 7)
      DEFINES += -DBL_PRVERSION=7
  endif
endif

ifeq ($(USE_VELOCITY), TRUE)
  DEFINES += -DBL_NOLINEVALUES -DBL_USE_VELOCITY -DBL_PARALLEL_IO
  #Bdirs   += pAmrvis
endif

ifeq ($(USE_FLCTS),TRUE)
  Bdirs   += iamrlib/InflowForce
  DEFINES += -DBL_DO_FLCT
endif


ifeq (${USE_SWIRL}, TRUE)
  DEFINES += -DSWIRL
endif

ifeq (${USE_SWIRL_DATA}, TRUE)
  DEFINES += -DSWIRL_DATA
endif

ifeq ($(USE_NEWMECH), TRUE)
  DEFINES += -DBL_NOLINEVALUES -DBL_USE_NEWMECH -DBL_PARALLEL_IO
  #Bdirs   += pAmrvis
endif

Bpack	:= $(foreach dir, $(Bdirs), $(TOP)/$(dir)/Make.package)
Blocs	:= $(foreach dir, $(Bdirs), $(TOP)/$(dir))

include $(Bpack)

ifdef FBOXLIB_HOME
  DEFINES += -DMG_USE_FBOXLIB

  include $(PBOXLIB_HOME)/MGT_Solver/Make.package
  INCLUDE_LOCATIONS += $(PBOXLIB_HOME)/MGT_Solver
  VPATH_LOCATIONS += $(PBOXLIB_HOME)/MGT_Solver

  include $(FBOXLIB_HOME)/mg/FParallelMG.mak
  INCLUDE_LOCATIONS += $(FBOXLIB_HOME)/mg/
  Fdirs   := boxlib mg extern/SPARSKIT extern/LAPACK
  Flocs   := $(foreach dir, $(Fdirs), $(FBOXLIB_HOME)/$(dir))
  VPATH_LOCATIONS += $(Flocs)
endif

INCLUDE_LOCATIONS += . $(Blocs)
VPATH_LOCATIONS   += . $(Blocs)

vpath %.c   $(VPATH_LOCATIONS)
vpath %.cpp $(VPATH_LOCATIONS)
vpath %.h   $(VPATH_LOCATIONS)
vpath %.H   $(VPATH_LOCATIONS)
vpath %.F   $(VPATH_LOCATIONS)
vpath %.f   $(VPATH_LOCATIONS)
vpath %.f90 $(VPATH_LOCATIONS)

all: $(executable)

o/3d.Linux.Intel.Intel.PROF.MPI.EXE/amr_real3d.2.o:
	ifort -module o/3d.Linux.Intel.Intel.PROF.MPI.EXE -O3 -hlo0 -ip -c ../../hgproj/amr_real3d.2.f -o o/3d.Linux.Intel.Intel.PROF.MPI.EXE/amr_real3d.2.o

o/3d.Linux.Intel.Intel.MPI.EXE/amr_real3d.2.o:
	ifort -module o/3d.Linux.Intel.Intel.MPI.EXE -O3 -hlo0 -ip -c ../../hgproj/amr_real3d.2.f -o o/3d.Linux.Intel.Intel.MPI.EXE/amr_real3d.2.o

o/3d.Linux.Intel.Intel.EXE/amr_real3d.2.o:
	ifort -module o/3d.Linux.Intel.Intel.EXE -O3 -hlo0 -ip -c ../../hgproj/amr_real3d.2.f -o o/3d.Linux.Intel.Intel.EXE/amr_real3d.2.o

include $(TOP)/mk/Make.rules
