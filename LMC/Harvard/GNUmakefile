PBOXLIB_HOME  = ../..
TOP           = $(PBOXLIB_HOME)
PRECISION      = DOUBLE
PROFILE        = FALSE
DIM    	       = 2
COMP           = Intel
FCOMP          = Intel
USE_RAD        = FALSE
BUILD_IN_PLACE = TRUE
EBASE          = amr
USE_HGPROJ_SERIAL = FALSE
USE_MPI           = FALSE
USE_MPI           = TRUE
DEBUG	          = FALSE

REACTION_MECHANISM=GRI30
REACTION_MECHANISM=GRI12
REACTION_MECHANISM=INERT30
REACTION_MECHANISM=GRI30NON
REACTION_MECHANISM=CH4-2STEP
REACTION_MECHANISM=PROPANE
REACTION_MECHANISM=CHEMH
REACTION_MECHANISM=GLARBORG
REACTION_MECHANISM=DRM19

include $(TOP)/mk/Make.defs


Bdirs := LMC LMC/Harvard BoxLib iamrlib amrlib tensorMG mglib bndrylib hgproj


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

include $(TOP)/mk/Make.rules
