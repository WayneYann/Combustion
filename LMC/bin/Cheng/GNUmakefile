TOP = ../../../..
BOXLIB_DIR = ${TOP}/BoxLib
IAMR_DIR = ${TOP}/IAMR
COMBUSTION_DIR = ${TOP}/Combustion

#
# Build configuration
#

PRECISION       = DOUBLE
DEBUG           = FALSE
DEBUG           = TRUE
DIM             = 3
COMP            = g++
FCOMP           = gfortran
USE_MPI         = TRUE
USE_MPI         = FALSE
USE_OMP         = FALSE
USE_F90_SOLVERS = FALSE
VERBOSE         = FALSE

# Enable flctin namelist
#USE_FLCTS       = TRUE
USE_FLCTS       = FALSE

#USE_PARTICLES   = TRUE
USE_PARTICLES   = FALSE

# Enable ht.velocity_plotfile option
USE_VELOCITY    = FALSE

# Enable time-dependent forcing term
USE_TDEPFORCING = TRUE
USE_TDEPFORCING = FALSE

# Choose model (from list below), and pmf file
#CHEMISTRY_MODEL = DRM19
#fEXE_sources += drm19Soln_seed_0.50.f

CHEMISTRY_MODEL = GRI30
fEXE_sources += gri30_070.f

EBASE = amr
Blocs   := .

CEXE_sources +=
FEXE_sources += PROB_$(DIM)D.F
CEXE_headers += 
FEXE_headers += probdata.H

include $(COMBUSTION_DIR)/LMC/tools/make/Make.LMC

# Available chemistry models:
# CanuH2 chem-CH4-2step chem-H dme glar gri Hai H-CW header inert Konnov
# LiDryer Marinov prf_ethanol Roehl sandiego smooke usc
