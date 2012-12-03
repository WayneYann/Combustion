TOP = ../../../..
BOXLIB_DIR = ${TOP}/BoxLib
IAMR_DIR = ${TOP}/IAMR
COMBUSTION_DIR = ${TOP}/Combustion

#
# Build configuration
#

PRECISION       = DOUBLE
DEBUG           = TRUE
DEBUG           = FALSE
DIM             = 3
DIM             = 2
COMP            = g++
FCOMP           = gfortran
COMP            = Intel
FCOMP           = Intel
USE_MPI         = FALSE
USE_MPI         = TRUE
USE_OMP         = FALSE
USE_F90_SOLVERS = TRUE
VERBOSE         = FALSE

# Enable flctin namelist
#USE_FLCTS       = TRUE

USE_LMC_SDC = FALSE

#USE_PARTICLES   = TRUE

# Enable ht.velocity_plotfile option
#USE_VELOCITY    = TRUE

# Enable the restart-from-plotfile option (ht.pltfile)
#  ... default this off because when on, code expects it to be used
# USE_NEWMECH    = TRUE

# Enable time-dependent forcing term
#USE_TDEPFORCING = TRUE

# Choose model (from list below), and pmf file
#CHEMISTRY_MODEL = DRM19
#fEXE_sources += drm19Soln_seed_0.50.f

#CHEMISTRY_MODEL = LIDRYER
#fEXE_sources += LiDryer_040_01ATM.f

CHEMISTRY_MODEL = LUDME
fEXE_sources += LuDME_0700.f

EBASE = amr
Blocs   := .

CEXE_sources +=
FEXE_sources += PROB_$(DIM)D.F PROB_F.F
CEXE_headers += 
FEXE_headers += probdata.H PROB_F.H

include $(COMBUSTION_DIR)/LMC/tools/make/Make.LMC
# Available chemistry models:
# CanuH2 chem-CH4-2step chem-H dme glar gri Hai H-CW header inert Konnov
# LiDryer Marinov prf_ethanol Roehl sandiego smooke usc
