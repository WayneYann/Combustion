TOP = ../../../..
BOXLIB_DIR = ${TOP}/BoxLib
IAMR_DIR = ${TOP}/IAMR
COMBUSTION_DIR = ${TOP}/Combustion

#
# Build configuration
#

PRECISION       = DOUBLE

DEBUG           = TRUE

DIM             = 2

COMP            = g++
FCOMP           = gfortran

USE_MPI         = TRUE
USE_MPI         = FALSE

USE_OMP         = TRUE
USE_OMP         = FALSE

#USE_F90_SOLVERS = TRUE
USE_F90_SOLVERS = FALSE

# Enable flctin namelist
USE_FLCTS       = TRUE

USE_PARTICLES   = TRUE

# Enable ht.velocity_plotfile option
USE_VELOCITY    = TRUE

# This allows GetForce to call a Fortran routine
#DEFINES += -DGENGETFORCE

# Choose model from list below
CHEMISTRY_MODEL = DRM19
fEXE_sources += drm19Soln_seed_0.50.f

EBASE = amr
Blocs   := .

CEXE_sources +=
FEXE_sources += PROB_$(DIM)D.F
CEXE_headers += ArrayViewEXT.H
FEXE_headers += probdata.H PROB_F.H


include $(COMBUSTION_DIR)/LMC/tools/make/Make.LMC

# Available chemistry models:
# CanuH2 chem-CH4-2step chem-H dme glar gri Hai H-CW header inert Konnov
# LiDryer Marinov prf_ethanol Roehl sandiego smooke usc
