TOP = ../../../..
BOXLIB_DIR = ${TOP}/BoxLib
IAMR_DIR = ${TOP}/IAMR
COMBUSTION_DIR = ${TOP}/Combustion

#
# User-settable variables.
#
PRECISION         = DOUBLE
DEBUG	          = TRUE
DEBUG	          = FALSE
DIM    	          = 2
COMP              = g++
FCOMP             = gfortran
USE_MPI           = TRUE
USE_OMP           = FALSE
VERBOSE           = FALSE
USE_LMC_SDC       = TRUE

# Enable flctin namelist
USE_FLCTS         = FALSE

USE_PARTICLES     = FALSE

# Enable ht.velocity_plotfile option
USE_VELOCITY      = FALSE

# Enable time-dependent forcing term
USE_TDEPFORCING = FALSE

CHEMISTRY_MODEL= BURKEMOD
fEXE_sources += pmf_BurkeMOD_040_10atm.f

Blocs   := .

CEXE_sources +=
FEXE_sources += PROB_$(DIM)D.F
CEXE_headers += 
FEXE_headers += probdata.H

include $(COMBUSTION_DIR)/LMC/tools/make/Make.LMC
