TOP = ../../../..
BOXLIB_DIR = ${TOP}/BoxLib
IAMR_DIR = ${TOP}/IAMR
COMBUSTION_DIR = ${TOP}/Combustion

#
# User-settable variables.
#
PRECISION         = DOUBLE
DEBUG	          = FALSE
DIM    	          = 2
COMP              = g++
FCOMP             = gfortran
USE_MPI           = TRUE
USE_OMP           = FALSE
USE_F90_SOLVERS   = FALSE
VERBOSE           = FALSE

# Enable flctin namelist
USE_FLCTS         = TRUE

USE_PARTICLES     = FALSE

# Enable ht.velocity_plotfile option
USE_VELOCITY      = FALSE

# Enable time-dependent forcing term
USE_TDEPFORCING = FALSE

CHEMISTRY_MODEL=LIDRYER
fEXE_sources += LiDryer_060_05ATM.f

EBASE = amr
Blocs   := .

CEXE_sources +=
FEXE_sources += PROB_$(DIM)D.F
CEXE_headers += 
FEXE_headers += probdata.H

include $(COMBUSTION_DIR)/LMC/tools/make/Make.LMC
