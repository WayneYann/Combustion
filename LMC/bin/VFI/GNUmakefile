TOP =../../../..
BOXLIB_DIR = ${TOP}/BoxLib
IAMR_DIR = ${TOP}/IAMR
COMBUSTION_DIR = ${TOP}/Combustion

#
# Build configuration
#

PRECISION       = DOUBLE
DEBUG           = FALSE
DEBUG           = TRUE
DIM             = 2
COMP            = g++
FCOMP           = gfortran
USE_MPI         = FALSE
USE_OMP         = FALSE
USE_MPI         = TRUE
VERBOSE         = FALSE
VERBOSE         = TRUE

# Enable SDC code - add .SDC to executable
# and compile in the #ifdef LMC_SDC
USE_LMC_SDC = FALSE
USE_LMC_SDC = TRUE

# Enable flctin namelist
#USE_FLCTS       = TRUE

USE_PARTICLES   = FALSE

# Enable ht.velocity_plotfile option
#USE_VELOCITY    = FALSE

# Enable time-dependent forcing term
#USE_TDEPFORCING = FALSE

# Choose model (from list below), and pmf file
CHEMISTRY_MODEL = CHEMH
fEXE_sources += H-flame_0370.f

Blocs   := .

CEXE_sources +=
FEXE_sources += PROB_$(DIM)D.F
CEXE_headers += 
FEXE_headers += probdata.H PROB_F.H

include $(COMBUSTION_DIR)/LMC/tools/make/Make.LMC

# Available chemistry models:
# CanuH2 chem-CH4-2step chem-H dme glar gri Hai H-CW header inert Konnov
# LiDryer Marinov prf_ethanol Roehl sandiego smooke usc

#DEFINES += -DDO_JBB_HACK
DEFINES += -DTRIGGER_NEW_J
#DEFINES += -DALWAYS_NEW_J
#DEFINES += -DSOLN_IS_1D

#DEFINES += -DDO_JBB_HACK_TEMP
#DEFINES += -DHACK_TEMP_MIN=273.d0
#DEFINES += -DHACK_TEMP_MAX=1900.d0
