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
DIM             = 2
COMP            = g++
FCOMP           = gfortran
USE_MPI         = TRUE
USE_OMP         = FALSE
VERBOSE         = FALSE

# Enable flctin namelist
#USE_FLCTS       = TRUE

USE_LMC_SDC = TRUE

#USE_PARTICLES   = TRUE

# Enable ht.velocity_plotfile option
#USE_VELOCITY    = TRUE

# Enable the restart-from-plotfile option (ht.pltfile)
#  ... default this off because when on, code expects it to be used
# USE_NEWMECH    = TRUE

# Enable time-dependent forcing term
#USE_TDEPFORCING = TRUE

# Choose model and pmf file
#CHEMISTRY_MODEL = DRM19
#fEXE_sources += drm19Soln_seed_0.50.f

#CHEMISTRY_MODEL = LIDRYER
#fEXE_sources += LiDryer_040_01ATM.f

#CHEMISTRY_MODEL = LUDME
#fEXE_sources += LuDME_0700.f

#CHEMISTRY_MODEL = GRI30
#fEXE_sources += gri30_070.f

CHEMISTRY_MODEL = ALZETA
f90EXE_sources += pmf_generic.f90

Blocs   := .

CEXE_sources +=
FEXE_sources += PROB_$(DIM)D.F PROB_F.F
CEXE_headers += 
FEXE_headers += probdata.H PROB_F.H

include $(COMBUSTION_DIR)/LMC/tools/make/Make.LMC

#-----------------------------------------------------------------------------
# for debugging.  To see the value of a Makefile variable,
# e.g. Fmlocs, simply do "make print-Fmlocs".  This will
# print out the value.
print-%: ; @echo $* is $($*)
