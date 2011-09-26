TOP = ../../../..
BOXLIB_DIR = ${TOP}/BoxLib
CHEMISTRY_DIR = ${TOP}/Combustion/Chemistry

#
# Build configuration
#

PRECISION       = DOUBLE
DEBUG           = TRUE
DEBUG           = FALSE
DIM             = 2
COMP            = g++
FCOMP           = gfortran
USE_MPI         = FALSE
USE_OMP         = FALSE

# Choose model (from list below), and pmf file
CHEMISTRY_MODEL = CHEMH
CHEMISTRY_MODEL = DRM19

EBASE = cdvode
Blocs := .

CEXE_sources += vodeDriver.cpp
FEXE_sources += 
CEXE_headers += 
FEXE_headers += 

include $(CHEMISTRY_DIR)/tools/make/Make.CHEM

# Available chemistry models:
# CanuH2 chem-CH4-2step chem-H dme glar gri Hai H-CW header inert Konnov
# LiDryer Marinov prf_ethanol Roehl sandiego smooke usc


