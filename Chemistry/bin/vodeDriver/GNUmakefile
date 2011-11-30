PRECISION       = DOUBLE
DEBUG           = TRUE
DEBUG           = FALSE
DIM             = 2
COMP            = g++
FCOMP           = gfortran
USE_MPI         = FALSE
USE_OMP         = FALSE

# Choose model (from list below), and pmf file
CHEMISTRY_MODEL = DRM19
CHEMISTRY_MODEL = CHEMH

# Set paths to BoxLib and Chemistry
TOP = ../../../..
BOXLIB_DIR = ${TOP}/BoxLib
CHEMISTRY_DIR = ${TOP}/Combustion/Chemistry

EBASE = cdvode
Blocs := .

CEXE_sources += vodeDriver.cpp

include $(CHEMISTRY_DIR)/tools/make/Make.CHEM
