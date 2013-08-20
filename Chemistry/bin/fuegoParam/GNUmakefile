PRECISION       = DOUBLE
DEBUG           = FALSE
DEBUG           = TRUE
DIM             = 2
COMP            = g++
FCOMP           = gfortran
USE_MPI         = FALSE
USE_OMP         = FALSE
USE_SDC         = FALSE
USE_SDC         = TRUE

# Choose model (from list below), and pmf file
#CHEMISTRY_MODEL = DRM19
#CHEMISTRY_MODEL = CHEMH
CHEMISTRY_MODEL = LUDME
#CHEMISTRY_MODEL = GRI30

ifeq (${USE_SDC}, TRUE)
  DEFINES += -DLMC_SDC
  DEFINES += -DDO_TEMPERATURE
  #DEFINES += -DDO_AJAC
endif
CFLAGS += -std=c99

# Set paths to BoxLib and Chemistry
TOP = ../../../..
BOXLIB_DIR = ${TOP}/BoxLib
CHEMISTRY_DIR = ${TOP}/Combustion/Chemistry

EBASE = driver
Blocs := . hacked_src data/Lu 
INCLUDE_LOCATIONS += hacked_src data/Lu

CEXE_headers += Observation.H
CEXE_sources += ${EBASE}.cpp Observation.cpp
f90EXE_sources += driver_f.f90

include $(CHEMISTRY_DIR)/tools/make/Make.CHEM
