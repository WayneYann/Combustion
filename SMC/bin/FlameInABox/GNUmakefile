BOXLIB_HOME ?= ../../../../BoxLib

NDEBUG := t
MPI    := t
OMP    := t
PROF   := t

COMP := gfortran

MKVERBOSE := t

#override F90FLAGS = -J$(mdir) -I $(mdir) -O3 -msse3 -fno-range-check -fopenmp


# define the location of the SMC top directory
SMC_TOP_DIR := ../..

# define the location of Chemistry directory
CHEMISTRY_DIR := ../../../Chemistry
CHEMISTRY_MODEL := LIDRYER


# include the SMC build stuff
include $(SMC_TOP_DIR)/GSMC.mak

