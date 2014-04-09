BOXLIB_HOME ?= /path/to/BoxLib
SDCLIB_HOME ?= /path/to/SDCLib

TOP = ../..

PRECISION  = DOUBLE
PROFILE    = FALSE

DEBUG      = FALSE
#DEBUG      = TRUE

DIM        = 1

COMP	   = gcc
FCOMP	   = gfortran

USE_MPI    = TRUE
USE_OMP    = FALSE

USE_SDCLIB = TRUE

CHEMISTRY_MODEL = LOCAL

# Is this a convergence test?  For convergence tests, T is computed to higher accuracy,
# and iterative routines in EGLIB are not used. 
CONVERGENCE = FALSE

Bpack   := ./Make.package
Blocs   := .

include ../Make.RNS
