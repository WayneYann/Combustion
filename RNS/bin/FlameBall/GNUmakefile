BOXLIB_HOME ?= /path/to/BoxLib
SDCLIB_HOME ?= /path/to/SDCLib

TOP = ../..

PRECISION  = DOUBLE
PROFILE    = FALSE

DEBUG      = FALSE
DEBUG      = TRUE

DIM        = 2

COMP	   = gcc
FCOMP	   = gfortran

USE_MPI    = FALSE
USE_OMP    = FALSE

USE_SDCLIB = TRUE

CHEMISTRY_MODEL = LIDRYER

# Is this a convergence test?  For convergence tests, T is computed to higher accuracy,
# and iterative routines in EGLIB are not used. 
CONVERGENCE = FALSE
CONVERGENCE = TRUE

Bpack   := ./Make.package
Blocs   := .

include ../Make.RNS
