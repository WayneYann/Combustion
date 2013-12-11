BOXLIB_HOME ?= /path/to/BoxLib

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

CHEMISTRY_MODEL = LIDRYER

# Is this a convergence test?  For convergence tests, T is computed to higher accuracy,
# and iterative routines in EGLIB are not used. 
CONVERGENCE = FALSE
CONVERGENCE = TRUE

Bpack   := ./Make.package
Blocs   := .

include ../Make.RNS
