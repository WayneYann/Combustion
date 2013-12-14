BOXLIB_HOME ?= /path/to/BoxLib
SDCLIB_HOME ?= ../../../../SDCLib

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

CHEMISTRY_MODEL = NULL

Bpack   := ./Make.package
Blocs   := .

#ZMQ = $(HOME)/anaconda
ifdef ZMQ
  cEXE_sources += zmq.c
endif

VERBOSE=FALSE

include ../Make.RNS

ifdef ZMQ
  CXXFLAGS += -DZMQ
  FFLAGS   += -DZMQ
  CFLAGS   += -DZMQ -I$(ZMQ)/include
  LDFLAGS  += -L$(ZMQ)/lib -lzmq
endif


