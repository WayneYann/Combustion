PBOXLIB_HOME= ../../..

FBOXLIB_HOME= ../../../../fParallel

TOP = ../..

PRECISION  = DOUBLE
PROFILE    = FALSE

DEBUG      = FALSE
DEBUG      = TRUE

DIM        = 3

COMP	   = gcc
FCOMP	   = gfortran

USE_MPI    = FALSE

USE_CHEM = TRUE

REACTION_MECHANISM=CHEMH

# These set the EOS and reaction network directories
EOS_dir     := $(FBOXLIB_HOME)/extern/EOS/gamma_law
Network_dir := $(FBOXLIB_HOME)/extern/networks/gammalaw

Bpack   := ./Make.package
Blocs   := .

include ../Make.CNSReact
