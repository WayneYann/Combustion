PBOXLIB_HOME= ../../..
FBOXLIB_HOME= ../../../../fParallel

PRECISION  = DOUBLE
PROFILE    = FALSE

DEBUG      = TRUE
DEBUG      = FALSE

#DIM        = 1
#DIM        = 2
DIM        = 3

COMP	   = gcc
FCOMP	   = gfortran

COMP	   = Intel
FCOMP	   = Intel

USE_MPI    = FALSE

USE_REACT      = FALSE
USE_REACT      = TRUE

USE_DIFFUSION  = TRUE
USE_DIFFUSION  = FALSE

# These set the EOS and reaction network directories
EOS_dir     := $(FBOXLIB_HOME)/extern/EOS/gamma_law
Network_dir := $(FBOXLIB_HOME)/extern/networks/gammalaw

Bpack   := ./Make.package
Blocs   := .

include ../Make.CNSReact
