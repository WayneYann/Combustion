PBOXLIB_HOME= ../../..
FBOXLIB_HOME= ../../../../fParallel

PRECISION  = DOUBLE
PROFILE    = FALSE

DEBUG      = FALSE
DEBUG      = TRUE

DIM        = 2

COMP	   = Intel
FCOMP	   = Intel

USE_MPI    = FALSE

USE_DIFFUSION = TRUE
USE_REACT     = TRUE

# These set the EOS and reaction network directories
EOS_dir     := $(FBOXLIB_HOME)/extern/EOS/gamma_law
Network_dir := $(FBOXLIB_HOME)/extern/networks/gammalaw

Bpack   := ./Make.package
Blocs   := .

include ../Make.CNSReact
