#--------------------------------------------------------------------------
#                       Michael A.G. Aivazis
#                California Institute of Technology
#                   (C) 1999 All Rights Reserved
#
#--------------------------------------------------------------------------

# $Log: Make.mm,v $
# Revision 1.2  2000/09/12 23:56:31  cummings
# Added LF77FLAGS to compile-and-link command.  Added empty
# definition of PROJ_DISTCLEAN, since we are not generating
# dependency files here.
#
# Revision 1.1  2000/08/14 21:32:05  aivazis
# Original source
#

PROJECT = test
PROJ_CLEAN = $(PROJ_BIN)
PROJ_DISTCLEAN = 

PROJ_BIN = hello

all: hello

hello:	hello.f
	$(F77) $(F77FLAGS) $(LF77FLAGS) hello.f -o $(PROJ_BIN)
	$(PROJ_BIN)


#
# End of file
