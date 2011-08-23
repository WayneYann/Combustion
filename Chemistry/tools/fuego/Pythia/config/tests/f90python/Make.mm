#--------------------------------------------------------------------------
#
#                       Michael A.G. Aivazis
#                California Institute of Technology
#                   (C) 1999 All Rights Reserved
#
#--------------------------------------------------------------------------

# $Log: Make.mm,v $
# Revision 1.4  2000/09/13 19:04:17  aivazis
# Cleaned up a bit
#
# Revision 1.3  2000/09/12 23:59:32  cummings
# Added BLD_ prefix to BINDIR and LIBDIR.  Added LCXX_SOFLAGS to command
# for creating shared object PROJ_DLL.
#
# Revision 1.2  2000/09/12 04:46:48  aivazis
# Temporarily pollute with rm3d calls.
#
# Revision 1.1  2000/09/11 02:15:13  aivazis
# Original source
#
# Revision 1.1  2000/09/06 17:40:09  aivazis
# Original source
#

#
# Local makefile

include local.def

#
PROJECT = pyf90
PACKAGE =

BLD_BINDIR = .
PROJ_DLL = $(BLD_BINDIR)/$(PROJECT).$(EXT_SO)
PROJ_LIB = $(BLD_LIBDIR)/$(PROJECT).$(EXT_LIB)
PROJ_CLEAN += $(PROJ_DLL)


PROJ_SRCS = \
    hello.$(EXT_F77)

MODINIT = $(PROJECT).cc

LIBRARIES = $(PROJ_LIB) $(DEV_LCXX_LIBRARIES) $(EXTERNAL_LIBS) $(COMPILER_LCXX_FORTRAN)

all: $(PROJ_DLL)

$(PROJ_DLL): $(MODINIT) $(PROJ_LIB)
	$(CXX) $(MODINIT) -o $(PROJ_DLL) $(CXXFLAGS) $(LCXX_SOFLAGS) $(LCXXFLAGS) $(LIBRARIES)

export:: all export-python-modules

#
# End of file
