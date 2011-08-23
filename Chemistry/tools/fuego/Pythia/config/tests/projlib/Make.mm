#--------------------------------------------------------------------------
#                       Michael A.G. Aivazis
#                California Institute of Technology
#                   (C) 1999 All Rights Reserved
#
#--------------------------------------------------------------------------

# $Log: Make.mm,v $
# Revision 1.10  2001/11/10 02:19:00  cummings
# Added modified build rule for compiling main test code when using
# the Tau target.  We parse and instrument the code first, then compile.
#
# Revision 1.9  2000/09/28 22:51:47  cummings
# Changed PROJ_LIB to PROJ_CXX_LIB in order to force template closure
# of C++ object files in library.
#
# Revision 1.8  2000/09/13 00:03:17  cummings
# Added BLD_ prefix to TMPDIR, LIBDIR and BINDIR.  Removed some extra
# definitions and dependencies that are no longer needed.  Added
# LCXXFLAGS to compile-and-link command.  Added execution of "hello"
# test program to target.
#
# Revision 1.7  2000/08/31 23:13:43  cummings
# Changed PROJ_LIB_DEPENDS back to PROJ_DEPENDENCIES here, since
# PROJ_LIB_DEPENDS is gone now.  Added "hello" target to "all" target.
#
# Revision 1.6  2000/08/31 21:52:34  cummings
# Added target "hello", which builds an executable to test the project
# library by calling Hello World! subroutines in each language.
#
# Revision 1.5  2000/08/31 04:45:38  aivazis
# Added include local.def. Renamed the FORTRAN test file with an uppercase
# extension in order to trigger the use of the c preprocessor
#
# Revision 1.4  2000/08/31 02:08:46  cummings
# Added f77_hello.f to PROJ_SRCS.  Changed PROJ_LIB to PROJ_LIBRARY in
# PROJ_CLEAN.  Added tag_projlib and product_dirs as prerequisites for
# default target "all".
#
# Revision 1.3  2000/08/30 02:12:53  cummings
# Modifications to illustrate new library archiving mechanism.
# Make.mm file defines PROJ_LIBRARY as its primary target and
# defines the output directory PROJ_TMPDIR where dependency
# and object files will be created.
#
# Revision 1.2  2000/08/08 18:49:57  aivazis
# Repaired the default target
#
# Revision 1.1  1999/11/27 00:46:42  aivazis
# Original source
#

include local.def

PROJECT = test
PROJ_CXX_LIB = $(BLD_LIBDIR)/$(PROJECT).$(EXT_LIB)
PROJ_TMPDIR = $(BLD_TMPDIR)/$(PROJECT)

PROJ_SRCS = \
    f77_hello.F \
    c_hello.c \
    cpp_hello.cc

all: hello

ifeq (tau, ${findstring tau, $(BLD_USER_TARGET)})
hello:	$(PROJ_CXX_LIB)
	$(PDTCXXPARSE) main.cc $(CXX_BUILD_DEFINES) $(CXX_BUILD_INCLUDES) \
	$(TAU_INCLUDE) $(TAU_DEFS)
	$(TAUINSTR) main.pdb main.cc -o main.inst.cc
	$(TAU_CXX) $(CXXFLAGS) -c main.inst.cc
	$(TAU_CXX) -o $(BLD_BINDIR)/hello main.inst.o \
	$(PROJ_CXX_LIB) $(LCXXFLAGS) $(LCXX_FORTRAN)
	$(RM) $(RMFLAGS) main.inst.o main.pdb main.inst.cc
	$(BLD_BINDIR)/hello
else
hello:	$(PROJ_CXX_LIB)
	$(CXX) -o $(BLD_BINDIR)/hello main.cc $(CXXFLAGS) \
	$(PROJ_CXX_LIB) $(LCXXFLAGS) $(LCXX_FORTRAN)
	$(BLD_BINDIR)/hello
endif

project: show-project


#
# End of file
