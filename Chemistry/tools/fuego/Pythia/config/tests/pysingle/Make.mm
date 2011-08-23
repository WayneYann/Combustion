# -*- Makefile -*-
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
#                               Michael A.G. Aivazis
#                        (C) 1998-2001  All Rights Reserved
#
# <LicenseText>
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#

include local.def

PROJECT = pytest
PACKAGE = 

LIBS = $(EXTERNAL_LIBS)
BLIBS = -L. $(EXTERNAL_LIBS)

PROJ_CLEAN = libsingle.$(EXT_SO) pya.$(EXT_SO) pyb.$(EXT_SO) *~ core


all: pya.$(EXT_SO) pyb.$(EXT_SO)


libsingle.$(EXT_SO): $(PROJ_CXX_LIB)
	$(CXX) Singleton.cc -o libsingle.$(EXT_SO) $(CXXFLAGS) $(LCXXFLAGS) $(LIBS)


pya.$(EXT_SO): libsingle.$(EXT_SO)
	$(CXX) pya.cc -o pya.$(EXT_SO) $(CXXFLAGS) $(LCXXFLAGS) $(BLIBS)


pyb.$(EXT_SO): libsingle.$(EXT_SO)
	$(CXX) pyb.cc -o pyb.$(EXT_SO) $(CXXFLAGS) $(LCXXFLAGS) $(BLIBS)


test: pya.$(EXT_SO) pyb.$(EXT_SO)
	@LD_LIBRARY_PATH=$$LD_LIBRARY_PATH:`pwd` python stest.py



# id
# $Id: Make.mm,v 1.3 2001/09/08 18:37:48 aivazis Exp $

# End of file
