#
# $Id: GNUmakefile,v 1.1 2003-08-08 20:09:41 lijewski Exp $
#
PBOXLIB_HOME = ..

TOP = $(PBOXLIB_HOME)
#
# Where libraries and include files will be installed.
#
INSTALL_ROOT = $(TOP)
#
# Variables for the user to set ...
#
PRECISION     = DOUBLE
DEBUG	      = TRUE
PROFILE       = FALSE
DIM	      = 3
COMP          = KCC
USE_TENSOR    = FALSE
USE_RAD       = FALSE
USE_WINDOWS   = FALSE
USE_MPI       = TRUE
#
# Use hgproj-serial -- only for testing.
#
# Touch Projection.cpp if you want to change the state of USE_HGPROJ_SERIAL.
#
USE_HGPROJ_SERIAL = TRUE
USE_HGPROJ_SERIAL = FALSE
#
# What stencil do you want to use?
#
PRVERSION = v9

ifeq ($(USE_HGPROJ_SERIAL),TRUE)
DEFINES += -DBL_USE_HGPROJ_SERIAL
endif

ifeq ($(USE_RAD),TRUE)
DEFINES += -DBL_USE_RADIATION
endif
#
# The base name of the library we're building.
#
LBASE = ht

include $(TOP)/mk/Make.defs ./Make.package

ifeq ($(USE_HGPROJ_SERIAL),FALSE)
INCLUDE_LOCATIONS += . $(TOP)/include
else
INCLUDE_LOCATIONS += . $(TOP)/hgproj-serial
INCLUDE_LOCATIONS += $(TOP)/hgproj-serial/include/$(DIM)d.$(PRVERSION)
INCLUDE_LOCATIONS += $(TOP)/include
endif


ifeq ($(KCC_VERSION),3.3)
CXXFLAGS += --one_instantiation_per_object
endif
#
# Turn on bounds checking for debug
#
FDEBF += -C

all: $(optionsLib)

#
# Libraries to close against.
#
ifeq ($(COMP),KCC)
LibsToCloseAgainst := $(TOP)/lib/$(machineSuffix)/libchem$(DIM)d.a
LibsToCloseAgainst += $(TOP)/lib/$(machineSuffix)/libiamr$(DIM)d.a
ifeq ($(USE_TENSOR),TRUE)
DEFINES += -DUSE_TENSOR
LibsToCloseAgainst += $(TOP)/lib/$(machineSuffix)/libmcmg$(DIM)d.a
LibsToCloseAgainst += $(TOP)/lib/$(machineSuffix)/libmg$(DIM)d.a
else
LibsToCloseAgainst += $(TOP)/lib/$(machineSuffix)/libmg$(DIM)d.a
endif
LibsToCloseAgainst += $(TOP)/lib/$(machineSuffix)/libamr$(DIM)d.a
LibsToCloseAgainst += $(TOP)/lib/$(machineSuffix)/libbndry$(DIM)d.a
ifeq ($(USE_HGPROJ_SERIAL),FALSE)
LibsToCloseAgainst += $(TOP)/lib/$(machineSuffix)/libproj$(DIM)d.a
else
LibsToCloseAgainst += $(TOP)/hgproj-serial/lib/$(machineSuffix)/libproj$(DIM)d.$(PRVERSION).a
endif
LibsToCloseAgainst += $(TOP)/lib/$(machineSuffix)/libbox$(DIM)d.a

$(optionsLib): $(LibsToCloseAgainst)
endif

include $(TOP)/mk/Make.rules
