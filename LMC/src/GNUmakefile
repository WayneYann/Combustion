#
# $Id: GNUmakefile,v 1.3 2010-06-20 01:35:02 almgren Exp $
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
USE_WINDOWS   = FALSE
USE_MPI       = TRUE
#
# What stencil do you want to use?
#
PRVERSION = v9

#
# The base name of the library we're building.
#
LBASE = ht

include $(TOP)/mk/Make.defs ./Make.package

INCLUDE_LOCATIONS += . $(TOP)/include

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
LibsToCloseAgainst += $(TOP)/lib/$(machineSuffix)/libproj$(DIM)d.a
LibsToCloseAgainst += $(TOP)/lib/$(machineSuffix)/libbox$(DIM)d.a

$(optionsLib): $(LibsToCloseAgainst)
endif

include $(TOP)/mk/Make.rules
