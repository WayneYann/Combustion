# A set of useful macros for putting together a SMC application.

# include the main Makefile stuff
include $(BOXLIB_HOME)/Tools/F_mk/GMakedefs.mak

# default target (make just takes the one that appears first)
ALL: main.$(suf).exe


#-----------------------------------------------------------------------------
# core BoxLib directories
BOXLIB_CORE := Src/F_BaseLib 


#-----------------------------------------------------------------------------
# core SMC directories
SMC_CORE := src


#-----------------------------------------------------------------------------
# Fmpack is the list of all the GPackage.mak files that we need to
# include into the build system to define the list of source files.
#
# Fmlocs is the list of all the directories that we want to search
# for the source files in -- this is usually going to be the
# same as the list of directories containing GPackage.mak defined
# above.
#
# Fincs is the list of directories that have include files that
# we need to tell the compiler about.

# SMC modules
Fmdirs += $(SMC_CORE)

Fmpack := $(foreach dir, $(Fmdirs), $(SMC_TOP_DIR)/$(dir)/GPackage.mak)
Fmlocs := $(foreach dir, $(Fmdirs), $(SMC_TOP_DIR)/$(dir))

# BoxLib
Fmpack += $(foreach dir, $(BOXLIB_CORE), $(BOXLIB_HOME)/$(dir)/GPackage.mak)
Fmlocs += $(foreach dir, $(BOXLIB_CORE), $(BOXLIB_HOME)/$(dir))


# any include directories
Fmincs := $(foreach dir, $(Fmincludes), $(SMC_TOP_DIR)/$(dir))


# include the necessary GPackage.mak files that define this setup
include $(Fmpack)

# vpath defines the directories to search for the source files

#  Note: GMakerules.mak will include '.' at the start of the
#  VPATH_LOCATIONS to first search in the problem directory
VPATH_LOCATIONS += $(Fmlocs)


# list of directories to put in the Fortran include path
FINCLUDE_LOCATIONS += $(Fmincs)


#-----------------------------------------------------------------------------
# define the build instructions for the executable
main.$(suf).exe: $(objects)
	$(HPCLINK) $(LINK.f90) -o main.$(suf).exe $(objects) $(libraries)
	@echo SUCCESS


#-----------------------------------------------------------------------------
# include the BoxLib Fortran Makefile rules
include $(BOXLIB_HOME)/Tools/F_mk/GMakerules.mak


#-----------------------------------------------------------------------------
# for debugging.  To see the value of a Makefile variable,
# e.g. Fmlocs, simply do "make print-Fmlocs".  This will
# print out the value.
print-%: ; @echo $* is $($*)

