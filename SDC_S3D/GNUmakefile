# BOXLIB_HOME defines the directory in which we will find all the BoxLib code
# If you set BOXLIB_HOME as an environment variable, this line will be ignored

BOXLIB_HOME ?= $(HOME)/Development/BoxLib
#PFASST_HOME ?= $(HOME)/Development/f90pfasst

NDEBUG    := t
MPI       := t
OMP       := t
PROF      :=
COMP      := gfortran
MKVERBOSE := 

include $(BOXLIB_HOME)/Tools/F_mk/GMakedefs.mak

# F90 = mpif90

include ./GPackage.mak
VPATH_LOCATIONS += .

include $(BOXLIB_HOME)/Src/F_BaseLib/GPackage.mak
VPATH_LOCATIONS += $(BOXLIB_HOME)/Src/F_BaseLib

# include $(PFASST_HOME)/mk/GPackage.mak
# VPATH_LOCATIONS += $(PFASST_HOME)/src

main.$(suf).exe: $(objects)
	$(LINK.f90) -o main.$(suf).exe $(objects) $(libraries)

include $(BOXLIB_HOME)/Tools/F_mk/GMakerules.mak
