NDEBUG := 
MPI    := 
OMP    :=

COMP := Intel
MKVERBOSE := T

FPARALLEL := ../../../../../VARDENS/LRvarden/VARDEN/fParallel

include $(FPARALLEL)/mk/GMakedefs.mak

#f90sources += strang_main.f90
#f90sources += subt_main.f90
f90sources += main.f90

include GPackage.mak

include $(FPARALLEL)/boxlib/GPackage.mak
VPATH_LOCATIONS += $(FPARALLEL)/boxlib

main.$(suf).exe: $(objects)
	$(LINK.f90) -o main.$(suf).exe $(objects) $(libraries)

include $(FPARALLEL)/mk/GMakerules.mak

