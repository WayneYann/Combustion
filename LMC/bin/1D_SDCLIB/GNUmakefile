USE_INTEL = FALSE
EXEname := lmc.exe

SDCLIB_HOME = ../../../../SDCLib

#
# Choose pmf file
#
#pmf_source := ../../../Chemistry/data/gri/PMFs/gri30_070.f
pmf_source := ../../../Chemistry/data/chem-H/PMFs/chemHSoln_seed_0.00.f
#pmf_source := ../../../Chemistry/data/Lu/PMFs/LuDME_0700.f
#
# Choose reaction mechanism
# Use GRI30 for pmf_source := gri30_070.f
# Use CHEMH for pmf_source := chemHSoln_seed_0.00.f
# Use LUDME for pmf_source := LuDME_0700.f
#
#REACTION_MECHANISM=GRI30
REACTION_MECHANISM=CHEMH
#REACTION_MECHANISM=LUDME

ifeq (${USE_INTEL},TRUE)
  CCOMP := icpc -Wno-deprecated -g -Wunused-variable -DBL_ICC_VERSION=10.1 -DBL_MAJOR_VERSION=10 -DBL_MINOR_VERSION=1 -DBL_FORT_USE_UNDERSCORE
  fCOMP := ifort -warn unused -g -u -fpconstant -fbounds-check
  # XTRALIBS :=  -L/opt/intel/fce/10.1.015/lib -lifcore -lm 
else
  fCOMP := gfortran -g -fbounds-check -O1 -Wall -ffixed-line-length-none
  CCOMP := g++ -g -O0 -fno-exceptions -Wall -Wno-sign-compare -DBL_FORT_USE_UNDERSCORE
  cCOMP := gcc -g -O0 -std=c99 -fno-exceptions -Wall -Wno-sign-compare -DBL_FORT_USE_UNDERSCORE -I$(SDCLIB_HOME)/include
  # XTRALIBS :=  -lm -lgfortran
endif

ifeq ($(REACTION_MECHANISM),GRI30)
  mech_source += ../../../Chemistry/data/gri/grimech30.c
endif
ifeq ($(REACTION_MECHANISM),CHEMH)
  mech_source += ../../../Chemistry/data/chem-H/chem-H.c
endif
ifeq ($(REACTION_MECHANISM),LUDME)
  mech_source += ../../../Chemistry/data/Lu/LuDME.c
endif

chem_sources := CD.f ../../../Chemistry/src/vode.f ../../../Chemistry/src/EGini.f ../../../Chemistry/src/EGSlib.f

pmf_object := $(pmf_source:%.f=%.o)
mech_object := $(mech_source:%.c=%.o)
chem_objects := $(chem_sources:%.f=%.o)

f_sources += lmc.f ${chem_sources} prob.f fill_ghost.f \
             calc_divu.f get_temp_visc_terms.f get_spec_visc_terms.f \
             get_vel_visc_terms.f \
             get_rhoh_visc_terms.f project.f strang_chem.f \
             pre_mac_predict.f compute_pthermo.f add_dpdt.f \
             mac_proj.f scal_aofs.f mkslopes.f \
             update_rho.f \
             rhoh_to_temp.f advance.f vel_edge_states.f update_vel.f \
             read_check.f write_check.f write_plt.f cn_solve.f ${pmf_source}

f_includes += spec.h
c_sources += main.c
f90sources += feval.f90
f90sources += sdc-interfaces.f90
f90sources += sdc-types.f90
f90sources += sdclib.f90

include $(SDCLIB_HOME)/Make.package
c_sources += $(cEXE_sources)

c_objects := $(c_sources:%.c=%.o) 
f_objects := $(f_sources:%.f=%.o)
f90_objects := $(f90sources:%.f90=%.o)


${EXEname}: ${f_objects} ${f90_objects} ${mech_object} ${c_objects} #${f_includes}
	${fCOMP} -o ${EXEname} $^

clean:
	rm -rf ${EXEname} ${f_objects} ${c_objects} Linkeg

%.o: %.cpp
	${CCOMP} -c $^ -o $*.o

%.o: %.c
	${cCOMP} -c $^ -o $*.o

%.o: %.f
	${fCOMP} -c $^ -o $*.o

%.o: %.f90
	${fCOMP} -c $^ -o $*.o

lmc.o: sdclib.o feval.o
sdclib.o: sdc-interfaces.o
sdc-interfaces.o: sdc-types.o


MySrcDirs += .
vpath %.cpp $(MySrcDirs)
vpath %.c   $(MySrcDirs) $(SDCLIB_HOME)/src
vpath %.h   $(MySrcDirs)
vpath %.f   $(MySrcDirs)
vpath %.f90 $(MySrcDirs) $(SDCLIB_HOME)/src 
