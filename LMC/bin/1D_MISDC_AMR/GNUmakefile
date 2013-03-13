USE_INTEL = FALSE
EXEname := lmc.exe

#
# Choose pmf file
#
#pmf_source := ../../../Chemistry/data/gri/PMFs/gri30_070.f
pmf_source := ../../../Chemistry/data/chem-H/PMFs/chemHSoln_seed_0.00.f

#
# Choose reaction mechanism
# Use GRI30 for pmf_source := gri30_070.f
# Use CHEMH for pmf_source := chemHSoln_seed_0.00.f
# Use CHEMH for pmf_source := H-flame_0370_MCDD.f
#
#REACTION_MECHANISM=GRI30
REACTION_MECHANISM=CHEMH

ifeq (${USE_INTEL},TRUE)
  CCOMP := icpc -Wno-deprecated -g -Wunused-variable -DBL_ICC_VERSION=10.1 -DBL_MAJOR_VERSION=10 -DBL_MINOR_VERSION=1 -DBL_FORT_USE_UNDERSCORE
  fCOMP := ifort -warn unused -g -u -fpconstant -fbounds-check
  XTRALIBS :=  -L/opt/intel/fce/10.1.015/lib -lifcore -lm 
else
  fCOMP := gfortran -g -fbounds-check -O1 -Wall -ffixed-line-length-none
  CCOMP := g++ -g -O0 -fno-exceptions -Wall -Wno-sign-compare -DBL_FORT_USE_UNDERSCORE
  XTRALIBS :=  -lm -lgfortran
endif

ifeq ($(REACTION_MECHANISM),GRI30)
  mech_source += ../../../Chemistry/data/gri/Trash/FROM_LMC/gri30.cpp
endif
ifeq ($(REACTION_MECHANISM),CHEMH)
  mech_source += ../../../Chemistry/data/chem-H/chem-H.cpp
endif

chem_sources := CD.f ../../../Chemistry/src/vode.f ../../../Chemistry/src/EGini.f ../../../Chemistry/src/EGSlib.f

pmf_object := $(pmf_source:%.f=%.o)
mech_object := $(mech_source:%.cpp=%.o)
chem_objects := $(chem_sources:%.f=%.o)

f_sources += lmc.f ${chem_sources} prob.f fill_ghost.f \
             calc_divu.f get_temp_visc_terms.f get_spec_visc_terms.f \
             get_vel_visc_terms.f \
             get_rhoh_visc_terms.f project.f strang_chem.f \
             pre_mac_predict.f compute_pthermo.f add_dpdt.f \
             mac_proj.f scal_aofs.f mkslopes.f \
             update_rho.f update_spec.f update_temp.f update_rhoh.f \
             rhoh_to_temp.f advance.f vel_edge_states.f update_vel.f \
             read_check.f write_check.f write_plt.f cn_solve.f ${pmf_source}

f_includes += spec.h
C_sources += main.cpp ${mech_source}

C_objects := $(C_sources:%.cpp=%.o)
f_objects := $(f_sources:%.f=%.o)

${EXEname}: ${f_objects} ${C_objects} ${f_includes}
	${CCOMP} -o ${EXEname} ${f_objects} ${C_objects} ${XTRALIBS}

test.exe: CDc.o driver.o prob.o ${pmf_object} ${mech_object} ${chem_objects}
	${CCOMP} CDc.o driver.o prob.o ${pmf_object} ${mech_object} ${chem_objects} -o test.exe ${XTRALIBS}

clean:
	\rm -rf ${EXEname} ${f_objects} ${C_objects} Linkeg

%.o: %.cpp
	${CCOMP} -c $^ -o $*.o

%.o: %.f
	${fCOMP} -c $^ -o $*.o


MySrcDirs += . .. ../PMFfiles
vpath %.cpp $(MySrcDirs)
vpath %.h   $(MySrcDirs)
vpath %.f   $(MySrcDirs)
