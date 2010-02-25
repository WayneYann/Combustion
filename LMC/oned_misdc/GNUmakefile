USE_INTEL = TRUE
EXEname := lmc.exe

#
# Choose pmf file
#
#pmf_source := H-flame_0370.f
#pmf_source := 2stepSoln0.8.f
pmf_source := chemHSoln_seed_0.00.f

#
# Choose reaction mechanism
#
#REACTION_MECHANISM=INERT30
#REACTION_MECHANISM=GRI12
#REACTION_MECHANISM=GRI30NON
#REACTION_MECHANISM=GRI30
#REACTION_MECHANISM=CH4-2STEP
#REACTION_MECHANISM=GRI30HNO
#REACTION_MECHANISM=DRM19
#REACTION_MECHANISM=GLARBORG
REACTION_MECHANISM=CHEMH
#REACTION_MECHANISM=LIDRYER

ifeq (${USE_INTEL},TRUE)
#  CCOMP := icpc -vec_report0 -DBL_FORT_USE_UNDERSCORE
  CCOMP := icpc -Wno-deprecated -g -Wunused-variable -DBL_ICC_VERSION=10.1 -DBL_MAJOR_VERSION=10 -DBL_MINOR_VERSION=1 -DBL_FORT_USE_UNDERSCORE
#  fCOMP := ifort -debug full -vec-report0
  fCOMP := ifort -warn unused -g -u -fpconstant -fbounds-check
#  fCOMP += -CB
  XTRALIBS :=  -L/opt/intel/fce/10.1.015/lib -lifcore -lm 
else
  fCOMP := g77 -g -fbounds-check -O -Wuninitialized 
  CCOMP := g++ -DBL_FORT_USE_UNDERSCORE
  XTRALIBS :=  -lm -lg2c
endif

ifeq ($(REACTION_MECHANISM),GLARBORG)
  mech_source += glarSkel.cpp
endif
ifeq ($(REACTION_MECHANISM),DRM19)
  mech_source += drm19.cpp
endif
ifeq ($(REACTION_MECHANISM),GRI12)
  mech_source += gri12.cpp
endif
ifeq ($(REACTION_MECHANISM),GRI30)
  mech_source += gri30.cpp
endif
ifeq ($(REACTION_MECHANISM),GRI30NON)
  mech_source += gri30-noN.cpp
endif
ifeq ($(REACTION_MECHANISM),GRI30HNO)
  mech_source += gri30_hno.cpp
endif
ifeq ($(REACTION_MECHANISM),GRI211HNO)
  mech_source += gri211-hno.cpp
endif
ifeq ($(REACTION_MECHANISM),CHEMH)
  mech_source += chem-H.cpp
endif
ifeq ($(REACTION_MECHANISM),INERT30)
  mech_source += inert30.cpp
endif
ifeq ($(REACTION_MECHANISM),CH4-2STEP)
  mech_source += chem-CH4-2step.cpp
endif
ifeq ($(REACTION_MECHANISM),PROPANE)
  mech_source += propane.cpp
endif
ifeq ($(REACTION_MECHANISM),LIDRYER)
  mech_source += LiDryer.cpp
endif

chem_sources := CD.f vode.f EGini.f EGSlib.f

pmf_object := $(pmf_source:%.f=%.o)
mech_object := $(mech_source:%.cpp=%.o)
chem_objects := $(chem_sources:%.f=%.o)

f_sources += lmc.f ${chem_sources} prob.f \
             calc_divu.f get_temp_visc_terms.f get_spec_visc_terms.f \
             get_vel_visc_terms.f divRhoDHgradY.f\
             get_rhoh_visc_terms.f est_dt.f project.f strang_chem.f \
             pre_mac_predict.f compute_pthermo.f add_dpdt.f \
             mac_proj.f scal_aofs.f mkslopes.f bdsslopes.f \
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
	\rm -rf ${EXEname} ${f_objects} ${C_objects}

%.o: %.cpp
	${CCOMP} -c $^ -o $*.o

%.o: %.f
	${fCOMP} -c $^ -o $*.o


MySrcDirs += . .. ../PMFfiles
vpath %.cpp $(MySrcDirs)
vpath %.h   $(MySrcDirs)
vpath %.f   $(MySrcDirs)
