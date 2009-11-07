
CCOMP := g++ -DBL_FORT_USE_UNDERSCORE
fCOMP := g77
XTRALIBS :=  -lm -lg2c
EXEname := lmc.exe

f_sources += lmc.f CD.f vode.f EGini.f EGSlib.f prob.f \
             calc_divu.f get_temp_visc_terms.f get_spec_visc_terms.f \
             get_rhoh_visc_terms.f est_dt.f project.f strang_chem.f \
             pre_mac_predict.f compute_pthermo.f add_dpdt.f \
             mac_proj.f scal_aofs.f mkslopes.f bdsslopes.f \
             update_rho.f update_spec.f update_temp.f update_rhoh.f \
             rhoh_to_temp.f advance.f vel_edge_states.f update_vel.f \
             read_check.f write_check.f cn_solve.f

f_includes += spec.h
C_sources += main.cpp

f_sources += H-flame_0370.f

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


ifeq ($(REACTION_MECHANISM),GLARBORG)
  C_sources += glarSkel.cpp
endif
ifeq ($(REACTION_MECHANISM),DRM19)
  C_sources += drm19.cpp
endif
ifeq ($(REACTION_MECHANISM),GRI12)
  C_sources += gri12.cpp
endif
ifeq ($(REACTION_MECHANISM),GRI30)
  C_sources += gri30.cpp
endif
ifeq ($(REACTION_MECHANISM),GRI30NON)
  C_sources += gri30-noN.cpp
endif
ifeq ($(REACTION_MECHANISM),GRI30HNO)
  C_sources += gri30_hno.cpp
endif
ifeq ($(REACTION_MECHANISM),GRI211HNO)
  C_sources += gri211-hno.cpp
endif
ifeq ($(REACTION_MECHANISM),CHEMH)
  C_sources += chem-H.cpp
endif
ifeq ($(REACTION_MECHANISM),INERT30)
  C_sources += inert30.cpp
endif
ifeq ($(REACTION_MECHANISM),CH4-2STEP)
  C_sources += chem-CH4-2step.cpp
endif
ifeq ($(REACTION_MECHANISM),PROPANE)
  C_sources += propane.cpp
endif
ifeq ($(REACTION_MECHANISM),LIDRYER)
  C_sources += LiDryer.cpp
endif

C_objects := $(C_sources:%.cpp=%.o)
f_objects := $(f_sources:%.f=%.o)


${EXEname}: ${f_objects} ${C_objects}
	${CCOMP} ${f_objects} ${C_objects} -o ${EXEname} ${XTRALIBS}

clean:
	\rm -rf ${EXEname} ${f_objects} ${C_objects}

%.o: %.cpp
	${CCOMP} -c -o $*.o $^

%.o: %.f
	${fCOMP} -c -o $*.o $^


MySrcDirs += . .. ../PMFfiles
vpath %.cpp $(MySrcDirs)
vpath %.h   $(MySrcDirs)
vpath %.f   $(MySrcDirs)
