USE_INTEL = FALSE
USE_INTEL = TRUE
EXEname := lmc.exe

#
# Choose reaction mechanism (note, pmf file corresponding to this must be set below....)
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
  CCOMP := icpc -Wno-deprecated -g -O3 -DBL_ICC_VERSION=10.1 -DBL_MAJOR_VERSION=10 -DBL_MINOR_VERSION=1 -DBL_FORT_USE_UNDERSCORE -diag-disable vec
  fCOMP := ifort -g -u -fpconstant -132 -O3 -check uninit -diag-disable vec
  XTRALIBS :=  -L/opt/intel/fce/10.1.015/lib -lifcore -lm 
else
  fCOMP := gfortran -g -fbounds-check -O -Wuninitialized -ffixed-line-length-132 -Warray-bounds -Waliasing -Wconversion -Wextra
  CCOMP := g++ -DBL_FORT_USE_UNDERSCORE
  XTRALIBS := $(shell gfortran -print-file-name=libgfortran.a) $(shell gfortran -print-libgcc-file-name) -lm -lg2c
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
  pmf_source := H-flame_0370.f
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
  pmf_source := pmf_LiDryer_060_01ATM.f
endif

chem_headers := spec.h eg.cmn
chem_sources := EGini.f EGSlib.f CD.f 

pmf_object := $(pmf_source:%.f=%.o)
mech_object := $(mech_source:%.cpp=%.o)
chem_objects := $(chem_sources:%.f=%.o)

f_sources += driver.f util.f coef.f prob.f ${pmf_source} ${chem_sources}
C_sources += main.cpp ${mech_source}

C_objects := $(C_sources:%.cpp=%.o)
f_objects := $(f_sources:%.f=%.o)

${EXEname}: ${f_objects} ${C_objects} ${f_includes}
	${CCOMP} -o ${EXEname} ${f_objects} ${C_objects} ${XTRALIBS}

clean:
	\rm -rf ${EXEname} ${f_objects} ${C_objects}

%.o: %.cpp
	${CCOMP} -c $^ -o $*.o

%.o: %.f
	${fCOMP} -c $^ -o $*.o


MySrcDirs += . ../.. ../../PMFfiles
vpath %.cpp $(MySrcDirs)
vpath %.h   $(MySrcDirs)
vpath %.f   $(MySrcDirs)
