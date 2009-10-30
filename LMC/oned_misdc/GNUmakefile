
CCOMP := g++ -DBL_FORT_USE_UNDERSCORE
fCOMP := g77
XTRALIBS :=  -lm -lg2c
EXEname := driver.exe

f_sources += CD.f vode.f EGini.f EGSlib.f
f_includes += spec.h
C_sources += CDc.cpp


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


MySrcDirs += . ..
vpath %.cpp $(MySrcDirs)
vpath %.h   $(MySrcDirs)
vpath %.f   $(MySrcDirs)
