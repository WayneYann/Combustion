SDCLIB_HOME = ../../../../SDCLib

#
# Choose PMF source
#
#pmf_source := ../../../Chemistry/data/gri/PMFs/gri30_070.f
pmf_source = chemHSoln_seed_0.00.f
#pmf_source := ../../../Chemistry/data/Lu/PMFs/LuDME_0700.f

#
# Choose reaction mechanism
#
# Use GRI30 for pmf_source := gri30_070.f
# Use CHEMH for pmf_source := chemHSoln_seed_0.00.f
# Use LUDME for pmf_source := LuDME_0700.f
#
#REACTION_MECHANISM=GRI30
REACTION_MECHANISM=CHEMH
#REACTION_MECHANISM=LUDME

CC = /home/memmett/gcc-4.7/bin/gcc
FC = /home/memmett/gcc-4.7/bin/gfortran

CFLAGS = -g -O0 -std=c99 -Wall -Wno-sign-compare -fno-exceptions -DBL_FORT_USE_UNDERSCORE
FFLAGS = -g -O1 -cpp -Wall -fbounds-check -Wno-unused-dummy-argument

FFLAGS  += -Ibuild -Jbuild
CFLAGS  += -I$(SDCLIB_HOME)/include

# FFLAGS  += -DZMQ
# CFLAGS  += -I$(HOME)/anaconda/include
# LDFLAGS += -L$(HOME)/anaconda/lib -lzmq
# csources   += zmq.c

ifeq ($(REACTION_MECHANISM),GRI30)
  mech_source += ../../../Chemistry/data/gri/grimech30.c
endif
ifeq ($(REACTION_MECHANISM),CHEMH)
  mech_source += chem-H.c
endif
ifeq ($(REACTION_MECHANISM),LUDME)
  mech_source += ../../../Chemistry/data/Lu/LuDME.c
endif

#chem_sources = CD.f vode.f EGini.f EGSlib.f
chem_sources = vode.f EGini.f EGSlib.f

f90sources += main.f90 feval.f90 encap.f90 debug.f90 probin.f90 lmc.f90 kernels.f90
f90sources += sdc-interfaces.f90 sdc-types.f90 sdclib.f90

csources += sdc-hooks.c
csources += sdc-imex-sweep.c
csources += sdc-imex.c
csources += sdc-logging.c
csources += sdc-mat.c
csources += sdc-mg-sweep.c
csources += sdc-mg.c
csources += sdc-mrex-sweep.c
csources += sdc-mrex.c
csources += sdc-nset.c
csources += sdc-poly.c
csources += sdc-print.c
csources += sdc-quadrature.c
csources += sdc-sweeper.c
csources += sdc-wrappers.c

f90sources += sdc-interfaces.f90
f90sources += sdclib.f90
f90sources += sdc-types.f90

csources += $(mech_source)
fsources += $(pmf_source) $(chem_sources)

OBJ = $(addprefix build/,$(csources:%.c=%.o) $(fsources:%.f=%.o) $(f90sources:%.f90=%.o))

lmc.exe: $(OBJ)
	$(FC) -o $@ $^ $(LDFLAGS)

clean:
	rm -rf *.exe build

# dependencies
build/lmc.o:     build/probin.o
build/main.o:    build/feval.o build/lmc.o build/probin.o
build/kernels.o: build/probin.o build/lmc.o
build/encap.o:   build/sdclib.o
build/feval.o:   build/encap.o build/debug.o build/lmc.o build/probin.o build/kernels.o
build/sdclib.o:         build/sdc-interfaces.o
build/sdc-interfaces.o: build/sdc-types.o

# rules
build/%.o: %.c
	@mkdir -p build
	$(CC) $(CFLAGS) -c $< -o $@

build/%.o: %.f
	@mkdir -p build
	$(FC) $(FFLAGS) -c $< -o $@

build/%.o: %.f90
	@mkdir -p build
	$(FC) $(FFLAGS) -c $< -o $@

vpath %.c   src $(SDCLIB_HOME)/src ../../../Chemistry/data/chem-H/
vpath %.f   src ../../../Chemistry/src ../../../Chemistry/data/chem-H/PMFs
vpath %.f90 src $(SDCLIB_HOME)/src
