
F90    ?= gfortran
FFLAGS  = -Ibuild -Jbuild -g -Wall -Wno-unused-dummy-argument

vpath %.f90 ../../../src_common
vpath %.f   ../../../src_common

all: coeffs.exe t1.exe

#
# rules
#

%.exe: %.f90 build/bdf.o build/LinAlg.o
	$(F90) $(FFLAGS) $^ -o $@

build/%.o: %.f
	@mkdir -p build
	$(F90) -c $(FFLAGS) $^ -o $@

build/%.o: %.f90
	@mkdir -p build
	$(F90) -c $(FFLAGS) $^ -o $@
