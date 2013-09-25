
F90 ?= gfortran

FFLAGS = -Ibuild -Jbuild -g -Wall -Wno-unused-dummy-argument

BDF = ../../../src_common/bdf.f90 ../../../src_common/LinAlg.f

%.exe: %.f90 $(BDF)
	@mkdir -p build
	$(F90) $(FFLAGS) $^ -o $@

all: coeffs.exe t1.exe