CHEMKIN_DIR := ../chemkin/src
CK_sources := ckinterp.f cklib_d.f cktime.f xerror.f
local_f77_sources =  tranfit_d.f 
local_f90_sources = strUtils.f90

CK_objects := $(CK_sources:%.f=%.o)
local_f77_objects := $(local_f77_sources:%.f=%.o)
local_f90_objects := $(local_f90_sources:%.f90=%.o)

vpath %.f . ${CHEMKIN_DIR}

all: convert.exe

convert.exe: ${CK_objects} ${local_f90_objects} ${local_f77_objects} convert.o
	gfortran ${CK_objects} ${local_f77_objects}  ${local_f90_objects} convert.o -o convert.exe

tranfit_d.o: tranfit_d.f string_utility.mod

string_utility.mod: ${local_f90_objects}

%.o: %.f
	gfortran -ffixed-line-length-100 -c $<

%.o: %.f90
	gfortran -c $<

clean:
	\rm -rf ${CK_objects} ${local_f77_objects} ${local_f90_objects} string_utility.mod *~ convert.o convert.exe
