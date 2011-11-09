BIN=../bin
FC=gfortran
#FC=ifort
PREMIX=${BIN}/premix.exe
OPPDIFF=${BIN}/oppdiff.exe
TRANFIT=${BIN}/tranfit.exe
CKINTERP=${BIN}/ckinterp.exe

SRC_OP=opdriv.f oppdif.f xerror.f math.f twopnt.f cklib_d.f tranlib_d.f mach.f refine.f cktime.f
SRC_PR=premix.f prdriv.f vode.f xerror.f math.f twopnt.f cklib_d.f tranlib_d.f mach.f refine.f cktime.f
SRC_TR=tranfit.f cklib_d.f xerror.f mcdriv.f cktime.f
SRC_CK=ckinterp.f cklib_d.f ckdriv.f cktime.f


OBJ_OP = $(SRC_OP:%.f=%.o)
OBJ_PR = $(SRC_PR:%.f=%.o)
OBJ_TR = $(SRC_TR:%.f=%.o)
OBJ_CK = $(SRC_CK:%.f=%.o)

%.o: %.f
	${FC} -c $*.f

all: ${OPPDIFF} ${PREMIX} ${TRANFIT} ${CKINTERP}

bindir:
	mkdir -p ${BIN}

${OPPDIFF}: ${OBJ_OP} bindir
	${FC} ${OBJ_OP} -o $@

${PREMIX}: ${OBJ_PR} bindir
	${FC} ${OBJ_PR} -o $@

${TRANFIT}: ${OBJ_TR} bindir
	${FC} ${OBJ_TR} -o $@

${CKINTERP}: ${OBJ_CK} bindir
	${FC} ${OBJ_CK} -o $@

clean:
	\rm -rf ${OBJ_OP} ${OBJ_PR} ${OBJ_TR} ${OBJ_CK} 


