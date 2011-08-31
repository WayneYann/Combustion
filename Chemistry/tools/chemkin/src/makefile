BIN=../bin
PREMIX=${BIN}/premix.exe
TRANFIT=${BIN}/tranfit.exe
CKINTERP=${BIN}/ckinterp.exe

SRC_PR=premix.f prdriv.f vode.f xerror.f math.f twopnt.f cklib_d.f tranlib_d.f mach.f refine.f cktime.f
SRC_TR=tranfit.f cklib_d.f xerror.f mcdriv.f cktime.f
SRC_CK=ckinterp.f cklib_d.f ckdriv.f cktime.f


OBJ_PR = $(SRC_PR:%.f=%.o)
OBJ_TR = $(SRC_TR:%.f=%.o)
OBJ_CK = $(SRC_CK:%.f=%.o)

%.o: %.f
	ifort -c $*.f

all: ${PREMIX} ${TRANFIT} ${CKINTERP}

${PREMIX}: ${OBJ_PR}
	ifort ${OBJ_PR} -o $@

${TRANFIT}: ${OBJ_TR}
	ifort ${OBJ_TR} -o $@

${CKINTERP}: ${OBJ_CK}
	ifort ${OBJ_CK} -o $@

clean:
	\rm -rf ${OBJ_PR} ${OBJ_TR} ${OBJ_CK} 


