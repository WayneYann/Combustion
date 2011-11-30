BIN = ../../bin
DATA = ../../../../data/chem-H

CHEMFILE=chem-H.inp
THERMFILE=thermo12.dat
TRANFILE=transport12.dat

CKINTERP=ckinterp.exe
TRANFIT=tranfit.exe
PREMIX=premix.exe

all: out

clean:
	rm -f *.asc *.bak *.bin *.dat *.out

out: premix.out

tranlink: ${DATA}/${TRANFILE}
	\rm -rf tran.dat
	ln -s ${DATA}/${TRANFILE} tran.dat

thermlink: ${DATA}/${THERMFILE}
	\rm -rf therm.dat
	ln -s ${DATA}/${THERMFILE} therm.dat

chem.out: ${BIN}/${CKINTERP} thermlink ${DATA}/${CHEMFILE}
	${BIN}/${CKINTERP} < ${DATA}/${CHEMFILE} > chem.out

tran.out: chem.out ${BIN}/${TRANFIT} tranlink
	${BIN}/${TRANFIT} > tran.out

premix.out: tran.out ${BIN}/${PREMIX} 
	${BIN}/${PREMIX} < premix.inp > premix.out
