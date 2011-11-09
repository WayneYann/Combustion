BIN = ../../../bin
DATA = ../../../../../data/chem-H

CHEMFILE=chem-H.inp
THERMFILE=thermo12.dat
TRANFILE=transport12.dat

CKINTERP=ckinterp.exe
TRANFIT=tranfit.exe
OPPDIFF=oppdiff.exe

all: out

clean:
	rm -f *.asc *.bak *.bin *.dat *.out

out: oppdif.out

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

oppdif.out: tran.out ${BIN}/${OPPDIFF} 
	${BIN}/${OPPDIFF} < oppdif.inp > oppdif.out
