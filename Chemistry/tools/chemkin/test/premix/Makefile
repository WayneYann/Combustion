BIN = ../../bin
DATA = ../../../../data/chem-H
DATA = ../../../../data/LiDryer

# Pick mechanism by setting relevant filenames

# -------- chem-H
CHEMFILE=chem-H.inp
THERMFILE=thermo12.dat
TRANFILE=transport12.dat

# -------- LiDryerMOD  (ie, LiDryer with the species order switched to match chem-H)
CHEMFILE=LiDryerMOD.mec
THERMFILE=LiDryer.therm
TRANFILE=LiDryer.trans


# Set names for executables
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
