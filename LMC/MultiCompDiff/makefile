AMRLIB = /home/sepp/src/Parallel/amrlib
BOXLIB = /home/sepp/src/Parallel/BoxLib
LMC = /home/sepp/src/Parallel/LMC

all:

clean:
	rm -f *~

links: \
	BC_TYPES.H \
	cdwrk.H \
	ChemDriver_F.H \
	REAL.H

ArrayLim.H:
	ln -s $(BOXLIB)/ArrayLim.H ArrayLim.H

BC_TYPES.H:
	ln -s $(AMRLIB)/BC_TYPES.H BC_TYPES.H

cdwrk.H:
	ln -s $(LMC)/cdwrk.H cdwrk.H

ChemDriver_F.H:
	ln -s $(LMC)/ChemDriver_F.H ChemDriver_F.H

#CONSTANTS.H:
#	ln -s $(BOXLIB)/CONSTANTS.H CONSTANTS.H

HEATTRANSFER_F.H:
	ln -s $(LMC)/HEATTRANSFER_F.H HEATTRANSFER_F.H

REAL.H:
	ln -s $(BOXLIB)/REAL.H REAL.H

realclean:
	rm -r d f o

sav:
	rm -r sav_old
	mv sav_new sav_old
	cp -r sav_old sav_new
	cp *.cpp *.F *.H sav_new

expanded.F: FORT_RK_MULTICOMPONENT.F ArrayLim.H BC_TYPES.H cdwrk.H ChemDriver_F.H HEATTRANSFER_F.H REAL.H
	rm -f expanded.F messages
	cpp FORT_RK_MULTICOMPONENT.F > expanded.F

messages: expanded.F
	ifort -c -o -132 expanded.F >& messages
