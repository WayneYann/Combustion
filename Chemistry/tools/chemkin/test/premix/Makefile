PREBIN =  ../../bin/premix.exe
CKBIN = ../../bin/ckinterp.exe
TRBIN = ../../bin/tranfit.exe

all: out

clean:
	rm -f *.asc *.bak *.bin *.out

out:
	$(CKBIN) < chem.inp > chem.out
	$(TRBIN) > tran.out
	$(PREBIN) < premix.inp | tee premix.out
#	rm -f recov.bin rest.bin save.bin
