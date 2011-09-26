CONVERT=../../tools/convert/convert.exe
ln -fs grimech30-noArN.dat grimech30-noArN.mec
ln -fs thermo30-noArN.dat grimech30-noArN.therm 
ln -fs transport30-noArN.dat grimech30-noArN.trans
python ../../tools/fuego/Pythia/products/bin/fmc.py -mechanism=grimech30-noArN.mec -thermo=grimech30.therm -name=mec.cpp
echo Compiling grimech30-noArN.cpp...
${CONVERT} model_files_grimech30-noArN.dat
cat mec.cpp grimech30-noArN-tran.cpp \
            ../header/header.start\
            ../header/header.mec   grimech30-noArN.mec\
            ../header/header.therm grimech30.therm\
            ../header/header.trans grimech30.trans\
            ../header/header.end > grimech30-noArN.cpp
rm -f mec.cpp grimech30-noArN-tran.cpp
