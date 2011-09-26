CONVERT=../../tools/convert/convert.exe
ln -fs grimech30.dat grimech30.mec
ln -fs thermo30.dat grimech30.therm 
ln -fs transport30.dat grimech30.trans
python ../../tools/fuego/Pythia/products/bin/fmc.py -mechanism=grimech30.mec -thermo=grimech30.therm -name=mec.cpp
echo Compiling grimech30.cpp...
${CONVERT} model_files_grimech30.dat
cat mec.cpp grimech30-tran.cpp \
            ../header/header.start\
            ../header/header.mec   grimech30.mec\
            ../header/header.therm grimech30.therm\
            ../header/header.trans grimech30.trans\
            ../header/header.end > grimech30.cpp
rm -f mec.cpp grimech30-tran.cpp
