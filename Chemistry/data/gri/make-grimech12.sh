CONVERT=../../tools/convert/convert.exe
ln -fs grimech12.dat grimech12.mec
ln -fs thermo12.dat grimech12.therm 
ln -fs transport12.dat grimech12.trans
python ../../tools/fuego/Pythia/products/bin/fmc.py -mechanism=grimech12.mec -thermo=grimech12.therm -name=mec.cpp
echo Compiling grimech12.cpp...
${CONVERT} model_files_grimech12.dat
cat mec.cpp grimech12-tran.cpp \
            ../header/header.start\
            ../header/header.mec   grimech12.mec\
            ../header/header.therm grimech12.therm\
            ../header/header.trans grimech12.trans\
            ../header/header.end > grimech12.cpp
rm -f mec.cpp grimech12-tran.cpp
