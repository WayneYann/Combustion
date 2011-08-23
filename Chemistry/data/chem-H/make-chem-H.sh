CONVERT=../../tools/convert/convert.exe
ln -fs chem-H.inp chem-H.mec
ln -fs thermo12.dat chem-H.therm 
ln -fs transport12.dat chem-H.trans
python ../../tools/fuego/Pythia/products/bin/fmc.py -mechanism=chem-H.mec -thermo=chem-H.therm -name=mec.cpp
echo Compiling chem-H.cpp...
${CONVERT} model_files_chem-H.dat
cat mec.cpp chem-H-tran.cpp \
            ../header/header.start\
            ../header/header.mec   chem-H.mec\
            ../header/header.therm chem-H.therm\
            ../header/header.trans chem-H.trans\
            ../header/header.end > chem-H.cpp
rm -f mec.cpp chem-H-tran.cpp
