CONVERT=../../tools/convert/convert.exe
python ../../tools/fuego/Pythia/products/bin/fmc.py -mechanism=chem.inp -thermo=therm.dat -name=mec.cpp
echo Compiling air.cpp...
${CONVERT} model_files_air.dat
cat mec.cpp air-tran.c \
            ../header/header.start\
            ../header/header.mec   chem.inp\
            ../header/header.therm therm.dat\
            ../header/header.trans tran.dat\
            ../header/header.end > air.c
rm -f mec.cpp air-tran.c
