CONVERT=../../tools/convert/convert.exe
python ../../tools/fuego/Pythia/products/bin/fmc.py -mechanism=chem.inp -thermo=therm.dat -name=mec.cpp
echo Compiling BurkeDryer.cpp...
${CONVERT} model_files_BurkeDryer.dat
cat mec.cpp BurkeDryer-tran.c \
            ../header/header.start\
            ../header/header.mec   chem.inp\
            ../header/header.therm therm.dat\
            ../header/header.trans tran.dat\
            ../header/header.end > BurkeDryer.c
rm -f mec.cpp BurkeDryer-tran.c
