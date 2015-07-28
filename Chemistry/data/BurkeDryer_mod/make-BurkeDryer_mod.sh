CONVERT=../../tools/convert/convert.exe
python ../../tools/fuego/Pythia/products/bin/fmc.py -mechanism=chem.inp -thermo=therm.dat -name=mec.cpp
echo Compiling BurkeDryer_mod.cpp...
${CONVERT} model_files_BurkeDryer_mod.dat
cat mec.cpp BurkeDryer_mod-tran.c \
            ../header/header.start\
            ../header/header.mec   chem.inp\
            ../header/header.therm therm.dat\
            ../header/header.trans tran.dat\
            ../header/header.end > BurkeDryer_mod.c
rm -f mec.cpp BurkeDryer_mod-tran.c
