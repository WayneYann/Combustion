CONVERT=../../tools/convert/convert.exe
python ../../tools/fuego/Pythia/products/bin/fmc.py -mechanism=chem.inp -thermo=therm.dat -name=mec.cpp
echo Compiling LuDME.cpp...
${CONVERT} model_files_LuDME.dat
cat mec.cpp LuDME-tran.cpp \
            ../header/header.start\
            ../header/header.mec   chem.inp\
            ../header/header.therm therm.dat\
            ../header/header.trans tran.dat\
            ../header/header.end > LuDME.cpp
rm -f mec.cpp LuDME-tran.cpp
