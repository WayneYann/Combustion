CONVERT=../../tools/convert/convert.exe
python ../../tools/fuego/Pythia/products/bin/fmc.py -mechanism=inert_3.dat -thermo=inert_3.therm -name=mec.cpp
echo Compiling inert_3.cpp...
${CONVERT} model_files_inert_3.dat
cat mec.cpp inert_3-tran.cpp \
            ../header/header.start\
            ../header/header.mec   inert_3.dat\
            ../header/header.therm inert_3.therm\
            ../header/header.trans inert_3.trans\
            ../header/header.end > inert_3.c
rm -f mec.cpp inert_3-tran.cpp
