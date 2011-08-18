python ../../tools/fuego/Fuego/Pythia/products/bin/fmc.py -mechanism=KonnovH2.mec-hack -thermo=KonnovH2.therm -name=mec.cpp
echo Compiling KonnovH2.cpp...
cat mec.cpp ../header/header.start\
            ../header/header.mec   KonnovH2.mec-hack\
            ../header/header.therm KonnovH2.therm\
            ../header/header.trans KonnovH2.trans\
            ../header/header.end > KonnovH2.cpp
rm -f mec.cpp
