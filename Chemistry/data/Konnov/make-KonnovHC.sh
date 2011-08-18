python ../../tools/fuego/Fuego/Pythia/products/bin/fmc.py -mechanism=KonnovHC.mec-hack -thermo=KonnovHC.therm-hack -name=mec.cpp
echo Compiling KonnovHC.cpp...
cat mec.cpp ../header/header.start\
            ../header/header.mec   KonnovHC.mec-hack\
            ../header/header.therm KonnovHC.therm-hack\
            ../header/header.trans KonnovHC.trans\
            ../header/header.end > KonnovHC.cpp
rm -f mec.cpp
