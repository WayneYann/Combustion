python ../../tools/fuego/Fuego/Pythia/products/bin/fmc.py -mechanism=Roehl.mec -thermo=Roehl.therm-hack -name=mec.cpp
echo Compiling Roehl.cpp...
cat mec.cpp ../header/header.start\
            ../header/header.mec   Roehl.mec\
            ../header/header.therm Roehl.therm-hack\
            ../header/header.trans Roehl.trans\
            ../header/header.end > Roehl.cpp
rm -f mec.cpp
