python ../../tools/fuego/Fuego/Pythia/products/bin/fmc.py -mechanism=prf_ethanol.mec -thermo=prf_ethanol.therm-hack -name=mec.cpp
echo Compiling prf_ethanol.cpp...
cat mec.cpp ../header/header.start\
            ../header/header.mec   prf_ethanol.mec\
            ../header/header.therm prf_ethanol.therm-hack\
            ../header/header.trans prf_ethanol.trans\
            ../header/header.end > prf_ethanol.cpp
rm -f mec.cpp
