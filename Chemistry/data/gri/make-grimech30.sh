ln -sf grimech30.dat grimech30.mec
ln -sf thermo30.dat grimech30.therm
ln -sf transport30.dat grimech30.trans
python ../../tools/fuego/Fuego/Pythia/products/bin/fmc.py -mechanism=grimech30.mec -thermo=grimech30.therm -name=mec.cpp
echo Compiling grimech30.cpp...
cat mec.cpp ../header/header.start\
            ../header/header.mec   grimech30.mec\
            ../header/header.therm grimech30.therm\
            ../header/header.trans grimech30.trans\
            ../header/header.end > grimech30.cpp
rm -f mec.cpp
