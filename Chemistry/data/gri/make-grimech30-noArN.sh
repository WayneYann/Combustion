ln -sf grimech30-noArN.dat grimech30-noArN.mec
ln -sf thermo30.dat grimech30-noArN.therm
ln -sf transport30.dat grimech30-noArN.trans
python ../../tools/fuego/Fuego/Pythia/products/bin/fmc.py -mechanism=grimech30-noArN.mec -thermo=grimech30-noArN.therm -name=mec.cpp
echo Compiling grimech30-noArN.cpp...
cat mec.cpp ../header/header.start\
            ../header/header.mec   grimech30-noArN.mec\
            ../header/header.therm grimech30-noArN.therm\
            ../header/header.trans grimech30-noArN.trans\
            ../header/header.end > grimech30-noArN.cpp
rm -f mec.cpp
