ln -sf grimech12.dat grimech12.mec
ln -sf thermo12.dat grimech12.therm
ln -sf transport12.dat grimech12.trans
python ../../tools/fuego/Fuego/Pythia/products/bin/fmc.py -mechanism=grimech12.mec -thermo=grimech12.therm -name=mec.cpp
echo Compiling grimech12.cpp...
cat mec.cpp ../header/header.start\
            ../header/header.mec   grimech12.mec\
            ../header/header.therm grimech12.therm\
            ../header/header.trans grimech12.trans\
            ../header/header.end > grimech12.cpp
rm -f mec.cpp
