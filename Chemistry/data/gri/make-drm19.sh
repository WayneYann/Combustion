ln -sf drm19.dat drm19.mec
ln -sf thermo12.dat drm19.therm
ln -sf transport12.dat drm19.trans
python ../../tools/fuego/Fuego/Pythia/products/bin/fmc.py -mechanism=drm19.mec -thermo=drm19.therm -name=mec.cpp
echo Compiling drm19.cpp...
cat mec.cpp ../header/header.start\
            ../header/header.mec   drm19.mec\
            ../header/header.therm drm19.therm\
            ../header/header.trans drm19.trans\
            ../header/header.end > drm19.cpp
rm -f mec.cpp
