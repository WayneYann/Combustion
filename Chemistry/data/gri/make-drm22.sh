ln -sf drm22.dat drm22.mec
ln -sf thermo12.dat drm22.therm
ln -sf transport12.dat drm22.trans
python ../../tools/fuego/Fuego/Pythia/products/bin/fmc.py -mechanism=drm22.mec -thermo=drm22.therm -name=mec.cpp
echo Compiling drm22.cpp...
cat mec.cpp ../header/header.start\
            ../header/header.mec   drm22.mec\
            ../header/header.therm drm22.therm\
            ../header/header.trans drm22.trans\
            ../header/header.end > drm22.cpp
rm -f mec.cpp
