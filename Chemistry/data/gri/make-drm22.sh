CONVERT=../../tools/convert/convert.exe
ln -fs drm22.dat drm22.mec
ln -fs thermo12.dat drm22.therm 
ln -fs transport12.dat drm22.trans
python ../../tools/fuego/Pythia/products/bin/fmc.py -mechanism=drm22.mec -thermo=drm22.therm -name=mec.cpp
echo Compiling drm22.cpp...
${CONVERT} model_files_drm22.dat
cat mec.cpp drm22-tran.cpp \
            ../header/header.start\
            ../header/header.mec   drm22.mec\
            ../header/header.therm drm22.therm\
            ../header/header.trans drm22.trans\
            ../header/header.end > drm22.cpp
rm -f mec.cpp drm22-tran.cpp
