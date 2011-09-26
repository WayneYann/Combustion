CONVERT=../../tools/convert/convert.exe
ln -fs drm19.dat drm19.mec
ln -fs thermo12.dat drm19.therm 
ln -fs transport12.dat drm19.trans
python ../../tools/fuego/Pythia/products/bin/fmc.py -mechanism=drm19.mec -thermo=drm19.therm -name=mec.cpp
echo Compiling drm19.cpp...
${CONVERT} model_files_drm19.dat
cat mec.cpp drm19-tran.cpp \
            ../header/header.start\
            ../header/header.mec   drm19.mec\
            ../header/header.therm drm19.therm\
            ../header/header.trans drm19.trans\
            ../header/header.end > drm19.cpp
rm -f mec.cpp drm19-tran.cpp
