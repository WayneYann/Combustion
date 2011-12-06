CONVERT=../../tools/convert/convert.exe
python ../../tools/fuego/Pythia/products/bin/fmc.py -mechanism=LiDryer.mec -thermo=LiDryer.therm -name=mec.cpp
echo Compiling LiDryer.cpp...
${CONVERT} model_files_LiDryer.dat
cat mec.cpp LiDryer-tran.cpp \
            ../header/header.start\
            ../header/header.mec   LiDryer.mec\
            ../header/header.therm LiDryer.therm\
            ../header/header.trans LiDryer.trans\
            ../header/header.end > LiDryer.cpp
rm -f mec.cpp LiDryer-tran.cpp
