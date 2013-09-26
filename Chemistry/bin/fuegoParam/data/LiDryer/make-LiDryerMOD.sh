CONVERT=../../tools/convert/convert.exe
python ../../tools/fuego/Pythia/products/bin/fmc.py -mechanism=LiDryerMOD.mec -thermo=LiDryer.therm -name=mec.cpp
echo Compiling LiDryerMOD.cpp...
${CONVERT} model_files_LiDryerMOD.dat
cat mec.cpp LiDryerMOD-tran.cpp \
            ../header/header.start\
            ../header/header.mec   LiDryerMOD.mec\
            ../header/header.therm LiDryer.therm\
            ../header/header.trans LiDryer.trans\
            ../header/header.end > LiDryerMOD.cpp
rm -f mec.cpp LiDryerMOD-tran.cpp
