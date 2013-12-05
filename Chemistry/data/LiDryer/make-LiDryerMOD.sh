CONVERT=../../tools/convert/convert.exe
python ../../tools/fuego/Pythia/products/bin/fmc.py -mechanism=LiDryerMOD.mec -thermo=LiDryer.therm -name=mec.cpp
echo Compiling LiDryerMOD.c...
${CONVERT} model_files_LiDryerMOD.dat
cat mec.cpp LiDryerMOD-tran.c \
            ../header/header.start\
            ../header/header.mec   LiDryerMOD.mec\
            ../header/header.therm LiDryer.therm\
            ../header/header.trans LiDryer.trans\
            ../header/header.end > LiDryerMOD.c
rm -f mec.cpp LiDryerMOD-tran.c
