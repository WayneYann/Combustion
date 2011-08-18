ln -fs chem-CH4-2step.inp chem-CH4-2step.mec
ln -fs thermo12.dat chem-CH4-2step.therm 
ln -fs transport12.dat chem-CH4-2step.trans
python ../../tools/fuego/Fuego/Pythia/products/bin/fmc.py -mechanism=chem-CH4-2step.mec -thermo=chem-CH4-2step.therm -name=mec.cpp
echo Compiling chem-CH4-2step.cpp...
cat mec.cpp ../header/header.start\
            ../header/header.mec   chem-CH4-2step.mec\
            ../header/header.therm chem-CH4-2step.therm\
            ../header/header.trans chem-CH4-2step.trans\
            ../header/header.end > chem-CH4-2step.cpp
rm -f mec.cpp
