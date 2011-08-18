ln -sf mech.txt Hai.mec
ln -sf Thermo.txt Hai.therm
ln -sf trandat.txt Hai.trans
python ../../tools/fuego/Fuego/Pythia/products/bin/fmc.py -mechanism=Hai.mec -thermo=Hai.therm-hack -name=mec.cpp
echo Compiling Hai.cpp...
cat mec.cpp ../header/header.start\
            ../header/header.mec   Hai.mec\
            ../header/header.therm Hai.therm-hack\
            ../header/header.trans Hai.trans\
            ../header/header.end > Hai.cpp
rm -f mec.cpp
