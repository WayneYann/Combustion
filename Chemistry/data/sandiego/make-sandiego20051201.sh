#ln -sf sandiego20051201.mec sandiego20051201.mec
ln -sf sandiego20050310.new.therm-hack sandiego20051201.therm
ln -sf sandiego20021001.new.trans sandiego20051201.trans
python ../../tools/fuego/Fuego/Pythia/products/bin/fmc.py -mechanism=sandiego20051201.mec -thermo=sandiego20051201.therm -name=mec.cpp
echo Compiling sandiego20051201.cpp...
cat mec.cpp ../header/header.start\
            ../header/header.mec   sandiego20051201.mec\
            ../header/header.therm sandiego20051201.therm\
            ../header/header.trans sandiego20051201.trans\
            ../header/header.end > sandiego20051201.cpp
rm -f mec.cpp
