#ln -sf sandiego20050310.mec sandiego20050310.mec
ln -sf sandiego20050310.new.therm-hack sandiego20050310.therm
ln -sf sandiego20021001.new.trans sandiego20050310.trans
python ../../tools/fuego/Fuego/Pythia/products/bin/fmc.py -mechanism=sandiego20050310.mec -thermo=sandiego20050310.therm -name=mec.cpp
echo Compiling sandiego20050310.cpp...
cat mec.cpp ../header/header.start\
            ../header/header.mec   sandiego20050310.mec\
            ../header/header.therm sandiego20050310.therm\
            ../header/header.trans sandiego20050310.trans\
            ../header/header.end > sandiego20050310.cpp
rm -f mec.cpp
