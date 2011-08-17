ln -sf sandiego20050310.mec sandiego.mec
ln -sf sandiego20050310.new.therm-hack sandiego.therm
ln -sf sandiego20021001.new.trans sandiego.trans
python ../fmc.py -mechanism=sandiego.mec -thermo=sandiego.therm -name=mec.cpp
echo Compiling sandiego.cpp...
cat mec.cpp ../header/header.start\
            ../header/header.mec   sandiego.mec\
            ../header/header.therm sandiego.therm\
            ../header/header.trans sandiego.trans\
            ../header/header.end > sandiego.cpp
rm -f mec.cpp
