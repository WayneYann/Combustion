ln -fs glarTherm.dat glarSkel.therm
ln -fs tranGlar.dat glarSkel.trans
python ../fmc.py -mechanism=glarSkel.mec -thermo=glarSkel.therm -name=mec.cpp
echo Compiling glarSkel.cpp...
cat mec.cpp ../header/header.start\
            ../header/header.mec   glarSkel.mec\
            ../header/header.therm glarSkel.therm\
            ../header/header.trans glarSkel.trans\
            ../header/header.end > glarSkel.cpp
rm -f mec.cpp
