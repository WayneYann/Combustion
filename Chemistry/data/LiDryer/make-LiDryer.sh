python ../fmc.py -mechanism=LiDryer.mec -thermo=LiDryer.therm -name=mec.cpp
echo Compiling LiDryer.cpp...
cat mec.cpp ../header/header.start\
            ../header/header.mec   LiDryer.mec\
            ../header/header.therm LiDryer.therm\
            ../header/header.trans LiDryer.trans\
            ../header/header.end > LiDryer.cpp
rm -f mec.cpp
