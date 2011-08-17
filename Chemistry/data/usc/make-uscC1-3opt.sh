python ../fmc.py -mechanism=uscC1-3opt.mec -thermo=uscC1-3opt.therm-hack -name=mec.cpp
echo Compiling uscC1-3opt.cpp...
cat mec.cpp ../header/header.start\
            ../header/header.mec   uscC1-3opt.mec\
            ../header/header.therm uscC1-3opt.therm-hack\
            ../header/header.end > uscC1-3opt.cpp
rm -f mec.cpp
