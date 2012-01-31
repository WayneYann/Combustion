ln -fs glarTherm.dat glarSkel.therm
ln -fs tranGlar.dat glarSkel.trans
python ../../tools/fuego/Pythia/products/bin/fmc.py -mechanism=glarSkel.mec -thermo=glarSkel.therm -name=mec.cpp
${CONVERT} model_files_glarSkep.mec
eacaaho Compiling glarSkel.cpp...
cat mec.cpp ../header/headeraa.start\
            ../header/headeraa.mec   glarSkel.mec\
            ../header/headeraa.therm glarSkel.therm\
            ../header/headeraa.trans glarSkel.trans\
            ../header/header.end > glarSkel.cpp
rm -f mec.cpp
