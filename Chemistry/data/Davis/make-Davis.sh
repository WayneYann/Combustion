CHEMTOOLSDIR=/scratch2/marc/src/CCSE/Combustion/Chemistry/tools
#CHEMTOOLSDIR=/Users/rgrout/Code/rev-controlled/LBNL/Combustion/Chemistry/tools

CHEMINP=mech.inp
THERMINP=therm.dat
TRANINP=tran.dat
FINALFILE=davis.c

CONVERT=${CHEMTOOLSDIR}/convert/convert.exe
FMC=${CHEMTOOLSDIR}/fuego/Pythia/products/bin/fmc.py

CHEMLK=chem.asc
LOG=chem.log
TRANC=tran.c
CHEMC=chem.c
TRANLOG=tran.log
HEADERDIR=${CHEMTOOLSDIR}/../data/header

python ${FMC} -mechanism=${CHEMINP} -thermo=${THERMINP} -name=${CHEMC}
echo Compiling ${FINALFILE}...
echo " &files"  > model_files.dat
echo "   CHEMKIN_input = \"$CHEMINP\"" >> model_files.dat
echo "   THERMO_input = \"$THERMINP\"" >> model_files.dat
echo "   TRANLIB_input = \"$TRANINP\"" >> model_files.dat
echo "   CHEMKIN_linking_file = \"$CHEMLK\"" >> model_files.dat
echo "   TRANLIB_c_file = \"$TRANC\"" >> model_files.dat
echo "   log_file = \"$TRANLOG\"" >> model_files.dat
echo " /" >> model_files.dat
${CONVERT} model_files.dat 2>&1 >> $TRANLOG
cat ${CHEMC} ${TRANC} \
          ${HEADERDIR}/header.start\
          ${HEADERDIR}/header.mec   ${CHEMINP}\
          ${HEADERDIR}/header.therm ${THERMINP}\
          ${HEADERDIR}/header.trans ${TRANINP}\
          ${HEADERDIR}/header.end > ${FINALFILE}
rm -f ${CHEMC} ${CHEMLK} ${LOG} ${TRANC} ${TRANLOG} model_files.dat
