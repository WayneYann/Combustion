#!/bin/sh

cd /home/wepazner/Documents/Code/Combustion/LMC/bin/1D_MISDC_WILL

touch probin

rm probin
cp probin256 probin
./lmc.exe

rm probin
cp probin512 probin
./lmc.exe

rm probin
cp probin1024 probin
./lmc.exe

#rm probin
#cp probin2048 probin
#./lmc.exe
