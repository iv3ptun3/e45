#!/bin/sh 
cd $HOME/private/E27_ichikawa_new1_1/for_TPC_new/git/k18geant4/shhwang_hypTPC_v26_new_KKpp
source /home/had/yudai/private/E27_ichikawa_new1_1/for_TPC_new/git/k18geant4/shhwang_hypTPC_v26_new_KKpp/KKpp_parameters_beamthrough.sh

export Generator=3102

for((i=0; $i<200; i++))
{
    export Out_ROOT_File_Name=./rootfiles/KKpp/BT/BT_geant_3102_$i.root

    bsub -q s -o "/dev/null" ./bin/Linux-g++/hyptpc1 run2.mac
}

