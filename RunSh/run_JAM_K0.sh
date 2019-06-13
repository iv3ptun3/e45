#!/bin/sh 
cd $HOME/private/E27_ichikawa_new1_1/for_TPC_new/git/k18geant4/shhwang_hypTPC_v26_new_KKpp

source /home/had/yudai/private/E27_ichikawa_new1_1/for_TPC_new/git/k18geant4/shhwang_hypTPC_v26_new_KKpp/KKpp_parameters_JAM.sh

export Generator=3101
export Input_JAM_File_Name=JAM_rootfile/rot_KKpp_ut0_K0.root

for((i=0; $i<848; i++))
{
    export Out_ROOT_File_Name=./rootfiles/KKpp/JAM/ut0/buf_K0/JAM_0_rot_geant_$i.root
 
    x=`expr 1 + 100000 \* $i`
    echo $x
    
    export Nbeam_first=$x
    
    
    bsub -q l -o "/dev/null" ./bin/Linux-g++/hyptpc1 runJAM_high.mac
}
