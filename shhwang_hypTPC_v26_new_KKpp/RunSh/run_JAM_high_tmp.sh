#!/bin/sh 
cd $HOME/private/E27_ichikawa_new1_1/for_TPC_new/shhwang_hypTPC_v26_new_KKpp

source /home/had/yudai/private/E27_ichikawa_new1_1/for_TPC_new/shhwang_hypTPC_v26_new_KKpp/KKpp_parameters_JAM.sh

export Generator=3101
export Input_JAM_File_Name=JAM_rootfile/rot_KKpp_ut0.root

#for((i=0; $i<6059; i++))
#for((i=6000; $i<6059; i++))
for((i=1002; $i<1003; i++))
{
    
    
    export Out_ROOT_File_Name=./rootfiles/KKpp/JAM/ut0/buf2/JAM_0_rot_geant_$i.root
    export Out_GEN_File_Name=./rootfiles/KKpp/JAM/ut0/buf2/JAM_0_rot_gen_$i.root
    x=`expr 1 + 10000 \* $i`
    echo $x
    
    export Nbeam_first=$x
    
    
#    bsub -q s -o "/dev/null" ./bin/Linux-g++/hyptpc1 runJAM.mac
    ./bin/Linux-g++/hyptpc1 runJAM.mac
}

#export Out_ROOT_File_Name=./rootfiles/KKpp/JAM/ut0/buf2/JAM_0_rot_geant_remain.root
#export Out_GEN_File_Name=./rootfiles/KKpp/JAM/ut0/buf2/JAM_0_rot_gen_remain.root
#export Nbeam_first=6059001

#bsub -q s -o "/dev/null" ./bin/Linux-g++/hyptpc1 runJAM_remain.mac
