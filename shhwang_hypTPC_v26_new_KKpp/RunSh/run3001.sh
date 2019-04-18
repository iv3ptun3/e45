source /home/had/yudai/private/E27_ichikawa_new1_1/for_TPC_new/git/k18geant4/shhwang_hypTPC_v26_new_KKpp/KKpp_parameters.sh

export Generator=3001

for((i=0; $i<100; i++))
{
    export Out_ROOT_File_Name=./rootfiles/KKpp/KKpp_LL1/back/test_woEM_3001_$i.root

    bsub -q s -o "/dev/null" ./bin/Linux-g++/hyptpc1 run2.mac
}

