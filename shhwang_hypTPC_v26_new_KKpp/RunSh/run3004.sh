source /home/had/yudai/private/E27_ichikawa_new1_1/for_TPC_new/shhwang_hypTPC_v26_new_KKpp/KKpp_parameters.sh

export Generator=3004

#for((i=0; $i<100; i++))
for((i=0; $i<10; i++))
{
    export Out_ROOT_File_Name=./rootfiles/KKpp/KKpp_LSpPim/back/test_woEM_3004_$i.root
    export Out_GEN_File_Name=./rootfiles/KKpp/KKpp_LSpPim/back/test_woEM_3004_gen_$i.root
#    bsub -q s -o "/dev/null" ./bin/Linux-g++/hyptpc1 run2.mac
    bsub -q s -o "/dev/null" ./bin/Linux-g++/hyptpc1 run1.mac
}

