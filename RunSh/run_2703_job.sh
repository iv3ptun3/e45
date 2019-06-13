source Param_E27/e27_parameters_2703.sh

export Generator=2703
for((i=0; $i<10; i++))
{
export Out_ROOT_File_Name=./rootfiles/Test_2/Kpp_LP/test_2703_$i.root
export Out_GEN_File_Name=./rootfiles/Test_2/Kpp_LP/test_2703_gen_$i.root
bsub -q s -o "/dev/null" ./bin/Linux-g++/hyptpc1 run.mac
}

