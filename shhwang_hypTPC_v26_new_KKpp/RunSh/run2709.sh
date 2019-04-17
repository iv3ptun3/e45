source Param_E27/e27_parameters_2703.sh

export Generator=2709
for((i=0; $i<100; i++))
{
export Out_ROOT_File_Name=./rootfiles/Kbeam/test_2709_$i.root
export Out_GEN_File_Name=./rootfiles/Kbeam/test_2709_gen_$i.root
bsub -q s -o "/dev/null" ./bin/Linux-g++/hyptpc1 run2.mac
}

