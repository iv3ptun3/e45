source Param_E27/e27_parameters_2703.sh

export Generator=2703
for((i=0; $i<20; i++))
{
export Out_ROOT_File_Name=./rootfiles/Test_2/Kpp_LP/test_2703_$i.root
export Out_GEN_File_Name=./rootfiles/Test_2/Kpp_LP/gen/test_2703_gen_$i.root
bsub -q s -o "/dev/null" ./bin/Linux-g++/hyptpc1 run2.mac
}

export Generator=2704
for((i=0; $i<20; i++))
{
export Out_ROOT_File_Name=./rootfiles/Test_2/Kpp_SzP/test_2704_$i.root
export Out_GEN_File_Name=./rootfiles/Test_2/Kpp_SzP/gen/test_2704_gen_$i.root
bsub -q s -o "/dev/null" ./bin/Linux-g++/hyptpc1 run2.mac
}

export Generator=2705
for((i=0; $i<20; i++))
{
export Out_ROOT_File_Name=./rootfiles/Test_2/Kpp_LPizP/test_2705_$i.root
export Out_GEN_File_Name=./rootfiles/Test_2/Kpp_LPizP/gen/test_2705_gen_$i.root
bsub -q s -o "/dev/null" ./bin/Linux-g++/hyptpc1 run2.mac
}

export Generator=2706
for((i=0; $i<20; i++))
{
export Out_ROOT_File_Name=./rootfiles/Test_2/Kpp_SzPizP/test_2706_$i.root
export Out_GEN_File_Name=./rootfiles/Test_2/Kpp_SzPizP/gen/test_2706_gen_$i.root
bsub -q s -o "/dev/null" ./bin/Linux-g++/hyptpc1 run2.mac
}

export Generator=2707
for((i=0; $i<20; i++))
{
export Out_ROOT_File_Name=./rootfiles/Test_2/Kpp_SpPimP/test_2707_$i.root
export Out_GEN_File_Name=./rootfiles/Test_2/Kpp_SpPimP/gen/test_2707_gen_$i.root
bsub -q s -o "/dev/null" ./bin/Linux-g++/hyptpc1 run2.mac
}

