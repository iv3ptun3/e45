source Param_E27/e27_parameters_2703.sh

export Generator=2705
export Out_ROOT_File_Name=./rootfiles/Test/Kpp_LPizP/test_2705_tmp.root
export Out_GEN_File_Name=./rootfiles/Test/Kpp_LPizP/test_2705_gen_tmp.root
bsub -q l -o "/dev/null" ./bin/Linux-g++/hyptpc1 run_tmp.mac

export Generator=2706
export Out_ROOT_File_Name=./rootfiles/Test/Kpp_SzPizP/test_2706_tmp.root
export Out_GEN_File_Name=./rootfiles/Test/Kpp_SzPizP/test_2706_gen_tmp.root
bsub -q l -o "/dev/null" ./bin/Linux-g++/hyptpc1 run_tmp.mac

export Generator=2707
export Out_ROOT_File_Name=./rootfiles/Test/Kpp_SpPimP/test_2707_tmp.root
export Out_GEN_File_Name=./rootfiles/Test/Kpp_SpPimP/test_2707_gen_tmp.root
bsub -q l -o "/dev/null" ./bin/Linux-g++/hyptpc1 run_tmp.mac


