#source ../env.csh
source e45_parameters.csh

setenv Generator 25
setenv Beam_mom 0.8
setenv Out_ROOT_File_Name ./grun_e45_elastic_08_pip.root
setenv Out_GEN_File_Name ./gen_e45_elastic_08_pip.root
./bin/Darwin-clang/hyptpc1 run.mac &

setenv Beam_mom 1.0
setenv Out_ROOT_File_Name ./grun_e45_elastic_10_pip.root
setenv Out_GEN_File_Name ./gen_e45_elastic_10_pip.root
./bin/Darwin-clang/hyptpc1 run.mac &

setenv Beam_mom 1.2
setenv Out_ROOT_File_Name ./grun_e45_elastic_12_pip.root
setenv Out_GEN_File_Name ./gen_e45_elastic_12_pip.root
./bin/Darwin-clang/hyptpc1 run.mac &


setenv Beam_mom 1.4
setenv Out_ROOT_File_Name ./grun_e45_elastic_14_pip.root
setenv Out_GEN_File_Name ./gen_e45_elastic_14_pip.root
./bin/Darwin-clang/hyptpc1 run.mac &

setenv Beam_mom 1.6
setenv Out_ROOT_File_Name ./grun_e45_elastic_16_pip.root
setenv Out_GEN_File_Name ./gen_e45_elastic_16_pip.root
./bin/Darwin-clang/hyptpc1 run.mac &


setenv Beam_mom 1.8
setenv Out_ROOT_File_Name ./grun_e45_elastic_18_pip.root
setenv Out_GEN_File_Name ./gen_e45_elastic_18_pip.root
./bin/Darwin-clang/hyptpc1 run.mac &

setenv Beam_mom 2.0
setenv Out_ROOT_File_Name ./grun_e45_elastic_20_pip.root
setenv Out_GEN_File_Name ./gen_e45_elastic_20_pip.root
./bin/Darwin-clang/hyptpc1 run.mac &



