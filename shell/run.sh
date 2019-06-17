#!/bin/sh

script_dir=$(dirname `readlink -f $0`)
bin_dir=$script_dir/../bin/Linux-g++
g4macro_dir=$script_dir/../g4macro

echo $script_dir

. $script_dir/e45_parameters.sh

export Generator=25
export Beam_mom=0.8
export Lambda_decay=1
export Ks_decay=1
export Out_ROOT_File_Name=./grun_e45_elastic_08_pip.root
export Out_GEN_File_Name=./gen_e45_elastic_08_pip.root
#$bin_dir/hyptpc1 $g4macro_dir/vis.mac
$bin_dir/hyptpc1 $g4macro_dir/run.mac

# export Beam_mom 1.0
# export Out_ROOT_File_Name ./grun_e45_elastic_10_pip.root
# export Out_GEN_File_Name ./gen_e45_elastic_10_pip.root
# $hyptpc run.mac &

# export Beam_mom 1.2
# export Out_ROOT_File_Name ./grun_e45_elastic_12_pip.root
# export Out_GEN_File_Name ./gen_e45_elastic_12_pip.root
# $hyptpc run.mac &


# export Beam_mom 1.4
# export Out_ROOT_File_Name ./grun_e45_elastic_14_pip.root
# export Out_GEN_File_Name ./gen_e45_elastic_14_pip.root
# $hyptpc run.mac &

# export Beam_mom 1.6
# export Out_ROOT_File_Name ./grun_e45_elastic_16_pip.root
# export Out_GEN_File_Name ./gen_e45_elastic_16_pip.root
# $hyptpc run.mac &


# export Beam_mom 1.8
# export Out_ROOT_File_Name ./grun_e45_elastic_18_pip.root
# export Out_GEN_File_Name ./gen_e45_elastic_18_pip.root
# $hyptpc run.mac &

# export Beam_mom 2.0
# export Out_ROOT_File_Name ./grun_e45_elastic_20_pip.root
# export Out_GEN_File_Name ./gen_e45_elastic_20_pip.root
# $hyptpc run.mac &
