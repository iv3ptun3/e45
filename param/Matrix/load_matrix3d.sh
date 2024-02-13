#!/bin/sh

set -e
ip_addr=192.168.11.13
work_dir=$(dirname $(readlink -f $0))
bin_dir=$work_dir/bin
param_dir=$work_dir/param/mtx3d
last_log=$work_dir/last_matrix3d.log
#
param_file=$param_dir/mtx3d_e42_Kaon_20210608.txt
# param_file=$param_dir/mtx3d_e42_data_Kbeam.txt
#param_file=$param_dir/mtx3d_e42_GEANT4.txt
#
# param_file=$param_dir/pattern.txt
# param_file=$param_dir/all_1.txt
# param_file=$param_dir/mtx3d_e03_GEANT4_normal.txt
# param_file=$param_dir/mtx3d_e03_GEANT4_tight.txt
# param_file=$param_dir/mtx3d_e03_GEANT4_wide.txt
#param_file=$param_dir/mtx3d_e03_GEANT4_normal_p3.txt
#param_file=$param_dir/mtx3d_e03_exdata_normal_plus3.txt
#
command="$bin_dir/load_matrix3d $ip_addr $param_file"
date | tee $last_log
echo $command | tee -a $last_log
exec $command | tee -a $last_log
