#!/bin/sh

set -e
ip_addr=192.168.11.13
work_dir=$(dirname $(readlink -f $0))
bin_dir=$work_dir/bin
param_dir=$work_dir/param/mtx2d2
last_log=$work_dir/last_matrix2d2.log
#
param_file=$param_dir/all_1.txt
#
# param_file=$param_dir/pattern.txt
# param_file=$param_dir/mtx2d2_e03_GEANT4_normal.txt
# param_file=$param_dir/mtx2d2_e03_GEANT4_tight.txt
# param_file=$param_dir/mtx2d2_e03_GEANT4_wide.txt
#param_file=$param_dir/mtx2d2_e03_GEANT4_wide_p3m6.txt
#param_file=$param_dir/mtx2d2_e03_exdata_normal.txt
#param_file=$param_dir/mtx2d2_e03_exdata_normal_minus1.txt
#param_file=$param_dir/mtx2d2_e03_exdata_normal_minus2.txt
#
#
command="$bin_dir/load_matrix2d $ip_addr $param_file 2"
date | tee $last_log
echo $command | tee -a $last_log
exec $command | tee -a $last_log
