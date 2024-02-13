#!/bin/sh

set -e
ip_addr=192.168.11.13
work_dir=$(dirname $(readlink -f $0))
bin_dir=$work_dir/bin
param_dir=$work_dir/param/dwg
last_log=$work_dir/last_dwg.log
#
param_file=$param_dir/dwg_register.txt
#
command="$bin_dir/load_register $ip_addr $param_file"
date | tee $last_log
echo $command | tee -a $last_log
exec $command | tee -a $last_log
