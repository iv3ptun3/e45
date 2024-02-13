#!/bin/sh

set -e
ip_addr=192.168.11.13
work_dir=$(dirname $(readlink -f $0))
bin_dir=$work_dir/bin
param_dir=$work_dir/param/nimo
last_log=$work_dir/last_nimo.log
#
param_file=$param_dir/nimo_register.txt
#
command="$bin_dir/load_nimo $ip_addr $param_file"
date | tee $last_log
echo $command | tee -a $last_log
exec $command | tee -a $last_log
