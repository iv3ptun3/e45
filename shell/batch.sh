#!/bin/sh

work_dir=$(dirname $(readlink -f $0))/..

#conf=$work_dir/param/conf/default.conf
#conf=$work_dir/param/conf/kpxi.conf
#conf=$work_dir/param/conf/e42_beam.conf
#conf=$work_dir/param/conf/e42_jam_kp.conf
conf=$work_dir/param/conf/e42_jam_proton.conf
g4macro_dir=$work_dir/g4macro
root_dir=$work_dir/rootfiles/root/all
log_dir=$work_dir/rootfiles/log

n_run=100

cd $work_dir
mkdir -p $root_dir
mkdir -p $log_dir

g4macro=$g4macro_dir/run.mac

#for i in $(seq 0 0); do
for i in $(seq 1 $n_run); do
    head=$(basename ${conf%.conf})
    log=$log_dir/${head}_`printf %05d $i`.log
    root=$root_dir/${head}_`printf %05d $i`.root
    rm -f $log
    bsub -q l -o $log hyptpc1 $conf $root $g4macro
done
