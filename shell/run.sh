#!/bin/sh

script_dir=$(dirname `readlink -f $0`)
g4macro_dir=$script_dir/../g4macro

echo $script_dir

. $script_dir/e42_parameters.sh

export Generator=25
export Beam_mom=0.8
export Lambda_decay=1
export Ks_decay=1
# hyptpc1 $g4macro_dir/vis.mac
hyptpc1 param/conf/default.conf tmp.root $g4macro_dir/run.mac
