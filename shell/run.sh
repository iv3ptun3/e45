#!/bin/sh

conf=param/conf/default.conf
root=tmp.root
macro=g4macro/run.mac
log=log/tmp.log
mkdir -p log

case $HOSTNAME in
    *.cc.kek.jp)
	bsub -q s -o $log hyptpc1 $conf $root $macro
	;;
    *)
	echo "hyptpc1 $conf $root $macro >$log 2>&1" | batch
	;;
esac
