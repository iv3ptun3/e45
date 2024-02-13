#!/bin/sh

disp_matrix()
{
  mtx_log=/data3/E03SubData/trigger_2020dec/mtx.log
  date | tee $mtx_log
  last_log=(
    # last_dwg.log
    last_matrix2d1.log
    last_matrix2d2.log
    last_matrix3d.log
    # last_nimo.log
  )
  for f in ${last_log[@]}
  do
    param=`head -2 $f | tail -1 | awk '{ print $3 }'`
    dir=$(basename $(dirname $param))
    echo -e "$dir:\t$(basename $param)" | tee -a $mtx_log
  done
}
export -f disp_matrix

watch -t disp_matrix
