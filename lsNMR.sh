#!/bin/bash 
set -o noclobber
#for i in {5187,5188,5191,5200,5207,5212,5219,5220,5221,5808,5814,5808,5814,5818,5821,5824,5828,5835,5837,5838,5840,5842,5844,5846,5847,5851,5851,5855,5856,5858,5860,5862,5864,5866}//Beamthrough
#for i in {5721,5764}
#for i in {6380..6383}
for i in {5123..6196}
#for i in {5177..5178}
do
	Kurama_F=`grep KURAMA_F /group/had/sks/E42/JPARC2021May/K18MagnetParam/K18MagParam_0$i | cut -c18-28`
	SHS_F=`grep SHS_F /group/had/sks/E42/JPARC2021May/K18MagnetParam/K18MagParam_0$i | cut -c18-28`
	echo "run0" "$i" "SHSNMR = $SHS_F , KURAMANMR = $Kurama_F"

done
