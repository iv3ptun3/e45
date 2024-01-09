#!/bin/bash 
set -o noclobber
#for i in {5187,5188,5191,5200,5207,5212,5219,5220,5221,5808,5814,5808,5814,5818,5821,5824,5828,5835,5837,5838,5840,5842,5844,5846,5847,5851,5851,5855,5856,5858,5860,5862,5864,5866}//Beamthrough
#for i in {5721,5764}
#for i in {6380..6383}
#for i in {5723..5756}
#for i in {5177..5178}
for i in {5113,5176,5177,5188,5219,5220,5221,5723,5724,5725,5726,5727,5728,5729,5730,5731,5732,5733,5734,5735,5736,5737,5738,5739,5740,5741,5742,5744,5745,5746,5747,5748,5749,5808,5831,5837,5838,5855}
do
#	file="./param/k18conf/run0${i}.conf"
	file="./param/k18conf_NoExtrap/run0${i}.conf"
#	file="./param/k18conf_NoScale/run0${i}.conf"
#	file="./param/k18conf_NoScaleNoExtrap/run0${i}.conf"
	echo "Generating $file"
	pK18=`grep pK18 /group/had/sks/E42/JPARC2021May/K18MagnetParam/K18MagParam_0$i | cut -c17-28`
	Kurama_F=`grep KURAMA_F /group/had/sks/E42/JPARC2021May/K18MagnetParam/K18MagParam_0$i | cut -c18-28`
	SHS_F=`grep SHS_F /group/had/sks/E42/JPARC2021May/K18MagnetParam/K18MagParam_0$i | cut -c18-28`
	echo "#">>$file
	echo "#  Generated Conf file for Run By Run PHC: run 0$i ">>$file
	echo "#">>$file
	echo "">>$file
	echo "BEAM		param/REALBEAM/run0${i}_K18HSTracking.root">>$file
	echo "DCGEO		param/DCGEO/DCGeomParam_E42_230919">>$file 
	echo "DSIZE		param/DSIZE/DetSize_E42_230919">>$file 
	echo "">>$file
	echo "KURAMAFLDCALC	0.749779">>$file 
	echo "KURAMAFLDNMR		$Kurama_F">>$file 
	echo "KURAMAFLDMAP		fieldmap/KuramaFieldMap_20210623">>$file 
	echo "">>$file
	echo "HSFLDCALIB	0.94838">>$file
#	echo "HSFLDCALIB	1.">>$file
	echo "SHSFLDCALC	0.90607">>$file
	echo "SHSFLDNMR	$SHS_F">>$file
#	echo "SHSFLDMAP	param/ShsFieldMap_20210526_Extrapolated">>$file
	echo "SHSFLDMAP	fieldmap/ShsFieldMap_20210526">>$file
	echo "ShsFieldOffset	0">>$file
	echo "##### General Parameter">>$file
	echo "EVDISP	0">>$file
	echo "Experiment	42">>$file
	echo "Generator	-135">>$file
	echo "TargetMaterial	Empty">>$file
	echo "BeamMom	1.8">>$file
	echo "KillStepInIron	0">>$file
	echo "##### KURAMA">>$file
	echo "ConstructKurama 1">>$file
	echo "KuramaFieldMap 1">>$file
	echo "KuramaField -0.7">>$file
	echo "SpectrometerAngle -0.0 # not supported">>$file

	echo "##### SHS">>$file
	echo "ShsFieldMap 1">>$file
	echo "ShsField -1.0">>$file

	echo "##### PHYSICS">>$file
	echo "Physics	USER">>$file
	echo "Decay	0">>$file
	echo "PolarizedDecay	0">>$file
	echo "EM	1">>$file
	echo "Hadron	0">>$file
	echo "HdibaryonMass	2.25 ">>$file
	echo "HdibaryonWidth	0.0 #GeV">>$file
	echo "HdibaryonLifeTime	10.0 #ns">>$file
	echo "##### HypTPC Parameter">>$file
	echo "TPCParam		param/TPCPRM/TPCParam_20230825">>$file 
	echo "##### Optional">>$file
	echo "TruncatedMeanCut 0.3">>$file
	echo "UseAC 0">>$file
	echo "UseNBar 0">>$file
	echo "NBeamFirst 0">>$file

done
