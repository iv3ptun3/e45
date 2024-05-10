#!/bin/bash 
set -o noclobber
#for i in {5187,5188,5191,5200,5207,5212,5219,5220,5221,5808,5814,5808,5814,5818,5821,5824,5828,5835,5837,5838,5840,5842,5844,5846,5847,5851,5851,5855,5856,5858,5860,5862,5864,5866}//Beamthrough
#for i in {5721,5764}
#for i in {6380..6383}
for i in {1,0.75,0.5,0.25,0,-0.25,-0.5,-0.75,-1}
#for i in {5177..5178}
do
	file="./param/conf/kpxi_data_P_${i}.conf"
	echo "Generating $file"
	echo "#">>$file
	echo "#  Generated Conf file data based Xi generator, PolXi = $i ">>$file
	echo "#">>$file
	echo "">>$file
	echo "BEAM		param/REALBEAM_WS/DstE42Summary.root">>$file
	echo "DCGEO		param/DCGEO/DCGeomParam_E42_230919">>$file 
	echo "DSIZE		param/DSIZE/DetSize_E42_230919">>$file 
	echo "TPCParam  param/TPCPRM/TPCParam_20230825">>$file
	echo "">>$file
	echo "KURAMAFLDMAP		fieldmap/KuramaFieldMap_20210623">>$file 
	echo "SHSFLDMAP		fieldmap/ShsFieldMap_20210623">>$file
	
	echo "KURAMAFLDCALC	0.749779">>$file 
	echo "KURAMAFLDNMR	0.760830712">>$file 
	echo "">>$file
	echo "HSFLDCALIB	1.">>$file
	echo "SHSFLDCALC	0.90607">>$file
	echo "SHSFLDNMR	0.8738">>$file
	echo "ShsFieldOffset	0">>$file
	echo "">>$file
	echo "##### General Parameter">>$file
	echo "">>$file
	echo "EVDISP	0">>$file
	echo "Experiment	42">>$file
	echo "Generator	181321 #KpXi RealData">>$file
	echo "TargetMaterial	CH2">>$file
	echo "BeamMom	1.8">>$file
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
	echo "Decay	1">>$file
	echo "PolarizedDecay	1">>$file
	echo "XiPolarization	${i}">>$file
	echo "EM	1">>$file
	echo "Hadron	1">>$file
	echo "KillStepInIron	1">>$file
	echo "HdibaryonMass	2.25 ">>$file
	echo "HdibaryonWidth	0.0 #GeV">>$file
	echo "HdibaryonLifeTime	10.0 #ns">>$file
	echo "##### HypTPC Parameter">>$file
	echo "##### Optional">>$file
	echo "TruncatedMeanCut 0.3">>$file
	echo "UseAC 0">>$file
	echo "UseNBar 0">>$file
	echo "NBeamFirst 0">>$file

done
