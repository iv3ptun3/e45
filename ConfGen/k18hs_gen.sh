#!/bin/bash 
set -o noclobber
run_start=5082
run_end=6372
for(( run=$run_start; run<=$run_end; run++));do
	
	beamname="param/KmBeam/run0${run}_K18HSTracking.root"
	if [[ -f "$beamname" ]]; then
		file="./param/k18hs/k18hs_run0${run}.conf"
		echo "Generating $file"
		echo "#">>$file
		echo "#  Generated Conf file for data based Beam simulation.">>$file
		echo "#">>$file
		echo "">>$file
		echo "BEAM		param/KmBeam/run0${run}_K18HSTracking.root">>$file
		echo "DCGEO		param/DCGEO/DCGeomParam_E42_230919">>$file 
		echo "DSIZE		param/DSIZE/DetSize_E42_230919">>$file 
		echo "TPCParam  param/TPCPRM/TPCParam_20230825">>$file
		echo "">>$file
		echo "KURAMAFLDMAP		fieldmap/KuramaFieldMap_20210526">>$file 
		echo "SHSFLDMAP		param/ShsFieldMap_20210526_Adjusted">>$file

		echo "KURAMAFLDCALC	0.749779">>$file 
		echo "KURAMAFLDNMR	0.760830712">>$file 
		echo "">>$file
		echo "HSFLDCALIB	0.992623">>$file
		echo "SHSFLDCALC	1.">>$file
		echo "SHSFLDNMR	1.0">>$file
		echo "ShsFieldOffset	0">>$file
		echo "">>$file
		echo "##### General Parameter">>$file
		echo "">>$file
		echo "EVDISP	0">>$file
		echo "DiscardData	1">>$file
		echo "Experiment	42">>$file
		echo "Generator	-493 #Km beam">>$file
		echo "TargetMaterial	Diamond">>$file
		echo "BeamMom	1.8">>$file
		echo "##### KURAMA">>$file
		echo "ConstructKurama 0">>$file
		echo "KuramaFieldMap 1">>$file
		echo "KuramaField -0.7">>$file

		echo "##### SHS">>$file
		echo "ShsFieldMap 1">>$file
		echo "ShsField -1.0">>$file

		echo "##### PHYSICS">>$file
		echo "Physics	USER">>$file
		echo "Decay	0">>$file
		echo "PolarizedDecay	0">>$file
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
	fi
done
