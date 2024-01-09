#for i in {5111,5113,5392,5393}
#for i in {5723..5732}
#for i in {5723..5756}
#for i in {5177,5723}
for i in {5113,5176,5177,5188,5219,5220,5221,5723,5724,5725,5726,5727,5728,5729,5730,5731,5732,5733,5734,5735,5736,5737,5738,5739,5740,5741,5742,5744,5745,5746,5747,5748,5749,5808,5831,5837,5838,5855}
	do
#		name1="bsub -q s ./bin/Linux-g++/hyptpc1 param/k18conf/run0${i}.conf rootfiles/run0${i}_Geant4ScaleExtrap.root g4macro/run1.mac"
		name2="bsub -q s ./bin/Linux-g++/hyptpc1 param/k18conf_NoExtrap/run0${i}.conf rootfiles/run0${i}_Geant4Scale.root g4macro/run1.mac"
#		name3="bsub -q s ./bin/Linux-g++/hyptpc1 param/k18conf_NoScale/run0${i}.conf rootfiles/run0${i}_Geant4NoScaleExtrap.root g4macro/run1.mac"
		name4="bsub -q s ./bin/Linux-g++/hyptpc1 param/k18conf_NoScaleNoExtrap/run0${i}.conf rootfiles/run0${i}_Geant4NoScale.root g4macro/run1.mac"
#		$name1
		$name2
#		$name3
		$name4
	done

