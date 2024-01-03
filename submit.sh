#for i in {5111,5113,5392,5393}
for i in {5723..5732}
#for i in {5723..5756}
#for i in {5177,5723}
	do
		name1="bsub -q s ./bin/Linux-g++/hyptpc1 param/k18conf/run0${i}.conf rootfiles/run0${i}_Geant4ScaleExtrap.root g4macro/run1.mac"
		name2="bsub -q s ./bin/Linux-g++/hyptpc1 param/k18conf_NoExtrap/run0${i}.conf rootfiles/run0${i}_Geant4Scale.root g4macro/run1.mac"
		name3="bsub -q s ./bin/Linux-g++/hyptpc1 param/k18conf_NoScale/run0${i}.conf rootfiles/run0${i}_Geant4NoScaleExtrap.root g4macro/run1.mac"
		name4="bsub -q s ./bin/Linux-g++/hyptpc1 param/k18conf_NoScaleNoExtrap/run0${i}.conf rootfiles/run0${i}_Geant4NoScale.root g4macro/run1.mac"
		$name1
		$name2
		$name3
		$name4
	done

