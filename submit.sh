#for i in {5723..5742}
for i in {5111,5113,5392,5393}
	do
		#name="bsub -q s ./bin/Linux-g++/hyptpc1 param/k18conf/run0$i.conf run0$i""_Geant4.root g4macro/run1.mac"
		name="bsub -q s ./bin/Linux-g++/hyptpc1 param/k18conf/run0${i}.conf run0${i}_Geant4.root g4macro/run1.mac"
		$name
	done

