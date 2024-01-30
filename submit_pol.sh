#for i in {5111,5113,5392,5393}
#for i in {5723..5732}
#for i in {5723..5756}
#for i in {5177,5723}
for i in {1,0.75,0.5,0.25,0,-0.25,-0.5,-0.75,-1}
	do
		echo ${i}
		#name="bsub -q l ./bin/Linux-g++/hyptpc1 param/conf/kpxi_data_P_${i}.conf rootfiles_ana/KpXi_Pol_${i}.root g4macro/run_1M.mac"
		name="bsub -q l ./bin/Linux-g++/hyptpc1 param/conf/kpxi_data_P_${i}.conf rootfiles_ana/KpXi_Pol_${i}.root g4macro/run_100k.mac"
		$name
#		$name1
#		$name2
#		$name3
#		$name4
	done

