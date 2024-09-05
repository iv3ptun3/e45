#for i in {5111,5113,5392,5393}
#for i in {5723..5732}
#for i in {5723..5756}
#for i in {5177,5723}
for i in {0..499}
do
  #name="bsub -q s ./bin/Linux-g++/hyptpc1 param/conf/KpUniformCH2.conf rootfiles_ana/KpHists/KpUniformHistCH2_0${i}_Geant4.root g4macro/run_1M.mac"
  name="bsub -q s ./bin/Linux-g++/hyptpc1 param/conf/KpUniformCarbon.conf rootfiles_ana/KpHists/KpUniformHistCarbon_0${i}_Geant4.root g4macro/run_1M.mac"
  $name
done

