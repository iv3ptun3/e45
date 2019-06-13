#!/bin/csh 
#### Generator
#setenv Generator 74 ####knp
#setenv Generator 64 ####lambda 1405 reso
#setenv Generator 1 ####inc
####setenv Generator 3 ####h-gen
#setenv Generator 10 ####beam study
setenv Generator 30 ####phase

#setenv Generator 25 ####ELASTIC SCATTERING pip
#setenv Generator 26 ####ELASTIC SCATTERING pin
#setenv Generator 21 ####p45 mode1:  pi- p --> pi+ pi- n
#setenv Generator 22 ####p45 mode2:  pi- p --> pi0 pi- p
#setenv Generator 23 ####p45 mode3:  pi+ p --> pi+ pi+ n
#setenv Generator 24 ####p45 mode4:  pi+ p --> pi0 pi+ p



#setenv Generator 99 ####single
#setenv Generator 78 ####dedx all
#setenv Generator 71

setenv Experiment_NUM 42 ####dedx test
setenv With_KURAMA 0 ### 1: w/, 0: w/o

setenv Trigger 0 ####dedx test
#trigger 0: no trigger
#trigger 1: TPC_TOF
#trigger 2: TPC_TOF & FTOF

## 1 : INC generator
## 2 : H generator
## 3 : H generator with forward K+ (15 deg)
## 4 : Test for the sako-san's code. It should be turned off magnetic field.
## 5 : Test for the sako-san's code. It should be turned off magnetic field.
## 6 : weak decay mode by using phase space decay channel. H -> Lppi-
## 7 : weak decay mode by using phase space decay channel. H -> S-p
## 9 : H generator by using PHSG. H -> LL. we should change the mass of H 
## 10 : study on K+ spectrometer, the generator produces the K+ and beam (K-). K+ kinematics uses H generator

## 11 : study on Kstar production, pi-p --> K0(892) L
## 12 : study on Kstar production, pi-p --> K0(892) S
## 13 : study on Kstar production, pi-p --> K0(892) L by using LL gen
## 14 : study on Kstar production, pi-p --> K0(892) S by using LL gen

## 70 : study on E07 --> kaon+ & beam
## 71 : study on E07 --> proton kaon+ pi+
## 72 : study on E07 --> k-p --> k-p
## 74 : study on E07 --> k-p --> k-p, with beam size(3cm x 1cm)
## 75 : study on E07 --> inc K+, with beam size(3cm x 1cm)
## 76 : study on E07 --> K+Xi-, with beam size(3cm x 1cm)
## 77 : study on E07 --> K+Xi-, with beam size(3cm x 1cm), only K+ gen
## 78 : study on E07 --> proton with momentum from 0.8 to 2.0 GeV
## 79 : study on E07 --> K+ with momentum from 0.4 to 1.5 GeV
## 80 : study on E07 --> K+Xi1530-, with beam size(3cm x 1cm)
###truncated mean cut
setenv Truncated_mean_cut 0.3

#setenv TOF_Segment 16
#setenv TOF_Segment 24
setenv TOF_Segment 32

setenv AC_Use 0. ##"1" is on. "0" is off
setenv N_Bar_Use 0. ##"1" is on. "0" is off

setenv AC_Skim 1. ##"1" is on. "0" is off

#### Target size for E42
setenv Target_Size_x 30.
setenv Target_Size_y 10. 
setenv Target_Size_z 15. 
#### Target size for E45
#setenv Target_Size_x 25. ### radius
#setenv Target_Size_y 0. 
#setenv Target_Size_z 25.  ### hight/2
#setenv Target_Size_x 300.
#setenv Target_Size_y 150.
#setenv Target_Size_z 100.

setenv Target_Pos_z -143.
#setenv Target_Pos_z -103.

setenv Out_ROOT_File_Name ./grun1.root
setenv Out_GEN_File_Name ./gen1.root
#It should be added some code.

#### Magnet field for Helmholtz coil 
setenv Helmholtz_fieldmap_on_off 0. # 0 : off, 1 : on
setenv Helmholtz_field -1. # only off case
#setenv Helmholtz_field 1. # only off case
#setenv Helmholtz_field 0.

#### Magnet field for Kurama dipole
setenv Kurama_fieldmap_on_off 0. # 0 : off, 1 : on
#setenv Kurama_field -1.25 # only off case Tm = 1. Tm
#setenv Kurama_gap 250.0 # only off case
#setenv Kurama_pos_z 1426.5 #minimum value : 1326.5
setenv Kurama_pos_z 1500. #minimum value : 1326.5. From Uguard to center of kurama : 820
setenv Kurama_gap 400.0 # only off case
setenv Kurama_move_x 5. # unit mm 
setenv Kurama_field -0.875 # only off case Tm = 0.7 Tm
#setenv Kurama_field 0.875 # only off case Tm = 0.7 Tm
#setenv Kurama_gap 800.0 # only off case
#setenv Helmholtz_field 0.

#### pad_design, 
#1: 9, 16, 10, 17
#2: 9, 13, 10, 20
setenv pad_configure 2. # 1: simple, 2:from final design, 3:right-left divide
setenv pad_length_in 9.
setenv pad_length_out 12.5
setenv pad_num_in 10.
setenv pad_num_out 22.
setenv pad_gap 0.5
setenv pad_width_in 2.5
setenv pad_width_out 2.5

### Target material
setenv Target_Material C ## LH2, Cu, C
## It should be added some code.

setenv H_decaytime 10.
setenv H_mass 2.250 #weak decay
setenv H_width 0.00 #weak decay
#setenv Weak_decay_mass 2.250 #for the generator #9

#### beam information
#setenv beam_position_x 0 #
#setenv beam_position_y 0 #
#setenv beam_position_z 0 #
setenv Beam_mom 1.8 #
#setenv beam_width 0.10 #

#### Spectrometer angle //plus left , minus right
setenv Spectrometer_ang 6.5
#setenv Spectrometer_ang 5.
#setenv Posz_ToF 3242.5 # note that target position is 242.3 mm
####setenv Angle_ToF -14.04
####setenv Angle_ToF 6.5 ####wo_kurama_case

#### GEM discharge on/off(1: design #1, 2: design #2, 3: design #3(now test on plane 3)0: off, 5:final design  )                                                               
setenv GEMDischarge 0.  ###design number
setenv DeadArea 50.  ### width of electrode
setenv GEMFixDead 1. ### on/off (1/0) when I fix the dead region
setenv GEMDeadPlane 3.  ###
setenv GEMDeadPlaneDivision 3.


## useage 
##char* OutFileName = getenv("Out_ROOT_File_Name");
##String to double
##G4String get_e = getenv("Gamma_Energy");
##Gamma_Energy = atof(get_e.c_str());

