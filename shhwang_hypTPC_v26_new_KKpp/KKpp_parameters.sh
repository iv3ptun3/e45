#!/bin/sh 
#### Generator
#export Generator 74 ####knp
#export Generator 64 ####lambda 1405 reso
#export Generator 1 ####inc
####export Generator 3 ####h-gen
#export Generator=10 ####beam study
#export Generator=30 ####phase

#export Generator 25 ####ELASTIC SCATTERING pip
#export Generator 26 ####ELASTIC SCATTERING pin
#export Generator 21 ####p45 mode1:  pi- p --> pi+ pi- n
#export Generator 22 ####p45 mode2:  pi- p --> pi0 pi- p
#export Generator 23 ####p45 mode3:  pi+ p --> pi+ pi+ n
#export Generator 24 ####p45 mode4:  pi+ p --> pi0 pi+ p
#export Generator 99 ####single
#export Generator 78 ####dedx all
#export Generator 71

##=====E27 Generator========
#export Generator=2701 #### pi+ beam through
#export Generator=2702 #### K+ scat gun
#export Generator=2703 #### K-pp -> Lp scat gun
#export Generator=2708 #### K-pp -> Lp scat gun
#export Generator=2709 #### K+ gun 2

##=====KKpp Generator========
#export Generator=3001 #### KKpp -> LL -> p pi p pi
export Generator=3002 #### KKpp -> LL 
#export Generator=3003 #### KKpp -> LSmPip
#export Generator=3004 #### KKpp -> LSpPim

export Experiment_NUM=45 ####dedx test
#export With_KURAMA=1 ### 1: w/, 0: w/o
export With_KURAMA=0 ### 1: w/, 0: w/o

export Trigger=0 ####dedx test
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
export Truncated_mean_cut=0.3

#export TOF_Segment 16
#export TOF_Segment 24
export TOF_Segment=32

export AC_Use=0. ##"1" is on. "0" is off
export N_Bar_Use=0. ##"1" is on. "0" is off

export AC_Skim=0. ##"1" is on. "0" is off

#### Target size for E42
#export Target_Size_x=30.
#export Target_Size_y=10. 
#export Target_Size_z=15. 
#### Target size for E45 and E27
export Target_Size_x=25. ### radius
export Target_Size_y=0. 
export Target_Size_z=25.  ### hight/2
#export Target_Size_x 300.
#export Target_Size_y 150.
#export Target_Size_z 100.

export Target_Pos_z=-143.
#export Target_Pos_z -103.

export Out_ROOT_File_Name=./grun1.root
export Out_GEN_File_Name=./gen1.root
#It should be added some code.

#### Magnet field for Helmholtz coil 
export Helmholtz_fieldmap_on_off=0 # 0 : off, 1 : on
#export Helmholtz_field=-1. # only off case
#export Helmholtz_field=1. # only off case
#export Helmholtz_field=0.
export Helmholtz_field=-1.
export sigma0x=0.204 #mm HIMAC analysis
export Dt=0.18 #mm/sqrt(cm) at 1T

#### Magnet field for Kurama dipole
##A/m = 1.25 uT
export Kurama_fieldmap_on_off=0 # 0 : off, 1 : on
export Kurama_field=-1.25 # only off case Tm = 1. Tm
#export Kurama_field=-10. # only off case Tm = 1. Tm
#export Kurama_field=0. # only off case Tm = 1. Tm
#export Kurama_gap 250.0 # only off case
#export Kurama_pos_z 1426.5 #minimum value : 1326.5
export Kurama_pos_z=1500. #minimum value : 1326.5. From Uguard to center of kurama : 820
export Kurama_gap=400.0 # only off case
export Kurama_move_x=5. # unit mm 
#export Kurama_field=-0.875 # only off case Tm = 0.7 Tm
#export Kurama_field 0.875 # only off case Tm = 0.7 Tm
#export Kurama_gap 800.0 # only off case
#export Helmholtz_field 0.

#### pad_design, 
#1: 9, 16, 10, 17
#2: 9, 13, 10, 20
export pad_configure=2. # 1: simple, 2:from final design, 3:right-left divide
export pad_length_in=9.
export pad_length_out=12.5
export pad_num_in=10.
export pad_num_out=22.
export pad_gap=0.5
export pad_width_in=2.5
export pad_width_out=2.5

### Target material
#export Target_Material=C ## LH2, Cu, C
export Target_Material=LD2 ## LH2, Cu, C
## It should be added some code.

export H_decaytime=10.
export H_mass=2.250 #weak decay
export H_width=0.00 #weak decay
#export Weak_decay_mass 2.250 #for the generator #9

#### beam information
#export beam_position_x 0 #
#export beam_position_y 0 #
#export beam_position_z 0 #
#export Beam_mom=1.8 #
export Beam_mom=2.0 #
#export Beam_mom=0.8 #test K+ beam
export Beam_width=0.021 # turtle 
export Beam_x0=0.0 # assumption 
export Beam_y0=0.0 # assumption
export Beam_dx=5.3 # turtle 
export Beam_dy=2.65 # turtle half of x
export Beam_u0=0.00 # assumption 
export Beam_v0=0.00 # assumption 
export Beam_du=0.0186 # turtle
export Beam_dv=0.0093 # turtle half of x


#### Spectrometer angle //plus left , minus right
export Spectrometer_ang=0.
#export Spectrometer_ang 5.
#export Posz_ToF 3242.5 # note that target position is 242.3 mm
####export Angle_ToF -14.04
####export Angle_ToF 6.5 ####wo_kurama_case

#### GEM discharge on/off(1: design #1, 2: design #2, 3: design #3(now test on plane 3)0: off, 5:final design  )                                                               
export GEMDischarge=0.  ###design number
export DeadArea=50.  ### width of electrode
export GEMFixDead=1. ### on/off (1/0) when I fix the dead region
export GEMDeadPlane=3.  ###
export GEMDeadPlaneDivision=3.


## useage 
##char* OutFileName = getenv("Out_ROOT_File_Name");
##String to double
##G4String get_e = getenv("Gamma_Energy");
##Gamma_Energy = atof(get_e.c_str());

