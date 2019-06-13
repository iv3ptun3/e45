source Param_E27/e27_parameters_2703.sh

export Generator=2708
export Target_Size_x=30.
export Target_Size_y=10. 
export Target_Size_z=15. 
export Target_Material=C ## LH2, Cu, C

for((i=0; $i<100; i++))
{
    export Out_ROOT_File_Name=./rootfiles/Kp_reaction/test_woEM_2708_$i.root
    export Out_GEN_File_Name=./rootfiles/Kp_reaction/test_woEM_2708_gen_$i.root
    bsub -q s -o "/dev/null" ./bin/Linux-g++/hyptpc1 run2.mac
}

