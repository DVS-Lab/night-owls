#!/bin/bash
#PBS -l walltime=6:00:00
#PBS -N tedana
#PBS -q normal
#PBS -m ae
#PBS -M matt.mattoni@temple.edu
#PBS -l nodes=1:ppn=28
cd $PBS_O_WORKDIR


all_cmds=(
    /gpfs/scratch/tug87422/smithlab-shared/night-owls/logs/cmd_tedana_sub-101_ses-10.txt
    /gpfs/scratch/tug87422/smithlab-shared/night-owls/logs/cmd_tedana_sub-104_ses-03.txt
    /gpfs/scratch/tug87422/smithlab-shared/night-owls/logs/cmd_tedana_sub-104_ses-06.txt
    /gpfs/scratch/tug87422/smithlab-shared/night-owls/logs/cmd_tedana_sub-105_ses-10.txt
)

for cmdfile in "${all_cmds[@]}"; do
    torque-launch "$cmdfile"
done
