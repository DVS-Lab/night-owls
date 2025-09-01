#!/bin/bash
#PBS -l walltime=12:00:00
#PBS -N imageSimilarity
#PBS -q normal
#PBS -m ae
#PBS -M matt.mattoni@temple.edu
#PBS -l nodes=1:ppn=28

projectdir=/gpfs/scratch/tug87422/smithlab-shared/night-owls
scriptdir=$projectdir/code

cd $PBS_O_WORKDIR
umask 0000

module load python

python imageSimilarity.py
