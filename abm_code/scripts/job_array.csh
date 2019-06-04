#!/bin/csh

#$ -pe smp 12     # Specify parallel environment and legal core size
#$ -q long        # Specify queue (use ‘short’ for development)
#$ -N dengue_sim  # Specify job name
#$ -t 1-1000      # Specify number of tasks in array

module load gcc/4.8.2
module load gsl/1.15
module load mpich2

mpiexec -n 10 ./DengueSim SimControl.csv