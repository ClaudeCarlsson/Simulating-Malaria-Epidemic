#!/bin/bash -l

#SBATCH -M snowy
#SBATCH -A uppmax2023-2-13
#SBATCH -p core -n 16
#SBATCH -t 0:30:00
#SBATCH -J simulate

module load gcc openmpi 

echo "Simulating malaria"


mpirun --bind-to none -n 16 ./simulate 3000000 output_histogram_3mil.csv output_timestat_3mil.txt
echo "Complete"