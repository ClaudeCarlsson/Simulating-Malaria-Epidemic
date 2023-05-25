#!/bin/bash -l

#SBATCH -M snowy
#SBATCH -A uppmax2023-2-13
#SBATCH -p core -n 16
#SBATCH -t 0:30:00
#SBATCH -J simulate

module load gcc openmpi 

echo "Simulating malaria"

mpirun --bind-to none -n 16 ./simulate 1000000 output_histogram_1mil.csv output_timestat_1mil.txt
mpirun --bind-to none -n 16 ./simulate 2000000 output_histogram_2mil.csv output_timestat_2mil.txt
mpirun --bind-to none -n 16 ./simulate 4000000 output_histogram_4mil.csv output_timestat_4mil.txt

echo "Complete"