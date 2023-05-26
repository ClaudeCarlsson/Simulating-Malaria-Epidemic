#!/bin/bash -l

#SBATCH -M snowy
#SBATCH -A uppmax2023-2-13
#SBATCH -t 01:15:00
#SBATCH -p node
#SBATCH -n 16
#SBATCH -o strong_scaling_%j.txt

module load gcc openmpi 

echo "Simulating malaria"

mpirun --bind-to none -n 8 ./simulate 50000 output_histogram514200.csv output_timestat501420.txt


echo "Complete"