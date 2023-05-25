#!/bin/bash -l

#SBATCH -M snowy
#SBATCH -A uppmax2023-2-13
#SBATCH -p core -n 16
#SBATCH -t 0:30:00
#SBATCH -J simulate

module load gcc openmpi 

echo "Simulating malaria"

mpirun --bind-to none -n 1 ./simulate 50000 output_histogram1.csv output_timestat1.txt
mpirun --bind-to none -n 2 ./simulate 100000 output_histogram2.csv output_timestat2.txt
mpirun --bind-to none -n 4 ./simulate 200000 output_histogram3.csv output_timestat3.txt
mpirun --bind-to none -n 8 ./simulate 400000 output_histogram4.csv output_timestat4.txt
mpirun --bind-to none -n 16 ./simulate 800000 output_histogram5.csv output_timestat5.txt

echo "Complete"