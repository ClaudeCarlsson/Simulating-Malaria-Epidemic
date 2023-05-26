# Simulating Malaria Epidemic


This project is a part of the course [Parallel and Distributed Programming](https://www.uu.se/en/admissions/freestanding-courses/course-syllabus/?kpid=48022&lasar=23%2F24&typ=1) at Uppsala Univeristy. 

The implementation is done in C and uses OpenMPI, the method for simulating are Monte Carlo in combination with Stochastic Simulation Algorithm.

# Usage

The command to use this program is:
mpirun --bind-to none -np <P> ./simulate <N> <output_histogram.csv> <output_timestat.txt>

The amount of processes <P> and the number of runs <N>

The histogram gets created in <output_histogram.csv> and the checkpoint timings gets created in <output_timestat.txt>
