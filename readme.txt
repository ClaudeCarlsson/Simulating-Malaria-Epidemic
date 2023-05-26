First, run the makefile:
makefile

The command to use this program is:
mpirun --bind-to none -np <P> ./simulate <N> <output_histogram.csv> <output_timestat.txt>

The program creates one .csv file and one .txt file.

The amount of processes <P> and the number of runs <N>

The histogram gets created in <output_histogram.csv> and the checkpoint timings gets created in <output_timestat.txt>

Insert your values and remove the <>