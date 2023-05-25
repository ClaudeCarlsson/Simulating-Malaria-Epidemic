/* Made by Claude Carlsson */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <mpi.h>
#include <limits.h>
#include "prop.h"

void accumulate_time(MPI_Win *window, double *time_run, int time_counter, int size, int rank, int driver_rank)
{
    MPI_Win_fence(0, *window);
    MPI_Accumulate(time_run, 1, MPI_DOUBLE, driver_rank, (rank + time_counter * size), 1, MPI_DOUBLE, MPI_SUM, *window);
    MPI_Win_fence(0, *window);
}

void generate_random_numbers(double *u1, double *u2)
{
    *u1 = (double)rand() / (double)RAND_MAX;
    *u2 = (double)rand() / (double)RAND_MAX;
}

void add_states(int dest[], int source[], int column)
{
    for (int i = 0; i < column; i++)
    {
        dest[i] += source[i];
    }
}

void copy_states(int *dest[], int source[], int column, int index)
{
    for (int j = 0; j < column; j++)
    {
        dest[j][index] = source[j];
    }
}

int main(int argc, char *argv[])
{
    // MPI Initialization
    MPI_Init(&argc, &argv);

    // Get rank and size
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    // Set driver rank
    int driver_rank = 0;

    // Check the input

    if (rank == driver_rank)
    {
        if (argc != 4)
        {
            printf("Usage: ./simulate <N> <output_file> <print_info>\n");
            return 1;
        }
    }

    /* Initialize all values*/

    // Set seed
    srand(time(NULL) + rank + 1998);

    // Input values
    int N = atoi(argv[1]);
    char *output_file = argv[2];
    int print_info = atoi(argv[3]);

    // Dimension/Amount values
    int P_row = 15;
    int P_col = 7;
    int w_size = P_row;
    int checkpoint_size = 4;
    int bins_amount = 20;
    int runs_per_process = N / size;
    int times_array_size = size * checkpoint_size;

    // Arrays
    int P[15][7] =
        {
            {1, 0, 0, 0, 0, 0, 0},
            {-1, 0, 0, 0, 0, 0, 0},
            {-1, 0, 1, 0, 0, 0, 0},
            {0, 1, 0, 0, 0, 0, 0},
            {0, -1, 0, 0, 0, 0, 0},
            {0, -1, 0, 1, 0, 0, 0},
            {0, 0, -1, 0, 0, 0, 0},
            {0, 0, -1, 0, 1, 0, 0},
            {0, 0, 0, -1, 0, 0, 0},
            {0, 0, 0, -1, 0, 1, 0},
            {0, 0, 0, 0, -1, 0, 0},
            {0, 0, 0, 0, -1, 0, 1},
            {0, 0, 0, 0, 0, -1, 0},
            {1, 0, 0, 0, 0, 0, -1},
            {0, 0, 0, 0, 0, 0, -1}};
    double w[w_size];
    int **simulation_states = calloc(P_col, sizeof(int *));
    for (int i = 0; i < P_col; i++)
    {
        simulation_states[i] = calloc(runs_per_process, sizeof(int));
    }

    int x_init[7] = {900, 900, 30, 330, 50, 270, 20};
    int x[7];
    int time_checkpoints[4] = {25, 50, 75, 100};
    double *time_array = (double *)calloc(times_array_size, sizeof(double));

    // Time variables
    int T = 100;
    double tau = 0.0;
    double t = 0.0;
    double start_time, end_time;
    int time_counter = 0;

    // Other values
    int r;
    double a0, u1, u2;

    // Additional MPI inits
    MPI_Win window;
    MPI_Win_create(time_array, times_array_size * sizeof(double), sizeof(double), MPI_INFO_NULL, MPI_COMM_WORLD, &window);
    MPI_Win_fence(0, window);

    /*Initialization complete*/

    start_time = MPI_Wtime();

    for (int i = 0; i < runs_per_process; i++)
    {
        t = 0.0;
        time_counter = 0;
        memcpy(x, x_init, sizeof(x_init));

        double time_run = 0;
        double time_elapse = MPI_Wtime();

        while (t < T)
        {
            if (time_counter < checkpoint_size && t > time_checkpoints[time_counter])
            {
                // Accumulate time_elapse into the time_array
                time_run = MPI_Wtime() - time_elapse;
                // Convert to milliseconds and average it
                time_run = time_run * 1000 / runs_per_process;
                accumulate_time(&window, &time_run, time_counter, size, rank, driver_rank);
                time_counter++;
            }

            // Compute reaction rates, w
            prop(x, w);

            // Compute a0
            a0 = 0.0;
            for (int j = 0; j < P_row; j++)
            {
                a0 += w[j];
            }

            // Generate random numbers
            generate_random_numbers(&u1, &u2);

            // Compute tau
            if (a0 != 0.0)
            {
                tau = -log(u1) / a0;
            }
            else
            {
                tau = 0.0;
            }

            // Compute r
            double cumsum = 0.0;
            for (r = 0; r < P_row; r++)
            {
                cumsum += w[r];
                if (cumsum / a0 > u2)
                {
                    break;
                }
            }

            // Update the state vector x
            add_states(x, P[r], P_col);

            // Increment time
            t += tau;
        }
        // Copy states
        copy_states(simulation_states, x, P_col, i);

        // Accumulate time_elapse into the time_array
        time_run = MPI_Wtime() - time_elapse;
        // Convert to milliseconds and average it
        time_run = time_run * 1000 / runs_per_process;
        accumulate_time(&window, &time_run, time_counter, size, rank, driver_rank);
    }

    /* Gather information about the run */
    int local_min = INT_MAX;
    int local_max = INT_MIN;
    int *local_susceptible = (int *)calloc(runs_per_process, sizeof(int));

    // Find local min and max for the bins
    for (int i = 0; i < runs_per_process; i++)
    {
        local_susceptible[i] = simulation_states[0][i];

        if (local_min > local_susceptible[i])
        {
            local_min = local_susceptible[i];
        }
        if (local_max < local_susceptible[i])
        {
            local_max = local_susceptible[i];
        }
    }

    /* Prepare and create a histogram */
    double interval;
    int *global_bins, *local_frequencies, *global_frequencies;
    int global_min, global_max, bin_index;

    // Find global min and max for the bins
    MPI_Allreduce(&local_min, &global_min, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
    MPI_Allreduce(&local_max, &global_max, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);

    // Create interval
    interval = (double)(global_max - global_min) / bins_amount;

    // Create local histograms
    global_bins = (int *)calloc(bins_amount + 1, sizeof(int));
    local_frequencies = (int *)calloc(bins_amount, sizeof(int));
    global_frequencies = (int *)calloc(bins_amount, sizeof(int));

    // Create the bin intervals
    for (int i = 0; i < bins_amount + 1; i++)
    {
        global_bins[i] = global_min + i * interval;
        if (i == bins_amount + 1)
        {
            global_bins[i] = global_max;
        }
    }

    for (int i = 0; i < runs_per_process; i++)
    {
        // Locate the index of the bin
        bin_index = (local_susceptible[i] - global_min) / interval;

        // Increment the bin
        local_frequencies[bin_index] += 1;
    }

    // Combine the local histograms
    MPI_Reduce(local_frequencies, global_frequencies, bins_amount, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);

    /* Prepare to print statistics */
    end_time = MPI_Wtime() - start_time;
    double max_time;
    MPI_Reduce(&end_time, &max_time, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

    // Print time statistics
    if (rank == driver_rank && print_info == 1)
    {
        printf("\n-------- Information --------\n\n");
        printf("Longest process time: %.2lf\n\n", max_time);
        printf("The average checkpoint times (milliseconds)\n");

        // Print header line with checkpoint labels
        printf("\t\t");
        for (int i = 0; i < checkpoint_size; ++i)
        {
            printf("chkp %d\t\t", time_checkpoints[i]);
        }
        printf("\n");
        // For each rank, print the rank and then the time for each checkpoint
        for (int rank_idx = 0; rank_idx < size; rank_idx++)
        {
            printf("Rank %d\t\t", rank_idx);

            for (int time_idx = 0; time_idx < checkpoint_size; time_idx++)
            {
                printf("%.3f\t\t", time_array[rank_idx + time_idx * size]);
            }

            printf("\n");
        }
    }
    // Print only time
    if (rank == driver_rank && print_info == 2)
    {
        printf("%lf\n", max_time);
    }

    // Write histogram to file
    if (rank == driver_rank)
    {
        FILE *file = fopen(output_file, "w");

        // Check if the file was opened successfully
        if (file == NULL)
        {
            printf("Error opening output file\n");
            return 1;
        }

        // Print histogram bins and frequencies
        fprintf(file, "bin,amount\n");
        for (int i = 0; i < bins_amount; i++)
        {
            fprintf(file, "%d,%d\n", global_bins[i], global_frequencies[i]);
        }

        // Close the file
        int result = fclose(file);
        if (result == EOF)
        {
            printf("Error closing output file\n");
            return 1;
        }
    }

    // Free all allocated memory
    for (int i = 0; i < P_col; i++)
    {
        free(simulation_states[i]);
    }

    free(simulation_states);
    free(local_susceptible);
    free(time_array);
    free(global_bins);
    free(local_frequencies);
    free(global_frequencies);

    // Free MPI window
    MPI_Win_free(&window);

    // Finalize the MPI environment.
    MPI_Finalize();

    return 0;
}