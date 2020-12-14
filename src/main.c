#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <omp.h>
#include <stdio.h>
#include <string.h>
#include <stdbool.h>
#include "../lib/pso.h"

#define for_i_swarm for (int i = 0; i < swarm_size; i++)
#define for_j_dimensions for (int j = 0; j < dimensions; j++)

const char m_max_flag[] = "-mig";
const char dimensions_flag[] = "-d";
const char swarm_size_flag[] = "-s";
const char max_value_flag[] = "-max";
const char verbose_flag[] = "-v";

// Square function
double objective_function(double *particle, int dimensions)
{
    double result = 0;

#pragma omp parallel for reduction(+ \
                                   : result)
    for_j_dimensions
    {
        result += pow(particle[j], 2);
    }
    return result;
}

double **generate_swarm(int swarm_size, int dimensions, double max_value)
{
    double **result;
    double random_value;
    unsigned int seed;

    result = malloc(swarm_size * dimensions * sizeof(double));

#pragma omp parallel for private(seed) shared(result)
    for_i_swarm
    {
        // trying to get unique seed
        seed = (time(NULL)) ^ (omp_get_thread_num() + 1 + 10 * i);
        result[i] = malloc(dimensions * sizeof(double));

        for_j_dimensions
        {
            // flush first value, since it is not really random
            rand_r(&seed);
            random_value = 2 * max_value * rand_r(&seed) / RAND_MAX - max_value;
            result[i][j] = random_value;
        }
    }
    return result;
}

int main(int argc, char *argv[])
{
    int dimensions = 30;
    double max_value = 50.0;
    int m_max = 100;
    int swarm_size = 30;
    double verbose = false;

    int i = 1;
    while (i < argc)
    {
        if ((i + 1) < argc)
        {
            if (strcmp(m_max_flag, argv[i]) == 0)
            {
                m_max = atoi(argv[i + 1]);
            }
            else if (strcmp(dimensions_flag, argv[i]) == 0)
            {
                dimensions = atoi(argv[i + 1]);
            }
            else if (strcmp(max_value_flag, argv[i]) == 0)
            {
                max_value = atoi(argv[i + 1]);
            }
            else if (strcmp(swarm_size_flag, argv[i]) == 0)
            {
                swarm_size = atoi(argv[i + 1]);
            }
            else if (strcmp(verbose_flag, argv[i]) == 0)
            {
                verbose = true;
                i--;
            }
            else
            {
                i--;
            }
            i += 2;
        }
        else
        {
            if (strcmp(verbose_flag, argv[i]) == 0)
            {
                verbose = true;
            }
            i++;
        }
    }

    if (verbose)
        printf("Generating swarm with %d individuals and %d dimensions\n", swarm_size, dimensions);
    double **swarm = generate_swarm(swarm_size, dimensions, max_value);

    printf("First iteration: %.2f\n", objective_function(swarm[0], dimensions));

    if (verbose)
        printf("Running pso algorithm with %d migrations\n", m_max);
    run_pso(swarm, &objective_function, swarm_size, dimensions, m_max, verbose);

    printf("Last iteration: %.2f\n", objective_function(swarm[0], dimensions));

    if (verbose)
        printf("Deallocating memory from main\n");

    for (int i = 0; i < swarm_size; i++)
    {
        free(swarm[i]);
    }
    free(swarm);

    return 0;
}
