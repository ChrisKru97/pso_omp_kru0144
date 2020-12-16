#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <stdio.h>
#include <omp.h>

#define for_i_swarm for (int i = 0; i < swarm_size; i++)
#define for_j_dimensions for (int j = 0; j < dimensions; j++)

// PSO constants
const double c1 = 2.0;
const double c2 = 2.0;
const double v_max = 2.0;
const double weight_s = 0.9;
const double weight_e = 0.4;

const int precalculated_inertia = weight_s * (weight_s - weight_e);

// Inertia weight for iteration
double get_inertia(int iteration, int m_max)
{
    return iteration * precalculated_inertia / m_max;
}

int find_global_best(double **swarm, int swarm_size, int particle_dimensions)
{
    int best_index = 0;
#pragma omp parallel for shared(best_index)
    for_i_swarm
    {
        if (swarm[i][particle_dimensions] < swarm[best_index][particle_dimensions])
        {
            best_index = i;
        }
    }
    return best_index;
}

void calculate_next_velocity(double *velocity, double inertia, int dimensions, double *global_best, double *personal_best, double *current_position)
{
    double next_velocity;
    for_j_dimensions
    {
        next_velocity = inertia * velocity[j] + c1 * rand() * (personal_best[j] - current_position[j]) / RAND_MAX + c2 * rand() * (global_best[j] - current_position[j]) / RAND_MAX;
        if (next_velocity > v_max)
        {
            next_velocity = v_max;
        }
        else if (next_velocity < -v_max)
        {
            next_velocity = -v_max;
        }
        velocity[j] = next_velocity;
    }
}

double ***prepare_data(double **swarm, double (*objective_function)(double *, int), int swarm_size, int dimensions)
{
    double ***result;

    int dimension_size = dimensions * sizeof(double);
    int values_size = swarm_size * (dimension_size + sizeof(double));
    int velocity_size = swarm_size * dimension_size;

    result = malloc(2 * values_size + velocity_size);
    // Current positions
    result[0] = malloc(values_size);
    // Current velocity, using calloc to make sure initial velocity is set to 0
    result[1] = calloc(swarm_size * dimensions, sizeof(double));
    // Personal best position in history
    result[2] = malloc(values_size);

#pragma omp parallel for default(shared)
    for_i_swarm
    {
        result[0][i] = malloc(dimension_size + sizeof(double));
        result[1][i] = calloc(dimensions, sizeof(double));
        result[2][i] = malloc(dimension_size + sizeof(double));

        memcpy(&result[0][i], &swarm[i], dimension_size);
        result[0][i][dimensions] = objective_function(result[0][i], dimensions);
        memcpy(&result[2][i], &result[0][i], dimension_size + sizeof(double));
    }

    return result;
}

void deallocate_data(double ***swarm, int swarm_size)
{
    for (int k = 0; k < 3; k++)
    {
        /*for_i_swarm
    {
        free(swarm[k][i]);
    }*/
        free(swarm[k]);
    }
    free(swarm);
}

void run_pso(double **swarm, double (*objective_function)(double *, int), int swarm_size, int dimensions, int m_max, bool verbose)
{
    if (verbose)
        printf("Preparing data for PSO\n");
    double ***swarm_data = prepare_data(swarm, objective_function, swarm_size, dimensions);
    double inertia;
    if (verbose)
        printf("Looking for best individual in first generation\n");
    int global_best_index = find_global_best(swarm_data[0], swarm_size, dimensions);
    if (verbose)
        printf("Best individual is %d with value %.2f\n", global_best_index, swarm_data[0][global_best_index][dimensions]);

    for (int m = 0; m < m_max; m++)
    {
        inertia = get_inertia(m, m_max);
#pragma omp parallel for default(shared)
        for_i_swarm
        {
            calculate_next_velocity(swarm_data[1][i], inertia, dimensions, swarm_data[0][global_best_index], swarm_data[2][i], swarm_data[0][i]);

            for_j_dimensions
            {
                swarm_data[0][i][j] = swarm_data[0][i][j] + swarm_data[1][i][j];
            }

            swarm_data[0][i][dimensions] = objective_function(swarm_data[0][i], dimensions);

            if (swarm_data[0][i][dimensions] < swarm_data[2][i][dimensions])
            {
                memcpy(&swarm[2][i], &swarm[0][i], (dimensions + 1) * sizeof(double));
            }

            if (swarm_data[0][i][dimensions] < swarm_data[0][global_best_index][dimensions])
            {
                global_best_index = i;
            }
        }
    }

    if (verbose)
        printf("Dealocatting data in pso\n");
    deallocate_data(swarm_data, swarm_size);
}
