#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <omp.h>
#include <math.h>

#define for_i_swarm for (int i = 0; i < swarm_size; i++)
#define for_j_dimensions for (int j = 0; j < dimensions; j++)

// Space constants
const int dimensions = 5000;
const double max_value = 50.0;

// PSO constants
const int m_max = 500;
const int swarm_size = 500;
const double c1 = 2.0;
const double c2 = 2.0;
const double v_max = 2.0;
const double weight_s = 0.9;
const double weight_e = 0.4;

const int precalculated_inertia = weight_s * (weight_s - weight_e) / m_max;

// Inertia weight for iteration
double get_inertia(int iteration)
{
    return iteration * precalculated_inertia;
}

// Make sure that particle is in bounds of max_value
void push_to_bounds(double *particle, int dimensions)
{
    for_j_dimensions
    {
        if (particle[j] > max_value)
        {
            particle[j] = max_value;
        }
        else if (particle[j] < -max_value)
        {
            particle[j] = -max_value;
        }
    }
}

// // Levy function
// double objective_function(double *particle, int dimensions, double **swarm, int global_best_index)
// {
//     double inertia = 1 + (particle[0] - 1) / 4;
//     double result = pow(sin(M_PI * inertia), 2);

//     // current_ofe += 1;
//     // if (current_ofe > max_ofe)
//     // {
//     //     printf("\nFinished with %d objective function evaluations\n", current_ofe);
//     //     printf("\nFirst: %.2f\n", swarm[0][dimensions]);
//     //     for_j_dimensions
//     //     {
//     //         printf("%.2f, ", j, swarm[global_best_index][j]);
//     //     }
//     //     printf("\nBest: %.2f\n", swarm[global_best_index][dimensions]);
//     //     exit(0);
//     // }

//     for (int j = 0; j < dimensions - 1; j++)
//     {
//         result += pow(inertia - 1, 2) * (1 + 10 * pow(sin(M_PI * inertia + 1), 2));
//         inertia = 1 + (particle[j + 1] - 1) / 4;
//     }

//     result += pow(inertia - 1, 2) * (1 + pow(sin(2 * M_PI * inertia), 2));
//     return result;
// }

// // Rosenbrock
// double objective_function(double *particle, int dimensions, double **swarm, int global_best_index)
// {
//     double result = 0;

//     current_ofe += 1;
//     if (current_ofe > max_ofe)
//     {
//         printf("\nFinished with %d objective function evaluations\n", current_ofe);
//         printf("\nFirst: %.2f\n", swarm[0][dimensions]);
//         for_j_dimensions
//         {
//             printf("%d: %.2f", j, swarm[global_best_index][j]);
//         }
//         printf("\nBest: %.2f\n", swarm[global_best_index][dimensions]);
//         exit(0);
//     }

//     for (int j = 0; j < dimensions - 1; j++)
//     {
//         result += 100 * pow((particle[j + 1] - pow(particle[j], 2)), 2) + pow((particle[j] - 1), 2);
//     }
//     return result;
// }

// Square function
double objective_function(double *particle, int dimensions, double **swarm, int global_best_index)
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

// // Zakharov function
// double objective_function(double *data, int size)
// {
//     double sum1 = 0, sum2 = 0;

//     current_ofe += 1;
//     if (current_ofe > max_ofe)
//     {
//         printf("\nFinished with %d objective function evaluations\n", current_ofe);
//         exit(0);
//     }

// #pragma omp parallel for reduction(+ \
//                                    : sum1, sum2) shared(data)
//     for (int i = 0; i < size; i++)
//     {
//         sum1 += pow(data[i], 2);
//         sum2 += 0.5 * (i + 1) * data[i];
//     }

//     return sum1 + pow(sum2, 2) + pow(sum2, 4);
// }

double ***generate_swarm(int swarm_size, int dimensions, double max_value)
{
    double ***result;
    double random_value;

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

    for_i_swarm
    {
        result[0][i] = malloc(dimension_size + sizeof(double));
        result[1][i] = calloc(dimensions, sizeof(double));
        result[2][i] = malloc(dimension_size + sizeof(double));
        for_j_dimensions
        {
            random_value = 2 * max_value * rand() / RAND_MAX - max_value;
            result[0][i][j] = result[2][i][j] = random_value;
        }
        result[0][i][dimensions] = result[2][i][dimensions] = objective_function(result[0][i], dimensions, result[0], 0);
    }
    return result;
}

int find_global_best(double **swarm, int swarm_size, int particle_dimensions)
{
    int best_index = 0;
    for_i_swarm
    {
        if (swarm[i][particle_dimensions] < swarm[best_index][particle_dimensions])
        {
            best_index = i;
        }
    }
    return best_index;
}

void copy(double *src, double *dest, int dimensions)
{
    for_j_dimensions
    {
        dest[j] = src[j];
    }
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

int main()
{
    srand(time(NULL));
    double ***swarm = generate_swarm(swarm_size, dimensions, max_value);
    double inertia;
    int global_best_index = find_global_best(swarm[0], swarm_size, dimensions);

    printf("\nFirst: %.2f\n", swarm[0][0][dimensions]);
    printf("\nBest: %.2f\n", swarm[0][global_best_index][dimensions]);

    for (int m = 0; m < m_max; m++)
    {
        inertia = get_inertia(m);
        for_i_swarm
        {
            calculate_next_velocity(swarm[1][i], inertia, dimensions, swarm[0][global_best_index], swarm[2][i], swarm[0][i]);
            for_j_dimensions
            {
                swarm[0][i][j] = swarm[0][i][j] + swarm[1][i][j];
            }
            push_to_bounds(swarm[0][i], dimensions);
            swarm[0][i][dimensions] = objective_function(swarm[0][i], dimensions, swarm[0], global_best_index);
            if (swarm[0][i][dimensions] < swarm[2][i][dimensions])
            {
                copy(swarm[0][i], swarm[2][i], dimensions);
            }
            if (swarm[0][i][dimensions] < swarm[0][global_best_index][dimensions])
            {
                global_best_index = i;
            }
        }
    }

    printf("\nFirst: %.2f\n", swarm[0][0][dimensions]);
    // for_j_dimensions
    // {
    //     printf("%d: %.2f", j, swarm[0][global_best_index][j]);
    // }
    printf("\nBest: %.2f\n", swarm[0][global_best_index][dimensions]);

    return 0;
}