#ifndef PSO
#define PSO

void deallocate_data(double ***swarm, int swarm_size);

void run_pso(double **swarm, double (*objective_function)(double *, int), int swarm_size, int dimensions, int m_max, bool verbose);

double ***prepare_data(double **swarm, double (*f)(double *, int), int swarm_size, int dimensions);

void calculate_next_velocity(double *velocity, double inertia, int dimensions, double *global_best, double *personal_best, double *current_position);

int find_global_best(double **swarm, int swarm_size, int particle_dimensions);

double get_inertia(int iteration, int m_max);

#endif
