#include "omp_utils.h"

/// @brief Application de la technique de pivot en OpenMP
/// @param linear_system le système linéaire
/// @param pivot_line la ligne du pivot choisi dans le système linéaire
void omp_apply_pivot(linear_system_t *linear_system, int pivot_line)
{
    int i, j;

    int nb_matrix_rows = linear_system->nb_unknowns;

    double pivot = linear_system->storage[pivot_line][pivot_line];

#pragma omp parallel private(i)
    {
// propagation du pivot sur les lignes en dessous
#pragma omp for schedule(dynamic, nb_matrix_rows)
        for (i = pivot_line + 1; i < nb_matrix_rows; i++)
        {
// propagation du pivot sur les coefficient de chaque lignes
#pragma omp simd
            for (j = pivot_line + 1; j <= nb_matrix_rows; j++)
            {
                linear_system->storage[i][j] -= ((linear_system->storage[i][pivot_line] / pivot) * linear_system->storage[pivot_line][j]);
            }
            linear_system->storage[i][pivot_line] = 0;
        }
    }
}

/* TODO: à finir, l'objectif est de paralléliser la recherche du maximum 
dans la colonne du système linéaire comme vu en classe avec une arborescence.


/// @brief Finds the pivot used in the currently provided linear system using a partial pivot strategy
/// @param linear_system the linear system
/// @param current_line
/// @param pivot_line
/// @return the pivot
double omp_select_current_pivot(linear_system_t *linear_system, int current_line, int *pivot_line)
{
    int nb_matrix_rows = linear_system->nb_unknowns;

    int step = 2;
    double global_max_pivot;
    double *max_arr = (double *)malloc(sizeof(double) * 0);
#pragma omp parallel shared(max_arr, global_max_pivot)
    {
        int thread_id = omp_get_thread_num();
        int y, start, stop;
        double local_max_pivot;
        int num_threads = omp_get_num_threads();

        start = thread_id * step;
        stop = (thread_id == (num_threads - 1)) ? nb_matrix_rows : start + step;

#pragma omp for schedule(dynamic)
        for (int y = current_line + 1; y < nb_matrix_rows; y++)
        {
            col_val = linear_system->storage[y][current_line];
            abs_val = fabs(col_val); // check if the absolute value of the pivot coefficient to be selected
            p_abs = fabs(p);
            pivot = MAX(p_abs, abs_val);
            if (p_abs != pivot)
                *pivot_line = y;                  // if the current pivot is updated, we update the line the pivot is on
            p = (pivot == abs_val) ? col_val : p; // determine whether to get the absolute value or not
        }
    }
}
*/