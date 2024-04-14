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
*/

/// @brief Finds the pivot used in the currently provided linear system using a partial pivot strategy
/// @param linear_system the linear system
/// @param current_line
/// @param pivot_line
/// @return the pivot
double omp_select_current_pivot(linear_system_t *linear_system, int current_line, int *pivot_line)
{
    int nb_matrix_rows = linear_system->nb_unknowns;

    double global_max_pivot = linear_system->storage[current_line][current_line];
    int global_pivot_line = *pivot_line;
    //int step = 2;
    //double *max_arr = (double *)malloc(sizeof(double) * 0);
    #pragma omp parallel shared(global_max_pivot, global_pivot_line)
    {
        int thread_id = omp_get_thread_num();
        int start, stop;
        double local_max_pivot = global_max_pivot;
        int local_pivot_line = global_pivot_line;
        int num_threads = omp_get_num_threads();
        int step = (nb_matrix_rows)/num_threads;

        start = (thread_id * step) + (current_line + 1);
        stop = (start + step > nb_matrix_rows) ? nb_matrix_rows : start + step;

        //printf(">>>>>>>>>>>>>>>thread %d: start=%d, stop=%d\n", thread_id, start, stop);

        double pivot = 0;
        double abs_val;
        double col_val;
        double p_abs;
        for (int y = start; y < stop; y++)
        {
            col_val = linear_system->storage[y][current_line];
            abs_val = fabs(col_val); // check if the absolute value of the pivot coefficient to be selected
            p_abs = fabs(local_max_pivot);
            pivot = MAX(p_abs, abs_val);
            if (p_abs != pivot)
                local_pivot_line = y;                  // if the current pivot is updated, we update the line the pivot is on
            local_max_pivot = (pivot == abs_val) ? col_val : local_max_pivot; // determine whether to get the absolute value or not
        }

        
        #pragma omp critical
        {
            double abs_max_pivot = fabs(local_max_pivot);
            if(abs_max_pivot > fabs(global_max_pivot)) {
                global_max_pivot = local_max_pivot;
                global_pivot_line = local_pivot_line;
            }
        }

        #pragma omp barrier
    }
    *pivot_line = global_pivot_line;
    return global_max_pivot;
}