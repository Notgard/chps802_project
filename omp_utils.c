#include "omp_utils.h"

/// @brief Application de la technique de pivot en OpenMP
/// @param linear_system le système linéaire
/// @param pivot_line la ligne du pivot choisi dans le système linéaire
void omp_apply_pivot(linear_system_t *linear_system, int pivot_line)
{
    int i, j;
    int size_page=4096/sizeof(double);

    int nb_matrix_rows = linear_system->nb_unknowns;

    double pivot = linear_system->storage[pivot_line][pivot_line];

    #pragma omp parallel private(i)
    {
        // propagation du pivot sur les lignes en dessous
        #pragma omp for schedule(static, size_page)
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

    int global_pivot_line = *pivot_line;
    double global_max_pivot = linear_system->storage[current_line][current_line];
    #pragma omp parallel shared(global_max_pivot, global_pivot_line)
    {
        int start, stop;
        int thread_id = omp_get_thread_num();
        double local_max_pivot = global_max_pivot;
        int local_pivot_line = global_pivot_line;

        int num_threads = omp_get_num_threads();
        int step = (nb_matrix_rows + num_threads - 1)/num_threads;

        start = (thread_id * step) + (current_line + 1);
        stop = start + step;
        if(stop > nb_matrix_rows) stop = nb_matrix_rows;

        //printf(">>>>>>>>>>>>>>>thread %d: start=%d, stop=%d %d\n", thread_id, start, stop, nb_matrix_rows);

        double abs_val;
        double col_val;
        double p_abs;
        //#pragma omp for schedule(static, step) this creates further problems if the step can't be divided evenly
        for (int y = start; y < stop; y++)
        {
            col_val = linear_system->storage[y][current_line];
            abs_val = fabs(col_val); // check if the absolute value of the pivot coefficient to be selected
            p_abs = fabs(local_max_pivot);
            if(abs_val > p_abs) {
                local_pivot_line = y;
                local_max_pivot = col_val;
            }
        }

        double abs_max_pivot = fabs(local_max_pivot);
        if(abs_max_pivot > fabs(global_max_pivot)){
            #pragma omp critical
            {
                abs_max_pivot = fabs(local_max_pivot);
                if(abs_max_pivot > fabs(global_max_pivot)) {
                    global_max_pivot = local_max_pivot;
                    global_pivot_line = local_pivot_line;
                }
            }
        }
    }
    *pivot_line = global_pivot_line;
    return global_max_pivot;
}