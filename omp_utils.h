#ifndef __OMP_UTILS_H__
#define __OMP_UTILS_H__

#include "config.h"

#include <stdlib.h>
#include <stdio.h>
#include <limits.h>
#include <sys/time.h>
#include <stdbool.h>
#include <time.h>

#include <math.h>

#include <omp.h>

/// @brief Propagation de chaque ligne du système linéaire en OpenMP
/// @param linear_system 
void omp_linear_system_propagation(linear_system_t * linear_system);

/// @brief Application de la technique de pivot en OpenMP
/// @param linear_system le système linéaire 
/// @param pivot_line la ligne du pivot choisi dans le système linéaire
void omp_apply_pivot(linear_system_t * linear_system, int pivot_line);

/// @brief Résolution du système linéaire avec OpenMP
/// @param linear_system the given linear system
/// @return the solutions to the linear system inside an array
double * omp_solve_linear_system(linear_system_t * linear_system);

/// @brief Finds the pivot used in the currently provided linear system using a partial pivot strategy
/// @param linear_system the linear system
/// @param current_line
/// @param pivot_line
/// @return the pivot
double omp_select_current_pivot(linear_system_t *linear_system, int current_line, int *pivot_line);

#endif