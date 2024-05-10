#ifndef __UTILS_H__
#define __UTILS_H__

#include <stdlib.h>
#include <stdio.h>
#include <limits.h>
#include <sys/time.h>
#include <stdbool.h>
#include <time.h>

#include <math.h>

#include "config.h"

/// @brief Stores the linear system read from the given file into a linear system structure
/// @param filename the given file's path
/// @param linear_system the linear system to create from the given file
void read_linear_system_from_file(char * input_filename, linear_system_t * linear_system);

/// @brief Frees the memory from the linear system structure members
/// @param linear_system 
void clean_linear_system_memory(linear_system_t * linear_system);

/// @brief Writes the result of the solved linear system into a file 
/// @param output_filename the ouput file's path
/// @param linear_system the linear system to write into the given output_file
/// @param solutions the solutions for the given linear system
void write_linear_system_to_file(char * output_filename, linear_system_t * linear_system, double * solutions);

/// @brief Finds the pivot used in the currently provided linear system using a partial pivot strategy
/// @param linear_system the linear system
/// @param current_line 
/// @param pivot_line
/// @return the pivot
double select_current_pivot(linear_system_t *linear_system, int current_line, int * pivot_line);

/// @brief Swaps row1 with row2 and vice versa from the given linear system
/// @param linear_system_t the linear system
/// @param row1 row index to swap with row2
/// @param row2 row index to swap with row1
void swap_linear_system_rows(linear_system_t * linear_system, int row1, int row2);

/// @brief 
/// @param linear_system 
void linear_system_propagation(linear_system_t * linear_system);

/// @brief Prints the contents of the linear system matrix
/// @param linear_system the given linear system
void print_linear_system_matrix(linear_system_t * linear_system);

/// @brief Application de la technique de pivot
/// @param linear_system le système linéaire 
/// @param pivot_line la ligne du pivot choisi dans le système linéaire
void apply_pivot(linear_system_t * linear_system, int pivot_line);

/// @brief Solves the linear system
/// @param linear_system the given linear system
/// @return the solutions to the linear system inside an array
double * solve_linear_system(linear_system_t * linear_system);

// Get the current time in seconds since the Epoch
double wtime(void);

double rand_double(unsigned int* seed, int min, int max);

//TODO : faire une fonction pour générer des systèmes linéaires résolvables 
//      - insertion directement dans un fichier en entrée
void generate_random_linear_system(char * output_filename, int size, bool solvable);

//prototype function for construct matrix pivot
void apply_matrix_pivot(linear_system_t *linear_system, int pivot_line, int pivot_col);

//cosntructs the elementary matrix used to pivot the linear system
void construct_matrix_pivot(linear_system_t *linear_system, double * e_matrix, int pivot_line, int pivot_col);

void linear_system_matrix_propagation(linear_system_t *linear_system);

#endif