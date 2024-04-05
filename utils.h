#ifndef __UTILS_H__
#define __UTILS_H__

#include <stdlib.h>
#include <stdio.h>

/// @brief linear system wrapper structure which contains the system matrix as a 1D array (data) 
///        and a pointer array (storage) to facilitate access to elements inside the linear system matrix
typedef struct linear_system_t {
    int nb_unknowns;
    double * data;
    double ** storage;
} linear_system_t;


/// @brief Stores the linear system read from the given file into a linear system structure
/// @param filename the given file's path
/// @param linear_system the linear system to create from the given file
void read_linear_system_from_file(char * input_filename, linear_system_t * linear_system);

/// @brief frees the memory from the linear system structure members
/// @param linear_system 
void clean_linear_system_memory(linear_system_t * linear_system);

/// @brief Writes the result of the solved linear system into a file 
/// @param output_filename the ouput file's path
/// @param linear_system the linear system to write into the given output_file
void write_linear_system_to_file(char * output_filename, linear_system_t * linear_system);

#endif