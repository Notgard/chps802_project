#ifndef __CUDA_UTILS_H__
#define __CUDA_UTILS_H__

#include <stdlib.h>
#include <stdio.h>
#include <limits.h>
#include <sys/time.h>
#include <stdbool.h>
#include <time.h>

#include "config.h"

/// @brief Stores the linear system read from the given file into a linear system structure
/// @param filename the given file's path
/// @param linear_system the linear system to create from the given file
void cuda_read_linear_system_from_file(char * input_filename, linear_system_t * linear_system);

/// @brief Frees the memory from the linear system structure members
/// @param linear_system 
void cuda_clean_linear_system_memory(linear_system_t * linear_system);

/// @brief Writes the result of the solved linear system into a file 
/// @param output_filename the ouput file's path
/// @param linear_system the linear system to write into the given output_file
/// @param solutions the solutions for the given linear system
void cuda_write_linear_system_to_file(char * output_filename, linear_system_t * linear_system, double * solutions);

#endif