#ifndef __MATMUL_H__
#define __MATMUL_H__

#include <stdlib.h>
#include <stdio.h>

#define TRUNC_THRESHOLD 10

void matrix_multiplication(double * matrix_a, 
                    double * matrix_b,
                    int a_height, int b_height, int b_width);

void ptr_matrix_multiplication(double * matrix_a, 
                    double *** matrix_b,
                    int a_height, int b_height, int b_width);

void matrix_vector_multiplication(double *matrix_a, double **vector_b,
                                  int a_height, int a_width, int b_height);

double * generate_identity_matrix(int size);

void print_flat_matrix(double *matrix, int size);

void print_matrix(double ** matrix, int width, int height);

#endif