#include <stdlib.h>
#include <stdio.h>

#include "matmul.h"

int main(int argc, char *argv[])
{
    printf("Starting %s...\n", argv[0]);

    double values[9] = {4, 5, 6, 0, 1.5, 1.0, 0, 0.5, 2};
    double ivaues[9] = {1, 0, 0, 0, 1, 0, 0, -(0.5/1.5), 1};

    double *matrix_a = values;
    double *matrix_b = ivaues;

    print_flat_matrix(matrix_b, 3);
    print_flat_matrix(matrix_a, 3);

    matrix_multiplication(matrix_b, matrix_a, 3, 3, 3);

    //double * c = dot_simple(matrix_b, 3, 3, matrix_a, 3, 3);

    print_flat_matrix(matrix_a, 3);

    return EXIT_SUCCESS;
}