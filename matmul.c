#include "matmul.h"

//optimize CPU cache hits for flat row major matrix
void matrix_multiplication(double * matrix_a, 
                    double * matrix_b,
                    int a_height, int b_height, int b_width) {
    unsigned int i, j, k;
    unsigned int A = a_height;
    unsigned int B = b_width;
    unsigned int K = b_height;
    double * product = (double *) malloc(A * B * sizeof(double));
    for (i = 0; i < A; i++) {
        for (k = 0; k < K; k++) {
            for (j = 0; j < B; j++) {
                product[i * A + j] += matrix_a[i * A + k] * matrix_b[k * K + j];
                //printf("product[%d] %f += %f + %f * %f\n", i * A + j, product[i * A + j], matrix_a[i * A + k] * matrix_b[k * K + j], matrix_a[i * A + k], matrix_b[k * K + j]);
            }
        }
    }

    for (i = 0; i < A; i++) {
        for (j = 0; j < B; j++) {
            matrix_b[i * A + j] = product[i * A + j];
        }
    }

    free(product);
}

double * generate_identity_matrix(int size) {
    double * matrix = (double *) malloc(size * size * sizeof(double));
    for (int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++) {
            matrix[i * size + j] = (i == j) ? 1 : 0;
        }
    }
    return matrix;
}

void print_flat_matrix(double *matrix, int size) {
    int i, j;
    printf("\n[\n");
    for(i = 0; i < size; i++) {
        if(i < TRUNC_THRESHOLD || i > size-TRUNC_THRESHOLD/2) {
            printf("\t[");
            for(j = 0; j < size; j++) {
                if(j > TRUNC_THRESHOLD) {
                    printf("...");
                    break;
                } else {
                    int index = i * size + j;
                    if(j != size-1)
                        printf("%.2f, ", matrix[index]);
                    else
                        printf("%.2f", matrix[index]);
                }
            }
            printf("],\n");
        }
        if(i == TRUNC_THRESHOLD) printf("\t...\n");
    }
    printf("\n];\n");
}

void matrix_vector_multiplication(double *matrix_a, double **vector_b,
                                  int a_height, int a_width, int b_height) {
    unsigned int i, j;
    unsigned int A = a_height;
    unsigned int B = a_width;

    if (a_width != b_height) {
        printf("Error: Incompatible dimensions for matrix-vector multiplication.\n");
        return;
    }

    double *product = (double *)malloc(A * sizeof(double));

    // Perform matrix-vector multiplication
    for (i = 0; i < A; i++) {
        product[i] = 0.0; // Initialize the result element to zero
        for (j = 0; j < B; j++) {
            product[i] += matrix_a[i * B + j] * (*vector_b[j]);
        }
    }

    for (i = 0; i < A; i++) {
        (*vector_b[i]) = product[i];
    }

    free(product);
}

void print_matrix(double ** matrix, int width, int height) {
    int i, j;
    printf("\n[\n");
    for(i = 0; i < height; i++) {
        if(i < TRUNC_THRESHOLD || i > height-TRUNC_THRESHOLD/2) {
            printf("\t[");
            for(j = 0; j < width; j++) {
                if(j > TRUNC_THRESHOLD) {
                    printf("...");
                    break;
                } else {
                    if(j != width-1)
                        printf("%.2f, ", matrix[i][j]);
                    else
                        printf("%.2f", matrix[i][j]);
                }
            }
            printf("],\n");
        }
        if(i == TRUNC_THRESHOLD) printf("\t...\n");
    }
    printf("\n];\n");
}

void ptr_matrix_multiplication(double * matrix_a, 
                    double *** matrix_b,
                    int a_height, int b_height, int b_width) {
    unsigned int i, j, k;
    unsigned int A = a_height;
    unsigned int B = b_width;
    unsigned int K = b_height;
    double * product = (double *) calloc(sizeof(double), A * B);
    for (i = 0; i < A; i++) {
        for (k = 0; k < K; k++) {
            for (j = 0; j < B; j++) {
                product[i * A + j] += matrix_a[i * A + k] * (*matrix_b[k][j]);
                //printf("product[%d] %f += %f + %f * %f\n", i * A + j, product[i * A + j], matrix_a[i * A + k] * (*matrix_b[k][j]), matrix_a[i * A + k], (*matrix_b[k][j]));
            }
        }
    }

    for (i = 0; i < A; i++) {
        for (j = 0; j < B; j++) {
            (*matrix_b[i][j]) = product[i * A + j];
        }
    }

    free(product);
}