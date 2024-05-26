#include <stdlib.h>
#include <stdio.h>

#include <cuda.h>
#include <cuda_runtime.h>
#include <cuda_runtime_api.h>

#include "config.h"

#define THREADS_PER_BLOCK 16

/// @brief Stores the linear system read from the given file into a linear system structure
/// @param input_filename the given file's path
/// @param linear_system the linear system to create from the given file
void cuda_read_linear_system_from_file(char *input_filename, linear_system_t *linear_system)
{
    FILE *file;

    int i;
    int status;
    int linear_system_unknowns, nb_matrix_rows, nb_matrix_cols;
    double file_value_pointer;

    // open the file in read-only mode
    if (!(file = fopen(input_filename, "r")))
    {
        fprintf(stderr, "Error reading file %s!\n", input_filename);
        exit(EXIT_FAILURE);
    }

    // read linear system configuration from given file
    if ((status = fscanf(file, "%d", &linear_system_unknowns)) == EOF)
    {
        fprintf(stderr, "[%d]Reached End Of File\n", status);
        exit(EXIT_FAILURE);
    }

    nb_matrix_rows = linear_system_unknowns;
    nb_matrix_cols = nb_matrix_rows + 1;

    int total_values = (nb_matrix_rows) * (nb_matrix_cols);

    // allocate memory for storing structures and check for out of memory allocations
    double *linear_system_matrix = (double *)malloc(sizeof(double) * total_values);
    if (linear_system_matrix == NULL)
    {
        fprintf(stderr, "Out of memory!\n");
        exit(EXIT_FAILURE);
    }

    double **pointer_storage_array = (double **)malloc(sizeof(double *) * nb_matrix_rows);
    if (pointer_storage_array == NULL)
    {
        fprintf(stderr, "Out of memory!\n");
        exit(EXIT_FAILURE);
    }

    // store the pointers into the pointer storage array
    for (i = 0; i < nb_matrix_rows; i++)
    {
        pointer_storage_array[i] = &linear_system_matrix[i * nb_matrix_cols];
    }

    double **vector_b_storage = (double **)malloc(sizeof(double *) * nb_matrix_rows);
    if (vector_b_storage == NULL)
    {
        fprintf(stderr, "Out of memory!\n");
        exit(EXIT_FAILURE);
    }

    // store the pointers into the pointer storage array
    for (i = 0; i < nb_matrix_rows; i++)
    {
        vector_b_storage[i] = &linear_system_matrix[i * nb_matrix_cols + nb_matrix_cols - 1];
    }

    double ***matrix_a = (double ***)malloc(sizeof(double **) * nb_matrix_rows);
    if (matrix_a == NULL)
    {
        fprintf(stderr, "Out of memory!\n");
        exit(EXIT_FAILURE);
    }

    // store the pointers into the pointer storage array
    for (i = 0; i < nb_matrix_rows; i++)
    {
        matrix_a[i] = (double **)malloc(sizeof(double *) * nb_matrix_rows);
    }

    // read the values of the linear system from the file and write them to the matrix
    for (i = 0; i < total_values; i++)
    {
        // read values from file as double
        if ((status = fscanf(file, "%lf", &file_value_pointer)) == EOF)
        {
            fprintf(stderr, "Reached End Of File\n");
            exit(EXIT_FAILURE);
        }

        linear_system_matrix[i] = file_value_pointer;
    }

    // copy content of linear system matrix to matrix_a
    for (i = 0; i < nb_matrix_rows; i++)
    {
        for (int j = 0; j < nb_matrix_rows; j++)
        {
            matrix_a[i][j] = &linear_system_matrix[i * nb_matrix_cols + j];
        }
    }

    // insert the extracted linear system matrix into the structure
    linear_system->data = linear_system_matrix; // matrice augmentée
    linear_system->storage = pointer_storage_array;

    linear_system->nb_unknowns = linear_system_unknowns;
    linear_system->vec_B = vector_b_storage;
    linear_system->matrix_A = matrix_a;

    if ((status = fclose(file)) == EOF)
    {
        fprintf(stderr, "Can't close file %s\n", input_filename);
        exit(EXIT_FAILURE);
    }
}

/// @brief Writes the result of the solved linear system into a file
/// @param output_filename the ouput file's path
/// @param linear_system the linear system to write into the given output_file
/// @param solutions the solutions for the given linear system
void cuda_write_linear_system_to_file(char *output_filename, linear_system_t *linear_system, double *solutions)
{
    FILE *file;

    int i;
    int status;

    int nb_matrix_rows = linear_system->nb_unknowns;

    // open the file in write mode, creates the file if doesn't exist
    if (!(file = fopen(output_filename, "w+")))
    {
        fprintf(stderr, "Error opening file %s\n", output_filename);
        exit(EXIT_FAILURE);
    }

    // write the number of unknowns in the linear system
    if ((status = fprintf(file, "%d\n", linear_system->nb_unknowns)) == EOF)
    {
        fprintf(stderr, "status: %d\n", status);
        perror("Can't write content to output file");
        exit(EXIT_FAILURE);
    }

    // loop over the 1D linear system matrix using the storage pointers

/*     int nb_matrix_cols = nb_matrix_rows + 1;
    for (i = 0; i < nb_matrix_rows; i++)
    {
        for (int j = 0; j < nb_matrix_cols; j++)
        {
            if (linear_system->storage[i][j] != 0.0f)
            {
                // write the contents of the linear system matrix into the output file
                if ((status = fprintf(file, "%.3lf ", linear_system->storage[i][j])) == EOF)
                {
                    perror("Can't write linear system to output file");
                    exit(EXIT_FAILURE);
                }
            }
            else if (linear_system->storage[i][j] != 0.0f && (linear_system->storage[i][j] < 1 && linear_system->storage[i][j] > -1))
            {                    // truncate the floating point value leading zero
                char buffer[20]; // Assuming a maximum length of the printed number
                double truncate_value = linear_system->storage[i][j];
                sprintf(buffer, "%.3lf", truncate_value);

                if ((status = fprintf(file, "%s ", buffer + 1)) == EOF)
                {
                    perror("Can't write linear system to output file");
                    exit(EXIT_FAILURE);
                }
            }
        }
        // add linebreak to file for each written row from the linear system matrix
        if ((status = fprintf(file, "%s", "\n")) == EOF)
        {
            perror("Can't write linebreak to output file");
            exit(EXIT_FAILURE);
        }
    } */

    if (solutions != NULL)
    {
        for (i = 0; i < nb_matrix_rows; i++)
        {
            // write solutions to the linear system to the file (if they exist)
            if ((status = fprintf(file, "%.3lf ", solutions[i])) == EOF)
            {
                perror("Can't write linear system to output file");
                exit(EXIT_FAILURE);
            }
        }
    }

    if ((status = fclose(file)) == EOF)
    {
        fprintf(stderr, "Can't close file %s\n", output_filename);
        exit(EXIT_FAILURE);
    }
}

/// @brief frees the memory from the linear system structure members
/// @param linear_system
void cuda_clean_linear_system_memory(linear_system_t *linear_system)
{
    free(linear_system->storage);
    free(linear_system->data);
    free(linear_system->vec_B);
    free(linear_system->matrix_A);
}

void cuda_print_linear_system_matrix(linear_system_t *linear_system)
{
    int i, j;

    int nb_matrix_rows = linear_system->nb_unknowns;
    int nb_matrix_cols = nb_matrix_rows + 1;

    for (i = 0; i < nb_matrix_rows; i++)
    {
        for (j = 0; j < nb_matrix_cols; j++)
        {
            printf("%-10.3lf ", linear_system->storage[i][j]);
        }
        printf("\n");
    }
}

double *cuda_solve_linear_system(linear_system_t *linear_system)
{
    int j = 0;
    double result = 0;

    int nb_matrix_rows = linear_system->nb_unknowns;
    int nb_matrix_cols = nb_matrix_rows + 1;

    // allocate memory for the solutions array, initialized with zeros
    double *solutions = (double *)calloc(linear_system->nb_unknowns, sizeof(double));
    if (solutions == NULL)
    {
        fprintf(stderr, "Out of memory!\n");
        exit(EXIT_FAILURE);
    }

    for (int i = nb_matrix_rows - 1; i >= 0; i--)
    { // commencer à la dernière ligne
        result = 0;
        for (j = i; j < nb_matrix_cols - 1; j++)
        {
            if (i != j)
            { // if the the two aren't the samen then an initial solution has been found
                // result += A[i][j] * R[j]
                result += linear_system->storage[i][j] * solutions[j];
            }
        }
        // R[i]=(A[i][dim-1]-result)/A[i][i]
        solutions[i] = (linear_system->storage[i][nb_matrix_cols - 1] - result) / linear_system->storage[i][i];
    }
    return solutions;
}

void cuda_swap_linear_system_rows(linear_system_t *linear_system, int row1, int row2)
{
    int i;
    int nb_matrix_cols = linear_system->nb_unknowns + 1;

    double temp_value;

    for (i = 0; i < nb_matrix_cols; i++)
    {
        temp_value = linear_system->storage[row1][i];
        linear_system->storage[row1][i] = linear_system->storage[row2][i];
        linear_system->storage[row2][i] = temp_value;
    }
}

// gauss elimination function taking the linear system matrix and loading the pivot line into shared memory
__global__ void find_pivot_and_swap(double *d_linear_system, int n_rows, int n_cols, int current_line, int *d_pivot_line)
{
    extern __shared__ double shared_data[];

    int tid = threadIdx.x;
    int col = current_line;
    double *local_max = shared_data;
    int *local_pivot = (int *)&local_max[blockDim.x];

    // Step 1: Identify the pivot row
    double pivot = 0;
    double p = d_linear_system[current_line * n_cols + col];
    double abs_val;
    double col_val;
    double p_abs;
    int pivot_row = current_line;
    for (int row = tid + current_line + 1; row < n_rows; row += blockDim.x)
    {
        col_val = d_linear_system[row * n_cols + col];
        abs_val = fabs(col_val); // check if the absolute value of the pivot coefficient to be selected
        p_abs = fabs(p);
        pivot = MAX(p_abs, abs_val);
        if (p_abs != pivot)
            pivot_row = row;
        p = (pivot == abs_val) ? col_val : p;
    }

    __syncthreads();

    local_max[tid] = p;
    local_pivot[tid] = pivot_row;

    __syncthreads();

    // Reduce to find the maximum value and corresponding row index
    for (int stride = blockDim.x / 2; stride > 0; stride >>= 1)
    {
        if (tid < stride)
        {
            if (fabs(local_max[tid]) < fabs(local_max[tid + stride]))
            {
                local_max[tid] = local_max[tid + stride];
                local_pivot[tid] = local_pivot[tid + stride];
            }
        }
        __syncthreads();
    }

    __syncthreads();

    if(tid == 0) {
        *d_pivot_line = local_pivot[0];
    }
/*     if (tid == 0)
    {
        if (current_line != local_pivot[0])
        {
            // Swap rows current_line and pivot_row
            for (int j = 0; j < n_cols; j++)
            {
                double temp = d_linear_system[current_line * n_cols + j];
                d_linear_system[current_line * n_cols + j] = d_linear_system[local_pivot[0] * n_cols + j];
                d_linear_system[local_pivot[0] * n_cols + j] = temp;
            }
        }
    }
    __syncthreads(); */
}

// paralell selection of the max absolute value gaussian pivot
template <int BLOCK_SIZE>
__global__ void gauss_elimination(double *d_linear_system, int n_rows, int n_cols, int pivot_line)
{
    // Thread ID
    int row = blockIdx.x * blockDim.x + threadIdx.x;
    int col = blockIdx.y * blockDim.y + threadIdx.y;

    // Allocate shared memory for the pivot row
    extern __shared__ double pivot_row[];

    if (row < n_rows && col < n_cols)
    {
        // Load the pivot row into shared memory
        if (row == pivot_line)
        {
            pivot_row[col] = d_linear_system[pivot_line * n_cols + col];
        }
        __syncthreads();

        // Perform Gaussian elimination for elements in rows below the pivot line
        if (row > pivot_line)
        {
            double pivot_value = pivot_row[pivot_line];
            double factor = d_linear_system[row * n_cols + pivot_line] / pivot_value;

            if (col >= pivot_line)
            {
                d_linear_system[row * n_cols + col] -= factor * d_linear_system[pivot_line * n_cols + col];
            }
        }
        __syncthreads();
    }
}

double GBPerSec(int bytes, double sec)
{
    return (double)(bytes) / (1024. * 1024. * 1024.) / sec;
}

void cudaMallocCheckError(void **addr, size_t size)
{
    cudaError_t code;
    code = cudaMalloc(addr, size);
    if (code != cudaSuccess)
    {
        fprintf(stderr, "GPUassert: %s\n", cudaGetErrorString(code));
        exit(EXIT_FAILURE);
    }
}

void printCudaInfo()
{
    int deviceCount = 0;
    cudaError_t err = cudaGetDeviceCount(&deviceCount);

    printf("---------------------------------------------------------\n");
    printf("Found %d CUDA devices\n", deviceCount);

    for (int i = 0; i < deviceCount; i++)
    {
        cudaDeviceProp deviceProps;
        cudaGetDeviceProperties(&deviceProps, i);
        printf("Device %d: %s\n", i, deviceProps.name);
        printf("   SMs:        %d\n", deviceProps.multiProcessorCount);
        printf("   Global mem: %.0f MB\n",
               static_cast<float>(deviceProps.totalGlobalMem) / (1024 * 1024));
        printf("   CUDA Cap:   %d.%d\n", deviceProps.major, deviceProps.minor);
    }
    printf("---------------------------------------------------------\n");
}

void pivot_de_gauss(linear_system_t *h_linear_system)
{
    double *d_linear_system = NULL;
    int* d_pivot_line = NULL;

    cudaError_t code;
    clock_t start, end;
    double elapsed;
    int n = h_linear_system->nb_unknowns;
    int h_system_size = n * (n + 1);
    size_t h_size = (n * (n + 1)) * sizeof(double);

    int h_pivot_line;

    cudaMallocCheckError((void **)&d_linear_system, h_size);
    cudaMallocCheckError((void **)&d_pivot_line, sizeof(int));

    code = cudaMemcpy(d_linear_system, h_linear_system->data, h_size, cudaMemcpyHostToDevice);

    if (code != cudaSuccess)
    {
        fprintf(stderr, "GPUassert: %s\n", cudaGetErrorString(code));
        exit(EXIT_FAILURE);
    }

    const int nb_blocks = (h_size + THREADS_PER_BLOCK - 1) / THREADS_PER_BLOCK;

    dim3 blockSize(32, 32);
    dim3 gridSize((n + blockSize.x - 1) / blockSize.x, ((n + 1) + blockSize.y - 1) / blockSize.y);
    int sharedMemSize = (n + 1) * sizeof(double);

    start = clock();
    for (int curr_line = 0; curr_line < n; curr_line++)
    {
        int pivot_line = curr_line;

        find_pivot_and_swap<<<nb_blocks, THREADS_PER_BLOCK, sharedMemSize>>>(d_linear_system, n, n + 1, pivot_line, d_pivot_line);
        code = cudaGetLastError();
        if (code != cudaSuccess)
        {
            fprintf(stderr, "GPUassert: %s\n", cudaGetErrorString(code));
            exit(EXIT_FAILURE);
        }
        cudaDeviceSynchronize();

        cudaMemcpy(&h_pivot_line, d_pivot_line, sizeof(int), cudaMemcpyDeviceToHost);
        //printf("Pivot line %d: %d\n", curr_line, h_pivot_line);

        cudaMemcpy(h_linear_system->data, d_linear_system, h_size, cudaMemcpyDeviceToHost);
        if(pivot_line != h_pivot_line) {
            cuda_swap_linear_system_rows(h_linear_system, pivot_line, h_pivot_line);
        }

/*         printf("After pivot line %d:\n", curr_line);
        // print the h_linear_system->data
        for (int i = 0; i < n; i++)
        {
            for (int j = 0; j < n + 1; j++)
            {
                printf("%-10.3lf ", h_linear_system->data[i * (n + 1) + j]);
            }
            printf("\n");
        }
        printf("---------------------------------------\n"); */

        cudaMemcpy(d_linear_system, h_linear_system->data, h_size, cudaMemcpyHostToDevice);
        
        gauss_elimination<1024><<<gridSize, blockSize, sharedMemSize>>>(d_linear_system, n, n + 1, pivot_line);
        cudaDeviceSynchronize();
/*         cudaMemcpy(h_linear_system->data, d_linear_system, h_size, cudaMemcpyDeviceToHost);
        // print the h_linear_system->data
        printf("After gauss elimination\n");
        for (int i = 0; i < n; i++)
        {
            for (int j = 0; j < n + 1; j++)
            {
                printf("%-10.3lf ", h_linear_system->data[i * (n + 1) + j]);
            }
            printf("\n");
        }
        printf("---------------------------------------\n"); */
        code = cudaGetLastError();
        if (code != cudaSuccess)
        {
            fprintf(stderr, "GPUassert: %s\n", cudaGetErrorString(code));
            exit(EXIT_FAILURE);
        }
    }
    end = clock();
    elapsed = (double)(end - start) / CLOCKS_PER_SEC;
    printf("Gaussian elimination time in CUDA: %.5lf\n", elapsed);

    cudaMemcpy(h_linear_system->data, d_linear_system, h_size, cudaMemcpyDeviceToHost);
    printf("Final linear system: \n");
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n + 1; j++)
        {
            printf("%f ", h_linear_system->data[i * (n + 1) + j]);
        }
        printf("\n");
    }
}

#define OUT_FILE "output.txt"

int main(int argc, char *argv[])
{
    printf("Starting %s...\n", argv[0]);

    if (argc < 2)
    {
        fprintf(stderr, "Incorrect arguments : No given filename!\n");
        exit(EXIT_FAILURE);
    }

    char *filename = argv[1];
    clock_t start, end;
    double elapsed;

    printCudaInfo();

    linear_system_t linear_system = {.nb_unknowns = 0, .data = NULL, .storage = NULL};

    cuda_read_linear_system_from_file(filename, &linear_system);

    start = clock();
    pivot_de_gauss(&linear_system);
    end = clock();

    printf("\n------------------------------------------\n");

    cuda_print_linear_system_matrix(&linear_system);

    double *solutions = cuda_solve_linear_system(&linear_system);

    cuda_write_linear_system_to_file(OUT_FILE, &linear_system, solutions);

    cuda_clean_linear_system_memory(&linear_system);

    free(solutions);

    elapsed = (double)(end - start) / CLOCKS_PER_SEC;
    double GBs = GBPerSec(linear_system.nb_unknowns * (linear_system.nb_unknowns + 1) * sizeof(double), elapsed);
    printf("Overall execution time for %f Gigabytes per second : %.5lf\n", GBs, elapsed);
    printf("----------------------------------------------------------------------------------\n");

    return EXIT_SUCCESS;
}