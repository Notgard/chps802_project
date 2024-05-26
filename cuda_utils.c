#include "cuda_utils.h"

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
        fprintf(stderr, "Reached End Of File\n");
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

    int i, j;
    int status;

    int nb_matrix_rows = linear_system->nb_unknowns;
    int nb_matrix_cols = nb_matrix_rows + 1;

    // open the file in write mode, creates the file if doesn't exist
    if (!(file = fopen(output_filename, "w+")))
    {
        fprintf(stderr, "Error opening file %s\n", output_filename);
        exit(EXIT_FAILURE);
    }

    // write the number of unknowns in the linear system
    if ((status = fprintf(file, "%d\n", linear_system->nb_unknowns)) == EOF)
    {
        perror("Can't write content to output file");
        exit(EXIT_FAILURE);
    }

    // loop over the 1D linear system matrix using the storage pointers
/*     for (i = 0; i < nb_matrix_rows; i++)
    {
        for (j = 0; j < nb_matrix_cols; j++)
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
