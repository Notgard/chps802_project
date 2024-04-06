#include "utils.h"

/// @brief Stores the linear system read from the given file into a linear system structure
/// @param input_filename the given file's path
/// @param linear_system the linear system to create from the given file
void read_linear_system_from_file(char *input_filename, linear_system_t *linear_system)
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

    // insert the extracted linear system matrix into the structure
    linear_system->data = linear_system_matrix; // matrice augmentée
    linear_system->storage = pointer_storage_array;

    linear_system->nb_unknowns = linear_system_unknowns;
}

/// @brief Writes the result of the solved linear system into a file
/// @param output_filename the ouput file's path
/// @param linear_system the linear system to write into the given output_file
void write_linear_system_to_file(char *output_filename, linear_system_t *linear_system)
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
    for (i = 0; i < nb_matrix_rows; i++)
    {
        for (j = 0; j < nb_matrix_cols; j++)
        {
            if (linear_system->storage[i][j] < 1 && linear_system->storage[i][j] > -1)
            {                    // truncate the floating point value leading zero
                char buffer[20]; // Assuming a maximum length of the printed number
                double truncate_value = linear_system->storage[i][j];

                sprintf(buffer, "%.2lf", truncate_value);

                if ((status = fprintf(file, "%s ", buffer + 1)) == EOF)
                {
                    perror("Can't write linear system to output file");
                    exit(EXIT_FAILURE);
                }
            }
            else
            {
                // write the contents of the linear system matrix into the output file
                if ((status = fprintf(file, "%.2lf ", linear_system->storage[i][j])) == EOF)
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
    }

    if ((status = fclose(file)) == EOF)
    {
        fprintf(stderr, "Can't close file %s\n", output_filename);
        exit(EXIT_FAILURE);
    }
}

/// @brief frees the memory from the linear system structure members
/// @param linear_system
void clean_linear_system_memory(linear_system_t *linear_system)
{
    free(linear_system->storage);
    free(linear_system->data);
}


/// @brief Finds the pivot used in the currently provided linear system using a partial pivot strategy
/// @param linear_system the linear system
/// @param current_line 
/// @param pivot_line
/// @return the pivot
double select_current_pivot(linear_system_t *linear_system, int current_line, int * pivot_line)
{
    /*     
        int x, y;
        double max = INT_MIN;
        double abs_max;

        int nb_matrix_cols = nb_matrix_rows + 1;

        // loop over the linear system matrix
        for (y = 0; y < nb_matrix_rows; y++)
        {
            for (x = 0; x < nb_matrix_cols; x++)
            {
                if (y == x)
                { // check diagonal of the matrix
                    // use the absolute value of the current max value for compare unless == INT_MIN
                    abs_max = (max == INT_MIN) ? max : fabs(max);
                    if (fabs(linear_system->storage[y][x]) > abs_max)
                    { // retrieve the maximum coefficient and position from the diagonal
                        r = y, c = x;
                        max = linear_system->storage[y][x];
                    }
                }
            }
        }

        *row = r;
        *col = c; */

    int r = -1, c = -1;
    int y;

    int nb_matrix_rows = linear_system->nb_unknowns;

    double pivot = 0;
    double p = linear_system->storage[current_line][current_line]; //pivot in diagonal
    for (y = current_line + 1; y < nb_matrix_rows; y++) //check for max value in same column
    {
        pivot = MAX(p, linear_system->storage[y][current_line]);
        p = pivot;
    }

    for(y = 0; y < nb_matrix_rows; y++) {
        if(linear_system->storage[y][current_line] == pivot)
            *pivot_line = y;
    }

    return linear_system->storage[r][c];
}

/// @brief Swaps row1 with row2 and vice versa from the given linear system
/// @param linear_system_t the linear system
/// @param row1 row index to swap with row2
/// @param row2 row index to swap with row1
void swap_linear_system_rows(linear_system_t *linear_system, int row1, int row2)
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

/// @brief
/// @param linear_system
void linear_system_propagation(linear_system_t *linear_system)
{
    int y;

    int nb_matrix_rows = linear_system->nb_unknowns;
    int nb_matrix_cols = nb_matrix_rows + 1;

    double pivot;
    int pivot_line;
    for(y = 0; nb_matrix_rows; y++) {
        //selection du pivot pour chaque ligne de la matrice augmentée du systeme lineaire
        pivot = select_current_pivot(linear_system, y, &pivot_line);

        //changement de ligne pour que le pivot soit sur la diagonale

        //pivotage de la matrice
        //A[i][j] = A[i][j] - ( (A[i][pi] / A[pi][pi]) * A[pi][j] )
    }
}

/// @brief Prints the contents of the linear system matrix
/// @param linear_system the given linear system
void print_linear_system_matrix(linear_system_t *linear_system)
{
    int i, j;

    int nb_matrix_rows = linear_system->nb_unknowns;
    int nb_matrix_cols = nb_matrix_rows + 1;

    for (i = 0; i < nb_matrix_rows; i++)
    {
        for (j = 0; j < nb_matrix_cols; j++)
        {
            printf("%lf ", linear_system->storage[i][j]);
        }
        printf("\n");
    }
}