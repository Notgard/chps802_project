#include "utils.h"
#include "omp_utils.h"

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
void write_linear_system_to_file(char *output_filename, linear_system_t *linear_system, double *solutions)
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
    }

    if (solutions != NULL)
    {
        for (i = 0; i < nb_matrix_rows; i++)
        {
            //write solutions to the linear system to the file (if they exist)
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
double select_current_pivot(linear_system_t *linear_system, int current_line, int *pivot_line)
{
    int y;

    int nb_matrix_rows = linear_system->nb_unknowns;

    double pivot = 0;
    double p = linear_system->storage[current_line][current_line]; // pivot in diagonal
    double abs_val;
    double col_val;
    double p_abs;
    for (y = current_line + 1; y < nb_matrix_rows; y++) // check for max value in same column
    {
        col_val = linear_system->storage[y][current_line];
        abs_val = fabs(col_val); // check if the absolute value of the pivot coefficient to be selected
        p_abs = fabs(p);
        pivot = MAX(p_abs, abs_val);
        if (p_abs != pivot)
            *pivot_line = y;                  // if the current pivot is updated, we update the line the pivot is on
        p = (pivot == abs_val) ? col_val : p; // determine whether to get the absolute value or not
    }

    return p;
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
    int curr_line;

    int nb_matrix_rows = linear_system->nb_unknowns;

    int pivot_line;

    // loop iterativly over each row of the linear system matrix
    for (curr_line = 0; curr_line < nb_matrix_rows; curr_line++)
    {
        pivot_line = curr_line;

        // selection du pivot pour chaque ligne de la matrice augmentée du systeme lineaire
#if _OMP_
        double pivot = omp_select_current_pivot(linear_system, curr_line, &pivot_line);
#elif _DEBUG_
        double pivot = select_current_pivot(linear_system, curr_line, &pivot_line);
#else
        select_current_pivot(linear_system, curr_line, &pivot_line);
#endif

#if _DEBUG_
        printf("\n[Ligne #%d] Valeur de pivot trouvée: %.3lf à la ligne %d\n",
               curr_line,
               pivot,
               pivot_line);
#endif
        // changement de ligne pour que le pivot soit sur la diagonale
        if (pivot_line != curr_line)
        {
            swap_linear_system_rows(linear_system, curr_line, pivot_line);
#if _DEBUG_
            printf("Permutation de la ligne %d avec la ligne %d:\n", curr_line, pivot_line);
            print_linear_system_matrix(linear_system);
#endif
        }

        // pivotage de la matrice
#if _OMP_
        omp_apply_pivot(linear_system, curr_line);
#else
        apply_pivot(linear_system, curr_line);
#endif

#if _DEBUG_
        printf("\naprès pivot:\n");
        print_linear_system_matrix(linear_system);
#endif
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
            printf("%-10.3lf ", linear_system->storage[i][j]);
        }
        printf("\n");
    }
}

/// @brief Application de la technique de pivot
/// @param linear_system le système linéaire
/// @param pivot_line la ligne du pivot choisi dans le système linéaire
void apply_pivot(linear_system_t *linear_system, int pivot_line)
{
    int i, j;

    int nb_matrix_rows = linear_system->nb_unknowns;

    double pivot = linear_system->storage[pivot_line][pivot_line];
    for (i = pivot_line + 1; i < nb_matrix_rows; i++) // propagation du pivot sur les lignes en dessous
    {
        for (j = pivot_line + 1; j <= nb_matrix_rows; j++) // propagation du pivot sur les coefficient de chaque lignes
        {
            // A[i][j] = A[i][j] - ( (A[i][pi] / A[pi][pi]) * A[pi][j] )
            linear_system->storage[i][j] = linear_system->storage[i][j] - ((linear_system->storage[i][pivot_line] / pivot) * linear_system->storage[pivot_line][j]);
        }
        linear_system->storage[i][pivot_line] = 0;
    }
}

/// @brief Solves the linear system
/// @param linear_system the given linear system
/// @return the solutions to the linear system inside an array
double *solve_linear_system(linear_system_t *linear_system)
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

// Get the current time in seconds since the Epoch
double wtime(void)
{
    struct timeval tv;
    gettimeofday(&tv, NULL);
    return tv.tv_sec + tv.tv_usec * 1e-6;
}

double rand_double(unsigned int *seed, int min, int max)
{
    return min + (double)rand_r(seed) / ((double)RAND_MAX / (max - min));
}

void generate_random_linear_system(char *output_filename, int size, bool solvable)
{
    int i, j;
    int status;
    double rand_val;
    unsigned int seed = (unsigned int)time(NULL);
    FILE *file;

    // open the file in write mode, creates the file if doesn't exist
    if (!(file = fopen(output_filename, "w+")))
    {
        fprintf(stderr, "Error opening file %s\n", output_filename);
        exit(EXIT_FAILURE);
    }

    if ((status = fprintf(file, "%d\n", size)) == EOF)
    {
        perror("Can't write content to output file");
        exit(EXIT_FAILURE);
    }

    if (solvable)
    {
        /*
        Whenever a linear system is considered to be solvable,
        i.e has at least one solution it means according to the Rouché–Capelli theorem
        that the linear system coefficient matrix and it's augmented matrix representation
        have the same rank. In that case would the linear system be solvable.
        Thus, we first have to store both of these matricies in memory and perform a
        matrix rank calculation in order to determine if they are the same
        (probably inside another seperate matrix rank function)
        */
        // TODO: -store both square matrix A and rectangular matrix A|b
        //-check the rank of both matrices
        //-call the function recursively again if the rank again
    }

    for (i = 0; i <= size; i++)
    {
        for (j = 0; j <= size; j++)
        {
            rand_val = rand_double(&seed, MIN_RAND_VAL, MAX_RAND_VAL);
            char *val = (j == size - 1) ? "%.3lf" : "%.3lf ";
            if ((status = fprintf(file, val, rand_val)) == EOF)
            {
                perror("Can't write content to output file");
                exit(EXIT_FAILURE);
            }
        }
        if ((status = fprintf(file, "\n")) == EOF)
        {
            perror("Can't write content to output file");
            exit(EXIT_FAILURE);
        }
    }

    if ((status = fclose(file)) == EOF)
    {
        fprintf(stderr, "Can't close file %s\n", output_filename);
        exit(EXIT_FAILURE);
    }
}