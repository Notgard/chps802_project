#include <stdlib.h>
#include <stdio.h>

#include "utils.h"

#include <omp.h>

#define OUT_FILE "output.txt"

int main(int argc, char *argv[])
{
    printf("Starting %s...\n", argv[0]);

    if (argc < 2)
    {
        fprintf(stderr, "Incorrect arguments : No given filename!\n");
        exit(EXIT_FAILURE);
    }

    int s;
    double start, end;
    char *filename = argv[1];

  int nb_threads = 0;
  #pragma omp parallel shared(nb_threads)
  #pragma omp master
    nb_threads = omp_get_num_threads();
  fprintf(stdout, "Nb threads: %d\n", nb_threads);

    // create a default empty linear system structure in which to store the contents read from the input file
    linear_system_t linear_system = {.nb_unknowns = 0, .data = NULL, .storage = NULL};

    read_linear_system_from_file(filename, &linear_system);

    // swap_linear_system_rows(&linear_system, 0, 1);
    printf("Matrice augmentée générée:\n");
    print_linear_system_matrix(&linear_system);

    printf("Start of the linear solver triangulation...\n");

    start = wtime();

    linear_system_propagation(&linear_system);

    double *solutions = solve_linear_system(&linear_system);

    end = wtime();
    
    for (s = 0; s < linear_system.nb_unknowns; s++)
    {
        printf("x%d = %.3lf\n", s + 1, solutions[s]);
    }

    // Output
    printf("-----------------------------------------------------\n");
    printf(" Total solver runtime: %lf seconds\n", end - start);
    printf("-----------------------------------------------------\n");

    write_linear_system_to_file(OUT_FILE, &linear_system);

    clean_linear_system_memory(&linear_system);

    free(solutions);

    //generate_random_linear_system("generated_out.txt", 32, false);

    return EXIT_SUCCESS;
}