#include <stdlib.h>
#include <stdio.h>

#include "utils.h"
#include "omp_utils.h"

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

    double start, end;
    char *filename = argv[1];
    
#if _OMP_
  int nb_threads = 0;
  #pragma omp parallel shared(nb_threads)
  #pragma omp master
    nb_threads = omp_get_num_threads();
  fprintf(stdout, "Nb threads: %d\n", nb_threads);
#endif

    // create a default empty linear system structure in which to store the contents read from the input file
    linear_system_t linear_system = {.nb_unknowns = 0, .data = NULL, .storage = NULL};

    read_linear_system_from_file(filename, &linear_system);

    start = wtime();

    int line = 0;
    double pivot;
    pivot = select_current_pivot(&linear_system, 0, &line);
    printf("found pivot line: %d (%lf)\n", line, pivot);
    end = wtime();

    // Output
    printf("-----------------------------------------------------\n");
    printf(" select pivot exec time : %lf seconds\n", end - start);
    printf("-----------------------------------------------------\n");

    start = wtime();

    pivot = omp_select_current_pivot(&linear_system, 0, &line);
    printf("found pivot line: %d (%lf)\n", line, pivot);
    end = wtime();

    // Output
    printf("-----------------------------------------------------\n");
    printf(" OMP select pivot exec time: %lf seconds\n", end - start);
    printf("-----------------------------------------------------\n");

    clean_linear_system_memory(&linear_system);

    //generate_random_linear_system("generated_out.txt", 1024, false);

    return EXIT_SUCCESS;
}