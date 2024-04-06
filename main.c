#include <stdlib.h>
#include <stdio.h>

#include "utils.h"

#define OUT_FILE "output.txt"

int main(int argc, char *argv[])
{
    printf("Starting %s...\n", argv[0]);

    if(argc < 2) {
        fprintf(stderr, "Incorrect arguments : No given filename!\n");
        exit(EXIT_FAILURE);
    }

    int s;
    char * filename = argv[1];

    //create a default empty linear system structure in which to store the contents read from the input file
    linear_system_t linear_system = {.nb_unknowns = 0, .data = NULL, .storage = NULL};

    read_linear_system_from_file(filename, &linear_system);

    //swap_linear_system_rows(&linear_system, 0, 1);
    print_linear_system_matrix(&linear_system);

    linear_system_propagation(&linear_system);

    double * solutions = solve_linear_system(&linear_system);

    for(s = 0; s < linear_system.nb_unknowns; s++) {
        printf("x%d = %.3lf\n", s+1, solutions[s]);
    }
    
    write_linear_system_to_file(OUT_FILE, &linear_system);
    
    clean_linear_system_memory(&linear_system);

    free(solutions);

    return EXIT_SUCCESS;
}