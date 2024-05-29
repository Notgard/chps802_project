#include <stdlib.h>
#include <stdio.h>

#include "utils.h"

#define OUT_FILE "generated_matrix.txt"

int main(int argc, char *argv[])
{
    printf("Starting %s...\n", argv[0]);

    if (argc < 2)
    {
        fprintf(stderr, "Incorrect arguments : No given number of unknowns of the linear system!\n");
        exit(EXIT_FAILURE);
    }

    int nb_unknowns = atoi(argv[1]);

    generate_random_linear_system(OUT_FILE, nb_unknowns, false);

    return EXIT_SUCCESS;
}