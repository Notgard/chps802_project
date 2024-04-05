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

    char * filename = argv[1];

    //create a default empty linear system structure in which to store the contents read from the input file
    linear_system_t linear_system = {.nb_unknowns = 0, .data = NULL, .storage = NULL};

    read_linear_system_from_file(filename, &linear_system);

    write_linear_system_to_file(OUT_FILE, &linear_system);

    clean_linear_system_memory(&linear_system);

    return 0;
}