#ifndef __STATS_H__
#define __STATS_H__

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>

#define COMMANDE_SIZE 4086
#define LINE_THICKNESS 2
#define LINE_TYPE "lines"
#define X_RANGE 5
#define Y_RANGE 10
#define FILE_SIZE 256

#define GET_OUTPUT false

/// @brief Plot solving statistics from given file
/// @param filename given file
void plot_stats(char *filename);

/// @brief Plot multiple solving statistics from given file using the program arguments
/// The argv parameters must follow a precise format where :
///     - first give the name of the file 
///     - second the caption of the graph
/// @param plot_title the title of the plot
/// @param argc number of arguments
/// @param argv arrar of arguments
void plot_statistics(char * plot_title, int argc, char *argv[]);

#endif