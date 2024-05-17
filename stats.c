#include "stats.h"

const char plot_colors[10][10] = {
    "black",
    "red",
    "green",
    "blue",
    "yellow",
    "purple",
    "cyan",
    "pink",
    "orange",
    "violet"
};

/// @brief Plot solving statistics from given file
/// @param filename given file
void plot_stats(char *filename)
{
    FILE *gnuplot = popen("gnuplot", "w");
    if (!gnuplot)
    {
        perror("popen");
        exit(EXIT_FAILURE);
    }

    int color = 0;

    fprintf(gnuplot, "set title \"%s\" font \"%s\"\n", "Evolution de la fonction", "Helvetica,18");
    fprintf(gnuplot, "set xlabel \"Redémarrage de l'algorithme\"\n");
    fprintf(gnuplot, "set ylabel \"Valeur de la fonction\"\n");
    fprintf(gnuplot, "set xtics 0, %d\n", X_RANGE);

    fprintf(gnuplot, "plot \"%s\" t 'Function' with %s linewidth %d linecolor \"%s\"\n",
            filename, LINE_TYPE, LINE_THICKNESS, plot_colors[color]);

    fflush(gnuplot);
    fprintf(stdout, "Click Ctrl+d to quit...\n");
    getchar();

    pclose(gnuplot);
    exit(EXIT_SUCCESS);
}

/// @brief Plot multiple solving statistics from given file using the program arguments
/// The argv parameters must follow a precise format where :
///     - first give the name of the file 
///     - second the caption of the graph
/// @param plot_title the title of the plot
/// @param argc number of arguments
/// @param argv arrar of arguments
void plot_statistics(char * plot_title, int argc, char *argv[])
{
    int color = 0;
    char command_buffer[COMMANDE_SIZE + 1];
    char command[COMMANDE_SIZE + 1] = "plot";

    FILE *gnuplot = popen("gnuplot", "w");
    if (!gnuplot)
    {
        perror("popen");
        exit(EXIT_FAILURE);
    }

    for (int i = 1, j = 2; i < argc; i+=2, j+=2)
    {
        snprintf(command_buffer,
                 FILE_SIZE + 1,
                 " \"%s\" t '%s' with %s linewidth %d linecolor \"%s\",",
                 argv[i], argv[j], LINE_TYPE, LINE_THICKNESS, plot_colors[color]);

        strcat(command, command_buffer);
        color++;
        printf("%s\n", argv[i]);
    }
    command[strlen(command) - 1] = '\0';
    printf("%s\n", command);

    fprintf(gnuplot, "set title \"%s\" font \"%s\"\n", plot_title, "Helvetica,18");
    //fprintf(gnuplot, "set xlabel \"Redémarrage de l'algorithme\"\n");
    fprintf(gnuplot, "set xlabel \"Nombre de processus\"\n");
    //fprintf(gnuplot, "set ylabel \"Valeur de la fonction de cout\"\n");
    fprintf(gnuplot, "set ylabel \"Temps d'execution de l'algorithme\"\n");
    //fprintf(gnuplot, "set xtics 0, %d\n", X_RANGE);
    // fprintf(gnuplot, "set ytics 0, %d\n", Y_RANGE);
    fprintf(gnuplot, "%s\n", command);
    fflush(gnuplot);
    fprintf(stdout, "Click Ctrl+d to quit...\n");
    getchar();

    pclose(gnuplot);
    exit(EXIT_SUCCESS);
}

int main(int argc, char *argv[])
{
    //char * file4 = "./data/0000183b305c-05-10-2023-(17-22-31).txt"; //with keep best strategy
    //char * file5 = "./data/0000183b305c-05-10-2023-(17-22-53).txt"; //without keep best strategy

    //char * title = "Evolution de la fonction de cout du recuit simmulé";
    char * title = "Evolution du temps d'execution par nombre de processus (OpenMP) en moyenne (sur 10)";

    plot_statistics(title, argc, argv);
    //sudoku_plot_multiple_stats("%s", file4, file5, NULL);
    //sudoku_plot_multiple_stats("%s", "./data/0000183b305c-07-10-2023-(15-09-25).txt",NULL);
    //sudoku_plot_multiple_stats("%s", "./data/0000183b305c-07-10-2023-(15-18-15).txt",NULL);

    return EXIT_SUCCESS;
}