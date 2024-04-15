# Résolution de système linéaire #

**Name:** Solveur de système linéaire
**Version:** 0.2

Algorithme de résolution de système linéaire écrit en C qui utilise une méthode de pivot de Gauss avec parallélisation en OpenMP.

## Installation ##

Utilisation du programme **gem_matrix** afin de générer des matrices de tailles souhaitées (compilation inclus dans le Makefile).

Une utilisation générique du programme est la suivante :
```console
foo@bar:~$ make clean && make

Use: ./gen_matrix nb_unknownks
Where :
  nb_unknownks   : Le nombre d'inconnues dans la matrice du système linéaire généré

foo@bar:~$ ./gen_matrix 1024
```

## Utilisation ##

Le Makefile possède différentes options de compilation de 'make', notamment celle pour la version parallèle :
* 'omp' : Permet de compiler le programme 'main', dans lequel l'algorithme de résolution du système linéaire est lancé en utilisant tous les threads renseignés.
* 'debug' : Permet de compiler le programme 'main' avec un affichage plus détaillées (non recommandé pour des matrices à plus de 10 inconnus).

Une utilisation générique du programme est la suivante :
```console
foo@bar:~$ make clean && make

Use: ./main matrix_fileame
Where :
  matrix_fileame   : Le fichier contenant le système linéaire à résoudre

foo@bar:~$ ./main generated_out.txt
```

Pour une analyse plus fine du temps d'exécution en moyenne de la version séquentielle et parallèle (OpenMP), le test bench_avg_runtime est mise en place, permettant de générer automatiquement des graphiques pour visualisation des données (Installation au préalable de la librairie gnuplot: 
```console
foo@bar:~$ sudo apt-get update -y
foo@bar:~$ sudo apt-get install -y gnuplot
```
).

## Configuration ##

Des options de configurations, comme l'intervalle des nombres aléatoires générés comme coefficient du système linéaire, sont disponibles dans le fichier config.h. 

Afin que les modifications faites dans ce fichier soient prises en compte, il est recommandé d'effectuer la commande **make clean** au préalable.

## Execution ROMEO ##

Pour une exécution sur le supercalculateur ROMEO, faire attention de charger la version la plus récente du module gcc:
```bash
module load gcc/10.2.0_spack2021_gcc-10.2.0-gdnn
```
Les scripts de soumissions omp_job.sh et benchmark_job.sh permettent de tester respectivement l'exécution du solveur de système linéaire pour un fichier donné, et l'autre de lancer un job pour le script bash **bench_avg_runtime** permettant d'évaluer de le temps d'exécution en moyenne de la version séquentielle et parallèle (avec affichage de graphe si gnuplot installée sur la machine, pas disponible sur ROMEO...).

## Changelog ##

### 0.2 ###
* Version parallèle avec OpenMP

### 0.1 ###
* Version séquentielle