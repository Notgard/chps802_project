# Résolution de système linéaire #

**Name:** Solveur de système linéaire
**Version:** 0.1

Algorithme de résolution de système linéaire écrit en CUDA qui utilise une méthode de pivot de Gauss.

## Utilisation ##

Les principaux programmes CUDA pour la résolution de système linéaire sont "latest_main.cu" et "dev_main_cuda.cu" qui sont respectivement la version la stable la plus récente et la version de développement, qui n'a pas été testé autant que la version stable et pour contenir des bugs pour d'autres instances de système linéaire.

Un Makefile est mis à votre disposition afin de compiler cette version CUDA par la commande suivante:

make clean && make

De plus, un nouveau script bash "bench_runtime.sh" vous est mis à disposition afin d'effectuer une analyse plus fine des performances en temps d'exécution du programme, pour des tailles de matrices augmentant de manière itérative.