#!/usr/bin/bash
#check for n_sudokus argument
if [[ $# -lt 2 || $# -gt 2 ]] ; then
	echo 'Wrong arguments provided, need filename and number of tries'
	exit 1
fi
exec_file="main"
input=$1
n_tries=$2
n=0
#make compile main solving C program
make clean && make
#set variables
stats_file="./benchmark/stats-seq-$$.txt"
for ((i = 1; i < n_tries; i++)); do
    result=$(./$exec_file "$input")
    #get the execution time
    exec_time=$(echo "$result" | grep "Total solver runtime:" | awk '{print $4}')
    printf "%d %.3f\n" "$i" "$exec_time" >> "$stats_file"
done
#create benchmark graph
./stats $stats_file "Sequentielle"