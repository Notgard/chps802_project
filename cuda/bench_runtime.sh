#!/usr/bin/bash
#check for n_sudokus argument
if [[ $# -lt 1 || $# -gt 1 ]] ; then
	echo 'Wrong arguments provided, need n_tries'
	exit 1
fi
exec_file="dev_main_cuda"
input="generated_matrix.txt"
n=0
n_avg=2
n_tries=$1
#make compile main solving C program
make clean && make
#set variables
stats_file="./benchmark/bench_cuda_dev.txt"
echo -e "Starting benchmark for file $input..."
avg_overall_exec=0.0
sum_overall_exec=0.0
for ((i = 2; i <= n_tries*2; i*=2)); do
    ./gen_matrix "$i"
    avg_seq_time=0.0
    sum_exec_time=0.0
    for ((k = 1; k < n_avg; k++)); do
        result=$(./$exec_file $input)
        #get the execution time
        exec_time=$(echo "$result" | grep "Gaussian elimination time in CUDA:" | awk '{print $6}')
        sum_exec_time=$(echo "$sum_exec_time + $exec_time" | bc -l)
    done
    avg_seq_time=$(echo "$sum_exec_time / $n_avg" | bc -l)
    sum_overall_exec=$(echo "$sum_overall_exec + $avg_seq_time" | bc -l)
    printf "%d %f\n" "$i" "$avg_seq_time" >> "$stats_file"
done
avg_overall_exec=$(echo "$sum_overall_exec / $n_tries" | bc -l)
printf "Average squential execution time: %f sec. (for %s)\n" "$avg_overall_exec" "$input"
./stats $stats_file "CUDA"