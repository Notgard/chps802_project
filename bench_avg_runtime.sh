#!/usr/bin/bash
if [[ $# -lt 2 || $# -gt 2 ]] ; then
	echo 'Wrong arguments provided, need filename and number of tries (nb threads to use in total for benchmarking)'
	exit 1
fi
exec_file="main"
input=$1
n_tries=$2
n=0
n_avg=3
#make compile main solving C program
make clean && make
#set variables
stats_file="./benchmark/stats-seq-$$.txt"
echo -e "Starting benchmark for file $input..."
avg_overall_exec=0.0
sum_overall_exec=0.0
for ((i = 1; i < n_tries; i++)); do
    avg_seq_time=0.0
    sum_exec_time=0.0
    for ((k = 1; k < n_avg; k++)); do
        result=$(./$exec_file "$input")
        #get the execution time
        exec_time=$(echo "$result" | grep "Total solver runtime:" | awk '{print $4}')
        sum_exec_time=$(echo "$sum_exec_time + $exec_time" | bc -l)
    done
    avg_seq_time=$(echo "$sum_exec_time / $n_avg" | bc -l)
    sum_overall_exec=$(echo "$sum_overall_exec + $avg_seq_time" | bc -l)
    printf "%d %.3f\n" "$i" "$avg_seq_time" >> "$stats_file"
done
avg_overall_exec=$(echo "$sum_overall_exec / $n_tries" | bc -l)
printf "Average squential execution time: %.3f sec. (for %s)\n" "$avg_overall_exec" "$input"
echo -e "Finished running sequential tests, now running parallel benchmark..."
#make compile with omp target directive
make clean && make omp
#parallel tests
stats_file_p="./benchmark/stats-omp-$$.txt"
threads=1
avg_overall_exec=0.0
sum_overall_exec=0.0
for ((i = 1; i < n_tries; i++)); do
    avg_omp_time=0.0
    sum_exec_time=0.0
    for ((k = 1; k < n_avg; k++)); do
        result=$(OMP_NUM_THREADS=$threads ./$exec_file "$input")
        #get the execution time
        exec_time=$(echo "$result" | grep "Total solver runtime:" | awk '{print $4}')
        sum_exec_time=$(echo "$sum_exec_time + $exec_time" | bc -l)
    done
    avg_omp_time=$(echo "$sum_exec_time / $n_avg" | bc -l)
    sum_overall_exec=$(echo "$sum_overall_exec + $avg_seq_time" | bc -l)
    printf "%d %.3f\n" "$threads" "$avg_omp_time" >> "$stats_file_p"
    ((threads++))
done
avg_overall_exec=$(echo "$sum_overall_exec / $n_tries" | bc -l)
printf "Average parallel execution time: %.3f sec. (for %s)\n" "$avg_overall_exec" "$input"
echo -e "Finished benchmarking."
#create benchmark graph
./stats $stats_file "Sequentielle" $stats_file_p "OpenMP"