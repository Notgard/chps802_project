#!/bin/bash
#SBATCH --partition=instant   ### Partition
#SBATCH --job-name=sys_solver ### Job Name
#SBATCH --time=00:10:00     ### WallTime
#SBATCH --nodes=1           ### Number of Nodes
#SBATCH --ntasks-per-node=1 ### Number of tasks (MPI processes)
#SBATCH --cpus-per-task=16 ### Number of threads per task (OMP threads)
#SBATCH --reservation=CHPS
#SBATCH --exclusive
#SBATCH --output solver_job.out
# Print some information about the job
echo "Running on host $(hostname)"
echo "Time is $(date)"
echo "Directory is $(pwd)"
echo "Slurm job ID is $SLURM_JOB_ID"
echo
echo "This job runs on the following machines:"
echo "$SLURM_JOB_NODELIST" | uniq
echo
filename="generated_out.txt"
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
#export OMP_NUM_THREADS=32
make clean && make omp
./main "$filename"