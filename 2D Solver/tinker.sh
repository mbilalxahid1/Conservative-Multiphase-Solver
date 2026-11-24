#!/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH -t 10:00:00
#SBATCH -p normal_q
#SBATCH -A m2clab

# Go to directory where job was submitted
cd "$SLURM_SUBMIT_DIR"

# Activate Python environment
source ~/venvs/2dsolver/bin/activate

# Give Numba the 32 allocated CPU cores
export NUMBA_NUM_THREADS=$SLURM_CPUS_PER_TASK

# Run the 2D Python solver
srun python -u run_case_2d.py 2>&1 | tee log.out

exit