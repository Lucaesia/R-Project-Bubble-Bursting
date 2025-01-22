#!/bin/bash
# Avon submission script for PX457 2024 Assignment 4
# Runs on the "PX457" queue which is a 2-node partition

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=3700mb
#SBATCH --time=00:02:00

# Some comments on resource requests
# ==================================
#
# There are two ways to request resources on avon:
# One can either use nodes=q ntasks-per-node=p and be 
# allocated p cores on q physical nodes, or one
# can just ask for ntasks=X, and get X processors
# spread over any physical nodes. However, in this
# assignment, we have a reservation of 2 nodes, so
# as long as you use X<=48 they will usually all be on
# one node, but it is possible they will be spread out
# over multiple nodes.
# 
# NB: if the system gets busy, you may experience
# some variability of timings for larger X since your
# job might be split over both nodes. If that appears
# to be the case, and you need accurate timings, try
# the nodes=1:ntasks-per-node=X option instead.

# Change into the directory where I submitted the job
cd $SLURM_SUBMIT_DIR 

export OMP_NUM_THREADS=4
make clean
make stability.tst
exit 0
