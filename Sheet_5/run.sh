#!/bin/bash

# OpenMP run
export OMP_NUM_THREADS=4
make droplet_impact.tst

# MPI run
#C99='mpicc -std=c99' qcc -Wall -O2 -L$(BASILISK)/gl -lglutils -lfb_osmesa -lGLU -lOSMesa -lm -D_MPI=1 droplet_impact.c -o droplet_impact 

#mpirun -np 8 ./droplet_impact