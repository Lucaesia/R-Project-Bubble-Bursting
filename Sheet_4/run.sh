#!/bin/bash

# COMMANDS:
#eval `ssh-agent -s`
#ssh-add ~/.ssh/id_rsa

export OMP_NUM_THREADS=4
rm stability/*
make clean
make

