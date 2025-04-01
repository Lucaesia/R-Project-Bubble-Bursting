#!/bin/bash

# Author: Lyllian Chanerley
# Based on code by: Radu Cimpeanu
# Date: 26/09/2025

# removes previous runs
#rm -r Water-IC-R0.002-Level9
rm -r Water-IC-R0.0011-Level10
# Additional velocities or resolution levels can be added below
I="IC"
for R in 0.0011 ; do #0.002 0.0018 0.0016 0.0014 0.0012 0.001 0.0008 0.0006 0.0004 0.0002 ; do
	for LEVEL in 10; do
    for I in IC; do

      # Copy all files to renamed folder based on key parameters
      cp -r MasterImpact/ Water-$I-R$R-Level$LEVEL/
      cp MasterImpact/DropImpact.c Water-$I-R$R-Level$LEVEL/DropImpact.c
      cd Water-$I-R$R-Level$LEVEL/
      
      echo "COMPILING R: $R  LEVEL: $LEVEL"
      # Compile code to create the executable (including visualisation)
      qcc -O2 -w -fopenmp -Wall DropImpact.c -lm -o DropImpact -L$BASILISK/gl -lglutils -lfb_tiny
      echo "RUNNING R: $R  LEVEL: $LEVEL"
      # Specify parallelisation features		
      export OMP_NUM_THREADS=2
      
      # parameters:
      # 1. rho liquid (dimensional)
      # 2. rho gas (dimensional)
      # 3. mu liquid (dimensional)
      # 4. mu gas (dimensional)
      # 5. sigma (surface tension coefficient, dimensional)
      # 6. g acceleration (dimensional)
      # 7. drop radius (dimensional)
      # 8. initial drop velocity (dimensional)
      # 9. simulation end time (dimensionless, tailored to contact and bounce duration)
      # 10. max level
      VAL=0
      THETA=2.3
      #0.002 0.0018 0.0016 0.0014 0.0012 0.001 0.0008 0.0006 0.0004 0.0002
      [ "$R" == "0.0018" ] && THETA=2.6
      [ "$R" == "0.0016" ] && THETA=3.0
      [ "$R" == "0.0014" ] && THETA=3.5
      [ "$R" == "0.0012" ] && THETA=4.2
      [ "$R" == "0.001" ] && THETA=5.1
      [ "$R" == "0.0008" ] && THETA=6.5
      [ "$R" == "0.0006" ] && THETA=8.9
      [ "$R" == "0.0004" ] && THETA=13.5
      [ "$R" == "0.0002" ] && THETA=27

      [ "$I" == "IC" ] && VAL=1

      
      
      # Run executable ./DropImpact 998.0 1.21 0.998e-3 1.81e-5 0.0722 9.81 0.35e-3 0.3855 6.0 5
      ./DropImpact 998.0 1.21 0.998e-3 1.81e-5 0.0722 9.81 $R 1 0.5 $LEVEL $THETA $VAL
      
      cd ..
    done
	done
done 
