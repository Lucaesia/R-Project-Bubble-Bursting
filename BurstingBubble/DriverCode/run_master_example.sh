#!/bin/bash

# Author: Radu Cimpeanu
# Date: 26/09/2022

# removes previous runs
rm -r Water*/
# Additional velocities or resolution levels can be added below
for V0 in 0.3855; do
	for LEVEL in 7; do

		# Copy all files to renamed folder based on key parameters
		cp -r MasterImpact/ Water-V$V0-Level$LEVEL
		cd Water-V$V0-Level$LEVEL/
		
		# Compile code to create the executable (including visualisation)
		qcc -O2 -w -fopenmp -Wall DropImpact.c -lm -o DropImpact -L$BASILISK/gl -lglutils -lfb_tiny

		# Specify parallelisation features		
		export OMP_NUM_THREADS=4
		
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

		# Run executable ./DropImpact 998.0 1.21 0.998e-3 1.81e-5 0.0722 9.81 0.35e-3 0.3855 6.0 5
		./DropImpact 998.0 1.21 0.998e-3 1.81e-5 0.0722 9.81 0.35e-3 $V0 2.0 $LEVEL
		
		cd ..
	done
done 
