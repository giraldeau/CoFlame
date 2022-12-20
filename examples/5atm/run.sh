#!/bin/bash -x

# All commands assume the current directory contains the case file.
SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
cd $SCRIPT_DIR

# Compile the chemistry kinetic database
# Input files chem.inp therm.dat
# Output file: chem.bin
chem

# Compile the tran database tran.bin
# Input files chem.bin tran.dat
# Output file: tran.bin
tran

# Run the case
# Input files: meshing.dat input.dat initialvalues.dat chem.bin tran.bin 
# Output file:
mpirun -np 16 CoFlame

