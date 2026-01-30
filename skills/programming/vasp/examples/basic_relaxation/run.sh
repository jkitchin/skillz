#!/bin/bash
# Basic VASP relaxation run script

# Check for required files
if [ ! -f POTCAR ]; then
    echo "ERROR: POTCAR not found!"
    echo "Copy Cu POTCAR: cp \$VASP_PP/potpaw_PBE/Cu/POTCAR ."
    exit 1
fi

# Number of cores (adjust for your system)
NCORES=16

# Run VASP
echo "Running VASP relaxation..."
mpirun -np $NCORES vasp_std > vasp.out 2>&1

# Check completion
if grep -q "reached required accuracy" OUTCAR; then
    echo "SUCCESS: Relaxation converged!"
    echo ""
    echo "Results:"
    echo "--------"
    grep "energy  without" OUTCAR | tail -1
    echo ""
    echo "Final lattice constant:"
    grep "VOLUME and BASIS" OUTCAR -A 5 | tail -5
else
    echo "WARNING: Relaxation may not have converged."
    echo "Check OSZICAR for details."
fi

echo ""
echo "Relaxed structure saved in CONTCAR"
