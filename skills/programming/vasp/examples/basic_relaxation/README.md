# Basic Relaxation Example

Complete VASP input files for a simple bulk Cu relaxation.

## Files

- `INCAR` - Calculation parameters
- `POSCAR` - Atomic structure (Cu FCC)
- `KPOINTS` - k-point mesh
- `POTCAR.info` - POTCAR information (actual file not included)
- `run.sh` - Example run script

## How to Run

1. Copy your Cu POTCAR:
   ```bash
   cp $VASP_PP/potpaw_PBE/Cu/POTCAR .
   ```

2. Run VASP:
   ```bash
   mpirun -np 16 vasp_std > vasp.out
   ```

3. Check convergence:
   ```bash
   grep "reached required accuracy" OUTCAR
   tail OSZICAR
   ```

4. View relaxed structure:
   ```bash
   # CONTCAR contains relaxed structure
   ```

## Expected Results

- Final energy: ~-3.73 eV/atom (PBE)
- Lattice constant: ~3.63 Å (PBE overestimates, exp: 3.61 Å)
- Forces: < 0.02 eV/Å (as specified in EDIFFG)

## Notes

- This is a simple bulk metal relaxation
- Uses standard parameters: ENCUT=520, 8×8×8 k-mesh
- Relaxes both ions and cell (ISIF=3)
- Takes ~1-2 minutes on 16 cores
