# Convergence Tests Example

Automated scripts for testing ENCUT and k-point convergence.

## Files

- `test_encut.py` - ENCUT convergence test script
- `test_kpoints.py` - k-point convergence test script
- `base_INCAR` - Template INCAR
- `base_POSCAR` - Cu FCC structure
- `base_POTCAR.info` - POTCAR information

## Usage

### ENCUT Convergence

```bash
# Setup
cp $VASP_PP/potpaw_PBE/Cu/POTCAR .

# Run test
python test_encut.py

# Results will be in encut_convergence/
# Creates plot: encut_convergence.png
```

### k-Point Convergence

```bash
# Setup (requires converged ENCUT)
cp $VASP_PP/potpaw_PBE/Cu/POTCAR .

# Run test
python test_kpoints.py

# Results will be in kpoint_convergence/
# Creates plot: kpoint_convergence.png
```

## What the Scripts Do

### test_encut.py

1. Creates directories for each ENCUT value (300, 350, 400, 450, 500, 550, 600 eV)
2. Copies input files with modified ENCUT
3. Runs VASP in each directory (or submits jobs)
4. Collects energies from OSZICAR
5. Plots E vs ENCUT and ΔE vs ENCUT
6. Reports convergence threshold

### test_kpoints.py

1. Creates directories for each k-mesh (2³, 4³, 6³, 8³, 10³, 12³)
2. Generates KPOINTS files
3. Runs VASP in each directory
4. Collects energies
5. Plots E vs k-density
6. Reports convergence threshold

## Customization

Edit the scripts to change:
- Test values (ENCUT range, k-mesh sizes)
- Convergence criteria (default: 1 meV/atom)
- VASP executable (default: vasp_std)
- Parallel execution (MPI ranks)

## Expected Results

For Cu FCC:
- **ENCUT:** Converges at ~500-520 eV
- **k-points:** Converges at 8×8×8 for 4 Å lattice

## Notes

- Run ENCUT test first, use result for k-point test
- Tests use static calculations (IBRION=-1) for speed
- Can adapt for your own system by changing POSCAR and POTCAR
- Runtime: ~30-60 minutes for both tests (depends on system)
