# Band Structure Example

Two-step workflow for calculating electronic band structure of Cu FCC.

## Files

### Step 1: Self-Consistent Calculation
- `INCAR.scf` - SCF calculation parameters
- `POSCAR` - Atomic structure
- `KPOINTS.scf` - Regular k-mesh for charge density
- `POTCAR.info` - POTCAR information

### Step 2: Non-Self-Consistent Band Structure
- `INCAR.bands` - Non-SCF band parameters
- `KPOINTS.bands` - High-symmetry k-path
- Uses CHGCAR from step 1

## Workflow

### Step 1: Self-Consistent Calculation

Generate charge density:

```bash
# Copy files
cp INCAR.scf INCAR
cp KPOINTS.scf KPOINTS
cp $VASP_PP/potpaw_PBE/Cu/POTCAR .

# Run VASP
mpirun -np 16 vasp_std > vasp.scf.out

# Save outputs
mv OUTCAR OUTCAR.scf
mv OSZICAR OSZICAR.scf
# Keep CHGCAR for step 2
```

### Step 2: Band Structure

Calculate bands along high-symmetry path:

```bash
# Copy files
cp INCAR.bands INCAR
cp KPOINTS.bands KPOINTS
# POSCAR, POTCAR, CHGCAR already present

# Run VASP
mpirun -np 16 vasp_std > vasp.bands.out

# Save outputs
mv OUTCAR OUTCAR.bands
mv EIGENVAL EIGENVAL.bands
```

### Plotting

Extract and plot bands:

```bash
# Option 1: Use pymatgen
python plot_bands.py

# Option 2: Use p4vasp
p4v vasprun.xml

# Option 3: Use sumo
sumo-bandplot --code vasp
```

## Expected Results

Cu is a metal, so:
- No band gap (Fermi level crosses bands)
- d-bands below Fermi level (~2-5 eV below)
- s-p bands cross Fermi level
- Typical bandwidth: ~10 eV

## k-Path Used

FCC high-symmetry path: Γ-X-W-K-Γ-L

```
Γ = (0, 0, 0)
X = (0.5, 0, 0.5)
W = (0.5, 0.25, 0.75)
K = (0.375, 0.375, 0.75)
L = (0.5, 0.5, 0.5)
```

## Notes

- Step 1 generates charge density (CHGCAR)
- Step 2 uses CHGCAR, no self-consistency (ICHARG=11)
- Typical runtime: Step 1: 2-3 min, Step 2: 1-2 min (16 cores)
