# Convergence Testing Subskill

Systematic strategies for converging VASP calculations.

## Overview

**All VASP calculations must be converged** with respect to computational parameters. Unconverged results are meaningless.

**Key parameters to converge:**
1. Energy cutoff (ENCUT)
2. k-point mesh (KPOINTS)
3. Electronic convergence (EDIFF)
4. Smearing width (SIGMA)
5. Supercell size (for surfaces, defects)

## General Convergence Strategy

**Order of testing:**

1. **ENCUT first** (with moderate k-points)
2. **k-points second** (with converged ENCUT)
3. **Other parameters** (SIGMA, etc.)

**Convergence criterion:**
- **Energy:** ΔE < 1 meV/atom between successive values
- **Forces:** ΔF < 0.01 eV/Å
- **Stress:** Δσ < 0.1 kBar

**For property X:**
- **Converged when:** ΔX < acceptable error

## Energy Cutoff (ENCUT) Convergence

### Basic Test

**Method:**
1. Fix k-point mesh (moderate density, e.g., 4×4×4)
2. Vary ENCUT
3. Plot E vs ENCUT

**Test sequence:**
```bash
ENCUT: 300, 350, 400, 450, 500, 550, 600 eV
```

**INCAR template:**
```bash
SYSTEM = ENCUT convergence test

# Test parameter
ENCUT = 400         # Vary this

# Fixed k-points
# Use KPOINTS file with fixed mesh

# Electronic
EDIFF = 1E-6
PREC = Accurate
ISMEAR = 1
SIGMA = 0.2

# Static calculation (fast)
IBRION = -1
NSW = 0

# Output
LWAVE = .FALSE.
LCHARG = .FALSE.
```

**Analysis:**
```python
# Pseudo-code
ENCUT_values = [300, 350, 400, 450, 500, 550, 600]
E_values = []  # From OUTCAR

for i in range(len(ENCUT_values)-1):
    dE = abs(E_values[i+1] - E_values[i]) / N_atoms
    print(f"ENCUT {ENCUT_values[i]} → {ENCUT_values[i+1]}: ΔE = {dE*1000:.2f} meV/atom")

# Choose ENCUT where dE < 1 meV/atom
```

**Typical results:**
- **Below 1.3 × ENMAX:** Energy not converged
- **1.3-1.5 × ENMAX:** Usually converged
- **Above 1.5 × ENMAX:** Converged, but expensive

**Recommendation:**
```bash
ENCUT = 520 eV      # Standard for most PAW potentials
```

### Forces Convergence

**For relaxations, forces may need higher ENCUT.**

**Test:**
1. Relax structure at each ENCUT
2. Compare final forces

```bash
ENCUT: 500, 550, 600, 650 eV
# Check max force in OUTCAR
```

**Criterion:** Max force changes < 0.01 eV/Å

## k-Point Convergence

### Basic Test

**Method:**
1. Fix ENCUT (use converged value)
2. Vary k-point mesh
3. Plot E vs k-mesh density

**Test sequence:**
```bash
KPOINTS: 2×2×2, 4×4×4, 6×6×6, 8×8×8, 10×10×10, 12×12×12
```

**KPOINTS template:**
```bash
Automatic mesh
0
Gamma              # Or Monkhorst-Pack
  4 4 4            # Vary this
  0 0 0
```

**INCAR:**
```bash
SYSTEM = k-point convergence test

ENCUT = 520        # Fixed (converged)
EDIFF = 1E-6
PREC = Accurate
ISMEAR = 1
SIGMA = 0.2
IBRION = -1
NSW = 0
```

**Analysis:**
```python
# Pseudo-code
k_mesh = ['2×2×2', '4×4×4', '6×6×6', '8×8×8', '10×10×10']
E_values = []

for i in range(len(k_mesh)-1):
    dE = abs(E_values[i+1] - E_values[i]) / N_atoms
    print(f"{k_mesh[i]} → {k_mesh[i+1]}: ΔE = {dE*1000:.2f} meV/atom")

# Choose k-mesh where dE < 1 meV/atom
```

### k-Point Density Rule

**Alternative to manual testing:**

Use k-point density rule:
```
k-density = N_k / (2π/a)
```

**Recommendation:**
- **Testing:** 20-30 k-points/Å⁻¹
- **Relaxation:** 30-40 k-points/Å⁻¹
- **Static/DOS:** 40-50 k-points/Å⁻¹

**Example:**
```python
# For cubic lattice a = 4 Å
# Target: 40 k-points/Å⁻¹

k_per_direction = 40 × (4 / (2π))
k ≈ 8

# Use 8×8×8 mesh
```

**Or use KSPACING:**
```bash
KSPACING = 0.5     # Å⁻¹ (smaller = denser)
# Typical: 0.3-0.5 for production
```

### Special Cases

**Metals (Fermi surface):**
- Need denser k-mesh
- Use ISMEAR = 1 or 2
- Test k-mesh carefully

**Insulators:**
- Can use sparser k-mesh
- ISMEAR = 0 or -5

**Surfaces/slabs:**
```bash
# Dense in xy, sparse in z
KPOINTS
0
Gamma
  8 8 1             # Example
  0 0 0
```

**2D materials:**
```bash
# Dense in xy, Γ-only in z
KPOINTS
0
Gamma
  12 12 1
  0 0 0
```

## Electronic Convergence (EDIFF)

**EDIFF controls SCF convergence.**

**Test:**
```bash
EDIFF: 1E-4, 1E-5, 1E-6, 1E-7, 1E-8
```

**Recommendations:**
| Calculation Type | EDIFF |
|------------------|-------|
| Testing | 1E-4 |
| Standard | 1E-6 |
| Forces (tight) | 1E-7 |
| Phonons | 1E-8 |
| GW | 1E-8 |

**Check:**
```bash
grep "EDIFF" OUTCAR
# Verify "reached required accuracy"
```

## Smearing Width (SIGMA)

**SIGMA affects DOS smearing near Fermi level.**

### Test Metals

**Method:**
1. Fix ENCUT, k-mesh
2. Vary SIGMA
3. Check entropy term

**Test sequence:**
```bash
SIGMA: 0.05, 0.1, 0.2, 0.3, 0.4
```

**Check entropy:**
```bash
grep "entropy T\*S" OUTCAR
```

**Criterion:**
```
T×S < 1 meV/atom    (should be small)
```

**Typical:**
```bash
ISMEAR = 1
SIGMA = 0.2         # Metals
```

### Test Semiconductors

**For semiconductors/insulators:**
```bash
ISMEAR = 0
SIGMA = 0.05        # Small
```

**Or use tetrahedron for static:**
```bash
ISMEAR = -5         # No SIGMA needed
```

## Complete Convergence Workflow

### Automated Testing Script

**Python script to automate:**

```python
#!/usr/bin/env python
"""
ENCUT and k-point convergence test.
"""
import os
import numpy as np
from ase.io import read

def test_encut():
    """Test ENCUT convergence."""
    encut_values = [300, 350, 400, 450, 500, 550, 600]
    energies = []

    for encut in encut_values:
        # Modify INCAR
        with open('INCAR', 'r') as f:
            lines = f.readlines()

        with open('INCAR', 'w') as f:
            for line in lines:
                if 'ENCUT' in line:
                    f.write(f'ENCUT = {encut}\n')
                else:
                    f.write(line)

        # Run VASP
        os.system('vasp_std > vasp.out')

        # Extract energy
        with open('OSZICAR', 'r') as f:
            energy = float(f.readlines()[-1].split()[2])
        energies.append(energy)

        print(f"ENCUT = {encut} eV: E = {energy:.6f} eV")

    # Analyze
    atoms = read('POSCAR')
    n_atoms = len(atoms)

    print("\nConvergence:")
    for i in range(len(encut_values)-1):
        dE = abs(energies[i+1] - energies[i]) / n_atoms * 1000
        print(f"{encut_values[i]} → {encut_values[i+1]}: ΔE = {dE:.2f} meV/atom")

def test_kpoints():
    """Test k-point convergence."""
    k_meshes = [2, 4, 6, 8, 10, 12]
    energies = []

    for k in k_meshes:
        # Modify KPOINTS
        with open('KPOINTS', 'w') as f:
            f.write(f"Automatic mesh\n0\nGamma\n  {k} {k} {k}\n  0 0 0\n")

        # Run VASP
        os.system('vasp_std > vasp.out')

        # Extract energy
        with open('OSZICAR', 'r') as f:
            energy = float(f.readlines()[-1].split()[2])
        energies.append(energy)

        print(f"k-mesh = {k}×{k}×{k}: E = {energy:.6f} eV")

    # Analyze
    atoms = read('POSCAR')
    n_atoms = len(atoms)

    print("\nConvergence:")
    for i in range(len(k_meshes)-1):
        dE = abs(energies[i+1] - energies[i]) / n_atoms * 1000
        print(f"{k_meshes[i]}³ → {k_meshes[i+1]}³: ΔE = {dE:.2f} meV/atom")

if __name__ == '__main__':
    print("=" * 60)
    print("ENCUT Convergence Test")
    print("=" * 60)
    test_encut()

    print("\n" + "=" * 60)
    print("k-Point Convergence Test")
    print("=" * 60)
    test_kpoints()
```

### Step-by-Step Manual Workflow

**Step 1: ENCUT convergence**

Create directories:
```bash
mkdir convergence_tests
cd convergence_tests
mkdir encut_300 encut_350 encut_400 encut_450 encut_500 encut_550 encut_600
```

For each directory:
```bash
cd encut_300
# Copy POSCAR, POTCAR, KPOINTS (fixed 4×4×4)
# Create INCAR with ENCUT=300
vasp_std > vasp.out
grep "energy  without" OUTCAR | tail -1
cd ..
```

Analyze:
```bash
# Compare energies, find convergence
```

**Step 2: k-point convergence**

```bash
mkdir kpoints_tests
cd kpoints_tests
mkdir k222 k444 k666 k888 k101010 k121212

# Use converged ENCUT from step 1
# Vary k-mesh in each directory
```

**Step 3: Document results**

```bash
# Create convergence_report.txt
echo "ENCUT Convergence:" > convergence_report.txt
echo "Converged at ENCUT = 520 eV (ΔE < 1 meV/atom)" >> convergence_report.txt
echo "" >> convergence_report.txt
echo "k-Point Convergence:" >> convergence_report.txt
echo "Converged at 8×8×8 (ΔE < 1 meV/atom)" >> convergence_report.txt
```

## Convergence for Specific Properties

### Total Energy

**Criterion:** ΔE < 1 meV/atom

**Recommended:**
```bash
ENCUT = 520 eV
k-mesh: 8×8×8 (for ~4 Å lattice constant)
EDIFF = 1E-6
```

### Formation Energy

**Criterion:** ΔE < 5 meV/atom (less stringent)

**Ensure:**
- Same ENCUT for all phases
- Same POTCARs
- Same k-point density

### Forces

**Criterion:** ΔF < 0.01 eV/Å

**Recommended:**
```bash
ENCUT = 550 eV      # Higher than energy
EDIFF = 1E-7        # Tighter
k-mesh: 8×8×8+
```

### Stress/Pressure

**Criterion:** Δσ < 0.1 kBar

**Recommended:**
```bash
ENCUT = 600 eV      # Higher
PREC = Accurate
EDIFF = 1E-6
k-mesh: Very dense
```

### Phonons

**Criterion:** Δω < 1 cm⁻¹

**Recommended:**
```bash
ENCUT = 550-600 eV
EDIFF = 1E-8        # Very tight
k-mesh: Dense
Tight relaxation: EDIFFG = -0.01
```

### Band Gap

**Criterion:** ΔE_gap < 0.01 eV

**Recommended:**
```bash
ENCUT = 520 eV
k-mesh: Very dense (12×12×12+)
EDIFF = 1E-6
ISMEAR = -5         # Tetrahedron
```

### Elastic Constants

**Criterion:** ΔC_ij < 1 GPa

**Recommended:**
```bash
ENCUT = 600 eV
EDIFF = 1E-7
k-mesh: Very dense
PREC = Accurate
```

## Convergence for Surfaces and Slabs

### Slab Thickness

**Test:**
```bash
# Number of layers: 3, 5, 7, 9, 11
# Check surface energy convergence
```

**Criterion:**
```
ΔE_surf < 0.01 eV/Å²
```

**Typical:** 5-7 layers (with bottom layers fixed)

### Vacuum Thickness

**Test:**
```bash
# Vacuum: 10 Å, 12 Å, 15 Å, 18 Å, 20 Å
# Check total energy
```

**Criterion:**
```
ΔE < 1 meV/atom
```

**Typical:** 15 Å vacuum

## Convergence Documentation Template

```markdown
# Convergence Tests for [System Name]

## Parameters

- **System:** Cu FCC bulk
- **Date:** 2024-01-15
- **VASP version:** 6.3.0
- **POTCARs:** PAW_PBE Cu 06Sep2000

## ENCUT Convergence

Tested: 300, 350, 400, 450, 500, 550, 600 eV
Fixed k-mesh: 4×4×4

| ENCUT (eV) | Energy (eV) | ΔE (meV/atom) |
|------------|-------------|---------------|
| 300 | -3.72534 | - |
| 350 | -3.72891 | 3.57 |
| 400 | -3.73102 | 2.11 |
| 450 | -3.73201 | 0.99 |
| 500 | -3.73245 | 0.44 |
| 550 | -3.73261 | 0.16 |
| 600 | -3.73267 | 0.06 |

**Converged at ENCUT = 520 eV** (ΔE < 1 meV/atom at 550 eV)

## k-Point Convergence

Tested: 2×2×2, 4×4×4, 6×6×6, 8×8×8, 10×10×10
Fixed ENCUT: 520 eV

| k-mesh | Energy (eV) | ΔE (meV/atom) |
|--------|-------------|---------------|
| 2³ | -3.71234 | - |
| 4³ | -3.72891 | 16.57 |
| 6³ | -3.73201 | 3.10 |
| 8³ | -3.73278 | 0.77 |
| 10³ | -3.73291 | 0.13 |

**Converged at 8×8×8** (ΔE < 1 meV/atom)

## Final Parameters

```bash
ENCUT = 520 eV
k-mesh: 8×8×8 (Γ-centered)
EDIFF = 1E-6
ISMEAR = 1
SIGMA = 0.2
```

## Notes

- Used for all subsequent calculations
- Formation energy calculations will use same parameters
```

## Troubleshooting

### Energy Not Converging with ENCUT

**Problem:** ΔE still large even at ENCUT = 600 eV

**Solutions:**
1. **Check POTCAR:**
   - Look at ENMAX in POTCAR
   - Ensure ENCUT > 1.3 × ENMAX

2. **Try PREC = Accurate:**
   ```bash
   PREC = Accurate
   ```

3. **Check for hard pseudopotentials (_h):**
   - These have higher ENMAX
   - May need ENCUT = 700-800 eV

### k-Points Not Converging

**Problem:** Energy oscillates with k-mesh

**Solutions:**
1. **Try Γ-centered instead of MP:**
   ```bash
   # In KPOINTS
   Gamma              # Instead of Monkhorst-Pack
   ```

2. **Increase SIGMA (metals):**
   ```bash
   SIGMA = 0.3        # Smoother DOS
   ```

3. **Check for Fermi surface:**
   - Metals may need very dense mesh
   - Try ISMEAR = 2 (MP order 2)

### Different Properties Need Different Convergence

**Problem:** Energy converged but forces not

**Solution:**
- This is normal
- Use tighter parameters for forces
- Document separate convergence for each property

## Best Practices

1. **Always test convergence:**
   - Never assume default parameters are sufficient
   - Test for your specific system

2. **Test in order:**
   - ENCUT first (cheaper)
   - k-points second (with converged ENCUT)

3. **Document everything:**
   - Save convergence tests
   - Create convergence report
   - Use for publication supplementary info

4. **Use consistent parameters:**
   - Same ENCUT, k-density for all related calculations
   - Same POTCARs

5. **Be conservative:**
   - Choose parameters slightly beyond convergence
   - Cost of over-convergence < cost of unconverged results

6. **Automate when possible:**
   - Use scripts to test systematically
   - Reduces human error

7. **Plot results:**
   - Visual inspection helps identify convergence
   - Look for plateaus in E vs parameter plots

## Quick Reference

| Property | ENCUT | k-mesh | EDIFF | Notes |
|----------|-------|--------|-------|-------|
| Energy | 520 eV | 8³ | 1E-6 | Standard |
| Forces | 550 eV | 8³ | 1E-7 | Higher ENCUT |
| Stress | 600 eV | 12³ | 1E-6 | Very dense k |
| Phonons | 550 eV | 8³ | 1E-8 | Tight all |
| Band gap | 520 eV | 12³ | 1E-6 | Dense k, ISMEAR=-5 |
| Elastic | 600 eV | 12³ | 1E-7 | Highest requirements |

## References

- VASP Wiki: https://vasp.at/wiki/index.php/Category:Accuracy_and_Convergence

## See Also

- All other subskills for calculation-specific convergence
- Main SKILL.md - Full parameter reference
