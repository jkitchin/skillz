# Phonons Subskill

Detailed guidance for phonon calculations in VASP using DFPT and finite differences.

## Overview

Phonon calculations determine vibrational modes and thermodynamic properties. VASP offers two methods:

1. **DFPT (Density Functional Perturbation Theory):** Direct calculation via IBRION=5,6,7,8
2. **Finite differences:** Numerical derivatives, use with Phonopy

## Method 1: DFPT (Built-in VASP)

DFPT calculates force constants directly using perturbation theory.

### Basic DFPT Phonons (IBRION=6)

**INCAR:**
```bash
SYSTEM = Phonon DFPT

# Electronic (TIGHT convergence required!)
ENCUT = 550         # Higher than relaxation
EDIFF = 1E-8        # Very tight!
PREC = Accurate
LREAL = .FALSE.     # Reciprocal space (more accurate)

# DFPT phonons
IBRION = 6          # Finite differences for phonons
NFREE = 2           # Central differences
POTIM = 0.015       # Displacement magnitude (Å)
NSW = 1             # Only one displacement per direction

# Smearing (small)
ISMEAR = 0
SIGMA = 0.01

# Symmetry
ISYM = 2            # Use symmetry

# Output
LWAVE = .FALSE.
LCHARG = .FALSE.
```

**KPOINTS (dense):**
```bash
Automatic mesh
0
Gamma
  8 8 8             # Denser than relaxation
  0 0 0
```

**Prerequisites:**
- **Relaxed structure** (forces < 0.01 eV/Å)
- Small primitive cell (not supercell)

**Output:**
- OUTCAR contains force constants, eigenvalues, eigenvectors
- vasprun.xml has complete phonon data

### DFPT for Specific q-Points (IBRION=7,8)

**IBRION=7:** Phonons at specific q-points

**INCAR:**
```bash
IBRION = 7
NFREE = 2
POTIM = 0.015
# ... same as IBRION=6
```

**KPOINTS (for q-point path):**
```bash
Phonon dispersion
10
Line-mode
Reciprocal
  0.0  0.0  0.0   !Γ
  0.5  0.0  0.5   !X

  0.5  0.0  0.5   !X
  0.5  0.5  0.5   !L
```

**IBRION=8:** Phonons on regular q-mesh

### Parameters

| Parameter | Purpose | Typical Value |
|-----------|---------|---------------|
| IBRION | Method | 6 (DFPT), 7 (q-path), 8 (q-mesh) |
| NFREE | Difference order | 2 (central), 4 (higher order) |
| POTIM | Displacement | 0.015 Å (test 0.01-0.02) |
| EDIFF | Convergence | 1E-8 eV (tight!) |

## Method 2: Finite Differences with Phonopy

**More flexible and commonly used.**

### Workflow

#### Step 1: Relax Structure

**INCAR:**
```bash
IBRION = 2
ISIF = 3
EDIFFG = -0.01      # Tight forces!
EDIFF = 1E-7
ENCUT = 550
```

**Must have tight convergence** or phonons will have errors.

#### Step 2: Create Supercell with Displacements

Use Phonopy to generate displaced structures:

```bash
# Create supercell (e.g., 2×2×2)
phonopy -d --dim="2 2 2"
```

This creates:
- `SPOSCAR`: Supercell
- `POSCAR-001`, `POSCAR-002`, ...: Displaced structures

**Parameters:**
- `--dim="2 2 2"`: Supercell size (adjust for system)

#### Step 3: Calculate Forces for Each Displacement

For each `POSCAR-XXX`:

```bash
# Copy POSCAR-001 to POSCAR
cp POSCAR-001 POSCAR
```

**INCAR (force calculation):**
```bash
SYSTEM = Phonopy force calculation

# Electronic (tight)
ENCUT = 550
EDIFF = 1E-8
PREC = Accurate
LREAL = .FALSE.

# Static calculation
IBRION = -1         # No relaxation!
NSW = 0

# Smearing
ISMEAR = 0
SIGMA = 0.01

# Symmetry
ISYM = 0            # Off for forces

# Output
LWAVE = .FALSE.
LCHARG = .FALSE.
```

**KPOINTS:**
```bash
Automatic mesh
0
Gamma
  4 4 4             # Sparse for supercell
  0 0 0
```

**Run VASP** for each displacement, save `vasprun.xml` as `vasprun-001.xml`, etc.

#### Step 4: Collect Forces

```bash
# Extract forces from all calculations
phonopy -f vasprun-{001..010}.xml
```

This creates `FORCE_SETS`.

#### Step 5: Post-Processing

**Calculate phonon dispersion:**
```bash
# Create band.conf
cat > band.conf <<EOF
DIM = 2 2 2
BAND = 0 0 0  0.5 0 0.5,  0.5 0 0.5  0.5 0.5 0.5
EOF

# Run phonopy
phonopy --dim="2 2 2" -p band.conf
```

**Calculate DOS:**
```bash
# Create mesh.conf
cat > mesh.conf <<EOF
DIM = 2 2 2
MP = 20 20 20
EOF

phonopy --dim="2 2 2" -p mesh.conf
```

**Calculate thermal properties:**
```bash
phonopy --dim="2 2 2" -t -p mesh.conf
```

Produces `thermal_properties.yaml` with:
- Free energy
- Entropy
- Heat capacity

## Supercell Size Selection

**Rule of thumb:**
```
Supercell size × k-mesh ≈ constant
```

**Examples:**
| Primitive | Supercell | k-mesh |
|-----------|-----------|--------|
| 1×1×1 | 1×1×1 | 8×8×8 |
| 1×1×1 | 2×2×2 | 4×4×4 |
| 1×1×1 | 3×3×3 | 3×3×3 |
| 1×1×1 | 4×4×4 | 2×2×2 |

**Choose supercell:**
- **Too small:** Short-range interactions only, inaccurate
- **Too large:** Very expensive
- **Typical:** 2×2×2 to 4×4×4

**For phonon DOS, use larger supercell.**

## Complete Phonopy Example

### Cu FCC Phonons

**Step 1: Relax**

POSCAR (primitive cell):
```bash
Cu primitive
1.0
  0.0  1.805  1.805
  1.805  0.0  1.805
  1.805  1.805  0.0
Cu
  1
Direct
  0.0 0.0 0.0
```

INCAR:
```bash
ENCUT = 550
EDIFF = 1E-7
IBRION = 2
ISIF = 3
EDIFFG = -0.01
```

Run VASP → get relaxed structure

**Step 2: Generate displacements**

```bash
phonopy -d --dim="3 3 3"
```

**Step 3: Calculate forces**

For each POSCAR-XXX:
```bash
cp POSCAR-001 POSCAR
# Run VASP (IBRION=-1, EDIFF=1E-8)
mv vasprun.xml vasprun-001.xml
```

**Step 4: Process**

```bash
phonopy -f vasprun-*.xml
```

**Step 5: Band structure**

band.conf:
```
DIM = 3 3 3
PRIMITIVE_AXES = AUTO
BAND = 0 0 0  0.5 0 0.5  0.5 0.5 0.5  0 0 0  0.5 0.5 0
BAND_LABELS = $\Gamma$ X L $\Gamma$ K
BAND_POINTS = 51
```

```bash
phonopy --dim="3 3 3" -p band.conf
```

Produces `band.yaml` and plot.

**Step 6: Thermal properties**

mesh.conf:
```
DIM = 3 3 3
MP = 20 20 20
TMIN = 0
TMAX = 1000
TSTEP = 10
```

```bash
phonopy --dim="3 3 3" -t -p mesh.conf
```

## Analyzing Results

### Check for Imaginary Frequencies

**Negative eigenvalues = imaginary frequencies**

**From Phonopy:**
```bash
phonopy --dim="2 2 2" -p band.conf
# Look at band.yaml for negative frequencies
```

**Causes:**
1. **Structure not relaxed:** Tighten EDIFFG
2. **Unstable structure:** System wants to distort
3. **Numerical noise:** Increase ENCUT, tighter EDIFF

**If imaginary modes at Γ:**
- Structure is unstable
- May indicate phase transition
- Relax with lower symmetry

### Thermal Properties

**From thermal_properties.yaml:**
- **Free energy F(T)**
- **Entropy S(T)**
- **Heat capacity C_v(T)**
- **Internal energy E(T)**

**Zero-point energy (ZPE):**
```
E_ZPE = F(T=0)
```

## Advanced: Raman and IR

### Raman-Active Modes

**INCAR:**
```bash
IBRION = 6
LEPSILON = .TRUE.   # Calculate dielectric tensor
LRAMAN = .TRUE.     # Raman tensors
NFREE = 2
POTIM = 0.015
```

**Requirements:**
- LEPSILON for Born effective charges
- LRAMAN for Raman susceptibility

**Output:**
- OUTCAR contains Raman tensors

### IR-Active Modes

**INCAR:**
```bash
IBRION = 6
LEPSILON = .TRUE.   # Born charges and dielectric
NFREE = 2
POTIM = 0.015
```

**Born effective charges in OUTCAR:**
```
BORN EFFECTIVE CHARGES
```

**Use for:**
- LO/TO splitting
- Infrared intensities

## Advanced: Non-Analytic Correction (NAC)

For polar materials, include LO/TO splitting.

### With Phonopy

**Step 1: Calculate Born charges**

**INCAR:**
```bash
IBRION = 6
LEPSILON = .TRUE.
NFREE = 2
```

**Step 2: Extract from OUTCAR**

```bash
phonopy-vasp-born > BORN
```

Creates BORN file with Born charges and dielectric tensor.

**Step 3: Use in Phonopy**

```bash
phonopy --dim="2 2 2" --nac -p band.conf
```

**band.conf must include:**
```
NAC = .TRUE.
```

## Convergence Testing

### Test POTIM

**DFPT:** Vary POTIM to ensure converged force constants.

```bash
# Test POTIM = 0.01, 0.015, 0.02
# Phonon frequencies should be converged
```

### Test Supercell Size

**Phonopy:** Increase supercell until frequencies converge.

```bash
# Test 2×2×2, 3×3×3, 4×4×4
# Check high-frequency modes converge
```

### Test ENCUT

```bash
# Test ENCUT = 500, 550, 600 eV
# Phonon frequencies should change < 1 cm⁻¹
```

### Test k-Points

```bash
# Test 2×2×2, 3×3×3, 4×4×4 k-mesh (for supercell)
# Converge to < 1 cm⁻¹
```

## Troubleshooting

### Large Imaginary Frequencies

**Problem:** Many negative modes, not just acoustic

**Solutions:**
1. **Tighten relaxation:**
   ```bash
   EDIFFG = -0.005     # Very tight
   EDIFF = 1E-8
   ```

2. **Check for instability:**
   - System may want to distort
   - Try lower symmetry

3. **Increase supercell:**
   - 2×2×2 → 3×3×3

### Small Imaginary Acoustic Modes at Γ

**Problem:** Acoustic modes at Γ should be zero, but small negative

**Cause:** Numerical noise (normal)

**Solution:**
- If frequencies < 10 cm⁻¹, ignore (artifact)
- Set to zero in post-processing

### Forces Not Consistent

**Problem:** Phonopy complains about forces

**Solutions:**
1. **Check ISYM:**
   ```bash
   ISYM = 0            # Off for force calculations
   ```

2. **Verify EDIFF:**
   ```bash
   EDIFF = 1E-8        # Very tight
   ```

3. **Check all calculations used same POTCAR**

## Performance Optimization

### DFPT

```bash
# DFPT is memory-intensive
# Use fewer cores per node
NCORE = 1

# Ensure LREAL = .FALSE. (accuracy)
# Use dense k-mesh
```

### Phonopy

```bash
# Force calculations are independent
# Run in parallel (job array)

# For supercell, use sparse k-mesh
KSPACING = 0.5      # Or manual sparse mesh

# Can use LREAL = Auto for large supercells
```

## Phonopy Configuration Files

### band.conf (Dispersion)

```
DIM = 2 2 2
PRIMITIVE_AXES = AUTO
BAND = 0 0 0  0.5 0 0.5,  0.5 0 0.5  0.5 0.5 0.5,  0.5 0.5 0.5  0 0 0
BAND_LABELS = $\Gamma$ X L $\Gamma$
BAND_POINTS = 51
BAND_CONNECTION = .TRUE.
NAC = .TRUE.          # If polar
```

### mesh.conf (DOS, Thermal)

```
DIM = 2 2 2
MP = 20 20 20
PRIMITIVE_AXES = AUTO
TMIN = 0
TMAX = 1000
TSTEP = 10
NAC = .TRUE.
```

### pdos.conf (Projected DOS)

```
DIM = 2 2 2
MP = 20 20 20
PRIMITIVE_AXES = AUTO
PDOS = 1 2, 3 4 5    # Atom indices to project
```

## Visualization

**Plot with Phonopy:**
```bash
phonopy -p band.yaml
```

**Export for external plotting:**
```bash
phonopy --writefc -p band.conf
# Creates band.yaml, parse with Python
```

**Python example:**
```python
import phonopy
from phonopy import load

# Load phonopy calculation
phonon = load(supercell_matrix=[[2,0,0],[0,2,0],[0,0,2]],
              primitive_matrix='auto',
              unitcell_filename='POSCAR')

# Plot band structure
phonon.auto_band_structure(plot=True)
```

## Comparison: DFPT vs Phonopy

| Aspect | DFPT (IBRION=6) | Phonopy |
|--------|-----------------|---------|
| Accuracy | High | High |
| Flexibility | Limited | High |
| Cost | Moderate | High (many calculations) |
| Supercell | Primitive only | Any size |
| NAC | Built-in (LEPSILON) | Manual (BORN file) |
| Parallelization | Limited | Embarrassingly parallel |
| Ease | Single calculation | Multi-step |
| Recommendation | Small systems | Most cases |

**Recommendation:** Use Phonopy for most cases (more flexible, easier to converge).

## Best Practices

1. **Relax structure first:**
   - Forces < 0.01 eV/Å (tight!)
   - Use EDIFFG = -0.01 or tighter

2. **Use tight convergence:**
   ```bash
   EDIFF = 1E-8
   ENCUT = 550
   ```

3. **Test convergence:**
   - Supercell size
   - k-points
   - ENCUT
   - POTIM (for DFPT)

4. **Check acoustic sum rule:**
   - Acoustic modes should be zero at Γ
   - Small imaginary modes (< 10 cm⁻¹) are OK

5. **Include NAC for polar systems:**
   - Calculate Born charges (LEPSILON)
   - Use in Phonopy (--nac)

6. **Document everything:**
   - Supercell size
   - Convergence parameters
   - Method (DFPT or Phonopy)

## References

- VASP DFPT: https://vasp.at/wiki/index.php/Phonons_from_finite_differences
- Phonopy: https://phonopy.github.io/phonopy/
- NAC: https://phonopy.github.io/phonopy/formulation.html#non-analytical-term-correction

## See Also

- `relaxation.md` - Structure optimization before phonons
- `convergence.md` - Systematic convergence testing
- Main SKILL.md - Full parameter reference
