# Electronic Structure Subskill

Detailed guidance for band structure and density of states calculations in VASP.

## Overview

Electronic structure calculations reveal the quantum mechanical properties of materials: band gaps, band dispersion, density of states, orbital character, and more.

## Workflow Overview

Typical electronic structure workflow:

1. **Structure optimization** (see `relaxation.md`)
2. **Static SCF calculation** (high accuracy)
3. **Band structure** (non-SCF along high-symmetry path)
4. **Density of states** (non-SCF with dense k-mesh)

## Static Calculation (Single-Point Energy)

After relaxation, calculate accurate energy with tight settings.

**INCAR:**
```bash
SYSTEM = Static calculation

# Electronic
ENCUT = 520
EDIFF = 1E-6        # Or 1E-7 for high accuracy
PREC = Accurate

# No ionic movement
IBRION = -1
NSW = 0

# Accurate DOS
ISMEAR = -5         # Tetrahedron method
# No SIGMA needed for ISMEAR=-5

# High precision
ALGO = Normal
LREAL = .FALSE.     # Reciprocal space (more accurate)

# Output
LWAVE = .FALSE.
LCHARG = .FALSE.
```

**KPOINTS (dense mesh):**
```bash
Automatic mesh
0
Gamma
  12 12 12          # Dense for accurate energy
  0 0 0
```

**Use static calculations for:**
- Final accurate total energy
- Formation energies
- Equation of state
- Reference for band structure

## Band Structure Calculation

Calculate electronic bands along high-symmetry path in reciprocal space.

### Two-Step Procedure

#### Step 1: Self-Consistent Calculation

Generate charge density with regular k-mesh.

**INCAR:**
```bash
SYSTEM = SCF for band structure

# Electronic
ENCUT = 520
EDIFF = 1E-6
PREC = Accurate

# No ionic movement
IBRION = -1
NSW = 0

# Charge generation
ICHARG = 2          # From atomic charge
LCHARG = .TRUE.     # Write CHGCAR

# Smearing
ISMEAR = 0          # Gaussian
SIGMA = 0.05

# Output
LWAVE = .FALSE.
```

**KPOINTS (regular mesh):**
```bash
Automatic mesh
0
Gamma
  8 8 8
  0 0 0
```

**Run VASP** → produces CHGCAR

#### Step 2: Non-Self-Consistent Band Structure

Read charge density, calculate eigenvalues along k-path.

**INCAR:**
```bash
SYSTEM = Band structure

# Electronic
ENCUT = 520
EDIFF = 1E-6
PREC = Accurate

# No ionic movement
IBRION = -1
NSW = 0

# Non-self-consistent
ICHARG = 11         # Read CHGCAR, no update

# Projections for orbital character
LORBIT = 11         # Write PROCAR

# Smearing (not critical for bands)
ISMEAR = 0
SIGMA = 0.05

# Number of bands
NBANDS = 50         # Include empty bands

# Output
LWAVE = .FALSE.
LCHARG = .FALSE.
```

**KPOINTS (line mode):**
```bash
k-points for band structure
20                  # Points between high-symmetry points
Line-mode
Reciprocal
  0.0  0.0  0.0   !Γ
  0.5  0.0  0.5   !X

  0.5  0.0  0.5   !X
  0.5  0.25 0.75  !W

  0.5  0.25 0.75  !W
  0.375 0.375 0.75 !K

  0.375 0.375 0.75 !K
  0.0  0.0  0.0   !Γ

  0.0  0.0  0.0   !Γ
  0.5  0.5  0.5   !L
```

**Run VASP** → produces EIGENVAL, PROCAR

### High-Symmetry Paths by Crystal System

#### FCC (e.g., Cu, Al, Ni)
```
Γ-X-W-K-Γ-L-U-W-L-K|U-X
```

**KPOINTS:**
```bash
k-path FCC
20
Line-mode
Reciprocal
  0.0  0.0  0.0   !Γ
  0.5  0.0  0.5   !X

  0.5  0.0  0.5   !X
  0.5  0.25 0.75  !W

  0.5  0.25 0.75  !W
  0.375 0.375 0.75 !K

  0.375 0.375 0.75 !K
  0.0  0.0  0.0   !Γ

  0.0  0.0  0.0   !Γ
  0.5  0.5  0.5   !L
```

#### BCC (e.g., Fe, Cr, W)
```
Γ-H-N-Γ-P-H|P-N
```

**KPOINTS:**
```bash
k-path BCC
20
Line-mode
Reciprocal
  0.0  0.0  0.0   !Γ
  0.5 -0.5  0.5   !H

  0.5 -0.5  0.5   !H
  0.0  0.0  0.5   !N

  0.0  0.0  0.5   !N
  0.0  0.0  0.0   !Γ

  0.0  0.0  0.0   !Γ
  0.25 0.25 0.25  !P
```

#### Simple Cubic (e.g., SC)
```
Γ-X-M-Γ-R-X|M-R
```

**KPOINTS:**
```bash
k-path SC
20
Line-mode
Reciprocal
  0.0  0.0  0.0   !Γ
  0.5  0.0  0.0   !X

  0.5  0.0  0.0   !X
  0.5  0.5  0.0   !M

  0.5  0.5  0.0   !M
  0.0  0.0  0.0   !Γ

  0.0  0.0  0.0   !Γ
  0.5  0.5  0.5   !R
```

#### Hexagonal (e.g., graphene, hBN)
```
Γ-M-K-Γ
```

**KPOINTS:**
```bash
k-path Hexagonal
20
Line-mode
Reciprocal
  0.0  0.0  0.0      !Γ
  0.5  0.0  0.0      !M

  0.5  0.0  0.0      !M
  0.3333 0.3333 0.0  !K

  0.3333 0.3333 0.0  !K
  0.0  0.0  0.0      !Γ
```

**Tool for automatic paths:**
Use SeeK-path (https://www.materialscloud.org/work/tools/seekpath) or pymatgen to generate high-symmetry paths automatically.

## Density of States (DOS)

Calculate DOS with dense k-mesh.

**INCAR:**
```bash
SYSTEM = DOS calculation

# Electronic
ENCUT = 520
EDIFF = 1E-6
PREC = Accurate

# No ionic movement
IBRION = -1
NSW = 0

# DOS parameters
ISMEAR = -5         # Tetrahedron (best for DOS)
LORBIT = 11         # Projected DOS (orbital character)
NEDOS = 3000        # DOS resolution

# Optional: read CHGCAR from SCF
ICHARG = 11

# Output
LWAVE = .FALSE.
LCHARG = .FALSE.
```

**KPOINTS (very dense):**
```bash
Automatic mesh
0
Gamma
  16 16 16          # Very dense for smooth DOS
  0 0 0
```

**Run VASP** → produces DOSCAR

### Projected DOS (PDOS)

**LORBIT parameter:**

| LORBIT | Output | Use |
|--------|--------|-----|
| 10 | PROCAR, DOSCAR (no RWIGS needed) | Projected DOS |
| 11 | PROCAR, DOSCAR (no RWIGS needed) | Projected DOS (recommended) |
| 12 | PROCAR, DOSCAR (no phase) | Projected DOS |

**Recommended:**
```bash
LORBIT = 11         # Works without RWIGS
```

**For specific site projections:**
```bash
LORBIT = 11
# VASP automatically decomposes by atom and orbital
```

DOSCAR contains:
- Total DOS
- Integrated DOS
- Projected DOS by atom and orbital (s, p, d, f)

## Band Gap Determination

### Direct vs Indirect Gap

**From OUTCAR:**
```bash
grep "E-fermi" OUTCAR         # Fermi energy
grep "HOMO-LUMO" OUTCAR       # For molecules
```

**From band structure:**
- **Direct gap:** VBM and CBM at same k-point
- **Indirect gap:** VBM and CBM at different k-points

**Accurately determine gap:**

1. Find VBM (valence band maximum):
   ```bash
   # Look at EIGENVAL or DOSCAR
   # Highest occupied state
   ```

2. Find CBM (conduction band minimum):
   ```bash
   # Lowest unoccupied state
   ```

3. Calculate gap:
   ```
   E_gap = E_CBM - E_VBM
   ```

**Note:** DFT systematically underestimates band gaps. For accurate gaps, use:
- Hybrid functionals (HSE06)
- GW approximation
- DFT+U for correlated systems

## Advanced: Hybrid Functionals (HSE06)

For accurate band gaps.

**INCAR:**
```bash
SYSTEM = HSE06 calculation

# Hybrid functional
LHFCALC = .TRUE.    # Activate hybrid
HFSCREEN = 0.2      # HSE screening parameter
AEXX = 0.25         # Exact exchange fraction (25%)

# Electronic
ENCUT = 520
EDIFF = 1E-6
PREC = Accurate

# Convergence (HSE harder to converge)
ALGO = All          # Or Damped
TIME = 0.4          # Damping parameter
NELM = 100

# Reduced k-mesh for HF
NKRED = 2           # Reduce k-mesh by factor 2 for HF part
# Or NKREDX, NKREDY, NKREDZ individually

# No ionic movement
IBRION = -1
NSW = 0

# Smearing
ISMEAR = 0
SIGMA = 0.05

# Performance
PRECFOCK = Fast     # Fast FFT for HF
LHFCALC = .TRUE.
```

**KPOINTS (can use coarser mesh):**
```bash
Automatic mesh
0
Gamma
  4 4 4             # HSE very expensive, use sparse mesh
  0 0 0
```

**Workflow:**
1. Relax with PBE (cheaper)
2. Static with HSE on relaxed structure
3. Optional: Band structure with HSE (very expensive)

**Cost:** HSE is ~10-100× more expensive than PBE.

## Advanced: GW Approximation

Most accurate band gaps (quasi-particle energies).

### G₀W₀ Workflow

#### Step 1: DFT (PBE)

**INCAR:**
```bash
SYSTEM = DFT for GW

# Electronic
ENCUT = 520
EDIFF = 1E-8        # Tight!
PREC = Accurate
ALGO = Exact        # Or Normal

# Many bands needed
NBANDS = 200        # Rule: 2-4× occupied bands

# No ionic movement
IBRION = -1
NSW = 0

# Smearing
ISMEAR = 0
SIGMA = 0.01

# Optics
LOPTICS = .TRUE.    # Frequency-dependent dielectric function

# Output
LWAVE = .TRUE.      # Need WAVECAR for GW
```

#### Step 2: G₀W₀

**INCAR:**
```bash
SYSTEM = G0W0

# Read from DFT
ALGO = GW0          # Or EVGW for self-consistent

# GW parameters
NOMEGA = 50         # Frequency grid points
NBANDS = 200        # Same as step 1

# Response function
LOPTICS = .FALSE.
ENCUTGW = 200       # Cutoff for response (lower than ENCUT)

# No ionic movement
IBRION = -1
NSW = 0
```

**Run on DFT WAVECAR** → produces quasi-particle energies in OUTCAR

**Cost:** GW is ~100-1000× more expensive than PBE.

## Orbital Character Analysis

Use PROCAR to analyze orbital contributions.

**INCAR:**
```bash
LORBIT = 11         # Write orbital projections
```

**PROCAR format:**
```
k-point 1:
band  1:  energy = -10.5 eV
  ion   s      py     pz     px    dxy    dyz    dz2    dxz    dx2
   1  0.05   0.10   0.10   0.10   0.01   0.01   0.01   0.01   0.01
   2  0.05   0.12   0.12   0.12   0.00   0.00   0.00   0.00   0.00
  tot 0.10   0.22   0.22   0.22   0.01   0.01   0.01   0.01   0.01
```

**Analysis:**
- Which orbitals contribute to VBM, CBM?
- Bonding vs antibonding character
- Hybridization

**Visualization:**
Use pymatgen, sumo, or vaspkit to plot "fat bands" (band structure with orbital weights).

## Effective Mass

Calculate effective mass from band curvature.

**Method 1: Finite differences**

1. Calculate bands along Γ-X direction
2. Fit parabola near band extremum:
   ```
   E(k) = E₀ + ℏ²k²/(2m*)
   ```
3. Extract m* from curvature

**Method 2: VASP automatic (experimental)**
```bash
IBRION = 5          # Finite differences
LEFG = .TRUE.       # Calculate effective mass tensor
```

## Fat Bands (Orbital-Weighted Bands)

Combine band structure with orbital character.

**Workflow:**
1. Calculate bands with `LORBIT = 11`
2. Use analysis tool (pymatgen, sumo) to plot
3. Color/thickness indicates orbital contribution

**Example with pymatgen:**
```python
from pymatgen.electronic_structure.plotter import BSPlotter
from pymatgen.io.vasp import BSVasprun

# Read band structure
vasprun = BSVasprun("vasprun.xml")
bs = vasprun.get_band_structure(line_mode=True)

# Plot fat bands
plotter = BSPlotter(bs)
plotter.get_plot(vbm_cbm_marker=True)
```

## Spin-Polarized Calculations

For magnetic systems.

**INCAR:**
```bash
ISPIN = 2           # Spin-polarized

# Initial magnetic moments
MAGMOM = 5.0 5.0 0.0    # One value per atom

# For DOS
LORBIT = 11
```

**Output:**
- DOSCAR contains spin-up and spin-down DOS
- EIGENVAL shows both spin channels
- Magnetic moments in OUTCAR

## Troubleshooting

### Band Gap Too Small

**Problem:** DFT underestimates gap

**Solutions:**
1. Use HSE06 (typically adds 1-2 eV)
2. Use GW (most accurate)
3. Use DFT+U for correlated systems
4. Check k-point convergence
5. Verify structure is relaxed

### Band Structure Looks Wrong

**Problem:** Discontinuities, strange features

**Solutions:**
1. **Check KPOINTS path:**
   - Must be continuous
   - Duplicate points at joints (see examples above)

2. **Increase number of bands:**
   ```bash
   NBANDS = 100      # Default may be too few
   ```

3. **Check CHGCAR is from same structure:**
   - CHGCAR from step 1 must match POSCAR in step 2

### DOS Not Smooth

**Problem:** Noisy DOS

**Solutions:**
1. **Increase k-points:**
   ```bash
   # Use denser mesh, e.g., 20×20×20
   ```

2. **Increase NEDOS:**
   ```bash
   NEDOS = 5000      # Higher resolution
   ```

3. **Use tetrahedron method:**
   ```bash
   ISMEAR = -5       # Best for DOS
   ```

## Visualization Tools

**Recommended tools:**
- **pymatgen** (Python): Comprehensive materials analysis
- **sumo** (Python): Pretty band structure and DOS plots
- **p4vasp** (GUI): Interactive VASP output viewer
- **vaspkit** (CLI): Post-processing toolkit
- **VESTA** (GUI): Structure visualization

**Example: Plot bands with sumo:**
```bash
# After VASP calculation
sumo-bandplot
sumo-dosplot
```

## Best Practices

1. **Always relax first:**
   - Band structure is meaningless for unrelaxed structure

2. **Use consistent k-mesh:**
   - Same mesh for SCF and DOS
   - Line mode only for band structure

3. **Check Fermi level:**
   - Verify it's in gap (insulators) or at DOS peak (metals)

4. **Converge empty bands:**
   - For band structure, include enough NBANDS to see conduction bands

5. **Document everything:**
   - Record k-path used
   - Save KPOINTS, EIGENVAL, DOSCAR

6. **Compare to experiment/literature:**
   - Verify gap, dispersion are reasonable
   - DFT gaps typically 30-50% too small

## References

- VASP Wiki: https://vasp.at/wiki/index.php/Band-structure_calculation
- High-symmetry paths: https://www.materialscloud.org/work/tools/seekpath
- pymatgen docs: https://pymatgen.org/

## See Also

- `relaxation.md` - Structure optimization before electronic structure
- `advanced-functionals.md` - HSE06, GW details
- Main SKILL.md - Full parameter reference
