# Structure Relaxation Subskill

Detailed guidance for structure optimization in VASP.

## Overview

Structure relaxation finds the minimum energy configuration by optimizing atomic positions and/or cell parameters. This is typically the first step in any VASP workflow.

## Types of Relaxation

### 1. Ions Only (Fixed Cell)

**Use when:** You know the lattice parameters (e.g., from experiment) and only want to optimize atomic positions.

**INCAR:**
```bash
IBRION = 2          # Conjugate gradient
ISIF = 2            # Relax ions, keep cell fixed
NSW = 100           # Max ionic steps
EDIFFG = -0.02      # Force convergence (eV/Å)
```

### 2. Ions + Cell Volume (Fixed Shape)

**Use when:** You want to find the equilibrium volume but keep the cell shape.

**INCAR:**
```bash
IBRION = 2
ISIF = 4            # Relax ions + volume
NSW = 100
EDIFFG = -0.02
```

### 3. Full Relaxation (Ions + Cell)

**Use when:** Starting from scratch or studying pressure-induced phase transitions.

**INCAR:**
```bash
IBRION = 2
ISIF = 3            # Relax everything
NSW = 100
EDIFFG = -0.02
```

### 4. Volume-Only Relaxation

**Use when:** Doing equation of state calculations.

**INCAR:**
```bash
IBRION = 2
ISIF = 7            # Volume only, fixed shape
NSW = 100
EDIFFG = -0.02
```

## ISIF Parameter Summary

| ISIF | Ions | Cell Shape | Cell Volume | Stress Tensor | Use Case |
|------|------|------------|-------------|---------------|----------|
| 2 | Yes | No | No | No | Fixed cell, optimize positions |
| 3 | Yes | Yes | Yes | Yes | Full relaxation |
| 4 | Yes | Yes | Yes | Yes | Fixed shape, relax volume |
| 5 | Yes | No | No | Yes | Calculate stress only |
| 6 | Yes | No | Yes | Yes | Change volume, keep shape |
| 7 | Yes | Yes | Yes | Yes | Volume-only relaxation |

## Optimization Algorithms (IBRION)

### IBRION = 2 (Conjugate Gradient)

**Default choice for most systems.**

```bash
IBRION = 2
POTIM = 0.5         # Scaling factor (default)
```

**Pros:**
- Robust and reliable
- Good for most systems
- Handles both ions and cell relaxation

**Cons:**
- Can be slow for difficult landscapes

### IBRION = 1 (RMM-DIIS, Quasi-Newton)

**Better for systems near equilibrium.**

```bash
IBRION = 1
POTIM = 0.5
```

**Pros:**
- Faster convergence near minimum
- Uses force history

**Cons:**
- Can fail if starting far from minimum
- Requires good initial structure

### IBRION = 3 (Damped MD)

**For very rough energy landscapes.**

```bash
IBRION = 3
POTIM = 0.1         # Small time step
SMASS = 0.4         # Damping
```

**Pros:**
- Very robust for bad starting structures
- Can escape local minima

**Cons:**
- Slow
- May need many steps

## Convergence Criteria

### Force-Based (Recommended)

```bash
EDIFFG = -0.02      # Forces < 0.02 eV/Å (standard)
EDIFFG = -0.01      # Forces < 0.01 eV/Å (tight)
EDIFFG = -0.05      # Forces < 0.05 eV/Å (loose, testing)
```

**When to use tighter criteria:**
- Phonon calculations (use -0.01 or tighter)
- Elastic constants
- High-accuracy thermodynamics

### Energy-Based (Less Common)

```bash
EDIFFG = 0.001      # Energy change < 0.001 eV (positive = energy)
```

**When to use:**
- Very flat energy surfaces
- When forces are unreliable

## Complete Relaxation Template

### Metals

```bash
# System
SYSTEM = Cu bulk relaxation

# Electronic
ENCUT = 520
EDIFF = 1E-6
PREC = Accurate
ISMEAR = 1          # Methfessel-Paxton
SIGMA = 0.2
ALGO = Fast
LREAL = Auto

# Ionic relaxation
IBRION = 2          # Conjugate gradient
ISIF = 3            # Full relaxation
NSW = 100
EDIFFG = -0.02

# Output
LWAVE = .FALSE.
LCHARG = .FALSE.
```

**KPOINTS:**
```bash
Automatic mesh
0
Gamma
  8 8 8
  0 0 0
```

### Semiconductors/Insulators

```bash
# System
SYSTEM = Si bulk relaxation

# Electronic
ENCUT = 520
EDIFF = 1E-6
PREC = Accurate
ISMEAR = 0          # Gaussian
SIGMA = 0.05
ALGO = Fast
LREAL = Auto

# Ionic relaxation
IBRION = 2
ISIF = 3
NSW = 100
EDIFFG = -0.02

# Output
LWAVE = .FALSE.
LCHARG = .FALSE.
```

**KPOINTS:**
```bash
Automatic mesh
0
Gamma
  6 6 6
  0 0 0
```

### Molecules (Gas Phase)

```bash
# System
SYSTEM = H2O molecule

# Electronic
ENCUT = 400
EDIFF = 1E-6
PREC = Accurate
ISMEAR = 0
SIGMA = 0.01        # Small smearing
ALGO = Normal

# Ionic relaxation
IBRION = 2
ISIF = 2            # Don't relax cell!
NSW = 100
EDIFFG = -0.01      # Tight for molecules

# Output
LWAVE = .FALSE.
LCHARG = .FALSE.
```

**POSCAR:**
```bash
H2O molecule
1.0
  12.0  0.0  0.0    # Large box (no PBC interaction)
  0.0  12.0  0.0
  0.0  0.0  12.0
H O
  2 1
Cartesian
  ...positions...
```

**KPOINTS:**
```bash
Gamma point only
0
Gamma
  1 1 1
  0 0 0
```

### Surfaces and Slabs

```bash
# System
SYSTEM = Cu(111) surface

# Electronic
ENCUT = 520
EDIFF = 1E-6
PREC = Accurate
ISMEAR = 1
SIGMA = 0.2
ALGO = Fast
LREAL = Auto

# Ionic relaxation
IBRION = 2
ISIF = 2            # Don't relax cell (slab)
NSW = 100
EDIFFG = -0.02

# Dipole correction for slabs
LDIPOL = .TRUE.
IDIPOL = 3          # Direction perpendicular to surface

# Selective dynamics (freeze bottom layers)
# Add to POSCAR: Selective dynamics, then T T T or F F F

# Output
LWAVE = .FALSE.
LCHARG = .FALSE.
```

**KPOINTS:**
```bash
Automatic mesh
0
Gamma
  8 8 1             # Dense in xy, sparse in z
  0 0 0
```

## Monitoring Convergence

### Check OSZICAR During Run

```bash
# Watch convergence in real-time
tail -f OSZICAR
```

Output shows:
```
  1 F= -.10234E+02 E0= -.10234E+02  d E =-.102340E+02  mag=     0.0000
  2 F= -.10456E+02 E0= -.10456E+02  d E =-.222000E+00  mag=     0.0000
```

- Column 1: Ionic step
- F: Free energy
- d E: Energy change
- Convergence when d E and forces small

### Check Final Convergence

```bash
# Check if relaxation converged
grep "reached required accuracy" OUTCAR

# Check final forces
grep "TOTAL-FORCE" OUTCAR -A 50

# Check final stress (if ISIF=3)
grep "in kB" OUTCAR | tail -1
```

## Troubleshooting

### Relaxation Not Converging

**Problem:** NSW steps exhausted, forces still large

**Solutions:**
1. **Increase NSW:**
   ```bash
   NSW = 200
   ```

2. **Check initial structure:**
   ```bash
   # Look for overlapping atoms
   # Verify reasonable bond lengths
   ```

3. **Try different optimizer:**
   ```bash
   IBRION = 1    # Or 3 for very rough surfaces
   POTIM = 0.3   # Smaller step size
   ```

4. **Tighten electronic convergence:**
   ```bash
   EDIFF = 1E-7
   ALGO = All
   ```

### Oscillating Forces

**Problem:** Forces oscillate, never converge

**Solutions:**
1. **Reduce step size:**
   ```bash
   POTIM = 0.2
   ```

2. **Switch to RMM-DIIS:**
   ```bash
   IBRION = 1
   ```

3. **Check for imaginary phonons** (unstable structure)

### Cell Becoming Unrealistic

**Problem:** ISIF=3 gives distorted cell

**Solutions:**
1. **Check symmetry:**
   ```bash
   SYMPREC = 1E-5  # Default
   # Or disable: ISYM = 0
   ```

2. **Start with ISIF=2, then ISIF=3:**
   - First relax ions with fixed cell
   - Then relax cell with optimized ions

3. **Check for stress:**
   ```bash
   # May be under pressure
   # Use PSTRESS to apply external pressure
   ```

### ZBRENT Error

**Problem:** "ZBRENT: fatal error in bracketing"

**Solutions:**
```bash
POTIM = 0.1         # Much smaller step
# Or
IBRION = 3          # Use damped MD
SMASS = 0.4
```

## Selective Dynamics

For surfaces, often want to freeze bottom layers.

**POSCAR with selective dynamics:**
```bash
Cu(111) surface
1.0
  ... lattice vectors ...
Cu
  12
Selective dynamics     # Add this line
Direct
  0.0  0.0  0.0   F F F    # Fixed
  0.5  0.5  0.0   F F F    # Fixed
  0.0  0.0  0.1   T T T    # Relaxed
  ...
```

**INCAR:**
```bash
# No special flags needed, VASP reads from POSCAR
IBRION = 2
ISIF = 2            # Don't relax cell
EDIFFG = -0.02
```

## Constrained Relaxation

### Fix Lattice Parameters

```bash
# Keep a, b, c fixed but allow angles to change (rare)
# Use external tools or fix in POSCAR

# Common: ISIF=2 (standard fixed cell)
```

### Fix Certain Atoms

Use selective dynamics in POSCAR (see above).

### Apply External Pressure

```bash
PSTRESS = 100       # External pressure in kBar
ISIF = 3            # Allow volume to change
```

## Multi-Step Relaxation Strategy

For difficult systems, use hierarchical relaxation:

### Step 1: Rough Relaxation
```bash
ENCUT = 400         # Lower cutoff
EDIFF = 1E-4        # Loose electronic
EDIFFG = -0.05      # Loose ionic
ISMEAR = 0
SIGMA = 0.1
# Coarse k-points
```

### Step 2: Fine Relaxation
```bash
ENCUT = 520         # Higher cutoff
EDIFF = 1E-6        # Tight electronic
EDIFFG = -0.02      # Standard ionic
ISMEAR = 1
SIGMA = 0.2
# Dense k-points
```

### Step 3: Final Static
```bash
IBRION = -1         # No relaxation
ISMEAR = -5         # Tetrahedron
# Very dense k-points
```

## Equation of State Workflow

To find equilibrium volume and bulk modulus:

```python
# Pseudo-code workflow
volumes = [0.95, 0.97, 0.99, 1.00, 1.01, 1.03, 1.05]  # Relative to initial

for V in volumes:
    # Scale POSCAR
    scale = V ** (1/3)

    # Run VASP with ISIF=2 (fixed cell, relax ions)
    # Or ISIF=4 (relax ions + volume from this starting point)

    # Collect E(V)

# Fit E(V) to Birch-Murnaghan equation
# Extract V0, B0, B0'
```

**INCAR for each volume:**
```bash
IBRION = 2
ISIF = 2            # Or 4 for automatic volume search
EDIFFG = -0.01      # Tight for EOS
```

## Best Practices

1. **Always start with clean relaxation:**
   - Don't reuse WAVECAR from very different structure
   - Delete WAVECAR if in doubt

2. **Check forces in all directions:**
   - Verify forces are small for ALL atoms
   - Not just RMS force

3. **Verify cell shape for ISIF=3:**
   - Compare initial and final lattice vectors
   - Check if physically reasonable

4. **Save CONTCAR:**
   - This is your relaxed structure
   - Rename to POSCAR for next calculation

5. **Document convergence:**
   - Record final energy, forces, stress
   - Note number of ionic steps

6. **Test k-points and ENCUT:**
   - Converge on unrelaxed structure first
   - Then use converged parameters for relaxation

## See Also

- `convergence.md` - Systematic convergence testing
- `electronic-structure.md` - Static calculations after relaxation
- Main SKILL.md - Full parameter reference
