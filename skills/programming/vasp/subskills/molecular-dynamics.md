# Molecular Dynamics Subskill

Detailed guidance for MD simulations in VASP.

## Overview

Molecular dynamics (MD) in VASP integrates Newton's equations of motion to simulate time evolution of atomic positions. Use for:
- Finite temperature properties
- Diffusion coefficients
- Phase transitions
- Liquid structure
- Thermodynamic sampling

## Basic MD Setup

**INCAR:**
```bash
SYSTEM = MD simulation

# Electronic
ENCUT = 400         # Can use lower than static (faster)
EDIFF = 1E-5        # Looser convergence OK
PREC = Normal       # Or Accurate

# MD parameters
IBRION = 0          # Molecular dynamics
NSW = 1000          # Number of MD steps
POTIM = 1.0         # Time step (fs)

# Temperature
TEBEG = 300         # Starting temperature (K)
TEEND = 300         # Final temperature (K)

# Smearing
ISMEAR = 0          # Gaussian
SIGMA = 0.1

# Performance
LREAL = Auto        # Faster for large systems
ALGO = Fast

# Output
LWAVE = .FALSE.
LCHARG = .FALSE.
```

**KPOINTS (can use sparse mesh):**
```bash
Automatic mesh
0
Gamma
  1 1 1             # Γ-point only for large supercells
  0 0 0
```

## Ensembles

### NVE (Microcanonical)

**Constant: N, V, E**
**Use:** Test energy conservation, fundamental dynamics

**INCAR:**
```bash
IBRION = 0
MDALGO = 0          # Standard MD (NVE)
SMASS = -3          # Microcanonical ensemble
NSW = 1000
POTIM = 1.0
TEBEG = 300         # Initial velocities
```

**Notes:**
- Total energy should be conserved
- Check energy drift in OSZICAR
- Good test of time step

### NVT (Canonical)

**Constant: N, V, T**
**Use:** Most common, constant temperature

#### Nose-Hoover Thermostat (Recommended)

**INCAR:**
```bash
IBRION = 0
MDALGO = 2          # Nose-Hoover
SMASS = 0           # Auto-determine (or set manually)
NSW = 1000
POTIM = 1.0
TEBEG = 300
TEEND = 300
```

**SMASS parameter:**
- `SMASS = 0`: Automatic (recommended)
- `SMASS > 0`: Thermostat mass (smaller = stronger coupling)
- Typical: `SMASS = 0.1` to `1.0` for manual control

#### Andersen Thermostat

**INCAR:**
```bash
IBRION = 0
MDALGO = 1          # Andersen
ANDERSEN_PROB = 0.05  # Collision probability
NSW = 1000
POTIM = 1.0
TEBEG = 300
TEEND = 300
```

**Use:** Simple stochastic thermostat, good for equilibration

### NPT (Isothermal-Isobaric)

**Constant: N, P, T**
**Use:** Variable cell, constant pressure

**INCAR:**
```bash
IBRION = 0
MDALGO = 3          # Parrinello-Rahman
SMASS = 0
ISIF = 3            # Allow cell changes
PSTRESS = 0         # External pressure (kBar)
NSW = 1000
POTIM = 1.0
TEBEG = 300
TEEND = 300
```

**PMASS parameter:**
```bash
PMASS = 500         # Barostat mass (optional)
```

**Note:** NPT is computationally expensive and can be unstable.

### NVT with Langevin Thermostat

**INCAR:**
```bash
IBRION = 0
MDALGO = 3          # Langevin (newer VASP versions)
LANGEVIN_GAMMA = 10.0  # Friction coefficient
LANGEVIN_GAMMA_L = 1.0 # Lattice friction (if ISIF=3)
ISIF = 2            # Or 3 for NPT-like
NSW = 1000
POTIM = 1.0
TEBEG = 300
TEEND = 300
```

## Time Step Selection

**Rule of thumb:**
- **POTIM = 1.0 fs:** Standard for most systems
- **POTIM = 0.5 fs:** Tight for light atoms (H)
- **POTIM = 2.0 fs:** Loose for heavy atoms only

**Test time step:**
1. Run short NVE simulation
2. Monitor energy conservation
3. Energy drift < 1 meV/atom/ps is good

**Criteria:**
```
Δt < 1/(10 × ω_max)
```
where ω_max is highest phonon frequency.

**For hydrogen:**
```bash
POTIM = 0.5         # H vibrations are fast
```

**For metals only:**
```bash
POTIM = 2.0         # Can use larger
```

## Initial Velocities

VASP generates random velocities from Maxwell-Boltzmann distribution at TEBEG.

**Set random seed (optional):**
```bash
RANDOM_SEED = 12345  # For reproducibility
```

**Restart from previous MD:**
```bash
# Copy CONTCAR to POSCAR
# VASP reads velocities from POSCAR if present
```

## Temperature Control

### Constant Temperature (NVT)

```bash
TEBEG = 300
TEEND = 300         # Same as TEBEG
MDALGO = 2          # Nose-Hoover
```

### Temperature Ramp (Heating/Cooling)

```bash
TEBEG = 100         # Start at 100 K
TEEND = 1000        # End at 1000 K
NSW = 5000          # Linear ramp over 5000 steps
MDALGO = 2
```

Temperature at step n:
```
T(n) = TEBEG + (TEEND - TEBEG) × n / NSW
```

### Simulated Annealing

**Heat → cool cycle:**

**Phase 1: Heat (0 → 1000 K)**
```bash
TEBEG = 0
TEEND = 1000
NSW = 2000
```

**Phase 2: Hold (1000 K)**
```bash
TEBEG = 1000
TEEND = 1000
NSW = 5000
# Copy CONTCAR to POSCAR
```

**Phase 3: Cool (1000 → 0 K)**
```bash
TEBEG = 1000
TEEND = 0
NSW = 2000
# Copy CONTCAR to POSCAR
```

## Supercell Size

MD requires large supercell to avoid self-interaction through PBC.

**Guidelines:**
- **Minimum:** 3×3×3 unit cell (100-200 atoms)
- **Liquids:** 200-500 atoms
- **Solids:** 100+ atoms
- **Ensure:** Interactions beyond cutoff length

**Check size:**
```python
# Pseudo-code
cutoff = 2 × interaction_range
L_min = cutoff
N_cells = ceil(L_min / lattice_constant)
```

## Equilibration

**Typical protocol:**

1. **Energy minimization:**
   ```bash
   IBRION = 2
   ISIF = 2
   NSW = 100
   ```

2. **Short NVT equilibration (1-5 ps):**
   ```bash
   IBRION = 0
   MDALGO = 1        # Andersen (fast equilibration)
   NSW = 1000
   POTIM = 1.0
   TEBEG = 300
   ```

3. **Production NVT (10-100 ps):**
   ```bash
   MDALGO = 2        # Nose-Hoover (better dynamics)
   NSW = 10000
   ```

**Check equilibration:**
- Monitor temperature (should fluctuate around target)
- Monitor total energy (should stabilize)
- Monitor volume (if NPT)

## Analysis

### Temperature

**From OSZICAR:**
```bash
grep "T=" OSZICAR | awk '{print $3}'
```

**Average temperature:**
```bash
grep "T=" OSZICAR | awk '{sum+=$3; n++} END {print sum/n}'
```

### Energy Conservation (NVE)

**From OSZICAR:**
```bash
grep "E0=" OSZICAR | awk '{print $5}'
```

**Energy drift:**
```
ΔE = (E_final - E_initial) / N_steps
```

**Acceptable:** < 1 meV/atom/ps

### Radial Distribution Function (RDF)

**Post-process XDATCAR with tools:**
- OVITO (GUI)
- MDAnalysis (Python)
- VMD (GUI)

**Example with Python:**
```python
from ase.io import read
import numpy as np

# Read trajectory
traj = read('XDATCAR', index=':')

# Compute RDF (using tool of choice)
```

### Mean Square Displacement (MSD)

**For diffusion coefficient:**
```
D = lim(t→∞) MSD(t) / (6t)
```

**Post-process XDATCAR:**
```python
# Example with MDAnalysis or custom script
```

### Velocity Autocorrelation Function (VACF)

**For vibrational density of states:**
```
VACF(t) = ⟨v(0) · v(t)⟩
DOS(ω) = FT[VACF(t)]
```

## Output Files

| File | Contains |
|------|----------|
| XDATCAR | Trajectory (positions over time) |
| OSZICAR | Energy, temperature per step |
| OUTCAR | Detailed output (forces, stress) |
| vasprun.xml | Complete trajectory + metadata |

**XDATCAR format:**
```bash
System name
   1.0
   lattice vectors
   Elements
   Numbers
Direct configuration=     1
   positions at step 1
Direct configuration=     2
   positions at step 2
...
```

## Complete Examples

### Liquid Water MD (NVT)

**INCAR:**
```bash
SYSTEM = Liquid water MD

# Electronic
ENCUT = 400
EDIFF = 1E-5
PREC = Normal
ISMEAR = 0
SIGMA = 0.1

# MD
IBRION = 0
MDALGO = 2          # Nose-Hoover
SMASS = 0
NSW = 10000         # 10 ps
POTIM = 0.5         # 0.5 fs (hydrogen!)
TEBEG = 300
TEEND = 300
ISIF = 2

# Performance
LREAL = Auto
ALGO = Fast
NCORE = 4

# Output
LWAVE = .FALSE.
LCHARG = .FALSE.
NBLOCK = 100        # Write XDATCAR every 100 steps
```

**POSCAR (supercell):**
```bash
# 64 H2O molecules in cubic box (~12 Å)
# Pre-equilibrated or random
```

**KPOINTS:**
```bash
Gamma-point only
0
Gamma
  1 1 1
  0 0 0
```

### Metal Melting (NPT)

**INCAR:**
```bash
SYSTEM = Al melting

# Electronic
ENCUT = 240
EDIFF = 1E-5
PREC = Normal
ISMEAR = 1
SIGMA = 0.2

# MD NPT
IBRION = 0
MDALGO = 3          # Parrinello-Rahman
SMASS = 0
PMASS = 500
NSW = 5000
POTIM = 2.0
TEBEG = 1000        # Above melting point
TEEND = 1000
ISIF = 3            # Cell relaxation
PSTRESS = 0         # 0 kBar

# Performance
LREAL = Auto
ALGO = Fast

# Output
LWAVE = .FALSE.
LCHARG = .FALSE.
```

### Solid Phonon Sampling (NVT)

**INCAR:**
```bash
SYSTEM = Crystal phonons via MD

# Electronic
ENCUT = 520
EDIFF = 1E-6
PREC = Accurate
ISMEAR = 0
SIGMA = 0.05

# MD
IBRION = 0
MDALGO = 2
SMASS = 0
NSW = 5000
POTIM = 1.0
TEBEG = 300
TEEND = 300
ISIF = 2            # Fixed cell

# Output
LWAVE = .FALSE.
LCHARG = .FALSE.
NBLOCK = 10         # Frequent output for VACF
```

**Post-process:** Compute VACF → FT → phonon DOS

## Advanced: Metadynamics

**Not directly in VASP, use PLUMED plugin:**

```bash
# In INCAR
LUSE_VDW = .TRUE.   # If using PLUMED for vdW
```

**External:** Link VASP with PLUMED for enhanced sampling.

## Advanced: Ab Initio Path Integral MD (PIMD)

**For quantum nuclear effects (not standard VASP).**

Use external tools:
- i-PI (interface for path integrals)
- Couple with VASP as force engine

## Troubleshooting

### Temperature Fluctuates Wildly

**Problem:** T not stable

**Solutions:**
1. **Increase equilibration time:**
   ```bash
   NSW = 5000        # Longer equilibration
   ```

2. **Adjust thermostat:**
   ```bash
   SMASS = 1.0       # Stronger coupling
   ```

3. **Check supercell size:**
   - Too small → large fluctuations
   - Use N > 100 atoms

### Energy Drift (NVE)

**Problem:** Energy not conserved

**Solutions:**
1. **Reduce time step:**
   ```bash
   POTIM = 0.5       # Smaller Δt
   ```

2. **Tighten electronic convergence:**
   ```bash
   EDIFF = 1E-6
   ```

3. **Check for errors:**
   - OSZICAR shows large forces?
   - Structure exploding?

### Simulation "Explodes"

**Problem:** Atoms flying apart

**Solutions:**
1. **Reduce POTIM:**
   ```bash
   POTIM = 0.2       # Much smaller
   ```

2. **Start from relaxed structure:**
   - Don't start MD with bad geometry

3. **Equilibrate first:**
   - Use Andersen thermostat for equilibration
   - Then switch to Nose-Hoover for production

### Cell Becoming Unrealistic (NPT)

**Problem:** NPT gives weird cell shapes

**Solutions:**
1. **Reduce PMASS:**
   ```bash
   PMASS = 100       # Smaller (faster response)
   ```

2. **Use NVT instead:**
   - Fix cell, relax ions only

3. **Check pressure:**
   ```bash
   PSTRESS = 0       # Verify target pressure
   ```

## Performance Optimization

**For large MD simulations:**

```bash
# Parallelization
NCORE = 4           # Cores per orbital
LPLANE = .TRUE.

# Real-space projection (faster)
LREAL = Auto

# Reduce accuracy slightly
PREC = Normal       # Not Accurate
EDIFF = 1E-5        # Looser
ENCUT = 400         # Lower cutoff (still converge!)

# Reduce k-points
# Use Γ-point for large supercells

# Reduce output
NBLOCK = 100        # Write less frequently
LWAVE = .FALSE.
LCHARG = .FALSE.
```

## Best Practices

1. **Always equilibrate first:**
   - 1-5 ps equilibration
   - Check energy/temperature stabilize

2. **Use appropriate supercell:**
   - Large enough to avoid PBC artifacts
   - Typically 100+ atoms

3. **Choose correct time step:**
   - Test energy conservation (NVE)
   - Smaller for light atoms

4. **Monitor during run:**
   ```bash
   tail -f OSZICAR   # Watch temperature, energy
   ```

5. **Save trajectory:**
   - XDATCAR for post-processing
   - Can be large (compress after)

6. **Convergence tests:**
   - Converge ENCUT, k-points like static calculations
   - Test time step (NVE energy drift)

7. **Document setup:**
   - Record ensemble, thermostat, temperature
   - Note equilibration protocol

## Production Workflow

```bash
# 1. Optimize structure
IBRION=2, ISIF=3, NSW=100

# 2. Create supercell (external, e.g., ASE)
# 3×3×3 or larger

# 3. Quick equilibration (Andersen, 1 ps)
IBRION=0, MDALGO=1, NSW=1000, POTIM=1.0

# 4. Production NVT (Nose-Hoover, 10-100 ps)
IBRION=0, MDALGO=2, NSW=10000-100000

# 5. Post-process
# Compute RDF, MSD, VACF, etc.
```

## References

- VASP Wiki: https://vasp.at/wiki/index.php/Molecular_Dynamics
- Thermostats: https://vasp.at/wiki/index.php/Category:Molecular_Dynamics

## See Also

- `relaxation.md` - Optimize before MD
- Main SKILL.md - Full parameter reference
