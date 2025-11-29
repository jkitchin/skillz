# VASP Quick Reference

One-page parameter lookup for common VASP calculations.

## Common INCAR Templates

### Relaxation (Metals)
```bash
SYSTEM = Metal relaxation
ENCUT = 520
PREC = Accurate
IBRION = 2          # Conjugate gradient
ISIF = 3            # Relax cell + ions
NSW = 100
EDIFFG = -0.02      # Force convergence (eV/Å)
EDIFF = 1E-6
ISMEAR = 1          # Methfessel-Paxton
SIGMA = 0.2
LREAL = Auto
ALGO = Fast
```

### Relaxation (Semiconductors/Insulators)
```bash
SYSTEM = Semiconductor relaxation
ENCUT = 520
PREC = Accurate
IBRION = 2
ISIF = 3
NSW = 100
EDIFFG = -0.02
EDIFF = 1E-6
ISMEAR = 0          # Gaussian
SIGMA = 0.05
LREAL = Auto
ALGO = Fast
```

### Static (Accurate Energy)
```bash
SYSTEM = Static calculation
ENCUT = 520
PREC = Accurate
IBRION = -1         # No relaxation
NSW = 0
EDIFF = 1E-6
ISMEAR = -5         # Tetrahedron
LREAL = .FALSE.
ALGO = Normal
```

### Band Structure (Step 1: SCF)
```bash
SYSTEM = SCF for bands
ENCUT = 520
PREC = Accurate
IBRION = -1
NSW = 0
EDIFF = 1E-6
ISMEAR = 0
SIGMA = 0.05
LCHARG = .TRUE.     # Write CHGCAR
```

### Band Structure (Step 2: Bands)
```bash
SYSTEM = Band structure
ENCUT = 520
PREC = Accurate
IBRION = -1
NSW = 0
EDIFF = 1E-6
ISMEAR = 0
SIGMA = 0.05
ICHARG = 11         # Read CHGCAR
LORBIT = 11         # Write projections
# Use line-mode KPOINTS
```

### DOS (Density of States)
```bash
SYSTEM = DOS calculation
ENCUT = 520
PREC = Accurate
IBRION = -1
NSW = 0
EDIFF = 1E-6
ISMEAR = -5         # Tetrahedron
LORBIT = 11
NEDOS = 3000
# Use dense k-mesh
```

### Molecular Dynamics (NVT)
```bash
SYSTEM = MD NVT
ENCUT = 400         # Lower for speed
PREC = Normal
IBRION = 0          # MD
NSW = 1000          # Number of steps
POTIM = 1.0         # Time step (fs)
TEBEG = 300
TEEND = 300
MDALGO = 2          # Nose-Hoover
SMASS = 0           # Thermostat mass
ISMEAR = 0
SIGMA = 0.1
```

### Phonons (DFPT)
```bash
SYSTEM = Phonon DFPT
ENCUT = 550         # Higher for forces
PREC = Accurate
IBRION = 6          # DFPT
NFREE = 2
POTIM = 0.015       # Displacement (Å)
EDIFF = 1E-8        # Tight!
ISMEAR = 0
SIGMA = 0.01
```

## Parameter Quick Lookup

| Parameter | Values | Meaning |
|-----------|--------|---------|
| **IBRION** | -1 | No ionic movement (static) |
| | 0 | Molecular dynamics |
| | 1 | RMM-DIIS (quasi-Newton) |
| | 2 | Conjugate gradient |
| | 5, 6 | Finite differences (forces, phonons) |
| **ISIF** | 2 | Relax ions only |
| | 3 | Relax ions + cell shape/volume |
| | 4 | Relax ions + cell shape (fixed volume) |
| **ISMEAR** | -5 | Tetrahedron (static, DOS) |
| | -1 | Fermi (metals) |
| | 0 | Gaussian |
| | 1, 2 | Methfessel-Paxton order 1, 2 |
| **ALGO** | Normal | Standard (blocked Davidson) |
| | Fast | Fast (RMM-DIIS) |
| | All | Try all algorithms |
| | Damped | Damped MD (HSE convergence) |
| **PREC** | Low | Testing only |
| | Normal | Standard |
| | Accurate | High precision (forces, phonons) |

## Convergence Criteria

| Quantity | Parameter | Typical Values |
|----------|-----------|----------------|
| Electronic (energy) | EDIFF | 1E-6 (standard), 1E-8 (tight) |
| Ionic (forces) | EDIFFG | -0.02 (standard), -0.01 (tight) |
| Max electronic steps | NELM | 60 (default), 100-200 (hard cases) |
| Max ionic steps | NSW | 100-200 (relaxation) |

## Recommended Settings by System

| System Type | ISMEAR | SIGMA | ENCUT | k-points |
|-------------|--------|-------|-------|----------|
| Metal | 1 | 0.2 | 520 | 8×8×8 |
| Semiconductor | 0 | 0.05 | 520 | 6×6×6 |
| Insulator | 0 | 0.05 | 520 | 4×4×4 |
| Molecule (gas phase) | 0 | 0.01 | 400 | 1×1×1 |
| Surface | 1 | 0.2 | 520 | 8×8×1 |

## k-Point Density Rules

| Calculation Type | Density | Example (cubic 4Å cell) |
|------------------|---------|-------------------------|
| Testing | 20-30 k-pts/Å⁻¹ | 4×4×4 |
| Relaxation | 30-40 k-pts/Å⁻¹ | 6×6×6 |
| Static/DOS | 40-50 k-pts/Å⁻¹ | 8×8×8 to 12×12×12 |
| Accurate | 50+ k-pts/Å⁻¹ | 12×12×12+ |

**Automatic k-spacing:**
```bash
KSPACING = 0.5    # Å⁻¹ (replaces KPOINTS file)
```

## ENCUT Guidelines

| Purpose | ENCUT (eV) |
|---------|-----------|
| Quick test | 400 |
| Standard | 520 (1.3 × ENMAX) |
| High accuracy | 550-600 |
| Forces/phonons | 550-600 |
| Elastic constants | 600+ |

**Always test convergence:** Vary ENCUT until ΔE < 1 meV/atom

## Common Output Files

| File | Contains |
|------|----------|
| OUTCAR | Main output (energies, forces, stress) |
| OSZICAR | Electronic + ionic convergence |
| CONTCAR | Final structure |
| DOSCAR | Density of states |
| EIGENVAL | Eigenvalues (band structure) |
| vasprun.xml | Complete output (structured) |
| WAVECAR | Wavefunctions (restart) |
| CHGCAR | Charge density |

## Checking Convergence

```bash
# Electronic convergence
grep "DAV:" OSZICAR | tail -20

# Final energy
grep "energy  without" OUTCAR | tail -1

# Forces
grep "TOTAL-FORCE" OUTCAR -A 50

# Check reached accuracy
grep "reached required accuracy" OUTCAR

# Total time
grep "Total CPU time" OUTCAR
```

## Common KPOINTS Files

**Gamma-centered mesh:**
```bash
Automatic mesh
0
Gamma
  8 8 8
  0 0 0
```

**Monkhorst-Pack mesh:**
```bash
Automatic mesh
0
Monkhorst-Pack
  8 8 8
  0 0 0
```

**Automatic (using KSPACING in INCAR instead):**
```bash
# No KPOINTS file needed
# Add to INCAR: KSPACING = 0.5
```

**Band structure path (FCC):**
```bash
k-path for band structure
10
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

## Error Quick Fixes

| Error | Quick Fix |
|-------|-----------|
| ZBRENT: fatal error | Reduce POTIM = 0.2 |
| Eigenvalues not converged | NELM = 200, ALGO = All |
| Sub-space matrix not hermitian | POTIM = 0.1, check structure |
| SCF not converging | ALGO = All, AMIX = 0.2 |
| BRMIX error | Lower BMIX = 0.0001 |

## HSE06 (Hybrid Functional)

```bash
LHFCALC = .TRUE.
HFSCREEN = 0.2
AEXX = 0.25
ALGO = All
TIME = 0.4
PRECFOCK = Fast
NKRED = 2           # Reduce k-mesh for HF
```

## DFT+U

```bash
LDAU = .TRUE.
LDAUTYPE = 2        # Dudarev
LDAUL = 2 -1        # d orbitals, then s/p
LDAUU = 5.0 0.0     # U values (eV)
LDAUJ = 0.0 0.0     # J values
```

## van der Waals

```bash
# DFT-D3 (Grimme)
IVDW = 11

# Or vdW-DF
GGA = MK
LUSE_VDW = .TRUE.
AGGAC = 0.0000
```

## Performance Tuning

```bash
# Parallelization
NCORE = 4           # Cores per orbital
KPAR = 4            # k-point parallelization
LPLANE = .TRUE.

# Large systems
LREAL = Auto        # Real-space projection

# Memory
NCORE = 4           # Reduce memory/core
```

## Decision Tree

```
What do you want to calculate?

├─ Structure optimization
│  ├─ Ions only → IBRION=2, ISIF=2
│  └─ Cell + ions → IBRION=2, ISIF=3
│
├─ Single-point energy
│  └─ Static → IBRION=-1, ISMEAR=-5
│
├─ Electronic structure
│  ├─ DOS → Dense k-mesh, ISMEAR=-5, LORBIT=11
│  └─ Bands → SCF + non-SCF (ICHARG=11, line-mode)
│
├─ Dynamics
│  └─ MD → IBRION=0, NSW=1000, MDALGO=2
│
└─ Properties
   ├─ Phonons → IBRION=6, EDIFF=1E-8
   └─ Elastic → IBRION=6, ISIF=3, NFREE=4
```

## Workflow Checklist

**Before running:**
- [ ] POSCAR has correct structure
- [ ] POTCAR matches POSCAR element order
- [ ] k-points appropriate for system size
- [ ] ENCUT ≥ 1.3 × ENMAX
- [ ] ISMEAR appropriate for system type

**After running:**
- [ ] Check "reached required accuracy" in OUTCAR
- [ ] Verify forces < EDIFFG threshold
- [ ] Check OSZICAR for convergence
- [ ] Save OUTCAR, vasprun.xml, CONTCAR

**For production:**
- [ ] Converge k-points (ΔE < 1 meV/atom)
- [ ] Converge ENCUT (ΔE < 1 meV/atom)
- [ ] Document all settings
- [ ] Save all input/output files
