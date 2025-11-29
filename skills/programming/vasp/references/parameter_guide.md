# VASP Parameter Reference Guide

Comprehensive reference for all commonly used VASP parameters.

## Table of Contents

1. [Electronic Structure](#electronic-structure)
2. [Ionic Relaxation](#ionic-relaxation)
3. [Molecular Dynamics](#molecular-dynamics)
4. [Exchange-Correlation](#exchange-correlation)
5. [Performance](#performance)
6. [Output Control](#output-control)
7. [Advanced](#advanced)

---

## Electronic Structure

### ENCUT
**Description:** Plane-wave energy cutoff
**Type:** Real
**Default:** Maximum ENMAX from POTCAR
**Units:** eV
**Recommendations:**
- Standard: 1.3 × ENMAX (typically 520 eV for PAW)
- High accuracy: 1.5 × ENMAX (550-600 eV)
- Forces/phonons: 550-600 eV

### EDIFF
**Description:** Electronic convergence criterion
**Type:** Real
**Default:** 1E-4
**Units:** eV
**Recommendations:**
- Testing: 1E-4
- Standard: 1E-6
- Tight (forces): 1E-7
- Phonons/GW: 1E-8

### NELM
**Description:** Maximum electronic SCF steps
**Type:** Integer
**Default:** 60
**Typical:** 60-100 (increase for difficult convergence)

### PREC
**Description:** Precision level
**Type:** String
**Values:**
- `Low`: Testing only
- `Normal`: Standard calculations
- `Accurate`: High accuracy (forces, phonons)
**Default:** `Normal`

### ALGO
**Description:** Electronic minimization algorithm
**Type:** String
**Values:**
- `Normal`: Blocked Davidson
- `Fast`: RMM-DIIS (faster)
- `All`: Try all algorithms sequentially
- `Damped`: Damped MD (for HSE)
- `Exact`: Exact diagonalization (GW)
**Default:** `Normal`
**Recommendations:**
- Standard: `Fast`
- Convergence problems: `All`
- Hybrids: `All` or `Damped`

### ISMEAR
**Description:** Smearing method
**Type:** Integer
**Values:**
| Value | Method | Use Case |
|-------|--------|----------|
| -5 | Tetrahedron | Static, DOS, accurate |
| -4 | Tetrahedron+Blöchl | Like -5 |
| -1 | Fermi | Metals |
| 0 | Gaussian | General |
| 1 | Methfessel-Paxton 1 | Metals, relaxation |
| 2 | Methfessel-Paxton 2 | Metals |
**Recommendations:**
- Metals (relax): 1
- Semiconductors (relax): 0
- Static/DOS: -5
- Large systems: -1

### SIGMA
**Description:** Smearing width
**Type:** Real
**Units:** eV
**Recommendations:**
- Metals: 0.1-0.2
- Semiconductors: 0.01-0.05
- Check entropy term (T×S < 1 meV/atom)

---

## Ionic Relaxation

### IBRION
**Description:** Ionic update algorithm
**Type:** Integer
**Values:**
| Value | Method | Use Case |
|-------|--------|----------|
| -1 | None | Static calculation |
| 0 | Molecular dynamics | MD |
| 1 | RMM-DIIS (quasi-Newton) | Near minimum |
| 2 | Conjugate gradient | Standard relaxation |
| 3 | Damped MD | Rough energy surface |
| 5, 6 | Finite differences | Phonons, elastic |
**Default:** -1
**Recommendation:** Use 2 for most relaxations

### ISIF
**Description:** Calculate stress and relax ions/cell
**Type:** Integer
**Values:**
| Value | Ions | Shape | Volume | Use Case |
|-------|------|-------|--------|----------|
| 2 | Yes | No | No | Fixed cell |
| 3 | Yes | Yes | Yes | Full relaxation |
| 4 | Yes | Yes | Yes | Fixed shape |
| 7 | Yes | Yes | Yes | Volume only |
**Default:** 2
**Recommendation:** Use 3 for full relaxation, 2 for fixed cell

### NSW
**Description:** Maximum ionic steps
**Type:** Integer
**Default:** 0
**Typical:** 100-200 for relaxation

### EDIFFG
**Description:** Ionic convergence criterion
**Type:** Real
**Units:**
- Negative: Force threshold (eV/Å)
- Positive: Energy threshold (eV)
**Recommendations:**
- Standard: -0.02 (forces < 0.02 eV/Å)
- Tight: -0.01
- Loose: -0.05

### POTIM
**Description:** Time step (MD) or scaling factor (relaxation)
**Type:** Real
**Units:** fs (for MD)
**Default:** 0.5
**Recommendations:**
- MD: 0.5-2.0 fs (0.5 for H-containing)
- Relaxation: 0.1-0.5 (reduce if unstable)

---

## Molecular Dynamics

### MDALGO
**Description:** MD algorithm/thermostat
**Type:** Integer
**Values:**
| Value | Method | Ensemble |
|-------|--------|----------|
| 0 | Standard | NVE |
| 1 | Andersen | NVT |
| 2 | Nose-Hoover | NVT |
| 3 | Parrinello-Rahman | NPT |
**Default:** 0

### TEBEG / TEEND
**Description:** Starting/ending temperature
**Type:** Real
**Units:** K
**Example:** TEBEG = 300, TEEND = 300 (constant T)

### SMASS
**Description:** Nose mass (thermostat coupling)
**Type:** Real
**Default:** -3
**Recommendations:**
- 0: Automatic (recommended)
- 0.1-1.0: Manual control (smaller = stronger coupling)
- -3: Microcanonical (NVE)

### PMASS
**Description:** Barostat mass (NPT)
**Type:** Real
**Typical:** 500

---

## Exchange-Correlation

### GGA
**Description:** Exchange-correlation functional
**Type:** String
**Values:**
- `PE`: PBE (default)
- `91`: PW91
- `RP`: revPBE
- `AM`: AM05
- `MK`: optPBE-vdW
**Default:** `PE` (PBE)

### LHFCALC
**Description:** Activate hybrid functionals
**Type:** Logical
**Default:** .FALSE.
**Set to .TRUE. for HSE, PBE0**

### HFSCREEN
**Description:** Screening parameter for HSE
**Type:** Real
**Units:** Å⁻¹
**Values:**
- 0.2: HSE06
- 0.3: HSE03
- Omit for PBE0 (unscreened)

### AEXX
**Description:** Exact exchange fraction
**Type:** Real
**Default:** 0.25 (for hybrids)

---

## Performance

### NCORE
**Description:** Cores per orbital (band parallelization)
**Type:** Integer
**Recommendations:**
- Typical: 4-8
- Rule: ~sqrt(total cores) or cores_per_node/2
- Memory: Larger NCORE = less memory per core

### KPAR
**Description:** Number of k-point groups
**Type:** Integer
**Recommendations:**
- Set to number of k-points (if enough cores)
- Good for many k-points
- NCORE and KPAR are related: total_cores = NCORE × KPAR × (bands parallelization)

### LREAL
**Description:** Real-space projection
**Type:** Logical or String
**Values:**
- .FALSE.: Reciprocal space (accurate)
- .TRUE. or `Auto`: Real space (faster for large systems)
**Recommendations:**
- Small systems (< 20 atoms): .FALSE.
- Large systems (> 100 atoms): `Auto`

### LPLANE
**Description:** Plane-wise FFT distribution
**Type:** Logical
**Default:** .FALSE.
**Set to .TRUE. for better memory distribution**

---

## Output Control

### LWAVE
**Description:** Write WAVECAR
**Type:** Logical
**Default:** .TRUE.
**Set to .FALSE. to save disk space**

### LCHARG
**Description:** Write CHGCAR
**Type:** Logical
**Default:** .TRUE.
**Set to .FALSE. unless needed for band structure**

### LORBIT
**Description:** Projected DOS
**Type:** Integer
**Values:**
- 0: No projected DOS
- 10, 11: Projected DOS (DOSCAR, PROCAR)
- 11: Recommended (no RWIGS needed)

### NEDOS
**Description:** DOS resolution
**Type:** Integer
**Default:** 301
**High resolution:** 3000-5000

### NWRITE
**Description:** Output verbosity
**Type:** Integer
**Default:** 2
**Values:** 0 (minimal), 1 (normal), 2 (verbose), 3 (debug)

---

## Advanced

### DFT+U

**LDAU**
- Type: Logical
- Description: Activate DFT+U
- Default: .FALSE.

**LDAUTYPE**
- Type: Integer
- Description: 1 (Liechtenstein), 2 (Dudarev)
- Default: 2

**LDAUL**
- Type: Integer array
- Description: l quantum number per atom type (2=d, 3=f, -1=no U)

**LDAUU**
- Type: Real array
- Description: U value (eV) per atom type

**LDAUJ**
- Type: Real array
- Description: J value (eV), ignored if LDAUTYPE=2

### van der Waals

**IVDW**
- Type: Integer
- Description: vdW correction method
- Values:
  - 0: No correction
  - 11: DFT-D3 (Grimme)
  - 12: DFT-D3 with BJ damping

**LUSE_VDW**
- Type: Logical
- Description: Activate vdW-DF functionals
- Use with GGA = MK or ML

### Spin

**ISPIN**
- Type: Integer
- Description: 1 (non-spin-polarized), 2 (spin-polarized)
- Default: 1

**MAGMOM**
- Type: Real array
- Description: Initial magnetic moments per atom (μ_B)

### Spin-Orbit Coupling

**LSORBIT**
- Type: Logical
- Description: Activate spin-orbit coupling
- Default: .FALSE.

**LNONCOLLINEAR**
- Type: Logical
- Description: Non-collinear magnetism
- Required for SOC

### Constraints

**PSTRESS**
- Type: Real
- Description: External pressure (kBar)
- Default: 0

---

## Quick Selection Guide

### Bulk Metal Relaxation
```bash
ENCUT = 520
PREC = Accurate
IBRION = 2
ISIF = 3
NSW = 100
EDIFFG = -0.02
ISMEAR = 1
SIGMA = 0.2
ALGO = Fast
LREAL = Auto
```

### Static Calculation
```bash
ENCUT = 520
PREC = Accurate
IBRION = -1
NSW = 0
ISMEAR = -5
ALGO = Normal
LREAL = .FALSE.
```

### Band Structure (Step 2)
```bash
ENCUT = 520
IBRION = -1
ICHARG = 11
LORBIT = 11
ISMEAR = 0
SIGMA = 0.05
```

### Phonons (DFPT)
```bash
ENCUT = 550
EDIFF = 1E-8
PREC = Accurate
IBRION = 6
NFREE = 2
POTIM = 0.015
LREAL = .FALSE.
```

---

## Complete INCAR Index

For exhaustive list, see VASP Wiki: https://vasp.at/wiki/index.php/Category:INCAR

Common omissions from this guide:
- Magnetic (detailed spin settings)
- Optics (LOPTICS, CSHIFT)
- Time-dependent DFT
- Exact exchange details
- Advanced parallelization (NPAR deprecated)

This guide covers ~90% of typical VASP usage.
