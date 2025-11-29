# Advanced Functionals Subskill

Detailed guidance for hybrid functionals (HSE06, PBE0), GW, and DFT+U in VASP.

## Overview

Standard DFT (PBE, LDA) has limitations:
- **Band gaps:** Underestimated by 30-50%
- **Localized electrons:** Poorly described (transition metals, f-elements)
- **van der Waals:** Missing

Advanced methods address these issues at higher computational cost.

## Hybrid Functionals

### HSE06 (Recommended)

HSE06 is the most popular hybrid functional for solids. It uses screened Hartree-Fock exchange.

**INCAR:**
```bash
SYSTEM = HSE06 calculation

# Hybrid functional
LHFCALC = .TRUE.    # Activate exact exchange
HFSCREEN = 0.2      # HSE screening (Å⁻¹)
AEXX = 0.25         # 25% exact exchange
ALGO = All          # Recommended for hybrids
TIME = 0.4          # Damping for convergence

# Electronic
ENCUT = 520
EDIFF = 1E-6
PREC = Accurate
NELM = 100          # May need more SCF steps

# No ionic movement (static calculation)
IBRION = -1
NSW = 0

# Smearing
ISMEAR = 0
SIGMA = 0.05

# Performance
PRECFOCK = Fast     # Fast FFT for Fock exchange
NKRED = 2           # Reduce k-mesh for HF part by factor 2
# Or individually: NKREDX = 2, NKREDY = 2, NKREDZ = 2

# Output
LWAVE = .FALSE.
LCHARG = .FALSE.
```

**KPOINTS:**
```bash
# HSE is expensive, use sparser mesh than PBE
Automatic mesh
0
Gamma
  4 4 4             # Converge this!
  0 0 0
```

**Computational cost:** ~10-50× more expensive than PBE.

**Workflow:**
1. Relax with PBE (cheap)
2. Single-point HSE06 on PBE-relaxed structure
3. Optional: Relax with HSE (very expensive)

### HSE Parameters

| Parameter | HSE06 | HSE03 | Notes |
|-----------|-------|-------|-------|
| HFSCREEN | 0.2 | 0.3 | Screening (Å⁻¹) |
| AEXX | 0.25 | 0.25 | Exact exchange fraction |

**To use HSE03:**
```bash
HFSCREEN = 0.3
```

### PBE0 (Unscreened)

PBE0 uses 25% unscreened exact exchange (more expensive than HSE).

**INCAR:**
```bash
LHFCALC = .TRUE.
AEXX = 0.25
# No HFSCREEN (unscreened)
ALGO = All
TIME = 0.4
PRECFOCK = Fast
```

**Cost:** ~100× more expensive than PBE (no screening).

**Use:** When you need unscreened exchange (rare for solids).

### B3LYP (Molecules)

**INCAR:**
```bash
GGA = B3          # B3LYP functional
AEXX = 0.2
AGGAX = 0.72
AGGAC = 0.81
ALDAC = 0.19
LHFCALC = .TRUE.
```

**Use:** Molecules, comparison with quantum chemistry codes.

### Custom Hybrid

**Control exact exchange fraction:**
```bash
LHFCALC = .TRUE.
AEXX = 0.40         # 40% exact exchange (custom)
HFSCREEN = 0.2
```

**Tune for specific systems** (e.g., to match experimental gap).

## GW Approximation

GW calculates quasi-particle energies (most accurate band structures).

### G₀W₀ (One-Shot GW)

#### Step 1: DFT Ground State

**INCAR:**
```bash
SYSTEM = DFT for GW

# Electronic
ENCUT = 520
EDIFF = 1E-8        # Tight convergence!
PREC = Accurate
ALGO = Exact        # Or Normal

# Many bands
NBANDS = 200        # Rule: 2-4× number of occupied bands
# Estimate: occupied = N_electrons / 2

# No ionic movement
IBRION = -1
NSW = 0

# Smearing
ISMEAR = 0
SIGMA = 0.01

# Optics (for GW)
LOPTICS = .TRUE.    # Compute frequency-dependent ε
CSHIFT = 0.1        # Avoid divergence

# Output
LWAVE = .TRUE.      # Essential for GW
LCHARG = .FALSE.
```

**KPOINTS (sparse mesh OK):**
```bash
Automatic mesh
0
Gamma
  4 4 4
  0 0 0
```

**Run VASP** → produces WAVECAR, WAVEDER

#### Step 2: G₀W₀

**INCAR:**
```bash
SYSTEM = G0W0

# Read wavefunctions
ALGO = GW0          # One-shot GW
# Or EVGW for eigenvalue-only self-consistent
# Or QPGW for quasi-particle self-consistent

# GW parameters
NOMEGA = 50         # Frequency grid points
NBANDS = 200        # Must match step 1

# Response function cutoff
ENCUTGW = 200       # Lower than ENCUT (convergence test!)

# No ionic movement
IBRION = -1
NSW = 0

# Output
LWAVE = .FALSE.
```

**Run on WAVECAR from step 1** → produces quasi-particle energies in OUTCAR

**Extract band gap:**
```bash
grep "QP energies" OUTCAR -A 50
```

**Computational cost:** ~100-1000× more expensive than PBE.

### scGW (Self-Consistent GW)

**INCAR (step 2):**
```bash
ALGO = QPGW         # Quasi-particle self-consistent
NELM = 10           # GW iterations
```

**More accurate but extremely expensive.**

### GW Parameters

| Parameter | Typical Value | Notes |
|-----------|---------------|-------|
| NOMEGA | 50-100 | Frequency grid, converge this |
| ENCUTGW | 150-300 eV | Cutoff for response, test convergence |
| NBANDS | 2-4 × occupied | More = more accurate |
| OMEGAMAX | 30-50 eV | Max frequency |

### GW Convergence

**Must converge:**
1. **NBANDS:** Test 100, 150, 200, ...
2. **ENCUTGW:** Test 150, 200, 250 eV
3. **NOMEGA:** Test 50, 75, 100
4. **k-points:** Same as DFT

**Typical convergence criteria:**
- ΔE_gap < 0.1 eV

## DFT+U

DFT+U adds on-site Coulomb repulsion for localized electrons (d, f orbitals).

### Dudarev Method (Recommended)

**INCAR:**
```bash
SYSTEM = DFT+U calculation

# Electronic
ENCUT = 520
EDIFF = 1E-6
PREC = Accurate

# DFT+U
LDAU = .TRUE.
LDAUTYPE = 2        # Dudarev (U-J formalism)
LDAUL = 2 -1        # l quantum number per atom type (2=d, 3=f, -1=no U)
LDAUU = 5.0 0.0     # U value (eV) per atom type
LDAUJ = 0.0 0.0     # J value (ignored for LDAUTYPE=2)
LDAUPRINT = 2       # Output level (0=minimal, 2=verbose)

# Initial magnetic moments (often needed)
ISPIN = 2
MAGMOM = 5.0 5.0 0.0  # Initial moments per atom

# No ionic movement
IBRION = -1
NSW = 0

# Smearing
ISMEAR = 0
SIGMA = 0.05
```

**Example for Fe₂O₃:**
- Fe: 3d, apply U
- O: 2p, no U

```bash
# POSCAR order: Fe Fe O O O
LDAUL = 2 2 -1 -1 -1   # d orbitals for Fe, none for O
LDAUU = 5.0 5.0 0.0 0.0 0.0
LDAUJ = 0.0 0.0 0.0 0.0 0.0
```

### Liechtenstein Method

**INCAR:**
```bash
LDAUTYPE = 1        # Liechtenstein (separate U and J)
LDAUL = 2 -1
LDAUU = 5.0 0.0     # U value
LDAUJ = 1.0 0.0     # J value (not ignored)
```

**Use:** When U and J are separately defined.

### Choosing U Values

**Sources:**
1. **Literature:** Search for your system
2. **Linear response:** Calculate U from first principles (expensive)
3. **Empirical:** Tune to match experiment

**Typical U values (eV):**
| Element | Orbital | U value |
|---------|---------|---------|
| Fe | 3d | 4.0-5.0 |
| Co | 3d | 5.0-6.0 |
| Ni | 3d | 6.0-7.0 |
| Cu | 3d | 7.0-8.0 |
| Mn | 3d | 4.0-5.0 |
| Ti | 3d | 3.0-5.0 |
| V | 3d | 3.0-4.0 |
| Ce | 4f | 5.0-6.0 |
| U | 5f | 4.0-5.0 |

**Always test sensitivity to U!**

### Occupation Matrix Control

**INCAR:**
```bash
LDAU = .TRUE.
LDAUTYPE = 2

# Occupation matrix method
LMAXMIX = 4         # For d electrons (6 for f)
```

### DFT+U+V

Add inter-site interactions (advanced).

**INCAR:**
```bash
LDAUTYPE = 4        # DFT+U+V
LDAUL = 2 -1
LDAUU = 5.0 0.0
LDAUJ = 0.0 0.0
```

**Requires:** V matrix defined in INCAR or separate file.

## van der Waals Corrections

### DFT-D3 (Grimme)

**Simplest vdW correction, add to PBE:**

**INCAR:**
```bash
IVDW = 11           # DFT-D3 (Grimme)
```

**Or with Becke-Jonson damping:**
```bash
IVDW = 12           # DFT-D3 with BJ damping
```

**Cost:** Negligible (empirical correction).

### vdW-DF Functionals

**Non-local correlation functional:**

#### optPBE-vdW

**INCAR:**
```bash
GGA = MK            # optPBE
LUSE_VDW = .TRUE.
AGGAC = 0.0000
```

#### vdW-DF2

**INCAR:**
```bash
GGA = ML            # revPBE for vdW-DF2
LUSE_VDW = .TRUE.
AGGAC = 0.0000
Zab_vdW = -1.8867
```

**Cost:** ~2× PBE (non-local functional).

### DFT-D2 (Older)

**INCAR:**
```bash
IVDW = 1            # DFT-D2 (older version)
# Or IVDW = 10 for DFT-D2 with different params
```

## Spin-Orbit Coupling (SOC)

For heavy elements, include relativistic effects.

**INCAR:**
```bash
LSORBIT = .TRUE.    # Activate SOC

# Need special POTCARs
# Use non-collinear calculation
LNONCOLLINEAR = .TRUE.

# Increase NBANDS
NBANDS = 100        # SOC splits bands

# Optional: set spin quantization axis
SAXIS = 0 0 1       # z-axis
```

**Requirements:**
- Non-collinear POTCARs
- More bands (2× for SOC)

**Use:** Heavy elements (Pb, Bi, Pt, topological materials)

## Comparison of Methods

| Method | Accuracy | Cost | Use Case |
|--------|----------|------|----------|
| PBE | Baseline | 1× | Structure, trends |
| HSE06 | Good gaps | 10-50× | Band gaps, defects |
| PBE0 | Better gaps | 100× | Small systems |
| G₀W₀ | Best gaps | 100-1000× | Accurate gaps |
| scGW | Most accurate | 1000-10000× | Reference |
| DFT+U | Good for TM | 1-2× | Correlated systems |
| vdW-DF | vdW interactions | 2× | Layered, molecular |
| DFT-D3 | Quick vdW | 1× | Cheap vdW correction |

## Complete Examples

### HSE06 Band Gap

**Step 1: Relax with PBE**
```bash
# INCAR (PBE relaxation)
ENCUT = 520
IBRION = 2
ISIF = 3
EDIFFG = -0.02
```

**Step 2: HSE06 on relaxed structure**
```bash
# Copy CONTCAR to POSCAR

# INCAR (HSE06)
LHFCALC = .TRUE.
HFSCREEN = 0.2
AEXX = 0.25
ALGO = All
TIME = 0.4
PRECFOCK = Fast
NKRED = 2
IBRION = -1
NSW = 0
ISMEAR = 0
SIGMA = 0.05
```

**Extract gap:**
```bash
grep "E-fermi" OUTCAR
# Analyze EIGENVAL or vasprun.xml
```

### G₀W₀ Workflow

**Step 1: DFT**
```bash
# INCAR
ALGO = Exact
NBANDS = 200
LOPTICS = .TRUE.
CSHIFT = 0.1
EDIFF = 1E-8
LWAVE = .TRUE.
```

**Step 2: GW**
```bash
# New INCAR (keep POSCAR, KPOINTS, POTCAR)
ALGO = GW0
NOMEGA = 50
NBANDS = 200
ENCUTGW = 200
LWAVE = .FALSE.
```

**Extract QP gap:**
```bash
grep "QP energies" OUTCAR
```

### DFT+U on NiO

**INCAR:**
```bash
# POSCAR: Ni O alternating
# Apply U to Ni 3d

LDAU = .TRUE.
LDAUTYPE = 2
LDAUL = 2 -1        # Ni: d, O: none
LDAUU = 6.0 0.0     # U = 6 eV for Ni
LDAUJ = 0.0 0.0

# Spin-polarized
ISPIN = 2
MAGMOM = 2.0 0.0    # Initial moment on Ni

# Electronic
ENCUT = 520
EDIFF = 1E-6
ISMEAR = 0
SIGMA = 0.05

# Relaxation
IBRION = 2
ISIF = 3
EDIFFG = -0.02
```

## Troubleshooting

### HSE Not Converging

**Problem:** SCF cycles not converging

**Solutions:**
1. **Try ALGO = All or Damped:**
   ```bash
   ALGO = Damped
   TIME = 0.5
   ```

2. **Increase NELM:**
   ```bash
   NELM = 200
   ```

3. **Adjust mixing:**
   ```bash
   AMIX = 0.2
   BMIX = 0.0001
   ```

4. **Start from PBE:**
   - Copy WAVECAR from PBE calculation
   - May help convergence

### GW Gives Unrealistic Results

**Problem:** QP energies too large/small

**Solutions:**
1. **Converge ENCUTGW:**
   - Test 150, 200, 250 eV
   - Should converge to < 0.1 eV

2. **Converge NBANDS:**
   - Try 2×, 3×, 4× occupied bands

3. **Check DFT starting point:**
   - GW depends on DFT orbitals
   - Use tight EDIFF = 1E-8 for DFT

4. **Verify LOPTICS:**
   - Must be .TRUE. in step 1

### DFT+U Results Depend on Initial MAGMOM

**Problem:** Different results with different initial moments

**Solutions:**
1. **This is expected:**
   - DFT+U can have multiple solutions
   - Try different MAGMOM

2. **Find lowest energy state:**
   - Test multiple initial configurations
   - Choose lowest energy

3. **Use symmetry:**
   - Set MAGMOM consistent with expected magnetic structure

## Performance Optimization

### HSE06

```bash
# Parallelization
NCORE = 4
KPAR = 2            # k-point parallelization

# Reduce HF k-mesh
NKRED = 2           # Or NKREDX/Y/Z individually

# Fast Fock exchange
PRECFOCK = Fast

# Sparse k-mesh (but converge!)
KSPACING = 0.3      # Coarser than PBE
```

### GW

```bash
# Use sparse k-mesh
# 4×4×4 often sufficient

# Reduce NOMEGA if converged
NOMEGA = 50         # Test 50, 75

# Reduce ENCUTGW
ENCUTGW = 150       # But converge!

# Parallelization
NCORE = 1           # Better for GW
```

### DFT+U

```bash
# Similar to PBE (no extra cost)
# Standard parallelization
NCORE = 4
KPAR = 4
```

## Best Practices

1. **Always start with PBE:**
   - Relax with PBE
   - Use for initial structure

2. **Test convergence:**
   - HSE: k-points, ENCUT (like PBE)
   - GW: NBANDS, ENCUTGW, NOMEGA

3. **Document U values:**
   - Record source of U (literature, fit, etc.)
   - Test sensitivity

4. **Compare methods:**
   - PBE → HSE06 → GW hierarchy
   - Justify choice based on accuracy/cost

5. **Save outputs:**
   - vasprun.xml has all results
   - Large files, but comprehensive

## References

- HSE: https://vasp.at/wiki/index.php/Hybrid_functionals
- GW: https://vasp.at/wiki/index.php/GW_calculations
- DFT+U: https://vasp.at/wiki/index.php/DFT%2BU
- vdW: https://vasp.at/wiki/index.php/Van_der_Waals_interactions

## See Also

- `electronic-structure.md` - Band structure workflow
- Main SKILL.md - Full parameter reference
