# Calculator Comparison Guide

## Overview

Different DFT codes and calculators have strengths and weaknesses. This guide helps choose appropriate tools and understand their differences.

## Calculator Types in ASE

### 1. EMT (Effective Medium Theory)

**Pros:**
- Extremely fast
- No external dependencies
- Good for testing workflows
- Qualitatively reasonable for FCC metals

**Cons:**
- Limited to specific elements (Cu, Au, Ag, Pt, Pd, Ni, Al, Pb, Fe, C, H, O, N)
- Not quantitatively accurate
- No magnetism, no bandgap

**Use Cases:**
- Algorithm development
- Testing scripts
- Teaching
- Initial structure screening

```python
from ase.calculators.emt import EMT

atoms.calc = EMT()
```

### 2. GPAW (Grid-based PAW)

**Pros:**
- Open-source, free
- Good ASE integration
- Plane-wave (PW) and LCAO modes
- TDDFT for optical properties
- Good documentation

**Cons:**
- Slower than VASP for PW
- Memory intensive for large systems
- Limited parallelization vs commercial codes

**Use Cases:**
- Academic research
- Method development
- Optical properties
- Free alternative to VASP

**Setup:**
```python
from gpaw import GPAW, PW

atoms.calc = GPAW(
    mode=PW(500),  # 500 eV cutoff
    xc='PBE',
    kpts=(4, 4, 4),
    txt='output.txt',
    convergence={'energy': 1e-6}
)
```

### 3. VASP (Vienna Ab initio Simulation Package)

**Pros:**
- Industry standard
- Highly optimized
- Excellent parallelization
- Extensive features (GW, BSE, DFPT, etc.)
- Large user community

**Cons:**
- Commercial license required
- Closed source
- Steep learning curve
- Complex input/output

**Use Cases:**
- Production calculations
- High-accuracy needs
- Complex materials (strongly correlated, etc.)
- Large-scale screening

**Setup (via ASE):**
```python
from ase.calculators.vasp import Vasp

atoms.calc = Vasp(
    xc='PBE',
    encut=500,
    kpts=(4, 4, 4),
    ismear=1,
    sigma=0.1,
    prec='Accurate',
    ediff=1e-6
)
```

### 4. Quantum ESPRESSO

**Pros:**
- Open-source
- Very flexible
- DFPT implementation (phonons, dielectric)
- Good for complex materials

**Cons:**
- Complex input format
- ASE interface less mature than GPAW/VASP
- Steeper learning curve

**Setup:**
```python
from ase.calculators.espresso import Espresso

pseudopotentials = {'Cu': 'Cu.pbe-n-rrkjus_psl.1.0.0.UPF'}

atoms.calc = Espresso(
    pseudopotentials=pseudopotentials,
    tstress=True,
    tprnfor=True,
    ecutwfc=50,  # Ry
    kpts=(4, 4, 4)
)
```

### 5. Classical Potentials (LAMMPS)

**Pros:**
- Very fast (millions of atoms)
- Long timescales (ns-μs)
- Good for finite-temperature dynamics
- Many potentials available

**Cons:**
- Limited accuracy
- No electronic structure
- Transferability issues
- Requires parameterized potentials

**Use Cases:**
- Large-scale MD
- Mechanical properties
- Diffusion
- Initial structure generation

```python
from ase.calculators.lammpslib import LAMMPSlib

cmds = ["pair_style eam/alloy",
        "pair_coeff * * Cu_u3.eam Cu"]

atoms.calc = LAMMPSlib(lmpcmds=cmds)
```

## Feature Comparison

| Feature | EMT | GPAW | VASP | QE | LAMMPS |
|---------|-----|------|------|----|----- --|
| Speed | ★★★★★ | ★★★ | ★★★★ | ★★★ | ★★★★★ |
| Accuracy | ★ | ★★★★ | ★★★★★ | ★★★★ | ★★ |
| Ease of Use | ★★★★★ | ★★★★ | ★★ | ★★ | ★★★ |
| Parallelization | N/A | ★★★ | ★★★★★ | ★★★★ | ★★★★★ |
| Open Source | Yes | Yes | No | Yes | Yes |
| Cost | Free | Free | $$$$ | Free | Free |

## Exchange-Correlation Functionals

### LDA (Local Density Approximation)

- **Accuracy:** Moderate
- **Cost:** Cheap
- **Best for:** Lattice constants often too small
- **Avoid for:** Molecular systems, weak bonds

### GGA (Generalized Gradient Approximation)

**PBE** (most common):
- Balanced accuracy
- Lattice constants ~1% too large
- Standard choice

**PBEsol**:
- Better lattice constants
- Good for solids

**RPBE**:
- Better for molecules
- Adsorption energies
- Barriers

### Meta-GGA

**SCAN**:
- More accurate than PBE
- Higher cost
- Better all-around

### Hybrid Functionals

**HSE06, PBE0:**
- Accurate band gaps
- Much more expensive (10-100×)
- Use for semiconductors, insulators

### van der Waals

**DFT-D3, vdW-DF**:
- Essential for physisorption
- Layered materials
- Molecular crystals

```python
# GPAW with vdW
atoms.calc = GPAW(
    xc='vdW-DF',
    # or
    # xc={'name': 'PBE', 'backend': 'libvdwxc'}
)

# VASP with DFT-D3
atoms.calc = Vasp(
    ivdw=11  # DFT-D3
)
```

## Performance Comparison

### Benchmark: Cu bulk (4-atom cell)

| Calculator | Time | Memory | Accuracy (vs VASP) |
|------------|------|--------|-------------------|
| EMT | 0.1s | 10 MB | N/A |
| GPAW-LCAO | 5s | 200 MB | ±0.1 eV |
| GPAW-PW | 30s | 500 MB | ±0.01 eV |
| VASP | 10s | 400 MB | Reference |
| QE | 15s | 300 MB | ±0.01 eV |

*(Approximate, system-dependent)*

## When to Use What

### Workflow Development
→ **EMT** or **GPAW-LCAO**

### High-Throughput Screening
→ **VASP** (with good parallelization) or **GPAW**

### Band Gaps
→ **HSE06** (VASP or QE)

### Surfaces/Catalysis
→ **GPAW** or **VASP** with **vdW corrections**

### Phonons
→ **VASP** (DFPT) or **Phonopy** + any calculator

### Large Systems (>200 atoms DFT)
→ **VASP** (best scaling) or **LAMMPS** (if accuracy allows)

### Teaching/Learning
→ **GPAW** (free, good docs) or **EMT** (instant)

### Magnetic Systems
→ **VASP** (noncollinear) or **GPAW**

### Strongly Correlated (DFT+U)
→ **VASP** (most tested) or **QE**

## Calculator-Specific Tips

### GPAW

```python
# LCAO for quick calculations
calc_lcao = GPAW(mode='lcao', basis='dzp', xc='PBE')

# PW for accuracy
calc_pw = GPAW(mode=PW(500), xc='PBE')

# Save/load for analysis
calc.write('gpaw_restart.gpw')
calc = GPAW('gpaw_restart.gpw', txt=None)
```

### VASP

```python
# Important tags
calc = Vasp(
    prec='Accurate',      # Precision
    ediff=1e-6,           # SCF tolerance
    encut=520,            # 1.3× ENMAX from POTCAR
    ismear=1, sigma=0.1,  # Methfessel-Paxton for metals
    algo='Fast',          # Algorithm
    lreal='Auto'          # Real-space projection (large systems)
)

# For surfaces
calc.set(
    ismear=0,      # Gaussian smearing
    sigma=0.05,    # Small smearing
    idipol=3,      # Dipole correction (z-direction)
    ldipol=True
)
```

### Quantum ESPRESSO

```python
# Typical settings
input_data = {
    'control': {
        'calculation': 'scf',
        'prefix': 'pwscf',
        'pseudo_dir': './pseudo/',
    },
    'system': {
        'ecutwfc': 50,  # Ry
        'ecutrho': 400,
        'occupations': 'smearing',
        'smearing': 'mv',
        'degauss': 0.01
    },
    'electrons': {
        'conv_thr': 1e-8
    }
}

calc = Espresso(input_data=input_data,
                pseudopotentials=pseudopotentials,
                kpts=(4,4,4))
```

## Accuracy Validation

Always validate against:
1. **Experimental data:** Lattice constants, formation energies
2. **Benchmark databases:** Materials Project, OQMD
3. **High-level theory:** Quantum chemistry (CCSD(T)) for molecules
4. **Code comparison:** Verify with 2+ codes for important results

## Cost-Accuracy Tradeoff

```
EMT → GPAW-LCAO → GPAW-PW/VASP-GGA → VASP-Hybrid → Quantum Chemistry

↑ Speed                                                         Accuracy ↑
↓ Accuracy                                                      Cost ↓
```

## Migration Between Calculators

```python
def convert_gpaw_to_vasp(atoms_gpaw):
    """Convert GPAW calculation to VASP."""

    # Extract converged structure
    atoms = atoms_gpaw.copy()

    # Setup equivalent VASP calculation
    atoms.calc = Vasp(
        xc='PBE',
        encut=500,
        kpts=atoms_gpaw.calc.parameters['kpts']
    )

    return atoms
```

## Troubleshooting

### GPAW Out of Memory
- Reduce k-points or cutoff temporarily
- Use LCAO mode
- Use domain decomposition: `parallel={'domain': N}`

### VASP Not Converging
- Adjust ALGO (Fast → Normal → All)
- Increase NELM (max electronic steps)
- Check MAGMOM (magnetic systems)
- Reduce POTIM (ionic steps)

### Quantum ESPRESSO Slow
- Optimize ecutrho (8-12× ecutwfc typically)
- Use optimized pseudopotentials
- Proper k-point sampling

## References

- **GPAW:** https://wiki.fysik.dtu.dk/gpaw/
- **VASP:** https://www.vasp.at/
- **Quantum ESPRESSO:** https://www.quantum-espresso.org/
- **ASE Calculators:** https://wiki.fysik.dtu.dk/ase/ase/calculators/calculators.html

## See Also

- `convergence_testing.md`: Testing calculation settings
- `best_practices.md`: General workflow recommendations
