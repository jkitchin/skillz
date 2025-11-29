# Convergence Testing Workflow

## Overview

Systematic convergence testing is essential for reliable computational materials science results. This guide provides standard procedures for all critical parameters.

## General Principles

1. **One parameter at a time:** Test each parameter independently
2. **Document everything:** Save convergence data for publication
3. **Define tolerance:** e.g., ΔE < 1 meV/atom
4. **Start coarse, refine:** Save computational time
5. **Use same settings:** Keep other parameters fixed while testing one

## 1. k-Point Convergence

### For Bulk Systems

```python
from ase.build import bulk
from ase.optimize import BFGS
import numpy as np

atoms = bulk('Cu', 'fcc', a=3.6)

kpoint_meshes = [
    (2, 2, 2),
    (4, 4, 4),
    (6, 6, 6),
    (8, 8, 8),
    (10, 10, 10),
    (12, 12, 12)
]

energies = []
for k in kpoint_meshes:
    atoms_test = atoms.copy()
    atoms_test.calc = calculator
    atoms_test.calc.set(kpts=k)

    opt = BFGS(atoms_test, logfile=None)
    opt.run(fmax=0.01)

    E = atoms_test.get_potential_energy() / len(atoms_test)
    energies.append(E)

    print(f"k-points {k}: E = {E:.6f} eV/atom")

# Check convergence
for i in range(1, len(energies)):
    diff = abs(energies[i] - energies[i-1]) * 1000  # meV
    print(f"Δ({kpoint_meshes[i-1]} → {kpoint_meshes[i]}): {diff:.2f} meV/atom")
```

**Convergence Criteria:**
- Total energies: < 1 meV/atom
- Forces: < 10 meV/Å
- Stresses: < 0.1 GPa

### For Surfaces/Slabs

In-plane k-points scale with supercell size:
- 3×3 slab: ~4×4×1 k-points
- 4×4 slab: ~3×3×1 k-points

Single k-point in vacuum direction usually sufficient.

## 2. Energy Cutoff Convergence

### Plane-Wave Codes

```python
cutoffs = [300, 350, 400, 450, 500, 550, 600]  # eV

energies = []
for ecut in cutoffs:
    atoms_test = atoms.copy()
    atoms_test.calc = calculator
    atoms_test.calc.set(encut=ecut, kpts=(8,8,8))

    opt = BFGS(atoms_test, logfile=None)
    opt.run(fmax=0.01)

    E = atoms_test.get_potential_energy() / len(atoms_test)
    energies.append(E)

    print(f"Cutoff {ecut} eV: E = {E:.6f} eV/atom")
```

**Guidelines:**
- Forces require higher cutoff than energies
- Check pseudopotential recommendations
- Convergence within 1-5 meV/atom typical

### For LCAO Codes (GPAW)

Test basis set size:
- dzp (double-zeta + polarization)
- tzp (triple-zeta + polarization)

## 3. Slab Thickness (Surfaces)

```python
from ase.build import fcc111

n_layers = [5, 7, 9, 11, 13]
surface_energies = []

for n in n_layers:
    slab = fcc111('Cu', size=(3, 3, n), vacuum=10)

    # Fix bottom 2 layers
    z_pos = slab.get_positions()[:, 2]
    mask = [z < z_pos.min() + 4.0 for z in z_pos]
    slab.set_constraint(FixAtoms(mask=mask))

    slab.calc = calculator
    opt = BFGS(slab, logfile=None)
    opt.run(fmax=0.02)

    E_slab = slab.get_potential_energy()
    N = len(slab)
    A = np.linalg.norm(np.cross(slab.cell[0], slab.cell[1]))

    gamma = (E_slab - N * E_bulk_per_atom) / (2 * A)
    surface_energies.append(gamma * 1000)  # meV/Å²

    print(f"{n} layers: γ = {gamma*1000:.2f} meV/Å²")
```

**Convergence:** |γ(n) - γ(n-2)| < 1 meV/Å²

## 4. Vacuum Thickness (Surfaces)

```python
vacuum_sizes = [8, 10, 12, 15, 20]  # Å

for vac in vacuum_sizes:
    slab = fcc111('Cu', size=(3, 3, 7), vacuum=vac)
    # ... calculate surface energy ...
```

**Convergence:** < 0.5 meV/Å²

Especially important for:
- Polar surfaces
- Charged systems
- Work function calculations

## 5. Supercell Size (Defects, Adsorption)

```python
sizes = [(2,2), (3,3), (4,4), (5,5)]

for nx, ny in sizes:
    slab = fcc111('Cu', size=(nx, ny, 4), vacuum=10)
    # Add adsorbate
    # Calculate adsorption energy
```

**Convergence:** ΔE_ads < 0.05 eV

Ensure:
- Adsorbates don't interact
- Defects isolated
- Strain fields contained

## 6. SCF Convergence

Density/wavefunction self-consistency:

```python
calc.set(
    convergence={'energy': 1e-6,  # eV
                 'density': 1e-5}
)
```

Tighter for:
- Forces/stresses
- Response properties
- Magnetic systems

## 7. Optimizer Convergence

```python
# Geometry
fmax_values = [0.1, 0.05, 0.02, 0.01, 0.005]

# Test impact on final properties
for fmax in fmax_values:
    opt = BFGS(atoms, logfile=None)
    opt.run(fmax=fmax)

    E = atoms.get_potential_energy()
    # Calculate property of interest
```

**Standard:**
- Structural properties: fmax = 0.02 eV/Å
- Vibrational analysis: fmax = 0.01 eV/Å
- High accuracy: fmax = 0.005 eV/Å

## 8. Smearing/Broadening

For metals:

```python
smearings = [0.05, 0.1, 0.2, 0.3]  # eV

for sigma in smearings:
    calc.set(occupations={'name': 'fermi-dirac',
                          'width': sigma})
    E = calculate_energy()
```

**Guidelines:**
- Smaller better for accuracy
- Larger helps SCF convergence
- Typical: 0.05-0.2 eV
- Extrapolate to T=0 if needed

## 9. Exchange-Correlation Functional

Test sensitivity:

```python
functionals = ['LDA', 'PBE', 'RPBE', 'PBEsol']

results = {}
for xc in functionals:
    calc = set_calculator(xc=xc)
    property = calculate_property()
    results[xc] = property
```

Different functionals for different properties:
- Lattice constants: PBEsol, LDA
- Molecules: PBE, RPBE
- Surfaces: BEEF-vdW (with vdW)

## Best Practices

### Documentation Template

```python
# Convergence Parameters:
# k-points: 8×8×8 (converged to < 1 meV/atom)
# Cutoff: 500 eV (converged to < 1 meV/atom)
# Smearing: 0.1 eV Fermi-Dirac
# SCF tolerance: 1e-6 eV
# fmax: 0.02 eV/Å
# XC: PBE
# Pseudopotential: PAW
```

### Convergence Plot

```python
import matplotlib.pyplot as plt

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 4))

# k-point convergence
ax1.plot(k_densities, energies, 'o-')
ax1.axhline(y=energies[-1], linestyle='--', color='r',
            label='Converged')
ax1.set_xlabel('k-point density')
ax1.set_ylabel('Energy (eV/atom)')
ax1.legend()
ax1.grid(True, alpha=0.3)

# Cutoff convergence
ax2.plot(cutoffs, energies_cutoff, 'o-')
ax2.set_xlabel('Cutoff (eV)')
ax2.set_ylabel('Energy (eV/atom)')
ax2.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig('convergence_tests.png', dpi=150)
```

### Automated Testing

```python
def test_convergence(atoms, parameter, values, tolerance=0.001):
    """
    Automated convergence testing.

    Args:
        atoms: ASE Atoms object
        parameter: 'kpts', 'encut', etc.
        values: list of values to test
        tolerance: convergence criterion (eV/atom)

    Returns:
        converged_value, results
    """
    results = []

    for val in values:
        atoms_test = atoms.copy()
        atoms_test.calc = setup_calculator(**{parameter: val})

        opt = BFGS(atoms_test, logfile=None)
        opt.run(fmax=0.01)

        E = atoms_test.get_potential_energy() / len(atoms_test)
        results.append((val, E))

        if len(results) >= 2:
            diff = abs(results[-1][1] - results[-2][1])
            if diff < tolerance:
                print(f"Converged at {parameter}={val}")
                return val, results

    print(f"Warning: Not converged to {tolerance} eV/atom")
    return values[-1], results
```

## Common Pitfalls

1. **Testing on different structures:** Always use same atomic positions
2. **Changing multiple parameters:** Only vary one at a time
3. **Insufficient range:** Test beyond convergence point
4. **Ignoring forces:** Energy may converge before forces
5. **System-dependent:** Results for one material don't transfer

## Quick Reference

| Property | k-points | Cutoff | fmax | Slab | Vacuum |
|----------|----------|--------|------|------|--------|
| Bulk | 8×8×8 | 500 eV | 0.02 | - | - |
| Surface | 4×4×1 | 500 eV | 0.02 | 7-9 layers | 10-15 Å |
| Adsorption | 4×4×1 | 500 eV | 0.02 | 4 layers | 10 Å |
| Defects | 2×2×2 | 500 eV | 0.02 | - | - |
| Vibrations | 8×8×8 | 500 eV | 0.001 | - | - |

*(Values are starting points, always test convergence!)*

## References

- Lejaeghere et al., "Reproducibility in DFT codes," Science 351, aad3000 (2016)
- Hjorth Larsen et al., "The atomic simulation environment," J. Phys.: Condens. Matter 29, 273002 (2017)
