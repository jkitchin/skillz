# Equation of State (EOS)

## Overview

The equation of state relates energy, volume, and pressure. Used to determine:
- Equilibrium lattice constants
- Bulk modulus
- Pressure derivatives
- Phase transitions

## EOS Forms

### Birch-Murnaghan (most common)

```
E(V) = E₀ + (B₀V₀)/(B₀') [(V₀/V)^B₀' / (B₀'-1) + 1] - B₀V₀/(B₀'-1)
```

Parameters: E₀, V₀, B₀, B₀' (4 parameters, 3rd order)

### Murnaghan

```
E(V) = E₀ + B₀V/B₀' [(V₀/V)^B₀' / (B₀'-1) + 1] - B₀V₀/(B₀'-1)
```

### Vinet

```
E(V) = E₀ + 2B₀V₀/(B₀'-1)² {2 - [5 + 3B₀'(x-1) - 3x] exp[-3(B₀'-1)(x-1)/2]}
```

where x = (V/V₀)^(1/3)

## Calculation Procedure

```python
from ase.eos import EquationOfState
from ase.build import bulk
import numpy as np

# Scan volumes
atoms = bulk('Cu', 'fcc', a=3.6)
cell0 = atoms.get_cell()

volumes = []
energies = []

for x in np.linspace(0.95, 1.05, 11):
    atoms_scaled = atoms.copy()
    atoms_scaled.set_cell(cell0 * x, scale_atoms=True)
    atoms_scaled.calc = calculator

    E = atoms_scaled.get_potential_energy()
    V = atoms_scaled.get_volume()

    volumes.append(V)
    energies.append(E)

# Fit EOS
eos = EquationOfState(volumes, energies, eos='birchmurnaghan')
v0, e0, B = eos.fit()

# Plot
eos.plot('eos.png')

print(f"V₀ = {v0:.3f} Å³")
print(f"E₀ = {e0:.4f} eV")
print(f"B = {B/1e9:.2f} GPa")
```

## Typical B₀' Values

- Most metals: 4-6
- Ionic crystals: 4-5
- Covalent: 3-4

## References

- Birch, "Finite Elastic Strain of Cubic Crystals," Phys. Rev. 71, 809 (1947)
- Vinet et al., "Universal equation of state," J. Phys.: Condens. Matter 1, 1941 (1989)
