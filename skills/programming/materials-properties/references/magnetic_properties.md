# Magnetic Properties

## Overview

Calculate magnetic moments, exchange interactions, and magnetic ordering from first principles.

## Magnetic Moment

```python
from ase.build import bulk

# Enable spin polarization
atoms = bulk('Fe', 'bcc', a=2.87)
atoms.set_initial_magnetic_moments([2.0])  # Initial guess

atoms.calc = calculator  # Must support spin polarization
atoms.calc.set(spinpol=True)  # Calculator-specific

opt = BFGS(atoms)
opt.run(fmax=0.01)

# Get final magnetic moment
magmom = atoms.get_magnetic_moment()
print(f"Magnetic moment: {magmom:.2f} μ_B/atom")
```

## Typical Moments (μ_B/atom)

- Fe (BCC): ~2.2
- Co (HCP): ~1.6
- Ni (FCC): ~0.6
- Gd: ~7.6 (4f electrons)

## Magnetic Configurations

### Ferromagnetic (FM)

All spins parallel:
```python
magmoms = [2.0] * len(atoms)  # All positive
atoms.set_initial_magnetic_moments(magmoms)
```

### Antiferromagnetic (AFM)

Alternating spins:
```python
# For checkerboard AFM
magmoms = [2.0 if i % 2 == 0 else -2.0
           for i in range(len(atoms))]
atoms.set_initial_magnetic_moments(magmoms)
```

### Noncollinear

Requires special DFT treatment (VASP: LNONCOLLINEAR=.TRUE.)

## Exchange Interactions

Heisenberg model:
```
H = -Σᵢⱼ Jᵢⱼ Sᵢ·Sⱼ
```

Calculate J from energy differences:
```python
E_FM = calculate_energy(ferromagnetic)
E_AFM = calculate_energy(antiferromagnetic)

# For nearest-neighbor J only
J = (E_AFM - E_FM) / (2 * z * S²)
```

Where z = coordination number, S = spin.

## Curie/Néel Temperature

Mean-field estimate:
```
T_C = (2/3) z J S(S+1) / k_B
```

More accurate: Monte Carlo with Heisenberg model.

## DFT+U for Correlated Systems

For transition metal oxides, f-electron systems:
```python
# VASP example
calc.set(
    ldau=True,
    ldautype=2,  # Dudarev
    ldauu={'Fe': 4.0},  # U value in eV
    ldaul={'Fe': 2}  # l quantum number (d orbitals)
)
```

## Spin-Orbit Coupling (SOC)

Important for:
- Heavy elements
- Magnetic anisotropy
- Topological materials

Enable in calculator (expensive!).

## Magnetic Anisotropy Energy (MAE)

```
MAE = E[100] - E[001]
```

Energy difference between magnetization directions.

## Applications

- Permanent magnets
- Magnetic storage
- Spintronics
- Multiferroics

## References

- Sandratskii, "Noncollinear magnetism in itinerant systems," Adv. Phys. 47, 91 (1998)
- Eriksson et al., "Atomistic Spin Dynamics," Oxford (2017)

## See Also

- `electronic_structure.md`: Density of states, band structures
- DFT+U: Hubbard corrections for localized electrons
