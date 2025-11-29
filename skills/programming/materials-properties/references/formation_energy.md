# Formation Energy Calculations

## Overview

Formation energy measures the thermodynamic stability of a compound relative to its constituent elements.

## Definition

```
ΔH_f = E(compound) - Σᵢ nᵢ μᵢ(element)
```

Where:
- E(compound): DFT energy per formula unit
- nᵢ: Stoichiometric coefficient
- μᵢ: Chemical potential (energy per atom in elemental reference)

## Calculation Procedure

```python
from ase.build import bulk

# Example: Cu₃Au formation energy

# 1. Compound
cu3au = create_L12_structure()  # Ordered Cu₃Au
cu3au.calc = calculator
opt = BFGS(cu3au)
opt.run(fmax=0.01)
E_compound = cu3au.get_potential_energy()
n_atoms = len(cu3au)

# 2. Elemental references
cu_bulk = bulk('Cu', 'fcc', a=3.6)
cu_bulk.calc = calculator
E_Cu = cu_bulk.get_potential_energy() / len(cu_bulk)

au_bulk = bulk('Au', 'fcc', a=4.08)
au_bulk.calc = calculator
E_Au = au_bulk.get_potential_energy() / len(au_bulk)

# 3. Formation energy per atom
n_Cu = 3
n_Au = 1
E_f_per_atom = (E_compound - n_Cu*E_Cu - n_Au*E_Au) / n_atoms

print(f"ΔH_f = {E_f_per_atom:.3f} eV/atom")
```

## Interpretation

- ΔH_f < 0: Compound stable relative to elements
- ΔH_f > 0: Metastable or unstable
- |ΔH_f| ~ 0.1 eV/atom: moderate stability

## Convex Hull Analysis

For alloy systems, determine stable phases:

```python
# Calculate formation energies at different compositions
compositions = [0, 0.25, 0.5, 0.75, 1.0]
formation_energies = []

for x in compositions:
    E_f = calculate_formation_energy(x)
    formation_energies.append(E_f)

# Find convex hull
# Phases on hull are stable
# Above hull: driving force for decomposition
```

## Applications

- Phase diagram prediction
- Alloy design
- Catalysis (oxide stability)
- Battery materials

## References

- Curtarolo et al., "AFLOW: Materials Design from First-Principles," Comp. Mat. Sci. 58, 218 (2012)
- Jain et al., "Formation enthalpies by mixing GGA and GGA+U," Phys. Rev. B 84, 045115 (2011)
