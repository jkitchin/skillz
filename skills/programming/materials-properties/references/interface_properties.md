# Interface Properties

## Overview

Interfaces between materials control properties in:
- Heterostructures
- Grain boundaries
- Multilayers
- Catalytic supports

## Interface Energy

```
γ_interface = (E_interface - E_slab1 - E_slab2) / A
```

Similar to surface energy but for two materials.

## Calculation Procedure

```python
from ase.build import fcc111

# Create interface
slab1 = fcc111('Cu', size=(3, 3, 6), vacuum=0)
slab2 = fcc111('Ni', size=(3, 3, 6), vacuum=0)

# Stack with minimal vacuum
interface = stack_slabs(slab1, slab2, gap=2.0)

# Fix bottom of each slab
z_pos = interface.get_positions()[:, 2]
n1 = len(slab1)
mask1 = [i < n1 and z_pos[i] < ... for i in range(len(interface))]
mask2 = [i >= n1 and z_pos[i] < ... for i in range(len(interface))]
interface.set_constraint(FixAtoms(mask=mask1+mask2))

# Optimize
interface.calc = calculator
opt = BFGS(interface)
opt.run(fmax=0.02)

E_interface = interface.get_potential_energy()
γ = (E_interface - E_slab1 - E_slab2) / A
```

## Lattice Mismatch

```
ε = (a₂ - a₁) / a₁
```

- Small ε: Coherent interface
- Large ε: Dislocations form

## Work of Adhesion

```
W_ad = γ_1 + γ_2 - γ_interface
```

Positive W_ad indicates stable interface.

## Grain Boundaries

Special case: same material, different orientations.

**Σ notation:** Ratio of grain boundary unit cell to bulk.
- Σ3: Twin boundary (low energy)
- Σ5, Σ7: Low-angle boundaries

## Applications

- Thin film growth
- Composite materials
- Semiconductor devices
- Battery interfaces (solid electrolyte)

## References

- Sutton & Balluffi, "Interfaces in Crystalline Materials," Oxford (1995)
- Schusteritsch et al., "First-principles study of interface properties," Phys. Rev. B 94, 035426 (2016)

## See Also

- `surface_energy.md`: Single surface calculations
- `adsorption_energy.md`: Molecule-surface interfaces
