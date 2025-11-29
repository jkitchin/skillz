# Surface Energy Calculations

## Overview

Surface energy (γ) quantifies the excess energy per unit area required to create a surface by cleaving a bulk crystal. It is fundamental for understanding:
- Crystal morphology and growth
- Catalytic activity
- Wetting and adhesion
- Nanoparticle stability

## Theoretical Background

### Definition

The surface energy is defined as the reversible work per unit area to create a surface:

```
γ = (E_slab - N × E_bulk) / (2 × A)
```

Where:
- `E_slab`: Total energy of the slab
- `N`: Number of atoms in the slab
- `E_bulk`: Energy per atom in bulk
- `A`: Surface area
- Factor of 2: accounts for two surfaces (top and bottom)

### Relationship to Surface Stress

For a relaxed surface, the surface energy is related to surface stress (f):

```
f_αβ = γ δ_αβ + ∂γ/∂ε_αβ
```

For liquids and many metals, surface energy ≈ surface stress.

## Computational Details

### 1. Slab Construction

**Thickness Requirements:**
- Converge with respect to number of layers
- Typically 5-15 atomic layers
- Must ensure bulk-like behavior in center

**Vacuum Requirements:**
- Typically 10-20 Å vacuum gap
- Prevents interaction between periodic images
- Converge carefully for charged or polar surfaces

**Symmetric vs Asymmetric Slabs:**
- Symmetric: same termination on both sides (preferred)
- Asymmetric: different terminations require dipole correction

### 2. Convergence Tests

**Slab Thickness:**
```python
n_layers = [5, 7, 9, 11, 13, 15]
for n in n_layers:
    slab = construct_slab(n_layers=n)
    gamma = calculate_surface_energy(slab)
```

Converged when |γ(n) - γ(n-2)| < 1 meV/Å²

**Vacuum Thickness:**
```python
vacuum_sizes = [8, 10, 12, 15, 20]
for vac in vacuum_sizes:
    slab = construct_slab(vacuum=vac)
    gamma = calculate_surface_energy(slab)
```

Converged when |γ(v) - γ(v-2)| < 0.5 meV/Å²

**k-point Sampling:**
- Dense in-plane sampling (surface directions)
- Single k-point in vacuum direction acceptable
- Test convergence as for bulk calculations

### 3. Relaxation Strategy

**Fixed Layers:**
- Fix bottom 1-2 layers to represent bulk
- Allows surface reconstruction

**Full Relaxation:**
- Relax all atoms for small slabs
- More accurate but computationally expensive

**Selective Optimization:**
```python
# Fix bottom 2 layers
z_positions = slab.get_positions()[:, 2]
z_min = z_positions.min()
mask = [z < z_min + 4.0 for z in z_positions]
slab.set_constraint(FixAtoms(mask=mask))
```

### 4. Surface Reconstruction

Some surfaces undergo reconstruction:
- Si(111): 7×7 reconstruction
- Au(110): 2×1 missing row
- Pt(100): hex reconstruction

Account for reconstruction in surface energy calculations.

## Typical Values

### FCC Metals (111) Surface:
- Cu: 1.2-1.8 J/m²
- Al: 0.9-1.2 J/m²
- Ag: 1.0-1.3 J/m²
- Au: 1.3-1.6 J/m²
- Pt: 2.1-2.5 J/m²

### Surface Energy Ordering (FCC):
```
γ(111) < γ(100) < γ(110)
```

Most close-packed surface (111) has lowest energy.

### BCC Metals (110) Surface:
- Fe: 2.2-2.5 J/m²
- W: 3.0-4.0 J/m²

## Advanced Topics

### 1. Chemical Potential Dependence

For compounds (e.g., oxides):
```
γ = [E_slab - Σᵢ nᵢ μᵢ] / (2A)
```

Surface energy depends on environmental conditions (O₂ pressure, temperature).

### 2. Temperature Effects

Include vibrational contributions via:
- Harmonic approximation (phonons)
- Ab initio molecular dynamics
- Quasiharmonic approximation

```
γ(T) = γ₀ + γ_vib(T)
```

### 3. Anisotropy

Wulff construction predicts equilibrium crystal shape:
- Minimize total surface energy
- Shape determined by γ-plot (polar plot of γ vs orientation)

### 4. Surface Phase Diagrams

For compounds, construct surface phase diagrams showing:
- Stable surface terminations vs chemical potential
- Surface reconstruction transitions
- Adsorbate coverage

## Practical Tips

1. **Always Check Convergence:**
   - Slab thickness (most critical)
   - Vacuum thickness
   - k-points
   - Energy cutoff

2. **Use Symmetric Slabs:**
   - Avoid spurious dipole moments
   - Easier convergence
   - More physical

3. **Monitor Surface Relaxation:**
   - Calculate interlayer spacings
   - Typical: top layer contracts 1-5%
   - Compare with experiments (LEED)

4. **Energy Cutoff:**
   - Use same cutoff as bulk calculations
   - Higher cutoff for accurate forces
   - Test convergence

5. **Calculator Settings:**
   - Spin polarization for magnetic materials
   - Dispersion corrections for layered materials
   - DFT+U for correlated systems

## Example Workflow

```python
from ase.build import bulk, fcc111
from ase.optimize import BFGS
from ase.constraints import FixAtoms
import numpy as np

# 1. Get bulk energy
bulk_atoms = bulk('Cu', 'fcc', a=3.6)
bulk_atoms.calc = calculator
E_bulk_per_atom = bulk_atoms.get_potential_energy() / len(bulk_atoms)

# 2. Create slab
slab = fcc111('Cu', size=(3, 3, 9), vacuum=12.0)

# 3. Fix bottom layers
z_pos = slab.get_positions()[:, 2]
z_min = z_pos.min()
mask = [z < z_min + 4.0 for z in z_pos]
slab.set_constraint(FixAtoms(mask=mask))

# 4. Relax
slab.calc = calculator
opt = BFGS(slab)
opt.run(fmax=0.02)

# 5. Calculate surface energy
E_slab = slab.get_potential_energy()
N = len(slab)
cell = slab.get_cell()
A = np.linalg.norm(np.cross(cell[0], cell[1]))

gamma = (E_slab - N * E_bulk_per_atom) / (2 * A)

# Convert units
gamma_mJ_m2 = gamma * 16.0217  # eV/Å² to mJ/m²

print(f"Surface energy: {gamma_mJ_m2:.2f} mJ/m²")
```

## Common Pitfalls

1. **Insufficient Slab Thickness:**
   - Center layers must behave like bulk
   - Check interlayer spacings

2. **Vacuum Too Small:**
   - Leads to artificial interaction
   - Especially problematic for polar surfaces

3. **Inconsistent Settings:**
   - Use exact same DFT settings as bulk
   - Same XC functional, k-points density, cutoff

4. **Neglecting Reconstruction:**
   - Some surfaces have lower-energy reconstructed phases
   - Test different surface models

5. **Dipole Corrections:**
   - Necessary for asymmetric slabs
   - Most codes have built-in corrections

## References

1. **Fundamental Theory:**
   - Fiorentini & Methfessel, "Extracting convergent surface energies from slab calculations," J. Phys.: Condens. Matter 8, 6525 (1996)
   - Boettger, "Nonconvergence of surface energies obtained from thin-film calculations," Phys. Rev. B 49, 16798 (1994)

2. **Surface Science:**
   - Vitos et al., "The surface energy of metals," Surf. Sci. 411, 186 (1998)
   - Tran et al., "Surface energies of elemental crystals," Sci. Data 3, 160080 (2016)

3. **Computational Methods:**
   - Kresse & Joubert, "From ultrasoft pseudopotentials to the projector augmented-wave method," Phys. Rev. B 59, 1758 (1999)

4. **Reviews:**
   - Mehl et al., "The AFLOW library of crystallographic prototypes," Comput. Mater. Sci. 136, S1 (2017)

## See Also

- `adsorption_energy.md`: Molecular adsorption on surfaces
- `interface_properties.md`: Interfaces between materials
- `examples/surface_energy.py`: Working examples with convergence tests
