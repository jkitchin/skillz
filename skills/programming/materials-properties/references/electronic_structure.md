# Electronic Structure Analysis

## Overview

Electronic structure provides insight into bonding, optical properties, and transport. Key quantities:
- Band structure
- Density of states (DOS)
- Band gap
- Charge density
- Work function

## Band Structure

Maps energy eigenvalues along high-symmetry paths in k-space.

```python
from ase.dft.kpoints import bandpath

# Optimize structure first
atoms.calc = calculator
opt = BFGS(atoms)
opt.run(fmax=0.01)

# Define k-point path
lat = atoms.cell.get_bravais_lattice()
path = lat.bandpath(npoints=100)

# Calculate band structure
atoms.calc.set(kpts=path.kpts)
energies = atoms.calc.get_eigenvalues(kpt=range(len(path.kpts)))

# Plot
path.plot(emin=-10, emax=10)
```

### High-Symmetry Points

**FCC:** Γ-X-W-K-Γ-L
**BCC:** Γ-H-N-Γ-P-H
**Hexagonal:** Γ-M-K-Γ-A-L-H-A

## Density of States (DOS)

Energy histogram of states:
```
DOS(E) = Σ_nk δ(E - ε_nk)
```

```python
# Calculate DOS
from ase.dft.dos import DOS

dos = DOS(atoms.calc, width=0.1)
energies = dos.get_energies()
dos_values = dos.get_dos()

# Plot
import matplotlib.pyplot as plt
plt.plot(energies, dos_values)
plt.xlabel('Energy (eV)')
plt.ylabel('DOS (states/eV)')
```

### Projected DOS (PDOS)

Project onto atomic orbitals:
- s, p, d, f character
- Atom-resolved
- Useful for bonding analysis

## Band Gap

```
E_g = E_CBM - E_VBM
```

**Direct:** CBM and VBM at same k-point
**Indirect:** Different k-points

```python
# Get band gap
calc = atoms.calc
homo, lumo = calc.get_homo_lumo()
gap = lumo - homo

print(f"Band gap: {gap:.3f} eV")
```

### DFT Band Gap Problem

GGA/LDA underestimate gaps significantly:
- Si: exp ~1.1 eV, DFT ~0.6 eV
- GaAs: exp ~1.4 eV, DFT ~0.5 eV

**Solutions:**
- Hybrid functionals (HSE, PBE0)
- GW approximation
- Scissor operator (shift CBM up)

## Charge Density

```
n(r) = Σ_nk |ψ_nk(r)|²
```

Visualization:
```python
# Get charge density
rho = atoms.calc.get_pseudo_density()

# Write to cube file for visualization
from ase.io.cube import write_cube
write_cube('density.cube', atoms, rho)
```

View with VESTA, VMD, or similar.

## Work Function

For surfaces:
```
Φ = V_vacuum - E_F
```

```python
# Calculate work function
from ase.dft.dos import DOS

# Get Fermi level
E_F = atoms.calc.get_fermi_level()

# Get electrostatic potential
V = atoms.calc.get_electrostatic_potential()

# Extract vacuum potential
z_pos = np.linspace(0, atoms.cell[2,2], len(V))
V_vacuum = V[z_pos > atoms.positions[:,2].max() + 5].mean()

Phi = V_vacuum - E_F
print(f"Work function: {Phi:.2f} eV")
```

## Effective Mass

From band curvature:
```
m* = ℏ² / (∂²E/∂k²)
```

Extract from fitted band structure near extrema.

## Bader Charge Analysis

Partition charge density into atomic contributions:

```python
# Requires bader executable
from ase.calculators.bader import Bader

bader = Bader(atoms, calc=atoms.calc)
charges = bader.get_charges()

for i, q in enumerate(charges):
    print(f"Atom {i}: {q:.3f} e")
```

## Applications

- Semiconductors: band engineering
- Metals: Fermi surface analysis
- Catalysis: d-band center theory
- Optical properties: absorption spectra

## Advanced Analysis

### d-Band Center

For transition metal catalysis:
```
ε_d = ∫ E × DOS_d(E) dE / ∫ DOS_d(E) dE
```

Correlates with adsorption energies.

### Band Alignment

At interfaces, determine:
- Valence band offset (VBO)
- Conduction band offset (CBO)
- Type I, II, or III alignment

## Tools Beyond ASE

- **VASPKIT:** Post-processing VASP
- **pymatgen:** Electronic structure analysis
- **sumo:** Publication-quality plots
- **Lobster:** Chemical bonding analysis

## References

- Ashcroft & Mermin, "Solid State Physics," Holt-Saunders (1976)
- Martin, "Electronic Structure: Basic Theory," Cambridge (2004)

## See Also

- `magnetic_properties.md`: Spin-polarized calculations
- `formation_energy.md`: Thermodynamic stability from energies
