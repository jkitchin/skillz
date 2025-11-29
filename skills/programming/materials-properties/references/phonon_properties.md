# Phonon Properties with Phonopy

## Overview

Phonon calculations using phonopy provide comprehensive lattice dynamics information for crystalline solids:
- Phonon band structures
- Phonon density of states (DOS)
- Thermal properties (heat capacity, entropy, free energy)
- Thermal expansion
- Thermodynamic stability

## Theoretical Background

### Lattice Dynamics

The phonon eigenvalue problem:

```
D(q) · e_λ(q) = ω²_λ(q) e_λ(q)
```

Where:
- `D(q)`: Dynamical matrix at wavevector q
- `ω_λ(q)`: Phonon frequency for branch λ at q
- `e_λ(q)`: Phonon eigenvector (polarization)

### Force Constant Method

**Interatomic Force Constants:**
```
Φ_αβ(ll') = ∂²V / ∂u_α(l) ∂u_β(l')
```

Calculate by finite displacements in supercell.

**Fourier Transform to q-space:**
```
D_αβ(q) = (1/√(m_α m_β)) Σ_l' Φ_αβ(ll') exp(iq·r_ll')
```

### Thermodynamic Properties

**Vibrational Free Energy:**
```
F_vib = kT Σ_qλ ln[2 sinh(ℏω_qλ / 2kT)]
```

**Entropy:**
```
S = k Σ_qλ [ℏω/(2kT) coth(ℏω/2kT) - ln(2 sinh(ℏω/2kT))]
```

**Heat Capacity:**
```
C_v = k Σ_qλ (ℏω/kT)² exp(ℏω/kT) / [exp(ℏω/kT) - 1]²
```

## Computational Workflow

### 1. Setup with Phonopy

**Installation:**
```bash
pip install phonopy
```

**Integration with ASE:**
```python
from phonopy import Phonopy
from phonopy.structure.atoms import PhonopyAtoms

def ase_to_phonopy(atoms):
    """Convert ASE atoms to Phonopy format."""
    return PhonopyAtoms(
        symbols=atoms.get_chemical_symbols(),
        cell=atoms.get_cell(),
        positions=atoms.get_positions()
    )
```

### 2. Supercell Generation

**Create Phonopy Object:**
```python
from ase.build import bulk

# Primitive cell
atoms = bulk('Si', 'diamond', a=5.43)
phonopy_atoms = ase_to_phonopy(atoms)

# Phonopy with supercell
phonon = Phonopy(
    phonopy_atoms,
    supercell_matrix=[[2, 0, 0],
                      [0, 2, 0],
                      [0, 0, 2]]
)

# Generate displacements
phonon.generate_displacements(distance=0.01)
supercells = phonon.supercells_with_displacements
```

**Supercell Size:**
- Must capture force constant range
- Typical: 2×2×2 to 4×4×4
- Test convergence

### 3. Force Calculations

**For each displaced supercell:**
```python
from ase import Atoms as ASEAtoms

forces = []
for scell in supercells:
    # Convert to ASE
    ase_scell = ASEAtoms(
        numbers=scell.numbers,
        positions=scell.positions,
        cell=scell.cell,
        pbc=True
    )
    ase_scell.calc = calculator

    # Get forces
    f = ase_scell.get_forces()
    forces.append(f)

# Set forces in phonopy
phonon.forces = forces
```

### 4. Force Constants

**Produce force constants:**
```python
phonon.produce_force_constants()
```

**Acoustic Sum Rule:**
Automatically enforced by phonopy to ensure:
- ω(q→0) = 0 for acoustic branches
- Crystal momentum conservation

## Property Calculations

### 1. Phonon Band Structure

**Define high-symmetry path:**
```python
# For FCC: Γ-X-W-K-Γ-L
path = [
    [[0.0, 0.0, 0.0], [0.5, 0.0, 0.5]],  # Γ-X
    [[0.5, 0.0, 0.5], [0.5, 0.25, 0.75]],  # X-W
    [[0.5, 0.25, 0.75], [0.375, 0.375, 0.75]],  # W-K
    [[0.375, 0.375, 0.75], [0.0, 0.0, 0.0]],  # K-Γ
    [[0.0, 0.0, 0.0], [0.5, 0.5, 0.5]],  # Γ-L
]

labels = ['$\\Gamma$', 'X', 'W', 'K', '$\\Gamma$', 'L']

# Calculate and plot
phonon.run_band_structure(path, labels=labels)
phonon.plot_band_structure()
```

**Analysis:**
- Identify acoustic vs optical branches
- Check for imaginary frequencies (instability)
- Compare with experimental data (INS, IXS)

### 2. Phonon Density of States

**Calculate DOS:**
```python
# Mesh sampling of Brillouin zone
phonon.run_mesh([20, 20, 20])
phonon.run_total_dos()

# Plot
phonon.plot_total_dos()
```

**Partial DOS:**
```python
phonon.run_mesh([20, 20, 20])
phonon.run_projected_dos()

# Plot contributions from each atom type
```

**Mesh Convergence:**
- Test DOS with increasing mesh density
- Typically 20×20×20 to 50×50×50
- Converged when DOS shape stable

### 3. Thermal Properties

**Calculate over temperature range:**
```python
phonon.run_mesh([20, 20, 20])
phonon.run_thermal_properties(
    t_min=0,
    t_max=1000,
    t_step=10
)

# Get results
tp_dict = phonon.get_thermal_properties_dict()

temperatures = tp_dict['temperatures']
free_energy = tp_dict['free_energy']
entropy = tp_dict['entropy']
heat_capacity = tp_dict['heat_capacity']
```

**Plotting:**
```python
import matplotlib.pyplot as plt

fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2)

ax1.plot(temperatures, free_energy)
ax1.set_ylabel('F (kJ/mol)')

ax2.plot(temperatures, entropy)
ax2.set_ylabel('S (J/K/mol)')

ax3.plot(temperatures, heat_capacity)
ax3.set_ylabel('C_v (J/K/mol)')
ax3.axhline(y=3*8.314, linestyle='--', label='Dulong-Petit')

plt.tight_layout()
```

## Convergence Testing

### 1. Supercell Size

```python
supercell_sizes = [[2,2,2], [3,3,3], [4,4,4]]

for size in supercell_sizes:
    phonon = Phonopy(structure, supercell_matrix=size)
    # ... calculate DOS ...
```

Converged when DOS and thermal properties stable.

### 2. Displacement Magnitude

```python
displacements = [0.005, 0.01, 0.015, 0.02]

for disp in displacements:
    phonon.generate_displacements(distance=disp)
    # ... calculate forces ...
```

Default 0.01 Å usually good. Check for:
- Anharmonic effects (large disp)
- Numerical noise (small disp)

### 3. DFT Convergence

**k-point Sampling:**
- Denser than bulk (larger supercell)
- Automatic scaled k-points work
- Check force convergence

**Energy Cutoff:**
- Same as bulk calculations
- Forces more sensitive than energies

**SCF Convergence:**
- Tight thresholds for clean forces
- Especially important for supercells

## Advanced Features

### 1. Non-Analytical Corrections (NAC)

For polar materials (LO-TO splitting):

```python
# Born effective charges and dielectric tensor
# (From DFPT or finite differences)

phonon.nac_params = {
    'born': born_charges,
    'dielectric': epsilon,
    'factor': 14.399652  # Conversion factor
}
```

Essential for:
- Ionic crystals
- Perovskites
- Polar semiconductors

### 2. Group Velocity

```python
phonon.run_mesh([20, 20, 20], is_mesh_symmetry=False)
group_velocity = phonon.get_group_velocity()
```

Applications:
- Thermal conductivity
- Sound velocity
- Phonon transport

### 3. Grüneisen Parameter

Measure of anharmonicity:

```
γ_qλ = -V/ω ∂ω/∂V
```

Calculate via:
- Multiple volumes
- Quasiharmonic approximation
- Direct from strain derivatives

### 4. Mode Grüneisen Parameters

```python
# Calculate at multiple volumes
phonon_V1 = calculate_phonopy(volume=V1)
phonon_V2 = calculate_phonopy(volume=V2)

# Extract mode-by-mode Grüneisen
gruneisen = phonon.get_mesh_gruneisen_parameters()
```

## Quasiharmonic Approximation (QHA)

**Thermal Expansion:**
Calculate phonons at multiple volumes:

```python
from phonopy import PhonopyQHA

# Multiple volumes
volumes = [V1, V2, V3, V4, V5]
free_energies = []
entropies = []
heat_capacities = []

for V in volumes:
    atoms_V = optimize_at_volume(V)
    phonon_V = run_phonopy(atoms_V)

    tp = phonon_V.get_thermal_properties()
    free_energies.append(tp['free_energy'])
    entropies.append(tp['entropy'])
    heat_capacities.append(tp['heat_capacity'])

# QHA analysis
phonopy_qha = PhonopyQHA(
    volumes,
    electronic_energies,
    temperatures,
    free_energies,
    entropies,
    heat_capacities
)

# Get thermal expansion
alpha = phonopy_qha.get_thermal_expansion()
```

## Typical Results

### Debye Temperature (K):
- Al: ~430
- Cu: ~343
- Si: ~645
- Diamond: ~2230

### Heat Capacity:
- Low T: C_v ∝ T³ (Debye model)
- High T: C_v → 3Nk (Dulong-Petit)

### Phonon Frequencies:
- Acoustic branches: 0 at Γ point
- Optical branches: finite at Γ
- Max frequencies: material-dependent (0-1000 cm⁻¹ typical)

## Practical Tips

1. **Optimize Geometry First:**
   - Tight convergence (fmax < 0.01 eV/Å)
   - Converged cell parameters
   - Symmetry preserved

2. **Test Supercell:**
   - Start with 2×2×2
   - Check DOS convergence
   - Increase if needed

3. **Displacement Size:**
   - 0.01 Å standard
   - Smaller for hard materials
   - Check linearity (harmonic regime)

4. **Symmetry:**
   - Phonopy uses crystal symmetry
   - Reduces number of calculations
   - Preserve symmetry in geometry

5. **Save Results:**
   - Force constants: phonon.save('phonopy.yaml')
   - Can reload for analysis
   - Avoid recalculation

## Validation

**Compare with Experiments:**
- INS (Inelastic Neutron Scattering)
- IXS (Inelastic X-ray Scattering)
- Raman/IR spectroscopy
- Heat capacity measurements

**Typical DFT Accuracy:**
- GGA: ±2-5% for frequencies
- Better for relative trends
- Mode assignment usually reliable

## Common Pitfalls

1. **Small Supercell:**
   - Artificial periodicity
   - Wrong force constants
   - Always test convergence

2. **Poor Geometry:**
   - Loose optimization → imaginary modes
   - Broken symmetry → too many calculations
   - Inconsistent settings

3. **Coarse Mesh:**
   - DOS not smooth
   - Thermal properties inaccurate
   - Increase mesh density

4. **Ignoring NAC:**
   - Wrong LO-TO splitting
   - Incorrect Γ-point frequencies
   - Critical for polar materials

5. **Inconsistent DFT:**
   - Different settings for forces
   - Changed pseudopotentials
   - Non-converged calculations

## Example Complete Workflow

```python
from ase.build import bulk
from ase.optimize import BFGS
from ase.constraints import ExpCellFilter
from phonopy import Phonopy
from phonopy.structure.atoms import PhonopyAtoms

# 1. Optimize structure
atoms = bulk('Cu', 'fcc', a=3.6)
atoms.calc = calculator

ecf = ExpCellFilter(atoms)
opt = BFGS(ecf)
opt.run(fmax=0.001)  # Tight!

# 2. Setup phonopy
pho_atoms = PhonopyAtoms(
    symbols=atoms.get_chemical_symbols(),
    cell=atoms.get_cell(),
    positions=atoms.get_positions()
)

phonon = Phonopy(pho_atoms,
                 supercell_matrix=[[2, 0, 0],
                                   [0, 2, 0],
                                   [0, 0, 2]])

phonon.generate_displacements(distance=0.01)

# 3. Calculate forces
forces = []
for scell in phonon.supercells_with_displacements:
    from ase import Atoms as ASEAtoms

    ase_scell = ASEAtoms(
        numbers=scell.numbers,
        positions=scell.positions,
        cell=scell.cell,
        pbc=True
    )
    ase_scell.calc = calculator
    forces.append(ase_scell.get_forces())

phonon.forces = forces

# 4. Produce force constants
phonon.produce_force_constants()

# 5. Band structure
path = [[[0, 0, 0], [0.5, 0, 0.5]]]
labels = ['$\\Gamma$', 'X']
phonon.run_band_structure(path, labels=labels)
phonon.plot_band_structure()

# 6. DOS
phonon.run_mesh([20, 20, 20])
phonon.run_total_dos()
phonon.plot_total_dos()

# 7. Thermal properties
phonon.run_thermal_properties(t_min=0, t_max=1000, t_step=10)
tp = phonon.get_thermal_properties_dict()

# Analysis
print(f"Heat capacity at 300 K: {tp['heat_capacity'][30]:.2f} J/K/mol")
```

## References

1. **Phonopy:**
   - Togo & Tanaka, "First principles phonon calculations in materials science," Scr. Mater. 108, 1 (2015)
   - Documentation: https://phonopy.github.io/phonopy/

2. **Lattice Dynamics:**
   - Born & Huang, "Dynamical Theory of Crystal Lattices," Oxford (1954)
   - Dove, "Introduction to Lattice Dynamics," Cambridge (1993)

3. **Computational Methods:**
   - Parlinski, Li & Kawazoe, "First-Principles Determination of the Soft Mode in Cubic ZrO₂," Phys. Rev. Lett. 78, 4063 (1997)

4. **Thermal Properties:**
   - Togo, Chaput & Tanaka, "Distributions of phonon lifetimes in Brillouin zones," Phys. Rev. B 91, 094306 (2015)

5. **QHA:**
   - Baroni, de Gironcoli & Dal Corso, "Phonons and related crystal properties from density-functional perturbation theory," Rev. Mod. Phys. 73, 515 (2001)

## See Also

- `vibrational_modes.md`: Molecular vibrations with ASE
- `thermal_transport.md`: Phonon transport and conductivity
- `elastic_constants.md`: Mechanical properties
- `examples/phonon_calculation.py`: Complete working examples
