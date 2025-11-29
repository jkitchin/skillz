# Materials Properties Quick Reference

## Basic Setup

```python
from ase import Atoms
from ase.build import bulk, surface, molecule
from ase.optimize import BFGS, FIRE
from ase.calculators.emt import EMT
from ase.constraints import FixAtoms, ExpCellFilter
import spglib
```

## Structure Creation

### Bulk Crystals
```python
# FCC structure
atoms = bulk('Cu', 'fcc', a=3.6)

# BCC structure
atoms = bulk('Fe', 'bcc', a=2.87)

# Custom structure
atoms = Atoms('NaCl',
              positions=[[0,0,0], [0.5,0.5,0.5]],
              cell=[5.64, 5.64, 5.64],
              pbc=True)
```

### Surfaces
```python
from ase.build import fcc111, fcc100, bcc110

# FCC (111)
slab = fcc111('Cu', size=(3,3,4), vacuum=10)

# FCC (100)
slab = fcc100('Pt', size=(2,2,5), vacuum=15)

# Generic surface
from ase.build import surface
slab = surface('Cu', (1,1,1), layers=7)
```

### Molecules
```python
mol = molecule('H2O')
mol = molecule('CO2')
mol = molecule('CH4')
```

## Calculators

### EMT (Fast Testing)
```python
atoms.calc = EMT()
```

### GPAW
```python
from gpaw import GPAW, PW

atoms.calc = GPAW(mode=PW(500),
                 xc='PBE',
                 kpts=(8,8,8),
                 txt='gpaw.txt')
```

### VASP
```python
from ase.calculators.vasp import Vasp

atoms.calc = Vasp(xc='PBE',
                 encut=500,
                 kpts=(8,8,8))
```

## Structure Relaxation

### Geometry Optimization
```python
opt = BFGS(atoms, trajectory='opt.traj')
opt.run(fmax=0.01)
```

### Cell Optimization
```python
ecf = ExpCellFilter(atoms)
opt = BFGS(ecf)
opt.run(fmax=0.01)
```

### With Constraints
```python
# Fix atoms
c = FixAtoms(indices=[0,1,2])
atoms.set_constraint(c)

# Fix bottom layers
mask = [atom.position[2] < 5.0 for atom in atoms]
c = FixAtoms(mask=mask)
atoms.set_constraint(c)
```

## Property Calculations

### Energy
```python
E = atoms.get_potential_energy()  # eV
```

### Forces
```python
F = atoms.get_forces()  # eV/Å
max_force = max(np.linalg.norm(F, axis=1))
```

### Stress
```python
stress = atoms.get_stress()  # eV/ų
```

### Lattice Parameters
```python
a, b, c, alpha, beta, gamma = atoms.cell.cellpar()
V = atoms.get_volume()  # ų
```

### Space Group
```python
import spglib

cell = (atoms.cell, atoms.get_scaled_positions(), atoms.get_atomic_numbers())
spacegroup = spglib.get_spacegroup(cell, symprec=1e-5)
dataset = spglib.get_symmetry_dataset(cell)
```

## Common Formulas

### Surface Energy
```python
gamma = (E_slab - N * E_bulk_per_atom) / (2 * A)
# Units: eV/ų
```

### Adsorption Energy
```python
E_ads = E_slab_plus_ads - E_slab - E_molecule
# Negative = exothermic
```

### Formation Energy
```python
E_form = E_compound - sum(n_i * mu_i)
# mu_i = chemical potentials (reference states)
```

### Cohesive Energy
```python
E_coh = E_atom - E_bulk_per_atom
# Positive for stable solids
```

## Equation of State

```python
from ase.eos import calculate_eos

eos = calculate_eos(atoms, trajectory='eos.traj')
v, e, B = eos.fit()  # Volume, energy, bulk modulus
eos.plot('eos.png')

# Manual
from ase.eos import EquationOfState
volumes = []
energies = []
for x in np.linspace(0.95, 1.05, 11):
    atoms_scaled = atoms.copy()
    atoms_scaled.set_cell(atoms.cell * x, scale_atoms=True)
    energies.append(atoms_scaled.get_potential_energy())
    volumes.append(atoms_scaled.get_volume())

eos = EquationOfState(volumes, energies)
v0, e0, B = eos.fit()
```

## Surface Energy

```python
from ase.build import bulk, fcc111

# Bulk
bulk_atoms = bulk('Cu', 'fcc', a=3.6)
bulk_atoms.calc = calc
E_bulk = bulk_atoms.get_potential_energy() / len(bulk_atoms)

# Slab
slab = fcc111('Cu', size=(3,3,7), vacuum=10)
slab.calc = calc
E_slab = slab.get_potential_energy()

# Surface energy
N = len(slab)
A = slab.get_cell()[0,0] * slab.get_cell()[1,1]
gamma = (E_slab - N * E_bulk) / (2 * A)
```

## Adsorption Energy

```python
from ase.build import fcc111, molecule, add_adsorbate

# Clean slab
slab = fcc111('Pt', size=(3,3,4), vacuum=10)
slab.calc = calc
E_slab = slab.get_potential_energy()

# Molecule
mol = molecule('CO')
mol.calc = calc
E_mol = mol.get_potential_energy()

# Adsorb
add_adsorbate(slab, mol, height=2.0, position='ontop')
slab.calc = calc
E_total = slab.get_potential_energy()

# Adsorption energy
E_ads = E_total - E_slab - E_mol
```

## Nudged Elastic Band

```python
from ase.neb import NEB
from ase.optimize import BFGS

# Create image chain
images = [initial]
images += [initial.copy() for i in range(5)]
images += [final]

# Interpolate
neb = NEB(images)
neb.interpolate()

# Set calculators
for img in images[1:-1]:
    img.calc = calc

# Optimize
opt = BFGS(neb, trajectory='neb.traj')
opt.run(fmax=0.05)

# Get barrier
from ase.neb import NEBTools
nebtools = NEBTools(images)
Ef, dE = nebtools.get_barrier()
print(f"Barrier: {Ef:.2f} eV")
```

## Vibrations

```python
from ase.vibrations import Vibrations

vib = Vibrations(atoms, indices=[0,1])  # Atoms to vibrate
vib.run()
vib.summary()

# Zero-point energy
zpe = vib.get_zero_point_energy()

# Frequencies
freqs = vib.get_frequencies()
```

## Phonopy

```python
from phonopy import Phonopy

# Create phonopy
phonon = Phonopy(atoms, [[2,0,0],[0,2,0],[0,0,2]])

# Displacements
phonon.generate_displacements(distance=0.01)
supercells = phonon.supercells_with_displacements

# Calculate forces (loop over supercells)
# ... set forces back to phonopy

# Properties
phonon.produce_force_constants()
phonon.auto_band_structure()
phonon.plot_band_structure()
```

## Convergence Testing

### K-points
```python
for k in [2, 4, 6, 8, 10, 12]:
    atoms.calc.set(kpts=(k,k,k))
    E = atoms.get_potential_energy()
    print(f"k={k}: E={E:.4f} eV")
```

### Cutoff
```python
for ecut in [300, 400, 500, 600]:
    atoms.calc.set(encut=ecut)
    E = atoms.get_potential_energy()
    print(f"Ecut={ecut}: E={E:.4f} eV")
```

### Slab Thickness
```python
for n_layers in [5, 7, 9, 11]:
    slab = fcc111('Cu', size=(3,3,n_layers), vacuum=10)
    slab.calc = calc
    # Calculate surface energy
```

## Trajectory Analysis

```python
from ase.io import read

# Read trajectory
traj = read('opt.traj', ':')

# Energies
energies = [a.get_potential_energy() for a in traj]

# Forces
forces = [a.get_forces() for a in traj]
max_forces = [np.max(np.linalg.norm(f, axis=1)) for f in forces]

# Plot
import matplotlib.pyplot as plt
plt.plot(energies)
plt.xlabel('Step')
plt.ylabel('Energy (eV)')
plt.show()
```

## File I/O

### Read Structures
```python
from ase.io import read

atoms = read('structure.cif')
atoms = read('POSCAR')
atoms = read('structure.xyz')
traj = read('opt.traj', ':')  # All images
```

### Write Structures
```python
from ase.io import write

write('structure.cif', atoms)
write('POSCAR', atoms, format='vasp')
write('structure.xyz', atoms)
write('trajectory.traj', traj)
```

## Visualization

### View Structure
```python
from ase.visualize import view
view(atoms)
```

### Plot
```python
import matplotlib.pyplot as plt

# Energy vs volume
plt.plot(volumes, energies, 'o-')
plt.xlabel('Volume (ų)')
plt.ylabel('Energy (eV)')

# Forces
plt.plot(max_forces)
plt.xlabel('Step')
plt.ylabel('Max Force (eV/Å)')
```

## Utility Functions

### Get Neighbors
```python
from ase.neighborlist import NeighborList

nl = NeighborList([1.5]*len(atoms), self_interaction=False, bothways=True)
nl.update(atoms)
indices, offsets = nl.get_neighbors(0)  # Neighbors of atom 0
```

### Distance
```python
dist = atoms.get_distance(0, 1, mic=True)  # mic = minimum image convention
```

### Angle
```python
angle = atoms.get_angle(0, 1, 2)  # degrees
```

### Dihedral
```python
dihedral = atoms.get_dihedral(0, 1, 2, 3)  # degrees
```

## Common Optimizers

```python
from ase.optimize import BFGS, LBFGS, FIRE, GPMin

opt = BFGS(atoms, trajectory='bfgs.traj')
opt = LBFGS(atoms)
opt = FIRE(atoms)
opt = GPMin(atoms)

# All use same interface
opt.run(fmax=0.01, steps=100)
```

## Parallel Calculations

```python
# GPAW parallel
from gpaw import GPAW
calc = GPAW(..., parallel={'domain': 1, 'band': 4})

# Run with mpirun
# mpirun -np 4 python script.py
```

## Key Units

- Energy: eV
- Length: Å (Angstrom)
- Force: eV/Å
- Stress: eV/ų
- Pressure: eV/ų (to GPa: multiply by 160.2176)
- Time: fs (ASE dynamics)

## Common Pitfalls

1. **Forgot calculator**: `atoms.calc = calc`
2. **Wrong units**: ASE uses eV and Å
3. **Periodic boundaries**: Set `pbc=[True, True, True]`
4. **Constraints not set**: Use `atoms.set_constraint()`
5. **Not converged**: Check `fmax` and forces
6. **No vacuum**: Surfaces need vacuum layer
7. **Overlapping atoms**: Check initial geometry
8. **Wrong cell**: Use `atoms.set_cell()` carefully

## Tips

- Start with EMT for testing
- Always check convergence
- Save trajectories
- Use constraints wisely
- Visualize before calculating
- Keep a log file
- Version control your scripts
