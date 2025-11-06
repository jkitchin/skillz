# ASE Quick Reference

## Structure Building Cheat Sheet

```python
from ase import Atoms
from ase.build import bulk, molecule, surface, fcc111, bcc110, hcp0001

# Bulk structures
atoms = bulk('Al', 'fcc', a=4.05)
atoms = bulk('Fe', 'bcc', a=2.87)
atoms = bulk('Mg', 'hcp', a=3.21, c=5.21)
atoms = bulk('Si', 'diamond', a=5.43)

# Molecules (built-in database)
h2o = molecule('H2O')
co2 = molecule('CO2')
ch4 = molecule('CH4')

# Surfaces (metal, indices, size, vacuum)
slab = fcc111('Pt', size=(4, 4, 4), vacuum=10.0)
slab = bcc110('Fe', size=(3, 3, 4), vacuum=10.0)
slab = hcp0001('Mg', size=(4, 4, 4), vacuum=10.0)

# Manual construction
atoms = Atoms(
    'CO2',
    positions=[
        [0, 0, 0],      # C
        [-1.16, 0, 0],  # O
        [1.16, 0, 0]    # O
    ]
)
atoms.center(vacuum=10.0)
atoms.pbc = [True, True, True]
```

## Calculator Setup Cheat Sheet

```python
# EMT (for testing)
from ase.calculators.emt import EMT
atoms.calc = EMT()

# VASP - minimal
from ase.calculators.vasp import Vasp
calc = Vasp(xc='PBE', encut=400, kpts=(4,4,4))

# VASP - geometry optimization
calc = Vasp(
    xc='PBE', encut=400, kpts=(4,4,4),
    ibrion=2, nsw=100, ediff=1e-5,
    ismear=1, sigma=0.1
)

# VASP - static calculation
calc = Vasp(
    xc='PBE', encut=400, kpts=(8,8,8),
    ibrion=-1, nsw=0, ismear=-5
)

# GPAW
from gpaw import GPAW, PW
calc = GPAW(mode=PW(400), xc='PBE', kpts=(4,4,4))
```

## Optimization Cheat Sheet

```python
from ase.optimize import BFGS, LBFGS, FIRE

# Basic optimization
opt = BFGS(atoms, trajectory='opt.traj', logfile='opt.log')
opt.run(fmax=0.05)  # forces < 0.05 eV/Å

# With fixed atoms
from ase.constraints import FixAtoms
c = FixAtoms(indices=[0, 1, 2])  # fix first 3 atoms
atoms.set_constraint(c)

# Fix by position (e.g., bottom 2 layers)
mask = [atom.z < 10.0 for atom in atoms]
c = FixAtoms(mask=mask)
atoms.set_constraint(c)
```

## I/O Cheat Sheet

```python
from ase.io import read, write

# Read files
atoms = read('POSCAR')           # VASP
atoms = read('structure.cif')    # CIF
atoms = read('mol.xyz')          # XYZ
atoms = read('opt.traj', -1)     # last frame from trajectory
traj = read('md.traj', ':')      # all frames

# Write files
write('POSCAR', atoms)           # VASP
write('structure.cif', atoms)    # CIF
write('out.xyz', atoms)          # XYZ
write('snapshot.png', atoms)     # PNG image

# Trajectory writing during optimization
from ase.io.trajectory import Trajectory
traj = Trajectory('opt.traj', 'w', atoms)
opt.attach(traj.write, interval=1)
```

## Property Calculation Cheat Sheet

```python
# Single point calculations
energy = atoms.get_potential_energy()  # eV
forces = atoms.get_forces()            # eV/Å, shape (n_atoms, 3)
stress = atoms.get_stress()            # eV/Å³

# Distances and angles
d = atoms.get_distance(0, 1)                    # distance between atoms 0,1
angle = atoms.get_angle(0, 1, 2)                # angle 0-1-2 in degrees
dihedral = atoms.get_dihedral(0, 1, 2, 3)       # dihedral angle

# Cell properties
volume = atoms.get_volume()            # Å³
cell = atoms.get_cell()                # 3x3 array
```

## Constraint Types

```python
from ase.constraints import (
    FixAtoms,           # Fix atoms completely
    FixCartesian,       # Fix specific Cartesian components
    FixBondLength,      # Fix distance between two atoms
    FixLinearTriatomic, # Fix angle in molecule
    FixedPlane,         # Constrain to plane
)

# Fix specific atoms
c = FixAtoms(indices=[0, 1, 2])

# Fix only z-coordinate
c = FixCartesian(0, mask=[False, False, True])

# Fix bond length
c = FixBondLength(0, 1)

# Multiple constraints
atoms.set_constraint([c1, c2, c3])
```

## Molecular Dynamics Cheat Sheet

```python
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution
from ase.md.langevin import Langevin
from ase.md.verlet import VelocityVerlet
from ase import units

# Initialize velocities
MaxwellBoltzmannDistribution(atoms, temperature_K=300)

# NVE dynamics
dyn = VelocityVerlet(atoms, timestep=1.0*units.fs)

# NVT with Langevin
dyn = Langevin(
    atoms,
    timestep=1.0*units.fs,
    temperature_K=300,
    friction=0.002
)

# Run with trajectory
from ase.io.trajectory import Trajectory
traj = Trajectory('md.traj', 'w', atoms)
dyn.attach(traj.write, interval=10)
dyn.run(1000)  # 1000 steps
```

## k-points Cheat Sheet

```python
# Gamma point only
kpts=(1, 1, 1)

# Regular grid
kpts=(4, 4, 4)      # for bulk
kpts=(4, 4, 1)      # for slabs (no periodicity in z)

# Monkhorst-Pack with offset
kpts={'size': (4, 4, 4), 'gamma': True}

# From density
kpts={'density': 3.5}  # 3.5 k-points per Å⁻¹

# Band structure path
from ase.dft.kpoints import bandpath
path = atoms.cell.bandpath('GXWKG', npoints=100)
```

## Common Patterns

### Convergence Test (k-points)
```python
energies = []
for k in range(2, 12, 2):
    atoms = bulk('Si', 'diamond', a=5.43)
    atoms.calc = EMT()  # replace with real calculator
    calc.set(kpts=(k, k, k))
    energies.append(atoms.get_potential_energy())
```

### Convergence Test (ENCUT)
```python
energies = []
for encut in range(300, 600, 50):
    atoms = bulk('Si', 'diamond', a=5.43)
    calc = Vasp(xc='PBE', encut=encut, kpts=(8,8,8))
    atoms.calc = calc
    energies.append(atoms.get_potential_energy())
```

### Equation of State
```python
from ase.eos import calculate_eos

atoms = bulk('Cu', 'fcc', a=3.6)
atoms.calc = EMT()
eos = calculate_eos(atoms)
v, e, B = eos.fit()  # optimal volume, energy, bulk modulus
print(f"Bulk modulus: {B / ase.units.GPa:.1f} GPa")
```

### Energy vs. Distance Scan
```python
from ase import Atoms
import numpy as np

distances = np.linspace(0.5, 3.0, 20)
energies = []
for d in distances:
    mol = Atoms('H2', positions=[(0,0,0), (d,0,0)])
    mol.center(vacuum=10.0)
    mol.calc = EMT()
    energies.append(mol.get_potential_energy())
```

## Units Conversion

```python
from ase import units

# Energy
1 * units.eV    # = 1.0 (default)
1 * units.Ry    # Rydberg
1 * units.Ha    # Hartree

# Length
1 * units.Ang   # Ångström (default)
1 * units.Bohr  # Bohr radius

# Time
1 * units.fs    # femtosecond

# Pressure
1 * units.GPa   # Gigapascal
1 * units.eV / units.Ang**3  # eV/Å³

# Temperature
units.kB        # Boltzmann constant in eV/K
```

## Typical Convergence Criteria

| Property | k-points | ENCUT | fmax | ediff |
|----------|----------|-------|------|-------|
| Geometry | 2-4 / Å⁻¹ | Medium | 0.05 eV/Å | 1e-5 eV |
| Energy (total) | 4-6 / Å⁻¹ | High | - | 1e-6 eV |
| DOS/Bands | 6-8 / Å⁻¹ | High | - | 1e-6 eV |
| Forces | 3-5 / Å⁻¹ | High | - | 1e-6 eV |

## Troubleshooting

```python
# Check if atoms are too close
from ase.geometry import get_distances
dists = get_distances(atoms.positions)[1]
min_dist = dists[dists > 0].min()
if min_dist < 1.0:
    print("Warning: atoms too close!")

# Check if vacuum is sufficient
cell_z = atoms.cell[2, 2]
max_z = atoms.positions[:, 2].max()
min_z = atoms.positions[:, 2].min()
vacuum = cell_z - (max_z - min_z)
print(f"Vacuum: {vacuum:.2f} Å")

# Verify PBC
print(f"PBC: {atoms.pbc}")  # should be [T, T, F] for slabs

# Check calculator is attached
if atoms.calc is None:
    raise RuntimeError("No calculator attached!")
```
