# Pymatgen Core Objects

Guide to fundamental pymatgen classes: Element, Lattice, Structure, Molecule, Composition.

## Element and Species

### Element Class

**Basic usage:**
```python
from pymatgen.core import Element

# Create element
si = Element("Si")
fe = Element("Fe")

# Properties
print(si.atomic_mass)      # 28.0855
print(si.atomic_number)    # 14
print(si.symbol)           # Si
print(si.melting_point)    # 1687.0 K
print(si.boiling_point)    # 3538.0 K
print(si.X)                # Electronegativity
print(si.atomic_radius)    # Atomic radius
```

**Common properties:**
- `atomic_number`, `atomic_mass`, `symbol`
- `group`, `row` (periodic table position)
- `electronic_structure` (electron configuration)
- `common_oxidation_states`
- `melting_point`, `boiling_point`
- `X` (Pauling electronegativity)
- `ionization_energies`

### Species Class

**For ions and oxidation states:**
```python
from pymatgen.core import Species

# Create species with oxidation state
fe2 = Species("Fe", 2)   # Fe²⁺
fe3 = Species("Fe", 3)   # Fe³⁺

print(fe2.oxi_state)     # 2
print(fe2.ionic_radius)  # Ionic radius for Fe²⁺

# From string
o2minus = Species.from_string("O2-")  # O²⁻
```

## Lattice Class

### Creating Lattices

**From lattice parameters:**
```python
from pymatgen.core import Lattice

# Cubic
lattice = Lattice.cubic(5.0)  # 5 Å cubic cell

# Hexagonal
lattice = Lattice.hexagonal(a=3.0, c=5.0)

# From parameters (a, b, c, alpha, beta, gamma)
lattice = Lattice.from_parameters(
    a=4.0, b=4.0, c=6.0,
    alpha=90, beta=90, gamma=120
)
```

**From lattice vectors:**
```python
# Direct specification of lattice vectors
lattice = Lattice([
    [4.0, 0.0, 0.0],  # a vector
    [0.0, 4.0, 0.0],  # b vector
    [0.0, 0.0, 6.0]   # c vector
])

# Or as matrix
import numpy as np
matrix = np.array([[4, 0, 0], [0, 4, 0], [0, 0, 6]])
lattice = Lattice(matrix)
```

### Lattice Properties

```python
print(lattice.a, lattice.b, lattice.c)           # Lengths
print(lattice.alpha, lattice.beta, lattice.gamma) # Angles
print(lattice.volume)                             # Volume
print(lattice.matrix)                             # 3x3 matrix
print(lattice.reciprocal_lattice)                 # Reciprocal lattice
```

### Lattice Operations

```python
# Get fractional coordinates from Cartesian
cart_coords = [2.0, 2.0, 3.0]
frac_coords = lattice.get_fractional_coords(cart_coords)

# Get Cartesian from fractional
cart_coords = lattice.get_cartesian_coords([0.5, 0.5, 0.5])

# Distance in PBC
d = lattice.get_distance_and_image([0, 0, 0], [0.6, 0.6, 0.6])
```

## Structure Class

### Creating Structures

**Basic creation:**
```python
from pymatgen.core import Structure, Lattice

# Simple cubic structure
lattice = Lattice.cubic(4.2)
species = ["Si", "Si"]
coords = [[0, 0, 0], [0.25, 0.25, 0.25]]  # Fractional coords

structure = Structure(lattice, species, coords)
```

**With oxidation states:**
```python
from pymatgen.core import Species

species = [Species("Fe", 2), Species("O", -2)]
coords = [[0, 0, 0], [0.5, 0.5, 0.5]]
structure = Structure(lattice, species, coords)
```

**Site properties:**
```python
# Add properties to sites
structure = Structure(
    lattice,
    ["Fe", "O"],
    [[0, 0, 0], [0.5, 0.5, 0.5]],
    site_properties={"magmom": [5.0, 0.0]}  # Magnetic moments
)
```

### Structure Properties

```python
# Composition
print(structure.formula)                        # Full formula
print(structure.composition)                    # Composition object
print(structure.composition.reduced_formula)    # Reduced formula
print(structure.composition.alphabetical_formula)

# Geometry
print(structure.volume)          # Volume in Ų
print(structure.density)         # Density in g/cm³
print(structure.num_sites)       # Number of sites
print(len(structure))            # Also number of sites

# Access sites
print(structure[0])              # First site
print(structure[0].species_string)  # Element at site
print(structure[0].frac_coords)  # Fractional coordinates
print(structure[0].coords)       # Cartesian coordinates
```

### Structure Operations

**Adding/removing sites:**
```python
# Add site
structure.append("O", [0.5, 0.5, 0.5])

# Remove site
structure.remove_sites([0])  # Remove first site

# Replace species
structure.replace_species({"Li": "Na"})  # All Li → Na
structure[0] = "F"  # Change first site to F
```

**Supercells:**
```python
# Simple scaling
supercell = structure * (2, 2, 1)  # 2x2x1 supercell

# Matrix specification
structure.make_supercell([[2, 0, 0], [0, 2, 0], [0, 0, 1]])
```

**Perturbations:**
```python
# Small random perturbations
structure.perturb(0.1)  # 0.1 Å random displacements
```

**Sorting:**
```python
# Sort by electronegativity
structure.sort()

# Or by custom key
structure.sort(key=lambda site: site.species_string)
```

## Molecule Class

**For non-periodic systems:**
```python
from pymatgen.core import Molecule

# Create molecule
species = ["C", "O"]
coords = [[0, 0, 0], [1.2, 0, 0]]  # Cartesian coords
molecule = Molecule(species, coords)

# Properties
print(molecule.formula)
print(molecule.center_of_mass)
print(molecule.principal_axes)

# From file
molecule = Molecule.from_file("molecule.xyz")
```

**Molecule operations:**
```python
# Translate
molecule.translate_sites([0, 1], [1.0, 0, 0])

# Rotate
molecule.rotate_sites([0, 1], theta=90, axis=[0, 0, 1])

# Get bonds
bonds = molecule.get_covalent_bonds()
```

## Composition Class

**Creating compositions:**
```python
from pymatgen.core import Composition

# From formula
comp = Composition("Fe2O3")

# From dict
comp = Composition({"Fe": 2, "O": 3})

# Properties
print(comp.reduced_formula)       # Fe2O3
print(comp.alphabetical_formula)  # Fe2 O3
print(comp.weight)                # Molecular weight
print(comp.num_atoms)             # Total atoms
print(comp.elements)              # List of Element objects
```

**Composition operations:**
```python
# Arithmetic
comp1 = Composition("Fe2O3")
comp2 = Composition("Al2O3")

# Addition
comp3 = comp1 + comp2  # Fe2Al2O6

# Get amounts
print(comp1["Fe"])  # 2.0
print(comp1["O"])   # 3.0

# Fractional composition
frac = comp1.fractional_composition
print(frac)  # Fe0.4O0.6

# Get element fraction
print(comp1.get_atomic_fraction("Fe"))  # 0.4
```

## Site Class

**Direct site creation (uncommon):**
```python
from pymatgen.core import Site, PeriodicSite, Lattice

# Non-periodic site
site = Site("Fe", [1.0, 2.0, 3.0])

# Periodic site (with lattice)
lattice = Lattice.cubic(5.0)
psite = PeriodicSite("Fe", [0.5, 0.5, 0.5], lattice)

# Properties
print(psite.species_string)  # "Fe"
print(psite.frac_coords)     # [0.5, 0.5, 0.5]
print(psite.coords)          # Cartesian coordinates
```

**Usually accessed from Structure:**
```python
structure = Structure.from_file("POSCAR")
for site in structure:
    print(site.species_string, site.frac_coords)
    if hasattr(site, "magmom"):
        print(f"  Magnetic moment: {site.magmom}")
```

## Common Workflows

### Building Complex Structures

```python
from pymatgen.core import Structure, Lattice

# Start with lattice
lattice = Lattice.from_parameters(
    a=5.0, b=5.0, c=5.0,
    alpha=90, beta=90, gamma=90
)

# Build structure incrementally
structure = Structure(lattice, [], [])
structure.append("Fe", [0, 0, 0])
structure.append("Fe", [0.5, 0.5, 0.5])
structure.append("O", [0.25, 0.25, 0.25])
structure.append("O", [0.75, 0.75, 0.75])

print(structure)
```

### From Dictionary (Serialization)

```python
# Structure to dict
structure_dict = structure.as_dict()

# Save to JSON
import json
with open("structure.json", "w") as f:
    json.dump(structure_dict, f)

# Load from JSON
with open("structure.json", "r") as f:
    loaded_dict = json.load(f)

structure = Structure.from_dict(loaded_dict)
```

### Converting Molecule ↔ Structure

```python
# Molecule to Structure (add lattice)
from pymatgen.core import Molecule, Structure, Lattice

molecule = Molecule.from_file("molecule.xyz")

# Put in box
lattice = Lattice.cubic(20.0)  # Large box
structure = Structure(
    lattice,
    molecule.species,
    molecule.cart_coords,
    coords_are_cartesian=True
)

# Structure to Molecule (remove periodicity)
structure = Structure.from_file("POSCAR")
molecule = structure.to_molecule()  # Only if structure is molecular
```

## Best Practices

**Element vs Species:**
- Use `Element` for neutral atoms
- Use `Species` when oxidation states matter
- Pymatgen often auto-converts strings

**Coordinates:**
- Structure uses fractional coordinates by default
- Use `coords_are_cartesian=True` if providing Cartesian
- Molecule always uses Cartesian coordinates

**Serialization:**
- Use `as_dict()` / `from_dict()` for persistence
- Prefer over Python pickle
- Survives code changes better

**Site Properties:**
- Add custom properties via `site_properties` dict
- Useful for magnetic moments, selective dynamics, etc.

**Immutability:**
- Most operations return new objects
- Use `.copy()` to create independent copies
- In-place operations are clearly named (`.make_supercell()`, etc.)
