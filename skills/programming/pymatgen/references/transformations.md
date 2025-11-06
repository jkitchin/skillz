# Pymatgen Structure Transformations

Comprehensive guide to transforming, modifying, and manipulating crystal structures.

## Overview

Pymatgen provides two main approaches for structure transformations:
1. **Direct methods** - Simple operations on Structure objects
2. **Transformation classes** - Reproducible, composable transformations with history tracking

## Direct Structure Modifications

### Supercells

**Simple scaling:**
```python
from pymatgen.core import Structure

structure = Structure.from_file("POSCAR")

# Create 2x2x1 supercell
supercell = structure * (2, 2, 1)

print(f"Original sites: {len(structure)}")
print(f"Supercell sites: {len(supercell)}")

# Save supercell
supercell.to(filename="POSCAR_supercell")
```

**Matrix specification:**
```python
# More complex supercell with transformation matrix
import numpy as np

# [[2, 0, 0], [0, 2, 0], [0, 0, 1]]
scaling_matrix = [[2, 0, 0],
                  [0, 2, 0],
                  [0, 0, 1]]

structure.make_supercell(scaling_matrix)  # In-place

# Or non-destructive
supercell = structure.copy()
supercell.make_supercell(scaling_matrix)
```

**Non-diagonal supercells:**
```python
# Create supercell with non-diagonal matrix
# Useful for specific k-point sampling or interfaces
scaling_matrix = [[2, 1, 0],
                  [0, 2, 0],
                  [0, 0, 1]]

structure.make_supercell(scaling_matrix)
```

### Adding and Removing Sites

**Add sites:**
```python
# Add site by fractional coordinates
structure.append("O", [0.5, 0.5, 0.5])

# Add with Cartesian coordinates
structure.append("O", [2.5, 2.5, 2.5], coords_are_cartesian=True)

# Add with properties
structure.append(
    "Fe",
    [0.25, 0.25, 0.25],
    properties={"magmom": 5.0}
)
```

**Remove sites:**
```python
# Remove by index
structure.remove_sites([0])

# Remove multiple sites
structure.remove_sites([0, 2, 5])

# Remove by species
indices_to_remove = [
    i for i, site in enumerate(structure)
    if site.species_string == "H"
]
structure.remove_sites(indices_to_remove)
```

**Replace species:**
```python
# Replace all Li with Na
structure.replace_species({"Li": "Na"})

# Replace with oxidation states
from pymatgen.core import Species
structure.replace_species({
    "Fe": Species("Fe", 2),
    "O": Species("O", -2)
})

# Replace specific site
structure[0] = "F"  # Change first site to F
```

### Perturbations

**Random atomic displacements:**
```python
# Small random perturbations (useful for breaking symmetry)
structure.perturb(distance=0.1)  # 0.1 Å random displacement

# Larger perturbations
structure.perturb(distance=0.5)
```

**Controlled perturbations:**
```python
# Displace specific sites
import numpy as np

for i in range(len(structure)):
    site = structure[i]
    # Add small displacement
    displacement = np.random.rand(3) * 0.1  # Å
    new_coords = site.coords + displacement
    structure.replace(i, site.species_string, new_coords, coords_are_cartesian=True)
```

### Translation and Rotation

**Translate structure:**
```python
# Translate all atoms
structure.translate_sites(
    indices=range(len(structure)),
    vector=[0.1, 0.1, 0.1]  # Fractional coordinates
)

# Translate specific atoms
structure.translate_sites(
    indices=[0, 1, 2],
    vector=[0.5, 0, 0]
)
```

**For molecules:**
```python
from pymatgen.core import Molecule

molecule = Molecule.from_file("molecule.xyz")

# Translate
molecule.translate_sites(range(len(molecule)), [1.0, 0, 0])

# Rotate
molecule.rotate_sites(
    indices=range(len(molecule)),
    theta=90,  # degrees
    axis=[0, 0, 1]  # z-axis
)
```

### Sorting

```python
# Sort by electronegativity (default)
structure.sort()

# Sort by custom key
structure.sort(key=lambda site: site.species_string)

# Sort by z-coordinate
structure.sort(key=lambda site: site.frac_coords[2])
```

## Transformation Classes

### Why Use Transformations?

Transformation classes provide:
- **Reproducibility** - Save and replay transformations
- **History tracking** - Know what was done
- **Composability** - Chain multiple transformations
- **Reversibility** - Some transformations can be undone

### Standard Transformations

**SupercellTransformation:**
```python
from pymatgen.transformations.standard_transformations import SupercellTransformation

# Create transformation
transformation = SupercellTransformation(scaling_matrix=[[2, 0, 0],
                                                          [0, 2, 0],
                                                          [0, 0, 1]])

# Apply to structure
supercell = transformation.apply_transformation(structure)

# Transformation is repeatable
another_supercell = transformation.apply_transformation(structure)
```

**SubstitutionTransformation:**
```python
from pymatgen.transformations.standard_transformations import SubstitutionTransformation

# Replace all Li with Na
transformation = SubstitutionTransformation({"Li": "Na"})
new_structure = transformation.apply_transformation(structure)

# Partial substitution
transformation = SubstitutionTransformation({"Li": {"Na": 0.5, "K": 0.5}})
doped_structure = transformation.apply_transformation(structure)
```

**RemoveSitesTransformation:**
```python
from pymatgen.transformations.standard_transformations import RemoveSitesTransformation

# Remove specific sites
transformation = RemoveSitesTransformation([0, 2, 5])
reduced_structure = transformation.apply_transformation(structure)
```

**PerturbStructureTransformation:**
```python
from pymatgen.transformations.standard_transformations import PerturbStructureTransformation

# Random perturbation
transformation = PerturbStructureTransformation(distance=0.1)
perturbed = transformation.apply_transformation(structure)
```

**RotationTransformation:**
```python
from pymatgen.transformations.standard_transformations import RotationTransformation

# Rotate structure
transformation = RotationTransformation(
    axis=[0, 0, 1],
    angle=45,  # degrees
    angle_in_radians=False
)
rotated = transformation.apply_transformation(structure)
```

### Advanced Transformations

**OrderDisorderedStructureTransformation:**
```python
from pymatgen.transformations.advanced_transformations import (
    OrderDisorderedStructureTransformation
)

# For structures with partial occupancies
# Creates ordered structure from disordered one
transformation = OrderDisorderedStructureTransformation()
ordered = transformation.apply_transformation(disordered_structure)
```

**EnumerateStructureTransformation:**
```python
from pymatgen.transformations.advanced_transformations import (
    EnumerateStructureTransformation
)

# Enumerate all possible orderings
transformation = EnumerateStructureTransformation(
    max_cell_size=4  # Maximum supercell size
)

# Returns multiple structures
structures = transformation.apply_transformation(structure)
```

**SlabTransformation:**
```python
from pymatgen.transformations.advanced_transformations import SlabTransformation

# Create surface slab
transformation = SlabTransformation(
    miller_index=[1, 1, 1],
    min_slab_size=10.0,  # Å
    min_vacuum_size=15.0,  # Å
    center_slab=True
)

slab = transformation.apply_transformation(structure)
```

**SubstitutionPredictorTransformation:**
```python
from pymatgen.transformations.advanced_transformations import (
    SubstitutionPredictorTransformation
)

# Suggest substitutions based on chemical similarity
transformation = SubstitutionPredictorTransformation(threshold=0.001)
suggestions = transformation.apply_transformation(structure)
```

### Chaining Transformations

```python
from pymatgen.transformations.standard_transformations import (
    SupercellTransformation,
    SubstitutionTransformation,
    PerturbStructureTransformation
)

# Create transformation sequence
transformations = [
    SupercellTransformation([[2, 0, 0], [0, 2, 0], [0, 0, 1]]),
    SubstitutionTransformation({"Li": "Na"}),
    PerturbStructureTransformation(distance=0.05)
]

# Apply in sequence
result = structure
for trans in transformations:
    result = trans.apply_transformation(result)
```

## Surface Slabs

### Creating Slabs

```python
from pymatgen.core.surface import SlabGenerator

structure = Structure.from_file("POSCAR")

# Create slab generator
slabgen = SlabGenerator(
    initial_structure=structure,
    miller_index=(1, 1, 1),
    min_slab_size=10.0,  # Å
    min_vacuum_size=15.0,  # Å
    center_slab=True,
    primitive=True
)

# Generate all possible slabs
all_slabs = slabgen.get_slabs()

print(f"Generated {len(all_slabs)} slab configurations")

# Get specific slab
slab = all_slabs[0]
slab.to(filename="POSCAR_slab")
```

### Slab Properties

```python
# Slab has additional properties
print(f"Miller index: {slab.miller_index}")
print(f"Slab thickness: {slab.get_slab_thickness():.2f} Å")
print(f"Surface area: {slab.surface_area:.2f} Ų")

# Identify surface sites
surface_sites = slab.get_surface_sites()
print(f"Number of surface atoms: {len(surface_sites)}")
```

### Symmetric and Asymmetric Slabs

```python
# Symmetric slab (same termination on both sides)
slabgen = SlabGenerator(
    structure,
    miller_index=(1, 0, 0),
    min_slab_size=10.0,
    min_vacuum_size=15.0,
    center_slab=True,
    primitive=True,
    lll_reduce=True,
    symmetrize=True  # Enforce symmetry
)

# Asymmetric slab
slabgen_asym = SlabGenerator(
    structure,
    miller_index=(1, 0, 0),
    min_slab_size=10.0,
    min_vacuum_size=15.0,
    symmetrize=False
)
```

### Adding Adsorbates

```python
from pymatgen.core.surface import Slab
from pymatgen.analysis.adsorption import AdsorbateSiteFinder

slab = all_slabs[0]

# Find adsorption sites
asf = AdsorbateSiteFinder(slab)

# Get sites
ads_sites = asf.find_adsorption_sites()

print("Adsorption sites:")
print(f"  Ontop: {len(ads_sites['ontop'])}")
print(f"  Bridge: {len(ads_sites['bridge'])}")
print(f"  Hollow: {len(ads_sites['hollow'])}")

# Add adsorbate at ontop site
from pymatgen.core import Molecule

# Create O adsorbate
adsorbate = Molecule(["O"], [[0, 0, 0]])

# Add to slab
slab_with_ads = asf.add_adsorbate(adsorbate, ads_sites['ontop'][0])
slab_with_ads.to(filename="POSCAR_with_adsorbate")
```

## Defects

### Vacancies

```python
from pymatgen.analysis.defects.generators import VacancyGenerator

structure = Structure.from_file("POSCAR")

# Generate all possible vacancies
vac_gen = VacancyGenerator()
vacancies = vac_gen.generate(structure)

print(f"Generated {len(vacancies)} vacancy configurations")

# Create structure with vacancy
vacancy_defect = vacancies[0]
vacancy_structure = vacancy_defect.get_supercell_structure(sc_size=[2, 2, 2])
vacancy_structure.to(filename="POSCAR_vacancy")
```

### Substitutional Defects

```python
from pymatgen.analysis.defects.generators import SubstitutionGenerator

# Generate substitutional defects
sub_gen = SubstitutionGenerator()
substitutions = sub_gen.generate(structure, {0: ["Na", "K"]})  # Substitute site 0

# Create structure with substitution
sub_defect = substitutions[0]
sub_structure = sub_defect.get_supercell_structure(sc_size=[2, 2, 2])
```

### Interstitials

```python
from pymatgen.analysis.defects.generators import InterstitialGenerator

# Generate interstitial sites
int_gen = InterstitialGenerator()
interstitials = int_gen.generate(structure, insertions=["Li"])

# Create structure with interstitial
int_defect = interstitials[0]
int_structure = int_defect.get_supercell_structure(sc_size=[2, 2, 2])
```

## Strain and Deformation

### Apply Strain

```python
from pymatgen.transformations.standard_transformations import DeformStructureTransformation

# Apply 2% tensile strain in a-direction
deformation = [[1.02, 0, 0],
               [0, 1.0, 0],
               [0, 0, 1.0]]

transformation = DeformStructureTransformation(deformation=deformation)
strained = transformation.apply_transformation(structure)

print(f"Original volume: {structure.volume:.2f} Ų")
print(f"Strained volume: {strained.volume:.2f} Ų")
```

### Uniaxial Strain

```python
# Uniaxial strain (one direction)
strain = 0.02  # 2%
deformation = [[1 + strain, 0, 0],
               [0, 1, 0],
               [0, 0, 1]]

transformation = DeformStructureTransformation(deformation=deformation)
strained = transformation.apply_transformation(structure)
```

### Hydrostatic Strain

```python
# Hydrostatic (volumetric) strain
strain = 0.02
factor = (1 + strain) ** (1/3)  # Preserve shape
deformation = [[factor, 0, 0],
               [0, factor, 0],
               [0, 0, factor]]

transformation = DeformStructureTransformation(deformation=deformation)
strained = transformation.apply_transformation(structure)
```

## Doping and Substitution

### Random Doping

```python
import random

structure = Structure.from_file("POSCAR")

# Find all sites of specific element
li_sites = [i for i, site in enumerate(structure) if site.species_string == "Li"]

# Replace random 10% with Na
n_to_replace = int(0.1 * len(li_sites))
sites_to_replace = random.sample(li_sites, n_to_replace)

for site_idx in sites_to_replace:
    structure[site_idx] = "Na"

structure.to(filename="POSCAR_doped")
```

### Systematic Doping

```python
# Create multiple doping concentrations
concentrations = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5]

for conc in concentrations:
    doped = structure.copy()

    li_sites = [i for i, site in enumerate(doped) if site.species_string == "Li"]
    n_to_replace = int(conc * len(li_sites))

    for i in range(n_to_replace):
        doped[li_sites[i]] = "Na"

    doped.to(filename=f"POSCAR_doped_{int(conc*100):02d}")
```

## Interpolation (NEB)

### Linear Interpolation

```python
from pymatgen.core import Structure

start = Structure.from_file("POSCAR_start")
end = Structure.from_file("POSCAR_end")

# Linear interpolation for NEB
n_images = 5

interpolated = start.interpolate(end, nimages=n_images)

# Save images
for i, image in enumerate(interpolated):
    image.to(filename=f"POSCAR_{i:02d}")
```

### With Autosort

```python
# Automatically match atoms between structures
interpolated = start.interpolate(
    end,
    nimages=5,
    autosort_tol=0.5  # Tolerance for matching atoms
)
```

## High-Throughput Workflows

### TransformedStructure

```python
from pymatgen.alchemy.materials import TransformedStructure
from pymatgen.transformations.standard_transformations import (
    SupercellTransformation,
    SubstitutionTransformation
)

# Track transformation history
ts = TransformedStructure(structure)

# Apply transformations
ts.append_transformation(SupercellTransformation([[2, 0, 0], [0, 2, 0], [0, 0, 1]]))
ts.append_transformation(SubstitutionTransformation({"Li": "Na"}))

# Access final structure
final_structure = ts.final_structure

# View history
print("Transformation history:")
for trans in ts.history:
    print(f"  - {trans}")

# Serialize with history
ts_dict = ts.as_dict()

# Reconstruct
ts_loaded = TransformedStructure.from_dict(ts_dict)
```

### Batch Transformations

```python
from pymatgen.transformations.standard_transformations import SubstitutionTransformation
import itertools

# Generate all binary combinations
elements_to_try = ["Li", "Na", "K"]
base_structure = Structure.from_file("POSCAR")

structures = []
for elem in elements_to_try:
    trans = SubstitutionTransformation({"Li": elem})
    new_struct = trans.apply_transformation(base_structure)
    structures.append(new_struct)
    new_struct.to(filename=f"POSCAR_{elem}")

print(f"Generated {len(structures)} structures")
```

## Grain Boundaries and Interfaces

### Grain Boundary Builder

```python
from pymatgen.analysis.gb.grain import GrainBoundaryGenerator

structure = Structure.from_file("POSCAR")

# Generate grain boundaries
gb_gen = GrainBoundaryGenerator(structure)

# Specific grain boundary
gb = gb_gen.gb_from_parameters(
    rotation_axis=[0, 0, 1],
    rotation_angle=36.87,  # degrees
    expand_times=2,
    vacuum_thickness=0.0,
    ab_shift=[0, 0],
    normal=True
)

gb.to(filename="POSCAR_GB")
```

### Interface Builder

```python
from pymatgen.core.interface import Interface

film = Structure.from_file("POSCAR_film")
substrate = Structure.from_file("POSCAR_substrate")

# Create interface
interface = Interface.from_slabs(
    slab1=film,
    slab2=substrate,
    in_layers=True,
    gap=2.0,  # Å
    vacuum_over_film=15.0
)

interface.to(filename="POSCAR_interface")
```

## Molecule Operations

### Align Molecule

```python
from pymatgen.core import Molecule

molecule = Molecule.from_file("molecule.xyz")

# Align principal axis
molecule.principal_axes  # Get principal axes

# Center at origin
center = molecule.center_of_mass
molecule.translate_sites(range(len(molecule)), -center)
```

### Merge Molecules

```python
# Combine two molecules
mol1 = Molecule.from_file("mol1.xyz")
mol2 = Molecule.from_file("mol2.xyz")

# Translate mol2
mol2.translate_sites(range(len(mol2)), [5.0, 0, 0])

# Merge
combined_species = list(mol1.species) + list(mol2.species)
combined_coords = mol1.cart_coords.tolist() + mol2.cart_coords.tolist()

combined = Molecule(combined_species, combined_coords)
combined.to(filename="combined.xyz")
```

## Best Practices

### Copying Structures

```python
# Always copy before modifying if you need original
original = Structure.from_file("POSCAR")
modified = original.copy()
modified.perturb(0.1)

# original is unchanged
```

### Validation

```python
# Check structure validity after transformations
if not structure.is_valid():
    print("Warning: Structure may have overlapping atoms")

# Check for reasonable volumes
if structure.volume <= 0:
    print("Error: Negative or zero volume")

# Check for reasonable distances
min_dist = min([
    structure.get_distance(i, j)
    for i in range(len(structure))
    for j in range(i+1, len(structure))
])

if min_dist < 0.5:
    print(f"Warning: Very close atoms (min distance: {min_dist:.2f} Å)")
```

### Reproducibility

```python
# Use Transformation classes for reproducibility
# Save transformation parameters
import json

transformation_params = {
    "type": "SupercellTransformation",
    "scaling_matrix": [[2, 0, 0], [0, 2, 0], [0, 0, 1]]
}

with open("transformation.json", "w") as f:
    json.dump(transformation_params, f)

# Can recreate transformation later
```

## Quick Reference

### Direct Methods

| Operation | Code |
|-----------|------|
| **Supercell** | `structure * (2, 2, 1)` |
| **Add site** | `structure.append("O", [0.5, 0.5, 0.5])` |
| **Remove site** | `structure.remove_sites([0])` |
| **Replace species** | `structure.replace_species({"Li": "Na"})` |
| **Perturb** | `structure.perturb(distance=0.1)` |
| **Sort** | `structure.sort()` |
| **Interpolate** | `start.interpolate(end, nimages=5)` |

### Transformation Classes

| Transformation | Purpose |
|---------------|---------|
| `SupercellTransformation` | Create supercells |
| `SubstitutionTransformation` | Replace elements |
| `PerturbStructureTransformation` | Random displacements |
| `RemoveSitesTransformation` | Remove specific sites |
| `DeformStructureTransformation` | Apply strain |
| `SlabTransformation` | Create surfaces |
| `OrderDisorderedStructureTransformation` | Order partial occupancies |

### Surface Creation

```python
from pymatgen.core.surface import SlabGenerator

slabgen = SlabGenerator(structure, (1,1,1), 10.0, 15.0)
slabs = slabgen.get_slabs()
```
