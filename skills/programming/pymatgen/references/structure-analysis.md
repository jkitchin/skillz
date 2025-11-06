# Pymatgen Structure Analysis

Comprehensive guide to analyzing crystal structures: symmetry, space groups, structure comparison, and neighbor analysis.

## Overview

Pymatgen provides powerful tools for analyzing crystal structures:
- **SpacegroupAnalyzer** - Symmetry and space group determination
- **StructureMatcher** - Comparing and matching structures
- **VoronoiAnalyzer** - Coordination environment analysis
- Distance and neighbor calculations with periodic boundary conditions

## Symmetry Analysis

### SpacegroupAnalyzer

**Basic usage:**
```python
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.core import Structure

structure = Structure.from_file("POSCAR")
sga = SpacegroupAnalyzer(structure)

# Space group information
print(f"Space group symbol: {sga.get_space_group_symbol()}")
print(f"Space group number: {sga.get_space_group_number()}")
print(f"Crystal system: {sga.get_crystal_system()}")
print(f"Point group: {sga.get_point_group_symbol()}")
```

**Common methods:**
```python
# Space group
sga.get_space_group_symbol()      # e.g., "Fm-3m"
sga.get_space_group_number()      # e.g., 225
sga.get_hall()                    # Hall symbol

# Crystal system
sga.get_crystal_system()          # e.g., "cubic", "hexagonal"
sga.get_lattice_type()            # e.g., "cubic", "tetragonal"
sga.get_point_group_symbol()      # e.g., "m-3m"

# Structure info
sga.get_space_group_info()        # (symbol, number)
```

### Symmetry Tolerance

```python
# Default tolerance: symprec=0.01, angle_tolerance=5
sga = SpacegroupAnalyzer(structure, symprec=0.01, angle_tolerance=5)

# Tighter tolerance for high-quality structures
sga_tight = SpacegroupAnalyzer(structure, symprec=0.001, angle_tolerance=0.1)

# Looser tolerance for noisy/distorted structures
sga_loose = SpacegroupAnalyzer(structure, symprec=0.1, angle_tolerance=10)
```

### Primitive and Conventional Cells

**Get primitive cell:**
```python
sga = SpacegroupAnalyzer(structure)

# Get primitive structure
primitive = sga.get_primitive_standard_structure()

print(f"Original sites: {len(structure)}")
print(f"Primitive sites: {len(primitive)}")
```

**Get conventional cell:**
```python
# Get conventional standard structure
conventional = sga.get_conventional_standard_structure()

print(f"Original sites: {len(structure)}")
print(f"Conventional sites: {len(conventional)}")
```

**Refined structure:**
```python
# Get symmetrized (idealized) structure
refined = sga.get_refined_structure()

# This applies the detected symmetry operations
# to get an idealized structure
```

### Symmetry Operations

**Get symmetry operations:**
```python
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

sga = SpacegroupAnalyzer(structure)

# Get symmetry operations
sym_ops = sga.get_symmetry_operations()

print(f"Number of symmetry operations: {len(sym_ops)}")

# Each operation is a SymmOp object
for i, op in enumerate(sym_ops[:3]):  # First 3
    print(f"\nOperation {i}:")
    print(f"Rotation:\n{op.rotation_matrix}")
    print(f"Translation: {op.translation_vector}")
```

**Apply symmetry operations:**
```python
# Get a symmetry operation
sym_op = sym_ops[1]

# Apply to coordinates
coords = [0.25, 0.25, 0.25]
new_coords = sym_op.operate(coords)
print(f"Original: {coords}")
print(f"Transformed: {new_coords}")
```

### Wyckoff Positions

```python
sga = SpacegroupAnalyzer(structure)

# Get Wyckoff symbols for each site
wyckoffs = sga.get_symmetry_dataset()["wyckoffs"]

for i, site in enumerate(structure):
    print(f"Site {i} ({site.species_string}): Wyckoff {wyckoffs[i]}")
```

## Structure Comparison

### StructureMatcher

**Basic usage:**
```python
from pymatgen.analysis.structure_matcher import StructureMatcher

# Create matcher
matcher = StructureMatcher()

# Compare two structures
structure1 = Structure.from_file("POSCAR1")
structure2 = Structure.from_file("POSCAR2")

if matcher.fit(structure1, structure2):
    print("Structures match!")
else:
    print("Structures are different")
```

**Get RMS distance:**
```python
matcher = StructureMatcher()

# Get RMS distance between structures
rms = matcher.get_rms_dist(structure1, structure2)

if rms is not None:
    print(f"RMS distance: {rms[0]:.4f} Å")
else:
    print("Structures don't match")
```

### Matcher Tolerance

```python
# Default tolerances
matcher = StructureMatcher(
    ltol=0.2,          # Lattice length tolerance (fractional)
    stol=0.3,          # Site position tolerance (Å)
    angle_tol=5.0,     # Angle tolerance (degrees)
    primitive_cell=True, # Compare primitive cells
    scale=True         # Allow volume scaling
)

# Strict matching
strict_matcher = StructureMatcher(
    ltol=0.01,
    stol=0.01,
    angle_tol=0.1
)

# Loose matching
loose_matcher = StructureMatcher(
    ltol=0.5,
    stol=0.5,
    angle_tol=10
)
```

### Advanced Matching

**Compare multiple structures:**
```python
from pymatgen.analysis.structure_matcher import StructureMatcher

structures = [
    Structure.from_file(f"POSCAR_{i}")
    for i in range(10)
]

matcher = StructureMatcher()

# Group similar structures
groups = matcher.group_structures(structures)

print(f"Found {len(groups)} unique structure types")
for i, group in enumerate(groups):
    print(f"Group {i}: {len(group)} structures")
```

**Get mapping:**
```python
# Get atom mapping between structures
matcher = StructureMatcher()

mapping = matcher.get_mapping(structure1, structure2)

if mapping:
    print("Atom mapping:")
    for idx1, idx2 in enumerate(mapping):
        print(f"Atom {idx1} in struct1 → Atom {idx2} in struct2")
```

### Element Comparison Options

```python
# Compare only structure, ignore elements
matcher = StructureMatcher(comparator=ElementComparator())

# Compare oxidation states
from pymatgen.analysis.structure_matcher import OccupancyComparator
matcher = StructureMatcher(comparator=OccupancyComparator())

# Frame-invariant comparison
from pymatgen.analysis.structure_matcher import FrameworkComparator
matcher = StructureMatcher(comparator=FrameworkComparator())
```

## Neighbor Analysis

### Get Neighbors

**Find neighbors within distance:**
```python
structure = Structure.from_file("POSCAR")

# Get neighbors of first site within 3 Å
neighbors = structure.get_neighbors(structure[0], r=3.0)

print(f"Found {len(neighbors)} neighbors")
for neighbor in neighbors:
    site = neighbor
    distance = structure.get_distance(0, neighbor.index)
    print(f"  {site.species_string} at {distance:.3f} Å")
```

**With periodic images:**
```python
# Get neighbors considering periodic boundary conditions
neighbors = structure.get_neighbors(
    structure[0],
    r=3.0,
    include_index=True,
    include_image=True
)

for neighbor, distance, index, image in neighbors:
    print(f"{neighbor.species_string}: {distance:.3f} Å, image: {image}")
```

### All Neighbors

```python
# Get all neighbors for all sites
site_index = 0
neighbors = structure.get_all_neighbors(r=3.0)

# neighbors[i] contains neighbors of site i
for neighbor_info in neighbors[site_index]:
    site = neighbor_info
    print(f"{site.species_string} at {site.nn_distance:.3f} Å")
```

### Nearest Neighbors

```python
# Get k nearest neighbors
k = 4  # Number of nearest neighbors
site_index = 0

neighbors = structure.get_neighbors(structure[site_index], r=10.0)
# Sort by distance
neighbors_sorted = sorted(neighbors, key=lambda x: x.nn_distance)
nearest_k = neighbors_sorted[:k]

print(f"{k} nearest neighbors:")
for n in nearest_k:
    print(f"  {n.species_string}: {n.nn_distance:.3f} Å")
```

## Distance Calculations

### Site-to-Site Distance

```python
structure = Structure.from_file("POSCAR")

# Distance between two sites
dist = structure.get_distance(0, 1)
print(f"Distance: {dist:.3f} Å")

# Distance with periodic boundary conditions
dist, image = structure.get_distance(0, 1, jimage=[0, 0, 0])
print(f"Distance: {dist:.3f} Å, image: {image}")
```

### All Pairwise Distances

```python
import numpy as np

# Get distance matrix
n_sites = len(structure)
dist_matrix = np.zeros((n_sites, n_sites))

for i in range(n_sites):
    for j in range(i+1, n_sites):
        dist = structure.get_distance(i, j)
        dist_matrix[i, j] = dist
        dist_matrix[j, i] = dist

print("Distance matrix:")
print(dist_matrix)
```

### Distance Between Fractional Coordinates

```python
# Get distance between fractional coordinates
frac_coords1 = [0.0, 0.0, 0.0]
frac_coords2 = [0.5, 0.5, 0.5]

# Convert to Cartesian and compute distance
lattice = structure.lattice
cart1 = lattice.get_cartesian_coords(frac_coords1)
cart2 = lattice.get_cartesian_coords(frac_coords2)

import numpy as np
dist = np.linalg.norm(cart2 - cart1)
print(f"Distance: {dist:.3f} Å")
```

### Minimum Image Convention

```python
# Get shortest distance considering periodic images
lattice = structure.lattice

# Distance with image information
dist, image = lattice.get_distance_and_image(
    [0, 0, 0],
    [0.6, 0.6, 0.6]
)

print(f"Distance: {dist:.3f} Å")
print(f"Nearest image: {image}")  # Which periodic image
```

## Coordination Analysis

### VoronoiAnalyzer

**Basic coordination analysis:**
```python
from pymatgen.analysis.structure_analyzer import VoronoiAnalyzer

structure = Structure.from_file("POSCAR")
va = VoronoiAnalyzer()

# Analyze coordination
coordination = va.analyze(structure)

for i, coord_info in enumerate(coordination):
    site = structure[i]
    print(f"Site {i} ({site.species_string}):")
    print(f"  Coordination: {len(coord_info)}")
    for neighbor_idx, info in coord_info.items():
        print(f"    Neighbor {neighbor_idx}: {info}")
```

### Coordination Number

```python
from pymatgen.analysis.local_env import CrystalNN

# Use neural network method
cnn = CrystalNN()

structure = Structure.from_file("POSCAR")

# Get coordination for each site
for i, site in enumerate(structure):
    cn_info = cnn.get_cn_dict(structure, i)
    cn = sum(cn_info.values())
    print(f"Site {i} ({site.species_string}): CN = {cn:.1f}")
```

**Other coordination methods:**
```python
from pymatgen.analysis.local_env import (
    VoronoiNN,           # Voronoi-based
    MinimumDistanceNN,   # Distance cutoff
    MinimumOKeeffeNN,    # O'Keeffe method
    BrunnerNN_real       # Brunner method
)

# Voronoi nearest neighbors
voronoi_nn = VoronoiNN()
neighbors = voronoi_nn.get_nn_info(structure, 0)

# Distance-based
min_dist_nn = MinimumDistanceNN(cutoff=3.0)
neighbors = min_dist_nn.get_nn_info(structure, 0)
```

## Bond Analysis

### Get Bonds

```python
# For molecules or structures
from pymatgen.core import Structure

structure = Structure.from_file("POSCAR")

# Simple bond detection by distance
bonds = []
for i, site1 in enumerate(structure):
    for j, site2 in enumerate(structure[i+1:], start=i+1):
        dist = structure.get_distance(i, j)
        # Use covalent radii sum as threshold
        if dist < 3.0:  # Å
            bonds.append((i, j, dist))

print(f"Found {len(bonds)} bonds")
for i, j, dist in bonds:
    s1 = structure[i].species_string
    s2 = structure[j].species_string
    print(f"{s1}-{s2}: {dist:.3f} Å")
```

### Bond Valence Sum

```python
from pymatgen.analysis.bond_valence import BVAnalyzer

structure = Structure.from_file("POSCAR")

# Analyze oxidation states via bond valence
bva = BVAnalyzer()

try:
    # Get structure with oxidation states
    structure_with_oxi = bva.get_oxi_state_decorated_structure(structure)

    for site in structure_with_oxi:
        print(f"{site.species_string}: oxidation state")
except ValueError as e:
    print(f"Could not determine oxidation states: {e}")
```

## Structure Comparison Metrics

### Root Mean Square Displacement (RMSD)

```python
from pymatgen.analysis.structure_matcher import StructureMatcher

matcher = StructureMatcher()
structure1 = Structure.from_file("POSCAR1")
structure2 = Structure.from_file("POSCAR2")

# Get RMSD
result = matcher.get_rms_dist(structure1, structure2)

if result:
    rmsd, max_dist = result
    print(f"RMSD: {rmsd:.4f} Å")
    print(f"Max displacement: {max_dist:.4f} Å")
```

### Lattice Similarity

```python
# Compare lattices
lattice1 = structure1.lattice
lattice2 = structure2.lattice

# Volume difference
vol_diff = abs(lattice1.volume - lattice2.volume) / lattice1.volume
print(f"Volume difference: {vol_diff*100:.2f}%")

# Lattice parameters
params1 = lattice1.parameters  # (a, b, c, alpha, beta, gamma)
params2 = lattice2.parameters

for i, label in enumerate(['a', 'b', 'c', 'α', 'β', 'γ']):
    diff = abs(params1[i] - params2[i])
    print(f"{label}: {params1[i]:.3f} vs {params2[i]:.3f} (Δ={diff:.3f})")
```

### Composition Similarity

```python
# Compare compositions
comp1 = structure1.composition
comp2 = structure2.composition

# Check if same elements
if set(comp1.elements) == set(comp2.elements):
    print("Same elements")

    # Compare ratios
    formula1 = comp1.reduced_formula
    formula2 = comp2.reduced_formula

    if formula1 == formula2:
        print(f"Same composition: {formula1}")
    else:
        print(f"Different ratios: {formula1} vs {formula2}")
```

## Structure Properties

### Geometric Properties

```python
structure = Structure.from_file("POSCAR")

# Basic properties
print(f"Volume: {structure.volume:.3f} Ų")
print(f"Density: {structure.density:.3f} g/cm³")
print(f"Number of sites: {structure.num_sites}")

# Lattice properties
lattice = structure.lattice
print(f"Lattice parameters: {lattice.parameters}")
print(f"Lattice type: {lattice.lattice_type}")
print(f"Is orthogonal: {lattice.is_orthogonal}")
```

### Packing Fraction

```python
from pymatgen.core import Element

# Estimate packing fraction (simple sphere model)
total_volume = 0
for site in structure:
    element = Element(site.species_string)
    # Use covalent radius or atomic radius
    if hasattr(element, 'atomic_radius'):
        r = element.atomic_radius
    else:
        r = 1.5  # Default
    vol_sphere = (4/3) * 3.14159 * r**3
    total_volume += vol_sphere

packing_fraction = total_volume / structure.volume
print(f"Packing fraction: {packing_fraction:.3f}")
```

## Analyzing Structure Collections

### Compare Multiple Structures

```python
from pymatgen.analysis.structure_matcher import StructureMatcher
import glob

structures = []
filenames = []

# Load all structures
for filename in glob.glob("structures/*.cif"):
    struct = Structure.from_file(filename)
    structures.append(struct)
    filenames.append(filename)

# Group by similarity
matcher = StructureMatcher()
groups = matcher.group_structures(structures)

print(f"Found {len(groups)} unique structure types:")
for i, group in enumerate(groups):
    print(f"\nGroup {i}:")
    for struct in group:
        idx = structures.index(struct)
        print(f"  {filenames[idx]}")
```

### Structure Fingerprints

```python
from pymatgen.analysis.fingerprints import (
    SiteStatsFingerprint,
    StructureFingerprint
)

# Create fingerprint
fingerprinter = SiteStatsFingerprint.from_preset("SOAP")
fp1 = fingerprinter.fingerprint(structure1)
fp2 = fingerprinter.fingerprint(structure2)

# Compare fingerprints
similarity = fingerprinter.distance(fp1, fp2)
print(f"Fingerprint distance: {similarity:.4f}")
```

## Symmetry Detection Tips

### Handling Noisy Structures

```python
# For MD snapshots or experimental structures with noise
sga = SpacegroupAnalyzer(structure, symprec=0.1, angle_tolerance=10)

# If no symmetry detected, try looser tolerance
if sga.get_space_group_number() == 1:  # P1 (no symmetry)
    sga_loose = SpacegroupAnalyzer(structure, symprec=0.5, angle_tolerance=15)
    print(f"With loose tolerance: {sga_loose.get_space_group_symbol()}")
```

### Symmetrizing Structures

```python
# "Fix" slightly distorted structures
sga = SpacegroupAnalyzer(structure)

# Get idealized structure
symmetrized = sga.get_refined_structure()

# Compare
print(f"Original space group: {sga.get_space_group_symbol()}")
sga_new = SpacegroupAnalyzer(symmetrized)
print(f"Refined space group: {sga_new.get_space_group_symbol()}")

# Check RMSD
matcher = StructureMatcher()
if matcher.fit(structure, symmetrized):
    rms = matcher.get_rms_dist(structure, symmetrized)
    print(f"RMSD from original: {rms[0]:.4f} Å")
```

## Quick Reference

### Common Analysis Tasks

| Task | Code |
|------|------|
| **Space group** | `SpacegroupAnalyzer(s).get_space_group_symbol()` |
| **Primitive cell** | `SpacegroupAnalyzer(s).get_primitive_standard_structure()` |
| **Conventional cell** | `SpacegroupAnalyzer(s).get_conventional_standard_structure()` |
| **Compare structures** | `StructureMatcher().fit(s1, s2)` |
| **Get neighbors** | `structure.get_neighbors(site, r=3.0)` |
| **Distance** | `structure.get_distance(i, j)` |
| **Coordination** | `CrystalNN().get_cn_dict(structure, i)` |

### Symmetry Tolerances

- **High quality** (DFT): `symprec=0.01, angle_tolerance=5`
- **Normal**: `symprec=0.01, angle_tolerance=5` (default)
- **Noisy** (MD, experimental): `symprec=0.1, angle_tolerance=10`
- **Very noisy**: `symprec=0.5, angle_tolerance=15`

### Structure Matching Tolerances

- **Strict**: `ltol=0.01, stol=0.01, angle_tol=0.1`
- **Normal**: `ltol=0.2, stol=0.3, angle_tol=5` (default)
- **Loose**: `ltol=0.5, stol=0.5, angle_tol=10`
