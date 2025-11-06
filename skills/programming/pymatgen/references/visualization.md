# Pymatgen Structure Visualization

Guide to visualizing crystal structures, molecules, and analysis results using pymatgen.

## Overview

Pymatgen provides several visualization options:
- **Built-in plotters** - Using matplotlib
- **External viewers** - VESTA, XCrySDen, Jmol, Avogadro
- **Interactive 3D** - Jupyter notebooks with nglview, crystaltoolkit
- **Publication figures** - High-quality structure images

## Quick Visualization

### View in External Program

**VESTA (recommended for crystals):**
```python
from pymatgen.core import Structure

structure = Structure.from_file("POSCAR")

# Write CIF for VESTA
structure.to(filename="structure.cif")

# Open in VESTA manually
# Or use subprocess
import subprocess
subprocess.run(["vesta", "structure.cif"])
```

**XCrySDen:**
```python
# Write XSF format
structure.to(filename="structure.xsf", fmt="xsf")

# Open in XCrySDen
subprocess.run(["xcrysden", "--xsf", "structure.xsf"])
```

**Avogadro (for molecules):**
```python
from pymatgen.core import Molecule

molecule = Molecule.from_file("molecule.xyz")

# Avogadro reads XYZ directly
subprocess.run(["avogadro", "molecule.xyz"])
```

## Matplotlib Visualization

### Basic Structure Plot

```python
from pymatgen.core import Structure
from pymatgen.vis.structure_vtk import StructureVis
import matplotlib.pyplot as plt

structure = Structure.from_file("POSCAR")

# Simple plot (limited 3D capability)
# Pymatgen's matplotlib support is basic
# Better to use external viewers for 3D
```

### Plotting Specific Views

**Unit cell:**
```python
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np

fig = plt.figure(figsize=(8, 8))
ax = fig.add_subplot(111, projection='3d')

# Plot lattice vectors
lattice = structure.lattice
origin = np.array([0, 0, 0])

vectors = [lattice.matrix[i] for i in range(3)]
colors = ['r', 'g', 'b']
labels = ['a', 'b', 'c']

for vec, color, label in zip(vectors, colors, labels):
    ax.quiver(*origin, *vec, color=color, arrow_length_ratio=0.1, label=label)

# Plot atoms
for site in structure:
    coords = site.coords
    ax.scatter(*coords, s=200, alpha=0.6)

ax.set_xlabel("x (Å)")
ax.set_ylabel("y (Å)")
ax.set_zlabel("z (Å)")
ax.legend()
plt.tight_layout()
plt.savefig("structure_3d.png", dpi=300)
plt.show()
```

## Jupyter Notebook Visualization

### Using nglview

**Interactive 3D viewer:**
```python
# In Jupyter notebook
import nglview as nv
from pymatgen.core import Structure

structure = Structure.from_file("POSCAR")

# Convert to ASE for nglview
from pymatgen.io.ase import AseAtomsAdaptor

atoms = AseAtomsAdaptor.get_atoms(structure)

# View interactively
view = nv.show_ase(atoms)
view
```

**Customizing view:**
```python
# Control representation
view = nv.show_ase(atoms)
view.clear_representations()
view.add_ball_and_stick()
view.add_unitcell()

# Color by element
view.add_representation('ball+stick', selection='Fe', color='orange')
view.add_representation('ball+stick', selection='O', color='red')

view
```

### Using Crystal Toolkit

```python
from crystal_toolkit import StructureMoleculeComponent
from pymatgen.core import Structure

structure = Structure.from_file("POSCAR")

# Create interactive component
struct_component = StructureMoleculeComponent(structure)

# Display in Dash app (requires running app)
# See Crystal Toolkit documentation for full examples
```

## Publication-Quality Figures

### Using pymatgen + matplotlib

**2D projections:**
```python
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
import numpy as np

structure = Structure.from_file("POSCAR")

fig, ax = plt.subplots(figsize=(8, 8))

# Project along c-axis (view from top)
for site in structure:
    coords = site.frac_coords
    x, y = coords[0], coords[1]

    # Map to Cartesian for plotting
    cart = structure.lattice.get_cartesian_coords([x, y, 0])

    # Color by element
    element = site.species_string
    colors = {"Li": "purple", "Fe": "orange", "P": "gray", "O": "red"}
    color = colors.get(element, "black")

    # Size by element
    sizes = {"Li": 200, "Fe": 300, "P": 250, "O": 150}
    size = sizes.get(element, 200)

    circle = Circle((cart[0], cart[1]), radius=0.3, color=color, alpha=0.7)
    ax.add_patch(circle)
    ax.text(cart[0], cart[1], element, ha='center', va='center', fontsize=8)

# Plot unit cell
a = structure.lattice.matrix[0][:2]
b = structure.lattice.matrix[1][:2]

vertices = np.array([[0, 0], a, a+b, b, [0, 0]])
ax.plot(vertices[:, 0], vertices[:, 1], 'k-', linewidth=2)

ax.set_aspect('equal')
ax.set_xlabel("x (Å)", fontsize=14)
ax.set_ylabel("y (Å)", fontsize=14)
ax.set_title("Structure viewed along c-axis", fontsize=16)
plt.tight_layout()
plt.savefig("structure_projection.png", dpi=300, bbox_inches='tight')
plt.show()
```

### Using VESTA for Publication

```python
# Export to CIF
structure.to(filename="structure.cif")

# Open in VESTA
# Adjust visualization settings:
# - Style > Atoms (ball and stick)
# - Objects > Boundary (show unit cell)
# - Objects > Polyhedra (show coordination)
# - File > Export > Raster Image (high DPI)

# For publication:
# - Use high DPI (600+)
# - Adjust bond thickness
# - Set background to white
# - Label atoms if needed
```

## Specialized Visualizations

### Slabs and Surfaces

```python
from pymatgen.core.surface import SlabGenerator
import matplotlib.pyplot as plt

# Create slab
bulk = Structure.from_file("POSCAR")
slabgen = SlabGenerator(bulk, (1,1,1), 10, 15)
slab = slabgen.get_slabs()[0]

# Save for visualization
slab.to(filename="slab.cif")

# View slab showing surface
# In VESTA: view along surface normal
```

### Supercells

```python
# Create and visualize supercell
supercell = structure * (2, 2, 1)

supercell.to(filename="supercell.cif")
# View in VESTA to see repeated pattern
```

### Defects

```python
# Visualize vacancy
defect_structure = structure.copy()
defect_structure.remove_sites([0])  # Create vacancy

defect_structure.to(filename="with_vacancy.cif")
```

### Bond Highlighting

```python
# Create figure highlighting specific bonds
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

fig = plt.figure(figsize=(10, 10))
ax = fig.add_subplot(111, projection='3d')

# Plot atoms
for i, site in enumerate(structure):
    ax.scatter(*site.coords, s=500, alpha=0.8, label=site.species_string)

# Highlight specific bonds
bond_pairs = [(0, 1), (0, 2), (1, 3)]  # Atom indices

for i, j in bond_pairs:
    coords_i = structure[i].coords
    coords_j = structure[j].coords

    ax.plot([coords_i[0], coords_j[0]],
            [coords_i[1], coords_j[1]],
            [coords_i[2], coords_j[2]],
            'k-', linewidth=2, alpha=0.5)

ax.set_xlabel("x (Å)")
ax.set_ylabel("y (Å)")
ax.set_zlabel("z (Å)")
plt.legend()
plt.savefig("bonds_highlighted.png", dpi=300)
```

## Coordination Polyhedra

### Visualize Coordination Environment

```python
# In VESTA:
# 1. Open structure.cif
# 2. Objects > Polyhedra
# 3. Select central atom type
# 4. Select coordinating atom types
# 5. Set distance cutoff
# 6. Apply

# Example: Show FeO6 octahedra in LiFePO4
# - Central: Fe
# - Coordinating: O
# - Distance: 2.5 Å
```

### Export Polyhedra View

```python
structure.to(filename="structure_for_polyhedra.cif")

# Open in VESTA
# Set up polyhedra as above
# File > Export > Raster Image
# Choose high DPI for publication
```

## Animating Structures

### Create Animation from MD Trajectory

```python
from pymatgen.io.vasp.outputs import Xdatcar
import matplotlib.pyplot as plt
import matplotlib.animation as animation

# Parse XDATCAR (MD trajectory)
xdatcar = Xdatcar("XDATCAR")
structures = xdatcar.structures

# Create animation
fig = plt.figure(figsize=(8, 8))
ax = fig.add_subplot(111, projection='3d')

def update(frame):
    ax.clear()
    structure = structures[frame]

    for site in structure:
        ax.scatter(*site.coords, s=200)

    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_zlabel("z")
    ax.set_title(f"Frame {frame}")

    return ax,

ani = animation.FuncAnimation(fig, update, frames=len(structures), interval=100)
ani.save("trajectory.gif", writer='pillow', fps=10)
```

### NEB Path Visualization

```python
# Visualize reaction path
from pymatgen.core import Structure

# Read NEB images
images = []
for i in range(7):  # 7 images
    structure = Structure.from_file(f"0{i}/CONTCAR")
    images.append(structure)

# Save each image
for i, img in enumerate(images):
    img.to(filename=f"neb_image_{i:02d}.cif")

# View sequence in VESTA or create animation
```

## Property Visualization

### Charge Density

```python
from pymatgen.io.vasp.outputs import Chgcar
import matplotlib.pyplot as plt

chgcar = Chgcar.from_file("CHGCAR")

# Plot 2D slice
charge_data = chgcar.data['total']

# Slice along z at half-height
z_slice = charge_data.shape[2] // 2
slice_data = charge_data[:, :, z_slice]

plt.figure(figsize=(8, 8))
plt.imshow(slice_data.T, origin='lower', cmap='viridis')
plt.colorbar(label="Charge density (e/Ų)")
plt.xlabel("x")
plt.ylabel("y")
plt.title(f"Charge density at z = {z_slice}")
plt.savefig("charge_density_slice.png", dpi=300)
```

### Potential Profile

```python
from pymatgen.io.vasp.outputs import Locpot
import matplotlib.pyplot as plt

locpot = Locpot.from_file("LOCPOT")

# Average along axis
avg_potential = locpot.get_average_along_axis(2)  # z-axis

# Distance along z
z_pos = np.linspace(0, locpot.structure.lattice.c, len(avg_potential))

plt.figure(figsize=(10, 6))
plt.plot(z_pos, avg_potential, linewidth=2)
plt.xlabel("z position (Å)", fontsize=14)
plt.ylabel("Potential (eV)", fontsize=14)
plt.title("Planar-averaged electrostatic potential", fontsize=16)
plt.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig("potential_profile.png", dpi=300)
```

### Phonon Dispersion

```python
# Requires phonopy integration
# Example visualization of phonon modes

import matplotlib.pyplot as plt

# (Phonopy data parsing would go here)

# Plot phonon dispersion
plt.figure(figsize=(10, 8))
# plt.plot(qpoints, frequencies)
plt.xlabel("Wave vector")
plt.ylabel("Frequency (THz)")
plt.title("Phonon dispersion")
plt.savefig("phonon_dispersion.png", dpi=300)
```

## Comparing Structures

### Side-by-Side Visualization

```python
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

struct1 = Structure.from_file("POSCAR_initial")
struct2 = Structure.from_file("POSCAR_final")

fig = plt.figure(figsize=(16, 8))

# Initial structure
ax1 = fig.add_subplot(121, projection='3d')
for site in struct1:
    ax1.scatter(*site.coords, s=200)
ax1.set_title("Initial")

# Final structure
ax2 = fig.add_subplot(122, projection='3d')
for site in struct2:
    ax2.scatter(*site.coords, s=200)
ax2.set_title("Final")

plt.tight_layout()
plt.savefig("structure_comparison.png", dpi=300)
```

### Displacement Vectors

```python
# Show atomic displacements
from pymatgen.analysis.structure_matcher import StructureMatcher

matcher = StructureMatcher()
mapping = matcher.get_mapping(struct1, struct2)

if mapping:
    fig = plt.figure(figsize=(10, 10))
    ax = fig.add_subplot(111, projection='3d')

    # Plot atoms and displacement vectors
    for i, site1 in enumerate(struct1):
        j = mapping[i]
        site2 = struct2[j]

        # Initial position
        ax.scatter(*site1.coords, s=200, alpha=0.5, color='blue')

        # Final position
        ax.scatter(*site2.coords, s=200, alpha=0.5, color='red')

        # Displacement vector
        displacement = site2.coords - site1.coords
        ax.quiver(*site1.coords, *displacement,
                  color='gray', arrow_length_ratio=0.1)

    ax.set_xlabel("x (Å)")
    ax.set_ylabel("y (Å)")
    ax.set_zlabel("z (Å)")
    ax.set_title("Atomic displacements")
    plt.savefig("displacements.png", dpi=300)
```

## Export Formats

### For Different Programs

```python
structure = Structure.from_file("POSCAR")

# VESTA (CIF)
structure.to(filename="structure.cif", fmt="cif")

# XCrySDen (XSF)
structure.to(filename="structure.xsf", fmt="xsf")

# Gaussian
structure.to(filename="structure.gjf", fmt="gaussian")

# LAMMPS
structure.to(filename="structure.lammps", fmt="lammps-data")

# PDB (for molecules)
molecule.to(filename="molecule.pdb", fmt="pdb")
```

## Best Practices

### Publication Figures

**High-quality structure images:**
1. Export to CIF or XSF
2. Open in VESTA or similar
3. Settings:
   - Style: Ball-and-stick or space-filling
   - Background: White for publications
   - Bonds: Appropriate cutoffs
   - Labels: Element symbols or site numbers
   - Unit cell: Show or hide as appropriate
4. Export as:
   - Vector (PDF, SVG) for scalability
   - Raster (PNG, TIFF) at ≥300 DPI

**For presentations:**
- Higher contrast colors
- Larger labels
- Simpler backgrounds
- Highlight key features

**For SI (supporting information):**
- Include all structural details
- Show multiple views
- Label important atoms/bonds
- Include scale bars/cell parameters

### Color Schemes

```python
# Standard element colors (Jmol)
element_colors = {
    "H": "#FFFFFF",
    "C": "#909090",
    "N": "#3050F8",
    "O": "#FF0D0D",
    "S": "#FFFF30",
    "Fe": "#E06633",
    "Cu": "#C88033",
    "Zn": "#7D80B0",
}

# Use consistent colors across figures
```

### Viewing Angles

```python
# Choose informative viewing angles
# - Along high-symmetry directions
# - Show interesting structural features
# - Avoid cluttered views

# For layered materials: perpendicular to layers
# For 1D materials: along chain direction
# For molecular crystals: show hydrogen bonding
```

## Troubleshooting

### Issue: Structure not visible in VESTA

**Solution:** Check file format and structure validity
```python
structure = Structure.from_file("POSCAR")

# Validate
if not structure.is_valid():
    print("Structure has issues")

# Try different format
structure.to(filename="structure.cif")
structure.to(filename="structure.xsf")
```

### Issue: Bonds not showing

**Solution:** Adjust bond distance cutoffs
```python
# In VESTA:
# Style > Bond
# Set appropriate distance cutoffs for each element pair

# In pymatgen, check distances:
for i in range(len(structure)):
    neighbors = structure.get_neighbors(structure[i], r=3.0)
    print(f"Site {i}: {len(neighbors)} neighbors within 3 Å")
```

### Issue: Slow rendering for large structures

**Solution:** Simplify or use different tools
```python
# For very large structures (>1000 atoms):
# - Use VMD or OVITO for molecular dynamics
# - Use ParaView for charge density visualization
# - Reduce representation complexity in VESTA
# - Show representative portion only
```

## Quick Reference

### Visualization Tools

| Tool | Best for | Format |
|------|----------|--------|
| **VESTA** | Crystals, polyhedra | CIF, VASP |
| **XCrySDen** | Crystals, charge density | XSF |
| **Avogadro** | Molecules | XYZ, PDB |
| **VMD** | MD trajectories | XYZ, DCD |
| **OVITO** | Large systems, MD | XYZ, LAMMPS |
| **nglview** | Jupyter interactive | All (via ASE) |

### Export Formats

```python
# CIF (most universal)
structure.to(filename="structure.cif")

# XYZ (for molecules)
molecule.to(filename="molecule.xyz")

# VASP
structure.to(filename="POSCAR")

# XSF (XCrySDen)
structure.to(filename="structure.xsf", fmt="xsf")
```

### Common Visualizations

- **Unit cell**: VESTA with cell boundaries
- **Coordination**: VESTA polyhedra tool
- **Bonding**: Ball-and-stick representation
- **Packing**: Space-filling representation
- **Surfaces**: View slab perpendicular to surface
- **Defects**: Compare with/without defect
- **Charge density**: 2D slices or isosurfaces
