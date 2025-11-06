# Pymatgen File Input/Output

Comprehensive guide to reading and writing crystal structures and molecules in various file formats.

## Overview

Pymatgen supports many file formats with automatic format detection. The universal `from_file()` method detects format from extension, or you can use format-specific parsers for more control.

## Universal File I/O

### Reading Structures (Any Format)

```python
from pymatgen.core import Structure

# Auto-detect format from extension
structure = Structure.from_file("structure.cif")
structure = Structure.from_file("POSCAR")
structure = Structure.from_file("molecule.xyz")

# Pymatgen automatically detects:
# - CIF (.cif)
# - POSCAR/CONTCAR (VASP)
# - XYZ (.xyz)
# - PDB (.pdb)
# - JSON (.json)
# - YAML (.yaml)
# - And many more
```

### Writing Structures (Any Format)

```python
# Write to specific format
structure.to(filename="output.cif")
structure.to(filename="POSCAR")
structure.to(filename="output.xyz")

# Or specify format explicitly
structure.to(filename="output", fmt="poscar")
structure.to(filename="output", fmt="cif")
```

## CIF Files (Crystallographic Information File)

### Reading CIF

**Basic usage:**
```python
from pymatgen.core import Structure

# Simple read
structure = Structure.from_file("structure.cif")

# For more control, use CifParser
from pymatgen.io.cif import CifParser

parser = CifParser("structure.cif")

# Get all structures (CIF can contain multiple)
structures = parser.get_structures()
structure = structures[0]

# Or get by identifier
structure = parser.get_structures(primitive=True)[0]
```

**CIF with multiple structures:**
```python
parser = CifParser("multi_structure.cif")

# Get dictionary of structures
all_structures = parser.as_dict()

# Iterate through structures
for key, structure in all_structures.items():
    print(f"{key}: {structure.composition.reduced_formula}")
```

**Options:**
```python
# Get primitive cell
parser = CifParser("structure.cif", occupancy_tolerance=1.0)
structure = parser.get_structures(primitive=True)[0]

# Get conventional cell
structure = parser.get_structures(primitive=False)[0]
```

### Writing CIF

```python
from pymatgen.io.cif import CifWriter

# Basic write
structure.to(filename="output.cif")

# More control with CifWriter
writer = CifWriter(structure, symprec=0.1)
writer.write_file("output.cif")

# Write multiple structures
from pymatgen.io.cif import CifWriter

structures = [struct1, struct2, struct3]
writer = CifWriter(structures)
writer.write_file("multi_structure.cif")
```

## VASP Files

### POSCAR/CONTCAR

**Reading:**
```python
from pymatgen.io.vasp import Poscar

# Read POSCAR
poscar = Poscar.from_file("POSCAR")
structure = poscar.structure

# Or directly
structure = Structure.from_file("POSCAR")

# CONTCAR (output from VASP run)
structure = Structure.from_file("CONTCAR")
```

**Writing:**
```python
from pymatgen.io.vasp import Poscar

# Create Poscar object
poscar = Poscar(structure)

# Write to file
poscar.write_file("POSCAR_new")

# With selective dynamics
from pymatgen.io.vasp import Poscar

# Add selective dynamics to structure
selective_dynamics = [[True, True, True], [False, False, False]]
structure.add_site_property("selective_dynamics", selective_dynamics)

poscar = Poscar(structure)
poscar.write_file("POSCAR")
```

**POSCAR with velocities:**
```python
# Read CONTCAR with velocities
poscar = Poscar.from_file("CONTCAR", check_for_POTCAR=False)
if poscar.velocities is not None:
    print("Velocities present")
    print(poscar.velocities)
```

### OUTCAR

**Reading VASP outputs:**
```python
from pymatgen.io.vasp.outputs import Outcar

# Parse OUTCAR
outcar = Outcar("OUTCAR")

# Extract information
print(f"Final energy: {outcar.final_energy} eV")
print(f"Forces:\n{outcar.final_fr}")
print(f"Stress: {outcar.final_stress}")

# Magnetization
if outcar.total_mag is not None:
    print(f"Total magnetization: {outcar.total_mag}")
```

### vasprun.xml

**Complete VASP output analysis:**
```python
from pymatgen.io.vasp.outputs import Vasprun

# Parse vasprun.xml
vasprun = Vasprun("vasprun.xml")

# Basic properties
final_structure = vasprun.final_structure
final_energy = vasprun.final_energy
is_converged = vasprun.converged

# Electronic structure
band_structure = vasprun.get_band_structure()
dos = vasprun.complete_dos

# Ionic steps
for i, step in enumerate(vasprun.ionic_steps):
    energy = step["e_wo_entrp"]
    structure = step["structure"]
    print(f"Step {i}: {energy} eV")
```

**Parse large files efficiently:**
```python
# For large files, parse only what you need
vasprun = Vasprun("vasprun.xml", parse_dos=False, parse_eigen=False)
final_energy = vasprun.final_energy
```

### INCAR

**Reading INCAR:**
```python
from pymatgen.io.vasp.inputs import Incar

incar = Incar.from_file("INCAR")

# Access parameters
print(incar["ENCUT"])
print(incar["EDIFF"])

# Check if parameter exists
if "ISPIN" in incar:
    print(f"ISPIN = {incar['ISPIN']}")
```

**Writing INCAR:**
```python
from pymatgen.io.vasp.inputs import Incar

# Create INCAR
incar = Incar({
    "PREC": "Accurate",
    "ENCUT": 520,
    "EDIFF": 1e-5,
    "NSW": 100,
    "IBRION": 2,
    "ISIF": 3
})

# Write to file
incar.write_file("INCAR")
```

### KPOINTS

**Reading KPOINTS:**
```python
from pymatgen.io.vasp.inputs import Kpoints

kpoints = Kpoints.from_file("KPOINTS")
print(kpoints.style)  # e.g., "Gamma" or "Monkhorst"
print(kpoints.kpts)   # K-point mesh
```

**Writing KPOINTS:**
```python
from pymatgen.io.vasp.inputs import Kpoints

# Automatic k-point mesh
kpoints = Kpoints.automatic(100)  # k-point density
kpoints.write_file("KPOINTS")

# Gamma-centered
kpoints = Kpoints.gamma_automatic([4, 4, 4])
kpoints.write_file("KPOINTS")

# Monkhorst-Pack
kpoints = Kpoints.monkhorst_automatic([4, 4, 4])
kpoints.write_file("KPOINTS")

# Manual k-points
kpoints = Kpoints(
    kpts=[[0, 0, 0], [0.5, 0, 0], [0.5, 0.5, 0]],
    kpts_weights=[1, 2, 1],
    style="Reciprocal"
)
```

### POTCAR

**Reading POTCAR:**
```python
from pymatgen.io.vasp.inputs import Potcar

potcar = Potcar.from_file("POTCAR")

# Get symbols
print(potcar.symbols)  # e.g., ['Fe', 'O']
```

**Writing POTCAR:**
```python
from pymatgen.io.vasp.inputs import Potcar

# Requires PMG_VASP_PSP_DIR environment variable or config
potcar = Potcar(["Fe", "O"], functional="PBE")
potcar.write_file("POTCAR")
```

## XYZ Files (Molecules)

### Reading XYZ

```python
from pymatgen.core import Molecule

# Read molecule
molecule = Molecule.from_file("molecule.xyz")

# Properties
print(molecule.formula)
print(molecule.cart_coords)
print(molecule.species)
```

### Writing XYZ

```python
# Write molecule
molecule.to(filename="output.xyz")

# Or for structure (removes periodicity)
structure = Structure.from_file("POSCAR")
molecule = structure.to(filename="molecule.xyz")
```

## JSON and YAML (Serialization)

### JSON Format

**Writing to JSON:**
```python
import json
from monty.json import MontyEncoder

# Structure to dict
structure_dict = structure.as_dict()

# Save to JSON
with open("structure.json", "w") as f:
    json.dump(structure_dict, f, cls=MontyEncoder, indent=2)

# Or use convenience method
structure.to(filename="structure.json")
```

**Reading from JSON:**
```python
import json
from monty.json import MontyDecoder

# Load from JSON
with open("structure.json", "r") as f:
    structure_dict = json.load(f, cls=MontyDecoder)

structure = Structure.from_dict(structure_dict)

# Or use convenience method
structure = Structure.from_file("structure.json")
```

### YAML Format

```python
# Write to YAML
structure.to(filename="structure.yaml")

# Read from YAML
structure = Structure.from_file("structure.yaml")
```

**Why JSON/YAML:**
- Complete serialization (preserves all metadata)
- Survives code changes better than pickle
- Human-readable
- Cross-platform compatible

## PDB Files (Protein Data Bank)

```python
from pymatgen.core import Molecule

# Read PDB
molecule = Molecule.from_file("protein.pdb")

# Write PDB
molecule.to(filename="output.pdb")
```

## Gaussian Input/Output

### Gaussian Input

```python
from pymatgen.io.gaussian import GaussianInput

# Create Gaussian input
mol = Molecule.from_file("molecule.xyz")

gin = GaussianInput(
    mol,
    route_parameters={"#": "B3LYP/6-31G(d)", "opt": None},
    title="Optimization job",
    charge=0,
    spin_multiplicity=1
)

gin.write_file("molecule.gjf")
```

### Gaussian Output

```python
from pymatgen.io.gaussian import GaussianOutput

# Parse Gaussian output
gout = GaussianOutput("molecule.log")

# Get final structure
final_structure = gout.final_structure

# Get energies
energies = gout.energies  # All SCF energies
final_energy = gout.final_energy

# Check if converged
is_converged = gout.properly_terminated
```

## Q-Chem Input/Output

```python
from pymatgen.io.qchem.inputs import QCInput
from pymatgen.io.qchem.outputs import QCOutput

# Write Q-Chem input
mol = Molecule.from_file("molecule.xyz")
qcinput = QCInput(mol, jobtype="opt", basis="6-31G*", method="B3LYP")
qcinput.write_file("molecule.qin")

# Parse Q-Chem output
qcout = QCOutput("molecule.qout")
final_structure = qcout.data["molecule_from_optimized_geometry"]
final_energy = qcout.data["final_energy"]
```

## LAMMPS Data Files

```python
from pymatgen.io.lammps.data import LammpsData

# Read LAMMPS data
lammps_data = LammpsData.from_file("lammps.data")
structure = lammps_data.structure

# Write LAMMPS data
lammps_data = LammpsData.from_structure(structure)
lammps_data.write_file("output.data")
```

## PWscf (Quantum ESPRESSO) Files

### PWscf Input

```python
from pymatgen.io.pwscf import PWInput

# Create PWscf input
pwinput = PWInput(
    structure,
    control={
        "calculation": "scf",
        "pseudo_dir": "./pseudo",
    },
    pseudo={"Fe": "Fe.pbe-nd-rrkjus.UPF", "O": "O.pbe-rrkjus.UPF"},
    system={"ecutwfc": 50, "ecutrho": 400}
)

pwinput.write_file("pwscf.in")
```

### PWscf Output

```python
from pymatgen.io.pwscf import PWOutput

# Parse PWscf output
pwout = PWOutput("pwscf.out")
final_structure = pwout.final_structure
final_energy = pwout.final_energy
```

## Format Conversion

### Common Conversions

**CIF to POSCAR:**
```python
structure = Structure.from_file("structure.cif")
structure.to(filename="POSCAR")
```

**POSCAR to CIF:**
```python
structure = Structure.from_file("POSCAR")
structure.to(filename="structure.cif")
```

**Structure to Molecule:**
```python
structure = Structure.from_file("POSCAR")
molecule = structure.to(filename="molecule.xyz")
```

**Molecule to Structure (in box):**
```python
from pymatgen.core import Molecule, Structure, Lattice

molecule = Molecule.from_file("molecule.xyz")

# Create large box
lattice = Lattice.cubic(20.0)
structure = Structure(
    lattice,
    molecule.species,
    molecule.cart_coords,
    coords_are_cartesian=True
)
structure.to(filename="POSCAR")
```

### Batch Conversion

```python
import glob
from pymatgen.core import Structure

# Convert all CIF files to POSCAR
for cif_file in glob.glob("*.cif"):
    structure = Structure.from_file(cif_file)
    basename = cif_file.replace(".cif", "")
    structure.to(filename=f"POSCAR_{basename}")
```

## File I/O Best Practices

### Error Handling

```python
from pymatgen.core import Structure

try:
    structure = Structure.from_file("file.cif")
except FileNotFoundError:
    print("File not found")
except Exception as e:
    print(f"Error reading file: {e}")
```

### Checking File Format

```python
from pathlib import Path

filepath = Path("structure.cif")

if filepath.suffix == ".cif":
    # CIF-specific handling
    pass
elif filepath.name == "POSCAR" or filepath.name == "CONTCAR":
    # VASP handling
    pass
```

### Validating Structures

```python
structure = Structure.from_file("input.cif")

# Check for reasonable structure
if structure.volume < 0:
    print("Error: Negative volume")

if structure.num_sites == 0:
    print("Error: No atoms in structure")

# Check for overlapping atoms
if not structure.is_valid():
    print("Warning: Structure may have overlapping atoms")
```

## Working with Multiple Files

### Reading Multiple Structures

```python
import glob
from pymatgen.core import Structure

structures = []
for filename in glob.glob("structures/*.cif"):
    struct = Structure.from_file(filename)
    structures.append(struct)

# Process all structures
for struct in structures:
    print(f"{struct.composition.reduced_formula}: {struct.volume:.2f} Ų")
```

### Writing Multiple Structures

```python
# Write to separate files
for i, structure in enumerate(structures):
    structure.to(filename=f"structure_{i}.cif")

# Or to single CIF with multiple structures
from pymatgen.io.cif import CifWriter

writer = CifWriter(structures)
writer.write_file("all_structures.cif")
```

## Custom File Parsers

### Creating Custom Parser

```python
from pymatgen.core import Structure, Lattice

def read_custom_format(filename):
    """Read custom structure file format."""
    with open(filename, "r") as f:
        lines = f.readlines()

    # Parse your format
    # ... custom parsing logic ...

    # Create structure
    lattice = Lattice(lattice_matrix)
    structure = Structure(lattice, species, coords)

    return structure

# Use custom parser
structure = read_custom_format("custom.dat")
```

## Troubleshooting

### Issue: Format not recognized

**Problem:** `ValueError: Unable to determine file type`

**Solution:**
```python
# Specify format explicitly
structure = Structure.from_file("file", fmt="poscar")
```

### Issue: Encoding errors

**Problem:** `UnicodeDecodeError` when reading files

**Solution:**
```python
# Try different encodings
structure = Structure.from_file("file.cif", encoding="latin-1")
```

### Issue: Multiple structures in CIF

**Problem:** Only getting first structure

**Solution:**
```python
from pymatgen.io.cif import CifParser

parser = CifParser("multi.cif")
all_structures = parser.get_structures()

for struct in all_structures:
    print(struct.composition)
```

### Issue: POTCAR not found

**Problem:** Cannot write POTCAR files

**Solution:**
```python
# Set environment variable
import os
os.environ["PMG_VASP_PSP_DIR"] = "/path/to/vasp/potentials"

# Or configure globally
# pmg config --add PMG_VASP_PSP_DIR /path/to/vasp/potentials
```

## Quick Reference

| Format | Read | Write |
|--------|------|-------|
| **Any** | `Structure.from_file("file")` | `structure.to("file")` |
| **CIF** | `Structure.from_file("file.cif")` | `structure.to("file.cif")` |
| **POSCAR** | `Structure.from_file("POSCAR")` | `structure.to("POSCAR")` |
| **XYZ** | `Molecule.from_file("file.xyz")` | `molecule.to("file.xyz")` |
| **JSON** | `Structure.from_file("file.json")` | `structure.to("file.json")` |
| **vasprun.xml** | `Vasprun("vasprun.xml")` | N/A (output only) |
| **OUTCAR** | `Outcar("OUTCAR")` | N/A (output only) |

All file I/O supports automatic format detection from filename extensions.
