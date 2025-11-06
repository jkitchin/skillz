# Pymatgen VASP Integration

Comprehensive guide to generating VASP input files and parsing VASP output files using pymatgen.

## Overview

Pymatgen provides comprehensive VASP integration:
- **Input generation** - Create INCAR, POSCAR, KPOINTS, POTCAR
- **Input sets** - Pre-configured sets for common calculations
- **Output parsing** - Parse vasprun.xml, OUTCAR, CONTCAR, etc.
- **Workflow management** - Chain multiple calculations

## VASP Input Files

### POSCAR (Structure)

**Write POSCAR:**
```python
from pymatgen.core import Structure
from pymatgen.io.vasp import Poscar

structure = Structure.from_file("structure.cif")

# Write POSCAR
poscar = Poscar(structure)
poscar.write_file("POSCAR")

# Or directly
structure.to(filename="POSCAR")
```

**With selective dynamics:**
```python
from pymatgen.io.vasp import Poscar

# Mark which atoms can move
selective_dynamics = [
    [True, True, True],    # Site 0: all directions
    [False, False, False], # Site 1: fixed
    [True, True, False],   # Site 2: only x,y
]

structure.add_site_property("selective_dynamics", selective_dynamics)

poscar = Poscar(structure)
poscar.write_file("POSCAR")
```

**Read POSCAR:**
```python
from pymatgen.io.vasp import Poscar

poscar = Poscar.from_file("POSCAR")
structure = poscar.structure

print(f"Formula: {structure.composition.reduced_formula}")
print(f"Lattice: {structure.lattice.parameters}")

# Check for selective dynamics
if poscar.selective_dynamics is not None:
    print("Selective dynamics enabled")
```

### INCAR (Parameters)

**Create INCAR:**
```python
from pymatgen.io.vasp.inputs import Incar

# Define parameters
incar = Incar({
    "PREC": "Accurate",
    "ENCUT": 520,
    "EDIFF": 1e-5,
    "NSW": 100,
    "IBRION": 2,
    "ISIF": 3,
    "ISMEAR": 0,
    "SIGMA": 0.05,
    "LREAL": False
})

# Write to file
incar.write_file("INCAR")
```

**Read INCAR:**
```python
incar = Incar.from_file("INCAR")

# Access parameters
print(f"ENCUT: {incar['ENCUT']}")
print(f"EDIFF: {incar['EDIFF']}")

# Check if parameter exists
if "ISPIN" in incar:
    print(f"ISPIN = {incar['ISPIN']}")
else:
    print("ISPIN not set")
```

**Modify INCAR:**
```python
# Read existing INCAR
incar = Incar.from_file("INCAR")

# Modify parameters
incar["ENCUT"] = 600
incar["NSW"] = 200

# Add new parameters
incar["LCHARG"] = False
incar["LWAVE"] = False

# Remove parameter
if "LREAL" in incar:
    del incar["LREAL"]

# Write modified INCAR
incar.write_file("INCAR_new")
```

### KPOINTS

**Automatic k-points:**
```python
from pymatgen.io.vasp.inputs import Kpoints

# Automatic mesh (density-based)
kpoints = Kpoints.automatic(100)  # 100 k-points per Ų
kpoints.write_file("KPOINTS")
```

**Gamma-centered:**
```python
# Gamma-centered mesh
kpoints = Kpoints.gamma_automatic([4, 4, 4])
kpoints.write_file("KPOINTS")
```

**Monkhorst-Pack:**
```python
# Monkhorst-Pack mesh
kpoints = Kpoints.monkhorst_automatic([6, 6, 6])
kpoints.write_file("KPOINTS")
```

**Automatic density:**
```python
from pymatgen.io.vasp.inputs import Kpoints

# Automatic based on structure
structure = Structure.from_file("POSCAR")
kpoints = Kpoints.automatic_density(structure, kppa=1000)  # k-points per atom
kpoints.write_file("KPOINTS")
```

**Line mode (band structure):**
```python
from pymatgen.io.vasp.inputs import Kpoints
from pymatgen.symmetry.bandstructure import HighSymmKpath

structure = Structure.from_file("POSCAR")

# Get high-symmetry path
kpath = HighSymmKpath(structure)

# Create line-mode KPOINTS
kpoints = Kpoints.automatic_linemode(20, kpath)  # 20 points per segment
kpoints.write_file("KPOINTS")
```

**Explicit k-points:**
```python
# Manual k-point specification
kpoints = Kpoints(
    kpts=[[0.0, 0.0, 0.0], [0.5, 0.0, 0.0], [0.5, 0.5, 0.0]],
    kpts_weights=[1, 2, 1],
    style=Kpoints.supported_modes.Reciprocal
)
kpoints.write_file("KPOINTS")
```

### POTCAR

**Generate POTCAR:**
```python
from pymatgen.io.vasp.inputs import Potcar

# Requires PMG_VASP_PSP_DIR configured
# Or set: os.environ["PMG_VASP_PSP_DIR"] = "/path/to/vasp/potcars"

# Create POTCAR
potcar = Potcar(["Fe", "O"], functional="PBE")
potcar.write_file("POTCAR")

# Or from structure
structure = Structure.from_file("POSCAR")
symbols = [site.species_string for site in structure]
unique_symbols = list(dict.fromkeys(symbols))  # Maintain order
potcar = Potcar(unique_symbols, functional="PBE")
potcar.write_file("POTCAR")
```

**Read POTCAR:**
```python
potcar = Potcar.from_file("POTCAR")

print(f"Elements: {potcar.symbols}")
print(f"Functional: {potcar.functional}")

# Get details for each element
for i, symbol in enumerate(potcar.symbols):
    single_potcar = potcar[i]
    print(f"{symbol}: ENMAX = {single_potcar.enmax} eV")
```

## VASP Input Sets

### Pre-configured Input Sets

**Materials Project relaxation:**
```python
from pymatgen.io.vasp.sets import MPRelaxSet

structure = Structure.from_file("structure.cif")

# Create input set
input_set = MPRelaxSet(structure)

# Write all input files
input_set.write_input("relax_calc")

# This creates:
# - relax_calc/POSCAR
# - relax_calc/INCAR
# - relax_calc/KPOINTS
# - relax_calc/POTCAR
```

**Static calculation:**
```python
from pymatgen.io.vasp.sets import MPStaticSet

# Use relaxed structure
structure = Structure.from_file("CONTCAR")

input_set = MPStaticSet(structure)
input_set.write_input("static_calc")
```

**Band structure calculation:**
```python
from pymatgen.io.vasp.sets import MPNonSCFSet

# Requires previous static calculation
prev_vasprun = "static_calc/vasprun.xml"

input_set = MPNonSCFSet.from_prev_calc(
    prev_calc_dir="static_calc",
    mode="line"  # Line mode for band structure
)

input_set.write_input("bands_calc")
```

**DOS calculation:**
```python
input_set = MPNonSCFSet.from_prev_calc(
    prev_calc_dir="static_calc",
    mode="uniform"  # Uniform for DOS
)

input_set.write_input("dos_calc")
```

### Customizing Input Sets

**Override INCAR settings:**
```python
from pymatgen.io.vasp.sets import MPRelaxSet

# Custom INCAR parameters
user_incar_settings = {
    "ENCUT": 600,
    "EDIFF": 1e-6,
    "ALGO": "Fast",
    "LCHARG": False
}

input_set = MPRelaxSet(
    structure,
    user_incar_settings=user_incar_settings
)

input_set.write_input("custom_relax")
```

**Custom k-point density:**
```python
input_set = MPRelaxSet(
    structure,
    user_kpoints_settings={"reciprocal_density": 200}  # k-points per Ų
)
```

**Magnetic calculations:**
```python
from pymatgen.io.vasp.sets import MPRelaxSet

# For magnetic materials
structure.add_site_property("magmom", [5.0, -5.0, 0.0, 0.0])  # Initial magnetic moments

input_set = MPRelaxSet(structure)
input_set.write_input("magnetic_relax")

# INCAR will include ISPIN=2 and MAGMOM
```

### Available Input Sets

```python
from pymatgen.io.vasp.sets import (
    MPRelaxSet,         # Structure relaxation
    MPStaticSet,        # Static calculation
    MPNonSCFSet,        # Non-SCF (bands, DOS)
    MPSOCSet,           # Spin-orbit coupling
    MPScanRelaxSet,     # SCAN functional
    MPHSERelaxSet,      # HSE06 functional
    MPMetalRelaxSet,    # For metals
    MPNMRSet,           # NMR calculations
    MVLRelaxSet,        # Legacy MVL
)
```

## VASP Output Files

### vasprun.xml

**Parse vasprun.xml:**
```python
from pymatgen.io.vasp.outputs import Vasprun

vasprun = Vasprun("vasprun.xml")

# Basic info
print(f"Final energy: {vasprun.final_energy:.6f} eV")
print(f"Converged: {vasprun.converged}")

# Structure
final_structure = vasprun.final_structure
initial_structure = vasprun.initial_structure

# Electronic structure
dos = vasprun.complete_dos
if vasprun.bands:
    bs = vasprun.get_band_structure()
```

**Parse options:**
```python
# For large files, skip some parsing
vasprun = Vasprun(
    "vasprun.xml",
    parse_dos=False,       # Skip DOS
    parse_eigen=False,     # Skip eigenvalues
    parse_projected_eigen=False  # Skip projections
)

# Just get final energy (fastest)
final_energy = vasprun.final_energy
```

**Ionic steps:**
```python
# Access all ionic steps
for i, step in enumerate(vasprun.ionic_steps):
    energy = step["e_wo_entrp"]  # Energy without entropy
    structure = step["structure"]
    forces = step["forces"]

    print(f"Step {i}: E = {energy:.4f} eV")

    # Max force
    max_force = max([sum([f**2 for f in force])**0.5 for force in forces])
    print(f"  Max force: {max_force:.4f} eV/Å")
```

**Parameters used:**
```python
# Get INCAR parameters used
incar = vasprun.incar

print(f"ENCUT: {incar['ENCUT']}")
print(f"EDIFF: {incar['EDIFF']}")

# Get KPOINTS
kpoints = vasprun.kpoints
print(f"K-points: {kpoints}")
```

### OUTCAR

**Parse OUTCAR:**
```python
from pymatgen.io.vasp.outputs import Outcar

outcar = Outcar("OUTCAR")

# Final energy
print(f"Final energy: {outcar.final_energy:.6f} eV")

# Forces
final_forces = outcar.final_fr
print("Final forces (eV/Å):")
for i, force in enumerate(final_forces):
    print(f"  Atom {i}: {force}")

# Stress
final_stress = outcar.final_stress
print(f"Stress tensor:\n{final_stress}")
```

**Magnetization:**
```python
if outcar.total_mag is not None:
    print(f"Total magnetization: {outcar.total_mag:.3f} μB")

# Site-specific magnetization
if outcar.magnetization is not None:
    for i, mag in enumerate(outcar.magnetization):
        print(f"Site {i}: {mag['tot']:.3f} μB")
```

**Timing:**
```python
# Get computation time
print(f"Elapsed time: {outcar.run_stats['elapsed_time']:.2f} s")
print(f"User time: {outcar.run_stats['user_time']:.2f} s")
```

### CONTCAR

```python
from pymatgen.core import Structure

# Read final structure
final_structure = Structure.from_file("CONTCAR")

print(f"Final formula: {final_structure.composition.reduced_formula}")
print(f"Final volume: {final_structure.volume:.3f} Ų")

# Use as input for next calculation
final_structure.to(filename="POSCAR_next")
```

### DOSCAR and EIGENVAL

**Parse DOSCAR:**
```python
from pymatgen.io.vasp.outputs import Doscar

doscar = Doscar("DOSCAR")
dos = doscar.completedos

print(f"Fermi energy: {dos.efermi:.3f} eV")
```

**Parse EIGENVAL:**
```python
from pymatgen.io.vasp.outputs import Eigenval

eigenval = Eigenval("EIGENVAL")

# Get band structure
bs = eigenval.get_band_structure()
```

### LOCPOT and CHGCAR

**Parse LOCPOT:**
```python
from pymatgen.io.vasp.outputs import Locpot

locpot = Locpot.from_file("LOCPOT")

# Get planar average
avg_potential = locpot.get_average_along_axis(2)  # Along z

# Plot
import matplotlib.pyplot as plt
plt.plot(avg_potential)
plt.xlabel("z (Å)")
plt.ylabel("Potential (eV)")
plt.savefig("potential_profile.png")
```

**Parse CHGCAR:**
```python
from pymatgen.io.vasp.outputs import Chgcar

chgcar = Chgcar.from_file("CHGCAR")

# Get charge density
charge_density = chgcar.data['total']

# Get structure
structure = chgcar.structure
```

## Calculation Workflows

### Sequential Calculations

**Relaxation → Static → Band Structure:**
```python
from pymatgen.io.vasp.sets import MPRelaxSet, MPStaticSet, MPNonSCFSet
from pymatgen.core import Structure
import os

structure = Structure.from_file("structure.cif")

# 1. Relaxation
relax_set = MPRelaxSet(structure)
relax_set.write_input("01_relax")

# Run VASP in 01_relax/
# ... (external VASP execution)

# 2. Static calculation with relaxed structure
relaxed_structure = Structure.from_file("01_relax/CONTCAR")
static_set = MPStaticSet(relaxed_structure)
static_set.write_input("02_static")

# Run VASP in 02_static/

# 3. Band structure
bs_set = MPNonSCFSet.from_prev_calc(
    prev_calc_dir="02_static",
    mode="line"
)
bs_set.write_input("03_bands")

# Run VASP in 03_bands/
```

### Convergence Testing

**ENCUT convergence:**
```python
from pymatgen.io.vasp.sets import MPRelaxSet

structure = Structure.from_file("POSCAR")

encut_values = [400, 450, 500, 550, 600, 650]

for encut in encut_values:
    user_settings = {"ENCUT": encut, "NSW": 0}  # Static

    input_set = MPRelaxSet(
        structure,
        user_incar_settings=user_settings
    )

    dirname = f"encut_{encut}"
    input_set.write_input(dirname)

# Run all calculations and compare energies
```

**K-point convergence:**
```python
kpoint_densities = [100, 200, 400, 600, 800, 1000]

for kppa in kpoint_densities:
    user_kpoints = {"reciprocal_density": kppa}

    input_set = MPStaticSet(
        structure,
        user_kpoints_settings=user_kpoints
    )

    dirname = f"kpoints_{kppa}"
    input_set.write_input(dirname)
```

## Error Analysis

### Check Convergence

```python
from pymatgen.io.vasp.outputs import Vasprun

vasprun = Vasprun("vasprun.xml")

if vasprun.converged:
    print("Calculation converged")
else:
    print("WARNING: Calculation did not converge")

# Check electronic convergence
if vasprun.converged_electronic:
    print("Electronic structure converged")

# Check ionic convergence
if vasprun.converged_ionic:
    print("Ionic structure converged")
```

### Parse Errors from OUTCAR

```python
from pymatgen.io.vasp.outputs import Outcar

outcar = Outcar("OUTCAR")

# Check for errors
if hasattr(outcar, 'errors'):
    if outcar.errors:
        print("Errors found:")
        for error in outcar.errors:
            print(f"  - {error}")
```

### Validate Results

```python
from pymatgen.io.vasp.outputs import Vasprun

vasprun = Vasprun("vasprun.xml")

# Check if calculation finished
if not vasprun.converged:
    print("Not converged")

# Check for imaginary frequencies (phonons)
# Check negative eigenvalues (WAVECAR)
# Check for charge sloshing (OSZICAR)

# Validate structure
final_structure = vasprun.final_structure
if final_structure.volume < 0:
    print("ERROR: Negative volume")

# Check forces
outcar = Outcar("OUTCAR")
max_force = max([sum([f**2 for f in force])**0.5
                 for force in outcar.final_fr])

if max_force > 0.05:
    print(f"Large residual forces: {max_force:.4f} eV/Å")
```

## Advanced Topics

### Restart Calculations

```python
from pymatgen.io.vasp.sets import MPRelaxSet
import shutil

# Copy necessary files
shutil.copy("old_calc/CONTCAR", "restart/POSCAR")
shutil.copy("old_calc/WAVECAR", "restart/WAVECAR")
shutil.copy("old_calc/CHGCAR", "restart/CHGCAR")

# Create INCAR, KPOINTS, POTCAR
structure = Structure.from_file("restart/POSCAR")
input_set = MPRelaxSet(structure)

# Write only INCAR, KPOINTS, POTCAR (POSCAR already there)
input_set.incar.write_file("restart/INCAR")
input_set.kpoints.write_file("restart/KPOINTS")
input_set.potcar.write_file("restart/POTCAR")
```

### Spin-Orbit Coupling

```python
from pymatgen.io.vasp.sets import MPSOCSet

structure = Structure.from_file("POSCAR")

# SOC calculation
soc_set = MPSOCSet(structure)
soc_set.write_input("soc_calc")

# INCAR will include LSORBIT=True
```

### Hybrid Functionals

```python
from pymatgen.io.vasp.sets import MPHSERelaxSet

structure = Structure.from_file("POSCAR")

# HSE06 calculation
hse_set = MPHSERelaxSet(structure)
hse_set.write_input("hse_calc")

# Note: HSE is expensive, use smaller k-point mesh
```

### NEB (Nudged Elastic Band)

```python
from pymatgen.io.vasp.sets import MPNEBSet

# Create images (interpolated structures)
start = Structure.from_file("POSCAR_start")
end = Structure.from_file("POSCAR_end")

images = start.interpolate(end, nimages=5)

# Create NEB input
neb_set = MPNEBSet(images)
neb_set.write_input("neb_calc")

# Creates 00/, 01/, 02/, ..., with POSCAR in each
```

## Best Practices

### Input File Organization

```python
# Organize calculations in directories
import os

base_dir = "calculations"
os.makedirs(base_dir, exist_ok=True)

# Structure
os.makedirs(f"{base_dir}/01_relax", exist_ok=True)
os.makedirs(f"{base_dir}/02_static", exist_ok=True)
os.makedirs(f"{base_dir}/03_bands", exist_ok=True)

# Write inputs
relax_set.write_input(f"{base_dir}/01_relax")
# etc.
```

### Validation

```python
# Always validate inputs before running
from pymatgen.io.vasp.sets import MPRelaxSet

input_set = MPRelaxSet(structure)

# Check INCAR
incar = input_set.incar
assert incar["ENCUT"] > 400, "ENCUT too low"
assert incar["EDIFF"] < 1e-4, "EDIFF not tight enough"

# Check KPOINTS
kpoints = input_set.kpoints
# Validate k-point density

# Check structure
assert structure.volume > 0, "Invalid structure"
```

### Documentation

```python
# Save metadata with calculations
import json
from datetime import datetime

metadata = {
    "date": datetime.now().isoformat(),
    "structure": "LiFePO4",
    "calculation_type": "relaxation",
    "notes": "Testing new pseudopotentials"
}

with open("01_relax/metadata.json", "w") as f:
    json.dump(metadata, f, indent=2)
```

## Quick Reference

### Input Generation

```python
from pymatgen.io.vasp.sets import MPRelaxSet, MPStaticSet, MPNonSCFSet

# Relaxation
MPRelaxSet(structure).write_input("relax")

# Static
MPStaticSet(structure).write_input("static")

# Band structure
MPNonSCFSet.from_prev_calc("static", mode="line").write_input("bands")
```

### Output Parsing

```python
from pymatgen.io.vasp.outputs import Vasprun, Outcar

# Energy and structure
vasprun = Vasprun("vasprun.xml")
energy = vasprun.final_energy
structure = vasprun.final_structure

# Forces
outcar = Outcar("OUTCAR")
forces = outcar.final_fr
```

### Common Input Sets

| Calculation Type | Input Set |
|-----------------|-----------|
| Structure relaxation | `MPRelaxSet` |
| Static (single point) | `MPStaticSet` |
| Band structure | `MPNonSCFSet(..., mode="line")` |
| DOS | `MPNonSCFSet(..., mode="uniform")` |
| HSE06 | `MPHSERelaxSet` |
| SOC | `MPSOCSet` |
| Magnetic | `MPRelaxSet` with magmoms |
| Metal | `MPMetalRelaxSet` |

### Output Files

| File | Parser | Purpose |
|------|--------|---------|
| vasprun.xml | `Vasprun` | Complete output (energy, structure, DOS, bands) |
| OUTCAR | `Outcar` | Forces, stress, magnetization |
| CONTCAR | `Structure.from_file` | Final structure |
| POSCAR | `Poscar` | Initial structure |
| DOSCAR | `Doscar` | Density of states |
| EIGENVAL | `Eigenval` | Eigenvalues |
| LOCPOT | `Locpot` | Local potential |
| CHGCAR | `Chgcar` | Charge density |
