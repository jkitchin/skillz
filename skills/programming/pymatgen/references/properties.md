# Pymatgen Materials Properties

Guide to calculating and retrieving materials properties using pymatgen.

## Overview

Pymatgen can calculate or retrieve various materials properties:
- **Structural** - Volume, density, distances
- **Chemical** - Composition, electronegativity, oxidation states
- **Electronic** - Band gap, work function
- **Mechanical** - Elastic constants, bulk/shear modulus
- **Magnetic** - Magnetic moments, ordering
- **Thermal** - (via interfaces to phonon codes)

## Structural Properties

### Basic Geometric Properties

```python
from pymatgen.core import Structure

structure = Structure.from_file("POSCAR")

# Lattice parameters
print(f"a = {structure.lattice.a:.3f} Å")
print(f"b = {structure.lattice.b:.3f} Å")
print(f"c = {structure.lattice.c:.3f} Å")
print(f"α = {structure.lattice.alpha:.2f}°")
print(f"β = {structure.lattice.beta:.2f}°")
print(f"γ = {structure.lattice.gamma:.2f}°")

# Volume and density
print(f"Volume: {structure.volume:.3f} Ų")
print(f"Density: {structure.density:.3f} g/cm³")
```

### Distances and Coordination

```python
# Nearest neighbor distances
from pymatgen.analysis.local_env import CrystalNN

cnn = CrystalNN()

for i, site in enumerate(structure):
    nn_info = cnn.get_nn_info(structure, i)
    coord_num = len(nn_info)

    print(f"Site {i} ({site.species_string}):")
    print(f"  Coordination number: {coord_num}")

    if nn_info:
        avg_dist = sum([n['weight'] * structure.get_distance(i, n['site_index'])
                        for n in nn_info]) / coord_num
        print(f"  Average bond length: {avg_dist:.3f} Å")
```

## Composition Properties

### Formula and Stoichiometry

```python
from pymatgen.core import Composition

comp = structure.composition

print(f"Formula: {comp.formula}")
print(f"Reduced formula: {comp.reduced_formula}")
print(f"Alphabetical formula: {comp.alphabetical_formula}")
print(f"Molecular weight: {comp.weight:.3f} g/mol")
print(f"Number of atoms: {comp.num_atoms}")
```

### Elemental Properties

```python
from pymatgen.core import Element

# Average properties
elements = comp.elements
avg_electronegativity = sum([comp.get_atomic_fraction(el) * el.X
                              for el in elements if el.X is not None])

print(f"Average electronegativity: {avg_electronegativity:.2f}")

# Element statistics
for el in elements:
    fraction = comp.get_atomic_fraction(el)
    print(f"{el.symbol}:")
    print(f"  Atomic fraction: {fraction:.3f}")
    print(f"  Electronegativity: {el.X}")
    print(f"  Atomic radius: {el.atomic_radius} Å")
```

### Oxidation States

```python
from pymatgen.analysis.bond_valence import BVAnalyzer

bva = BVAnalyzer()

try:
    # Automatically determine oxidation states
    structure_with_oxi = bva.get_oxi_state_decorated_structure(structure)

    for site in structure_with_oxi:
        species = site.species
        for el, occ in species.items():
            print(f"{el}: {occ}")

except ValueError as e:
    print(f"Could not determine oxidation states: {e}")
```

## Symmetry Properties

### Space Group

```python
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

sga = SpacegroupAnalyzer(structure)

print(f"Space group: {sga.get_space_group_symbol()}")
print(f"Space group number: {sga.get_space_group_number()}")
print(f"Crystal system: {sga.get_crystal_system()}")
print(f"Point group: {sga.get_point_group_symbol()}")
print(f"Hall symbol: {sga.get_hall()}")
```

### Symmetry Operations

```python
# Number of symmetry operations
sym_ops = sga.get_symmetry_operations()
print(f"Number of symmetry operations: {len(sym_ops)}")

# Wyckoff positions
wyckoffs = sga.get_symmetry_dataset()["wyckoffs"]
for i, site in enumerate(structure):
    print(f"Site {i} ({site.species_string}): {wyckoffs[i]}")
```

## Electronic Properties

### Band Gap (from Materials Project)

```python
from mp_api.client import MPRester

with MPRester() as mpr:
    doc = mpr.materials.summary.get_data_by_id("mp-149")

    print(f"Band gap: {doc.band_gap:.3f} eV")
    print(f"Direct gap: {doc.is_gap_direct}")

    # Band structure and DOS
    bs = mpr.get_bandstructure_by_material_id("mp-149")
    if bs:
        gap_info = bs.get_band_gap()
        print(f"VBM: {bs.get_vbm()['energy']:.3f} eV")
        print(f"CBM: {bs.get_cbm()['energy']:.3f} eV")
```

### Work Function (from calculations)

```python
from pymatgen.io.vasp.outputs import Locpot, Vasprun

# Requires slab calculation with LOCPOT
locpot = Locpot.from_file("LOCPOT")
vasprun = Vasprun("vasprun.xml")

# Average potential along surface normal
avg_pot = locpot.get_average_along_axis(2)  # z-axis

# Work function = vacuum_level - fermi_level
vacuum_level = max(avg_pot)
fermi_level = vasprun.efermi

work_function = vacuum_level - fermi_level
print(f"Work function: {work_function:.2f} eV")
```

## Mechanical Properties

### Elastic Properties (from Materials Project)

```python
from mp_api.client import MPRester

with MPRester() as mpr:
    doc = mpr.materials.elasticity.get_data_by_id("mp-149")

    if doc:
        print(f"Bulk modulus (Voigt): {doc.k_voigt:.1f} GPa")
        print(f"Bulk modulus (Reuss): {doc.k_reuss:.1f} GPa")
        print(f"Bulk modulus (VRH): {doc.k_vrh:.1f} GPa")

        print(f"Shear modulus (Voigt): {doc.g_voigt:.1f} GPa")
        print(f"Shear modulus (Reuss): {doc.g_reuss:.1f} GPa")
        print(f"Shear modulus (VRH): {doc.g_vrh:.1f} GPa")

        print(f"Young's modulus: {doc.elastic_anisotropy:.1f} GPa")
        print(f"Poisson ratio: {doc.poisson_ratio:.3f}")
```

### From Elastic Tensor

```python
from pymatgen.analysis.elasticity import ElasticTensor
import numpy as np

# Example elastic tensor (Voigt notation: 6x6 matrix)
# Units: GPa
elastic_matrix = np.array([
    [500, 100, 100, 0, 0, 0],
    [100, 500, 100, 0, 0, 0],
    [100, 100, 500, 0, 0, 0],
    [0, 0, 0, 150, 0, 0],
    [0, 0, 0, 0, 150, 0],
    [0, 0, 0, 0, 0, 150]
])

et = ElasticTensor.from_voigt(elastic_matrix)

# Derived properties
print(f"Bulk modulus (Voigt): {et.k_voigt:.1f} GPa")
print(f"Shear modulus (Voigt): {et.g_voigt:.1f} GPa")
print(f"Young's modulus: {et.y_mod:.1f} GPa")
```

## Magnetic Properties

### From Materials Project

```python
from mp_api.client import MPRester

with MPRester() as mpr:
    docs = mpr.materials.summary.search(
        formula="Fe2O3",
        fields=["material_id", "formula_pretty", "ordering", "total_magnetization"]
    )

    for doc in docs:
        if hasattr(doc, 'ordering'):
            print(f"{doc.formula_pretty}:")
            print(f"  Magnetic ordering: {doc.ordering}")
            if hasattr(doc, 'total_magnetization'):
                print(f"  Total magnetization: {doc.total_magnetization:.2f} μB")
```

### From VASP Calculations

```python
from pymatgen.io.vasp.outputs import Outcar

outcar = Outcar("OUTCAR")

if outcar.total_mag is not None:
    print(f"Total magnetization: {outcar.total_mag:.3f} μB")

    # Site-specific magnetization
    if outcar.magnetization:
        for i, mag in enumerate(outcar.magnetization):
            print(f"Site {i}: {mag['tot']:.3f} μB")
```

## Thermodynamic Properties

### Formation Energy

```python
from mp_api.client import MPRester
from pymatgen.analysis.phase_diagram import PhaseDiagram

with MPRester() as mpr:
    # Get entries
    entries = mpr.get_entries_in_chemsys(["Fe", "O"])

    # Build phase diagram
    pd = PhaseDiagram(entries)

    # Formation energy for each stable phase
    for entry in pd.stable_entries:
        form_e = pd.get_form_energy_per_atom(entry)
        formula = entry.composition.reduced_formula

        print(f"{formula}: ΔHf = {form_e:.3f} eV/atom")
```

### Energy Above Hull (Stability)

```python
from mp_api.client import MPRester

with MPRester() as mpr:
    doc = mpr.materials.summary.get_data_by_id("mp-19017")  # LiFePO4

    print(f"Energy above hull: {doc.energy_above_hull*1000:.2f} meV/atom")

    if doc.energy_above_hull < 0.001:
        print("Stable phase")
    elif doc.energy_above_hull < 0.050:
        print("Metastable phase")
    else:
        print("Unstable phase")
```

## Surface Properties

### Surface Energy

```python
from pymatgen.core.surface import SlabGenerator, WorkFunctionAnalyzer

# Create slab
bulk = Structure.from_file("POSCAR_bulk")
slabgen = SlabGenerator(bulk, miller_index=(1,1,1), min_slab_size=10, min_vacuum_size=15)
slab = slabgen.get_slabs()[0]

# Surface energy = (E_slab - n * E_bulk) / (2 * A)
# Requires energies from DFT calculations
E_slab = -100.0  # eV (from calculation)
E_bulk = -10.0   # eV/atom (from calculation)
n_atoms = len(slab)

surface_energy = (E_slab - n_atoms * E_bulk) / (2 * slab.surface_area)
print(f"Surface energy: {surface_energy:.3f} eV/Ų")

# Convert to J/m²
surface_energy_SI = surface_energy * 16.02176  # eV/Ų to J/m²
print(f"Surface energy: {surface_energy_SI:.3f} J/m²")
```

## Optical Properties

### Dielectric Constant (from VASP)

```python
from pymatgen.io.vasp.outputs import Vasprun

vasprun = Vasprun("vasprun.xml")

# Electronic contribution (high-frequency)
if hasattr(vasprun, 'epsilon_electronic'):
    eps_electronic = vasprun.epsilon_electronic
    print("Electronic dielectric tensor:")
    print(eps_electronic)

    # Average
    eps_avg = (eps_electronic[0][0] + eps_electronic[1][1] + eps_electronic[2][2]) / 3
    print(f"Average electronic dielectric constant: {eps_avg:.2f}")

# Ionic contribution (from DFPT)
if hasattr(vasprun, 'epsilon_ionic'):
    eps_ionic = vasprun.epsilon_ionic
    eps_total = vasprun.epsilon_static

    print(f"Average ionic contribution: {sum(eps_ionic.diagonal())/3:.2f}")
    print(f"Average total (static) dielectric constant: {sum(eps_total.diagonal())/3:.2f}")
```

### Refractive Index

```python
import numpy as np

# From electronic dielectric constant
eps_avg = 10.0  # Average electronic dielectric constant

# Refractive index n = sqrt(ε)
n = np.sqrt(eps_avg)
print(f"Refractive index: {n:.2f}")
```

## Transport Properties

### Electrical Conductivity (qualitative)

```python
from mp_api.client import MPRester

with MPRester() as mpr:
    doc = mpr.materials.summary.get_data_by_id("mp-149")

    # Qualitative estimate from band gap
    if doc.band_gap == 0:
        print("Metallic conductor")
    elif doc.band_gap < 0.5:
        print("Semimetal or narrow-gap semiconductor")
    elif doc.band_gap < 3.0:
        print("Semiconductor")
    else:
        print("Insulator")
```

## Property Trends

### Analyze Property Trends

```python
from mp_api.client import MPRester
import matplotlib.pyplot as plt

with MPRester() as mpr:
    # Get all materials in chemical system
    docs = mpr.materials.summary.search(
        chemsys="Li-Fe-O",
        fields=["formula_pretty", "volume", "energy_per_atom"]
    )

    # Extract data
    formulas = [doc.formula_pretty for doc in docs]
    volumes = [doc.volume for doc in docs]
    energies = [doc.energy_per_atom for doc in docs]

    # Plot
    plt.figure(figsize=(10, 6))
    plt.scatter(volumes, energies)
    for i, formula in enumerate(formulas):
        plt.annotate(formula, (volumes[i], energies[i]))

    plt.xlabel("Volume per atom (Ų)")
    plt.ylabel("Energy per atom (eV)")
    plt.title("Volume-Energy relationship in Li-Fe-O system")
    plt.tight_layout()
    plt.savefig("volume_energy.png", dpi=300)
```

### Compare Properties

```python
# Compare properties across related materials
materials = {
    "Si": "mp-149",
    "Ge": "mp-32",
    "GaAs": "mp-2534"
}

properties = []

with MPRester() as mpr:
    for name, mp_id in materials.items():
        doc = mpr.materials.summary.get_data_by_id(mp_id)

        properties.append({
            'material': name,
            'band_gap': doc.band_gap,
            'volume': doc.volume / len(doc.structure),  # per atom
            'density': doc.density
        })

# Print comparison
print(f"{'Material':<10} {'Gap (eV)':<10} {'Vol/atom':<12} {'Density':<10}")
print("-" * 50)
for prop in properties:
    print(f"{prop['material']:<10} {prop['band_gap']:<10.2f} "
          f"{prop['volume']:<12.3f} {prop['density']:<10.3f}")
```

## Custom Property Calculations

### Packing Fraction

```python
from pymatgen.core import Element
import numpy as np

structure = Structure.from_file("POSCAR")

# Calculate packing fraction (sphere packing)
total_atomic_volume = 0
for site in structure:
    el = Element(site.species_string)
    # Use covalent radius or atomic radius
    if hasattr(el, 'atomic_radius') and el.atomic_radius:
        r = el.atomic_radius
    else:
        r = 1.5  # Default
    vol_sphere = (4/3) * np.pi * r**3
    total_atomic_volume += vol_sphere

packing_fraction = total_atomic_volume / structure.volume
print(f"Packing fraction: {packing_fraction:.3f}")
```

### Cohesive Energy

```python
# Cohesive energy = E_atoms - E_solid
# Requires total energies from calculations

E_solid = -50.0  # eV (per atom, from bulk calculation)

# Elemental reference energies (from isolated atom calculations)
E_atom = {"Li": -1.0, "O": -5.0}  # eV (example values)

# Calculate cohesive energy
comp = structure.composition
cohesive_energy = 0
for el, amt in comp.items():
    cohesive_energy += amt * (E_atom[str(el)] - E_solid)

cohesive_energy /= comp.num_atoms
print(f"Cohesive energy: {cohesive_energy:.3f} eV/atom")
```

## Property Databases

### Query Multiple Properties

```python
from mp_api.client import MPRester

with MPRester() as mpr:
    docs = mpr.materials.summary.search(
        elements=["Li", "Co", "O"],
        fields=[
            "material_id",
            "formula_pretty",
            "band_gap",
            "formation_energy_per_atom",
            "energy_above_hull",
            "density",
            "volume",
            "symmetry"
        ]
    )

    print(f"{'Formula':<15} {'Gap(eV)':<10} {'ΔHf':<10} {'E_hull':<10} {'ρ(g/cm³)':<10}")
    print("-" * 65)

    for doc in docs[:10]:  # First 10
        print(f"{doc.formula_pretty:<15} {doc.band_gap:<10.2f} "
              f"{doc.formation_energy_per_atom:<10.3f} "
              f"{doc.energy_above_hull*1000:<10.1f} {doc.density:<10.3f}")
```

## Best Practices

### Units

Pymatgen typically uses:
- **Energy**: eV
- **Length**: Å (Angstrom)
- **Volume**: Ų
- **Density**: g/cm³
- **Pressure**: GPa
- **Magnetic moment**: μB (Bohr magneton)

Always verify units when comparing with other sources.

### Accuracy Considerations

```python
# DFT accuracy limitations
print("Typical DFT accuracies:")
print("  Lattice parameters: ±1-2%")
print("  Formation energies: ±0.1-0.2 eV/atom")
print("  Band gaps: Underestimated by 30-50% (PBE)")
print("  Elastic constants: ±5-10%")
print("  Magnetic moments: ±10-20%")
```

### Property Validation

```python
# Always validate computed properties
def validate_structure_properties(structure):
    """Check if structure properties are reasonable."""
    issues = []

    # Volume check
    vol_per_atom = structure.volume / len(structure)
    if vol_per_atom < 5 or vol_per_atom > 100:
        issues.append(f"Unusual volume per atom: {vol_per_atom:.1f} Ų")

    # Density check
    if structure.density < 0.1 or structure.density > 25:
        issues.append(f"Unusual density: {structure.density:.1f} g/cm³")

    # Bond length check
    min_dist = min([structure.get_distance(i, j)
                    for i in range(len(structure))
                    for j in range(i+1, len(structure))])
    if min_dist < 0.8:
        issues.append(f"Very short bond: {min_dist:.2f} Å")

    return issues

# Use
issues = validate_structure_properties(structure)
if issues:
    print("Warnings:")
    for issue in issues:
        print(f"  - {issue}")
```

## Quick Reference

### Common Properties

| Property | Source | Method |
|----------|--------|--------|
| Volume | Structure | `structure.volume` |
| Density | Structure | `structure.density` |
| Band gap | Materials Project | `doc.band_gap` |
| Formation energy | Materials Project | `doc.formation_energy_per_atom` |
| Stability | Materials Project | `doc.energy_above_hull` |
| Space group | SpacegroupAnalyzer | `sga.get_space_group_symbol()` |
| Elastic moduli | Materials Project | `doc.k_vrh`, `doc.g_vrh` |
| Magnetization | OUTCAR | `outcar.total_mag` |
| Dielectric | Vasprun | `vasprun.epsilon_electronic` |

### Typical Values

**Band gaps:**
- Metals: 0 eV
- Semiconductors: 0.5-3 eV
- Insulators: > 3 eV

**Densities:**
- Light elements (Li, C): 0.5-3 g/cm³
- Common oxides: 3-6 g/cm³
- Transition metals: 7-12 g/cm³
- Heavy elements (Au, Pt): 15-25 g/cm³

**Elastic moduli:**
- Soft materials: < 50 GPa
- Common ceramics: 100-300 GPa
- Hard materials (SiC, diamond): > 400 GPa
