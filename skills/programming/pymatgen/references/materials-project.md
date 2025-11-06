# Materials Project API Integration

Comprehensive guide to accessing the Materials Project database using pymatgen's MPRester.

## Overview

The Materials Project (MP) is a database of computed materials properties from DFT calculations. Pymatgen's `MPRester` class provides programmatic access to:
- Crystal structures (>140,000 materials)
- Computed properties (formation energies, band gaps, elasticity, etc.)
- Phase diagrams
- Calculated data (DOS, band structures)
- Experimental data

## Setup

### Get API Key

1. Register at [materialsproject.org](https://materialsproject.org)
2. Get your API key from your dashboard
3. Configure pymatgen:

```bash
pmg config --add PMG_MAPI_KEY your_api_key_here
```

Or set environment variable:
```bash
export PMG_MAPI_KEY="your_api_key_here"
```

### Initialize MPRester

```python
from mp_api.client import MPRester

# With configured API key
with MPRester() as mpr:
    # Your queries here
    pass

# Or pass API key directly
with MPRester(api_key="your_key_here") as mpr:
    # Your queries here
    pass
```

**Always use context manager (`with`) to ensure proper resource cleanup.**

## Basic Queries

### Get Structure by Material ID

```python
from mp_api.client import MPRester

with MPRester() as mpr:
    # Get structure for specific material
    structure = mpr.get_structure_by_material_id("mp-149")  # Silicon

    print(f"Formula: {structure.composition.reduced_formula}")
    print(f"Space group: {structure.get_space_group_info()}")
    structure.to(filename="POSCAR_Si")
```

### Search by Formula

```python
with MPRester() as mpr:
    # Get all entries for a formula
    docs = mpr.materials.summary.search(formula="Fe2O3")

    print(f"Found {len(docs)} entries for Fe2O3")

    for doc in docs:
        print(f"mp-id: {doc.material_id}")
        print(f"  Formula: {doc.formula_pretty}")
        print(f"  Space group: {doc.symmetry.symbol}")
        print(f"  Energy: {doc.energy_per_atom:.3f} eV/atom")
```

### Get Most Stable Structure

```python
with MPRester() as mpr:
    # Get most stable polymorph
    docs = mpr.materials.summary.search(
        formula="TiO2",
        fields=["material_id", "formula_pretty", "energy_per_atom", "structure"]
    )

    # Sort by energy
    docs_sorted = sorted(docs, key=lambda x: x.energy_per_atom)
    most_stable = docs_sorted[0]

    print(f"Most stable TiO2: {most_stable.material_id}")
    structure = most_stable.structure
    structure.to(filename="POSCAR_TiO2_stable")
```

## Advanced Searches

### Search by Criteria

```python
with MPRester() as mpr:
    # Search with multiple criteria
    docs = mpr.materials.summary.search(
        elements=["Li", "Fe", "P", "O"],  # Must contain these elements
        num_elements=4,                    # Exactly 4 elements
        fields=["material_id", "formula_pretty", "energy_per_atom"]
    )

    print(f"Found {len(docs)} LiFePO4-type materials")
    for doc in docs:
        print(f"{doc.material_id}: {doc.formula_pretty}")
```

### Search by Properties

```python
with MPRester() as mpr:
    # Find materials with specific properties
    docs = mpr.materials.summary.search(
        band_gap=(1.0, 3.0),               # Band gap between 1-3 eV
        energy_above_hull=(0, 0.01),       # Stable or near-stable
        elements=["Ga", "N"],
        fields=["material_id", "formula_pretty", "band_gap", "energy_above_hull"]
    )

    for doc in docs:
        print(f"{doc.formula_pretty}: Eg = {doc.band_gap:.2f} eV, "
              f"E_hull = {doc.energy_above_hull*1000:.1f} meV/atom")
```

### Search by Composition

```python
with MPRester() as mpr:
    # Chemical system (any composition in Li-Fe-P-O)
    docs = mpr.materials.summary.search(
        chemsys="Li-Fe-P-O",
        fields=["material_id", "formula_pretty", "energy_per_atom"]
    )

    print(f"Found {len(docs)} materials in Li-Fe-P-O system")
```

### Search by Space Group

```python
with MPRester() as mpr:
    # Find materials in specific space group
    docs = mpr.materials.summary.search(
        formula="Fe2O3",
        spacegroup_symbol="R-3c",  # Corundum structure
        fields=["material_id", "formula_pretty", "symmetry"]
    )

    for doc in docs:
        print(f"{doc.material_id}: {doc.symmetry.symbol}")
```

## Retrieving Properties

### Formation Energy

```python
with MPRester() as mpr:
    doc = mpr.materials.summary.get_data_by_id("mp-149")

    print(f"Formation energy: {doc.formation_energy_per_atom:.3f} eV/atom")
    print(f"Energy above hull: {doc.energy_above_hull:.3f} eV/atom")

    # Is it stable?
    if doc.is_stable:
        print("This is a stable phase")
```

### Band Gap

```python
with MPRester() as mpr:
    # Search semiconductors
    docs = mpr.materials.summary.search(
        elements=["Ga", "As"],
        band_gap=(0.1, None),  # Non-zero band gap
        fields=["material_id", "formula_pretty", "band_gap", "is_gap_direct"]
    )

    for doc in docs:
        gap_type = "direct" if doc.is_gap_direct else "indirect"
        print(f"{doc.formula_pretty}: {doc.band_gap:.2f} eV ({gap_type})")
```

### Elasticity

```python
with MPRester() as mpr:
    # Get elastic properties
    doc = mpr.materials.elasticity.get_data_by_id("mp-149")

    if doc:
        print(f"Bulk modulus: {doc.k_voigt:.1f} GPa")
        print(f"Shear modulus: {doc.g_voigt:.1f} GPa")
        print(f"Young's modulus: {doc.universal_anisotropy:.1f}")
```

### Magnetic Properties

```python
with MPRester() as mpr:
    # Search magnetic materials
    docs = mpr.materials.summary.search(
        elements=["Fe", "O"],
        ordering="FM",  # Ferromagnetic
        fields=["material_id", "formula_pretty", "ordering", "total_magnetization"]
    )

    for doc in docs:
        print(f"{doc.formula_pretty}: {doc.total_magnetization:.2f} μB/cell")
```

## Electronic Structure

### Get Band Structure

```python
with MPRester() as mpr:
    # Get band structure
    bs = mpr.get_bandstructure_by_material_id("mp-149")

    if bs:
        print(f"Band gap: {bs.get_band_gap()['energy']:.3f} eV")
        print(f"Direct gap: {bs.get_band_gap()['direct']}")
        print(f"VBM: {bs.get_vbm()['energy']:.3f} eV")
        print(f"CBM: {bs.get_cbm()['energy']:.3f} eV")

        # Plot band structure
        from pymatgen.electronic_structure.plotter import BSPlotter
        plotter = BSPlotter(bs)
        plotter.save_plot("band_structure.png")
```

### Get Density of States

```python
with MPRester() as mpr:
    # Get DOS
    dos = mpr.get_dos_by_material_id("mp-149")

    if dos:
        print(f"Fermi level: {dos.efermi:.3f} eV")

        # Plot DOS
        from pymatgen.electronic_structure.plotter import DosPlotter
        plotter = DosPlotter()
        plotter.add_dos("Total", dos)
        plotter.save_plot("dos.png")
```

## Thermodynamics

### Phase Diagrams

```python
with MPRester() as mpr:
    # Get entries for phase diagram
    entries = mpr.get_entries_in_chemsys(["Li", "Fe", "P", "O"])

    # Build phase diagram
    from pymatgen.analysis.phase_diagram import PhaseDiagram
    pd = PhaseDiagram(entries)

    # Stability of specific composition
    from pymatgen.core import Composition
    comp = Composition("LiFePO4")
    decomp, ehull = pd.get_decomp_and_e_above_hull(entries[0])

    print(f"Energy above hull: {ehull*1000:.1f} meV/atom")

    if ehull < 0.001:
        print("This composition is stable")
```

### Pourbaix Diagrams

```python
with MPRester() as mpr:
    # Get entries for Pourbaix diagram
    entries = mpr.get_pourbaix_entries(["Fe", "O"])

    # Build Pourbaix diagram
    from pymatgen.analysis.pourbaix_diagram import PourbaixDiagram
    pb = PourbaixDiagram(entries)

    # Plot
    from pymatgen.analysis.pourbaix_diagram import PourbaixPlotter
    plotter = PourbaixPlotter(pb)
    plotter.save_plot("pourbaix.png")
```

## XAS Spectra

```python
with MPRester() as mpr:
    # Get XAS spectrum
    xas = mpr.xas.get_data_by_id("mp-149")

    if xas:
        # Access spectrum data
        for doc in xas:
            print(f"Edge: {doc.edge}")
            print(f"Absorbing element: {doc.absorbing_element}")
            # Plot spectrum
            # doc.spectrum contains the data
```

## Similarity Search

### Find Similar Structures

```python
with MPRester() as mpr:
    # Find structures similar to query structure
    query_structure = Structure.from_file("POSCAR")

    similar = mpr.materials.similarity.search(
        structure=query_structure,
        ltol=0.2,
        stol=0.3,
        angle_tol=5
    )

    print(f"Found {len(similar)} similar structures")
    for match in similar[:5]:  # Top 5
        print(f"{match.material_id}: {match.formula_pretty}")
```

## Provenance and Citations

### Get Calculation Details

```python
with MPRester() as mpr:
    # Get detailed calculation info
    doc = mpr.materials.summary.get_data_by_id("mp-149")

    print(f"Database version: {doc.database_version}")
    print(f"Last updated: {doc.last_updated}")
    print(f"Energy: {doc.energy_per_atom:.4f} eV/atom")
```

### Get Citation

```python
with MPRester() as mpr:
    # Get proper citation
    doc = mpr.materials.provenance.get_data_by_id("mp-149")

    if doc:
        print("Please cite:")
        print(doc.references)
```

## Batch Queries

### Query Multiple Materials

```python
with MPRester() as mpr:
    mp_ids = ["mp-149", "mp-66", "mp-1143"]

    structures = {}
    for mp_id in mp_ids:
        structure = mpr.get_structure_by_material_id(mp_id)
        structures[mp_id] = structure

    for mp_id, struct in structures.items():
        print(f"{mp_id}: {struct.composition.reduced_formula}")
```

### Systematic Screening

```python
with MPRester() as mpr:
    # Screen for battery cathodes
    docs = mpr.materials.summary.search(
        elements=["Li"],                   # Must contain Li
        energy_above_hull=(0, 0.02),       # Stable
        band_gap=(0, 0.5),                 # Metallic/small gap
        num_elements=(2, 4),               # Not too complex
        fields=["material_id", "formula_pretty", "energy_above_hull", "band_gap"]
    )

    print(f"Found {len(docs)} candidate cathodes")

    # Filter further
    candidates = []
    for doc in docs:
        # Add your criteria
        if "O" in doc.formula_pretty or "S" in doc.formula_pretty:
            candidates.append(doc)

    print(f"Filtered to {len(candidates)} candidates")
```

## Working with Entries

### ComputedEntry vs ComputedStructureEntry

```python
with MPRester() as mpr:
    # Get entries (includes energy and composition)
    entries = mpr.get_entries_in_chemsys(["Fe", "O"])

    for entry in entries[:5]:
        print(f"{entry.composition.reduced_formula}:")
        print(f"  Energy: {entry.energy_per_atom:.3f} eV/atom")
        print(f"  mp-id: {entry.entry_id}")

        # Some entries have structures
        if hasattr(entry, 'structure'):
            print(f"  Space group: {entry.structure.get_space_group_info()}")
```

### Correction Schemes

```python
from pymatgen.entries.compatibility import MaterialsProject2020Compatibility

with MPRester() as mpr:
    entries = mpr.get_entries_in_chemsys(["Fe", "O"])

    # Apply MP2020 corrections
    compat = MaterialsProject2020Compatibility()
    corrected_entries = compat.process_entries(entries)

    print(f"Original entries: {len(entries)}")
    print(f"After corrections: {len(corrected_entries)}")
```

## Data Downloads

### Download All Structures in System

```python
import os

with MPRester() as mpr:
    # Get all stable materials in system
    docs = mpr.materials.summary.search(
        chemsys="Li-Fe-P-O",
        energy_above_hull=(0, 0.01),
        fields=["material_id", "structure"]
    )

    # Download all structures
    os.makedirs("structures", exist_ok=True)
    for doc in docs:
        filename = f"structures/{doc.material_id}.cif"
        doc.structure.to(filename=filename)

    print(f"Downloaded {len(docs)} structures")
```

## Error Handling and Rate Limits

### Handle API Errors

```python
from mp_api.client import MPRester
from mp_api.client.core import MPRestError

try:
    with MPRester() as mpr:
        structure = mpr.get_structure_by_material_id("mp-invalid")
except MPRestError as e:
    print(f"API Error: {e}")
except Exception as e:
    print(f"Other error: {e}")
```

### Rate Limiting

```python
import time

with MPRester() as mpr:
    mp_ids = ["mp-" + str(i) for i in range(100, 200)]

    structures = []
    for i, mp_id in enumerate(mp_ids):
        try:
            structure = mpr.get_structure_by_material_id(mp_id)
            structures.append(structure)
        except:
            print(f"Failed to get {mp_id}")

        # Be nice to the API
        if i % 10 == 0:
            time.sleep(1)

    print(f"Retrieved {len(structures)} structures")
```

## Caching Results

### Save Results Locally

```python
import json
from monty.json import MontyEncoder, MontyDecoder

with MPRester() as mpr:
    # Query once
    docs = mpr.materials.summary.search(
        chemsys="Li-Fe-O",
        fields=["material_id", "formula_pretty", "energy_per_atom"]
    )

    # Save to file
    data = [{"material_id": doc.material_id,
             "formula": doc.formula_pretty,
             "energy": doc.energy_per_atom}
            for doc in docs]

    with open("mp_data.json", "w") as f:
        json.dump(data, f, cls=MontyEncoder, indent=2)

# Later, load from file (no API call)
with open("mp_data.json", "r") as f:
    cached_data = json.load(f)

print(f"Loaded {len(cached_data)} entries from cache")
```

## Advanced Usage

### Custom Queries

```python
with MPRester() as mpr:
    # More complex filtering
    docs = mpr.materials.summary.search(
        nelements=(2, 3),                  # Binary or ternary
        band_gap=(1.5, 3.0),               # Semiconductor
        energy_above_hull=(0, 0.01),       # Stable
        volume=(None, 50),                 # Small unit cell
        fields=["material_id", "formula_pretty", "band_gap", "volume"]
    )

    print(f"Found {len(docs)} materials matching criteria")
```

### Combining with Local Calculations

```python
from pymatgen.entries.computed_entries import ComputedEntry

with MPRester() as mpr:
    # Get MP entries
    mp_entries = mpr.get_entries_in_chemsys(["Li", "Fe", "O"])

    # Add your own calculated entry
    my_structure = Structure.from_file("POSCAR")
    my_energy = -123.45  # From your calculation

    my_entry = ComputedEntry(
        composition=my_structure.composition,
        energy=my_energy,
        entry_id="my-calc-001"
    )

    # Combine
    all_entries = mp_entries + [my_entry]

    # Build phase diagram with your data
    from pymatgen.analysis.phase_diagram import PhaseDiagram
    pd = PhaseDiagram(all_entries)

    # Check stability of your structure
    decomp, ehull = pd.get_decomp_and_e_above_hull(my_entry)
    print(f"Your structure's energy above hull: {ehull*1000:.1f} meV/atom")
```

## Best Practices

### Use Context Manager

```python
# Good
with MPRester() as mpr:
    structure = mpr.get_structure_by_material_id("mp-149")

# Bad - doesn't clean up resources
mpr = MPRester()
structure = mpr.get_structure_by_material_id("mp-149")
```

### Limit Fields Retrieved

```python
# Good - only get what you need
with MPRester() as mpr:
    docs = mpr.materials.summary.search(
        formula="Fe2O3",
        fields=["material_id", "energy_per_atom"]
    )

# Bad - retrieves everything (slow, wasteful)
with MPRester() as mpr:
    docs = mpr.materials.summary.search(formula="Fe2O3")
```

### Cache Results

```python
# Cache expensive queries
import pickle

cache_file = "mp_cache.pkl"

try:
    with open(cache_file, "rb") as f:
        entries = pickle.load(f)
    print("Loaded from cache")
except FileNotFoundError:
    with MPRester() as mpr:
        entries = mpr.get_entries_in_chemsys(["Li", "Fe", "P", "O"])

    with open(cache_file, "wb") as f:
        pickle.dump(entries, f)
    print("Downloaded and cached")
```

### Error Handling

```python
with MPRester() as mpr:
    mp_ids = ["mp-149", "mp-invalid", "mp-66"]

    for mp_id in mp_ids:
        try:
            structure = mpr.get_structure_by_material_id(mp_id)
            print(f"Got {mp_id}")
        except Exception as e:
            print(f"Failed {mp_id}: {e}")
            continue
```

## Quick Reference

### Common Queries

| Task | Code |
|------|------|
| **Get structure by ID** | `mpr.get_structure_by_material_id("mp-149")` |
| **Search by formula** | `mpr.materials.summary.search(formula="Fe2O3")` |
| **Search by elements** | `mpr.materials.summary.search(elements=["Li","Fe","O"])` |
| **Search by chemsys** | `mpr.materials.summary.search(chemsys="Li-Fe-O")` |
| **Get band structure** | `mpr.get_bandstructure_by_material_id("mp-149")` |
| **Get DOS** | `mpr.get_dos_by_material_id("mp-149")` |
| **Get entries** | `mpr.get_entries_in_chemsys(["Fe", "O"])` |

### Common Fields

- `material_id` - MP identifier
- `formula_pretty` - Formatted formula
- `structure` - Structure object
- `energy_per_atom` - Energy per atom (eV)
- `energy_above_hull` - Stability (eV/atom)
- `band_gap` - Band gap (eV)
- `is_gap_direct` - Direct/indirect gap
- `symmetry` - Space group info
- `volume` - Unit cell volume (ų)
- `density` - Density (g/cm³)

### Search Filters

- `elements` - Must contain these elements
- `formula` - Specific formula
- `chemsys` - Chemical system
- `num_elements` - Number of element types
- `band_gap` - (min, max) tuple
- `energy_above_hull` - (min, max) stability
- `spacegroup_symbol` - Space group
- `nelements` - (min, max) elements
