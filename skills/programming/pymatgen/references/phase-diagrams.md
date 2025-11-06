# Pymatgen Phase Diagrams

Comprehensive guide to constructing and analyzing phase diagrams using pymatgen.

## Overview

Pymatgen supports several types of phase diagrams:
- **PhaseDiagram** - Standard thermodynamic phase diagrams
- **CompoundPhaseDiagram** - For specific compositions
- **GrandPotentialPhaseDiagram** - Open chemical potential
- **PourbaixDiagram** - Electrochemical stability (pH-potential)

## Prerequisites

Phase diagrams require ComputedEntry objects with:
- Chemical composition
- Total energy (usually from DFT)

These typically come from Materials Project or your own calculations.

## Basic Phase Diagrams

### Creating a Phase Diagram

```python
from mp_api.client import MPRester
from pymatgen.analysis.phase_diagram import PhaseDiagram

# Get entries from Materials Project
with MPRester() as mpr:
    entries = mpr.get_entries_in_chemsys(["Li", "Fe", "O"])

# Create phase diagram
pd = PhaseDiagram(entries)

print(f"Phase diagram with {len(pd.stable_entries)} stable phases")
```

### Plotting Phase Diagrams

**Binary system:**
```python
from mp_api.client import MPRester
from pymatgen.analysis.phase_diagram import PhaseDiagram, PDPlotter

with MPRester() as mpr:
    entries = mpr.get_entries_in_chemsys(["Fe", "O"])

pd = PhaseDiagram(entries)

# Plot
plotter = PDPlotter(pd)
plotter.show()  # Interactive plot
plotter.write_image("fe_o_phase_diagram.png")  # Save to file
```

**Ternary system:**
```python
with MPRester() as mpr:
    entries = mpr.get_entries_in_chemsys(["Li", "Fe", "O"])

pd = PhaseDiagram(entries)
plotter = PDPlotter(pd, show_unstable=True)
plotter.show()
```

**Quaternary system:**
```python
with MPRester() as mpr:
    entries = mpr.get_entries_in_chemsys(["Li", "Fe", "P", "O"])

pd = PhaseDiagram(entries)

# For quaternary, need to specify projection
plotter = PDPlotter(pd)
plotter.show()  # Will show 3D projection
```

## Stability Analysis

### Energy Above Hull

**Check if composition is stable:**
```python
from pymatgen.core import Composition

comp = Composition("LiFePO4")

# Find entry for this composition
target_entry = None
for entry in entries:
    if entry.composition.reduced_formula == comp.reduced_formula:
        target_entry = entry
        break

if target_entry:
    decomp, ehull = pd.get_decomp_and_e_above_hull(target_entry)

    print(f"Energy above hull: {ehull*1000:.2f} meV/atom")

    if ehull < 0.001:  # Within 1 meV/atom
        print("This composition is stable!")
    else:
        print("This composition is unstable")
        print(f"Decomposition products: {decomp}")
```

### Decomposition Products

```python
# Get decomposition for unstable compound
entry = entries[5]  # Some entry

decomp, ehull = pd.get_decomp_and_e_above_hull(entry)

print(f"{entry.composition.reduced_formula}:")
print(f"  Energy above hull: {ehull*1000:.2f} meV/atom")
print(f"  Decomposes to:")
for product, amount in decomp.items():
    formula = product.composition.reduced_formula
    print(f"    {amount:.3f} × {formula}")
```

### Get All Stable Phases

```python
# Get all stable entries
stable_entries = pd.stable_entries

print(f"Found {len(stable_entries)} stable phases:")
for entry in stable_entries:
    formula = entry.composition.reduced_formula
    energy = entry.energy_per_atom
    print(f"  {formula}: {energy:.3f} eV/atom")
```

### Metastability

```python
# Find entries within certain energy above hull
threshold = 0.05  # 50 meV/atom

metastable = []
for entry in entries:
    decomp, ehull = pd.get_decomp_and_e_above_hull(entry)
    if 0 < ehull < threshold:
        metastable.append((entry, ehull))

print(f"Found {len(metastable)} metastable phases (<{threshold*1000} meV/atom):")
for entry, ehull in sorted(metastable, key=lambda x: x[1]):
    print(f"  {entry.composition.reduced_formula}: {ehull*1000:.2f} meV/atom")
```

## Formation Energy

### Calculate Formation Energy

```python
# Formation energy from elements
for entry in pd.stable_entries:
    # Get formation energy per atom
    formation_energy = pd.get_form_energy_per_atom(entry)

    formula = entry.composition.reduced_formula
    print(f"{formula}: ΔHf = {formation_energy:.3f} eV/atom")
```

### Reaction Energy

```python
from pymatgen.analysis.reaction_calculator import ComputedReaction

# Define reaction
reactants = ["Li", "FePO4"]
products = ["LiFePO4"]

# Get entries
reactant_entries = []
product_entries = []

for entry in entries:
    formula = entry.composition.reduced_formula
    if formula in reactants:
        reactant_entries.append(entry)
    elif formula in products:
        product_entries.append(entry)

# Calculate reaction energy
# This requires the ComputedReaction class
# reaction = ComputedReaction(reactant_entries, product_entries)
# print(f"Reaction energy: {reaction.calculated_reaction_energy:.3f} eV")
```

## Compound Phase Diagrams

### For Specific Composition

```python
from pymatgen.analysis.phase_diagram import CompoundPhaseDiagram

# Focus on LiFePO4 composition
comp = Composition("LiFePO4")

# Create compound phase diagram
cpd = CompoundPhaseDiagram(entries, comp)

print(f"Number of tie-line end points: {len(cpd.get_tie_line_end_points())}")

# Energy above convex hull
ehull = cpd.get_energy_above_hull()
print(f"Energy above hull: {ehull*1000:.2f} meV/atom")
```

## Grand Potential Phase Diagrams

### Open System (Variable Chemical Potential)

```python
from pymatgen.analysis.phase_diagram import GrandPotentialPhaseDiagram
from pymatgen.core import Element

with MPRester() as mpr:
    entries = mpr.get_entries_in_chemsys(["Li", "Fe", "O"])

# Fix oxygen chemical potential
chempots = {Element("O"): -5.0}  # eV/atom

# Create grand potential phase diagram
gpd = GrandPotentialPhaseDiagram(entries, chempots)

print(f"Stable phases at μ_O = -5.0 eV:")
for entry in gpd.stable_entries:
    print(f"  {entry.composition.reduced_formula}")
```

### Chemical Potential Range

```python
# Scan over chemical potential range
import numpy as np

mu_o_range = np.linspace(-8, -2, 20)

for mu_o in mu_o_range:
    chempots = {Element("O"): mu_o}
    gpd = GrandPotentialPhaseDiagram(entries, chempots)

    n_stable = len(gpd.stable_entries)
    print(f"μ_O = {mu_o:.2f} eV: {n_stable} stable phases")
```

## Pourbaix Diagrams

### Electrochemical Stability

```python
from mp_api.client import MPRester
from pymatgen.analysis.pourbaix_diagram import PourbaixDiagram, PourbaixPlotter

# Get Pourbaix entries
with MPRester() as mpr:
    entries = mpr.get_pourbaix_entries(["Fe"])

# Create Pourbaix diagram
pb = PourbaixDiagram(entries)

# Plot
plotter = PourbaixPlotter(pb)
plotter.show()
plotter.save_plot("fe_pourbaix.png")
```

### Multi-element Pourbaix

```python
with MPRester() as mpr:
    # Iron-chromium system
    entries = mpr.get_pourbaix_entries(["Fe", "Cr"])

pb = PourbaixDiagram(entries)
plotter = PourbaixPlotter(pb)
plotter.show()
```

### Get Stable Species at Conditions

```python
# Get stable species at specific pH and potential
pH = 7.0
V = 0.0  # vs SHE

stable_entry = pb.get_stable_entry(pH, V)
print(f"Stable at pH={pH}, V={V}: {stable_entry.name}")

# Scan conditions
import numpy as np

pH_range = np.linspace(0, 14, 50)
V_range = np.linspace(-2, 2, 50)

stability_map = {}
for pH in pH_range:
    for V in V_range:
        entry = pb.get_stable_entry(pH, V)
        key = (round(pH, 1), round(V, 1))
        stability_map[key] = entry.name
```

## Phase Diagram Analysis

### Get Composition Range

```python
# Find composition range for stable phase
from pymatgen.core import Composition

target_formula = "LiFePO4"

for entry in pd.stable_entries:
    if entry.composition.reduced_formula == target_formula:
        # Get neighboring phases
        neighbors = pd.get_composition_neighbors(entry)

        print(f"Neighbors of {target_formula}:")
        for neighbor in neighbors:
            print(f"  {neighbor.composition.reduced_formula}")
```

### Convex Hull

```python
# Get convex hull data
hull_data = pd.get_hull_data()

# Each point on convex hull
for point in hull_data:
    comp = point['composition']
    energy = point['energy']
    print(f"{comp}: {energy:.3f} eV/atom")
```

### Chemical Potential Limits

```python
from pymatgen.analysis.phase_diagram import PhaseDiagram

# Get chemical potential limits
for element in pd.elements:
    # Maximum chemical potential (pure element phase)
    mu_max = pd.get_equilibrium_reaction_energy(
        pd.el_refs[element]
    )
    print(f"μ_{element} max: {mu_max:.3f} eV/atom")
```

## Custom Entry Data

### Add Your Own Calculations

```python
from pymatgen.entries.computed_entries import ComputedEntry
from pymatgen.core import Composition

# Your calculated data
my_compositions = ["Li2FeO2", "Li3FeO3"]
my_energies = [-45.67, -68.90]  # Total energies from your DFT

my_entries = []
for comp_str, energy in zip(my_compositions, my_energies):
    comp = Composition(comp_str)
    entry = ComputedEntry(comp, energy)
    my_entries.append(entry)

# Combine with MP data
with MPRester() as mpr:
    mp_entries = mpr.get_entries_in_chemsys(["Li", "Fe", "O"])

all_entries = mp_entries + my_entries

# Build phase diagram
pd = PhaseDiagram(all_entries)

# Check if your phases are stable
for entry in my_entries:
    decomp, ehull = pd.get_decomp_and_e_above_hull(entry)
    print(f"{entry.composition.reduced_formula}: {ehull*1000:.2f} meV/atom")
```

### With Corrections

```python
from pymatgen.entries.computed_entries import ComputedEntry

# Add energy corrections (e.g., from experiments)
corrections = {"Li2FeO2": -0.1, "Li3FeO3": -0.05}  # eV/atom

corrected_entries = []
for entry in my_entries:
    formula = entry.composition.reduced_formula
    if formula in corrections:
        correction = corrections[formula] * entry.composition.num_atoms
        corrected_entry = ComputedEntry(
            entry.composition,
            entry.energy + correction,
            entry_id=entry.entry_id
        )
        corrected_entries.append(corrected_entry)
```

## Advanced Analysis

### Interface Reactions

```python
from pymatgen.analysis.interface_reactions import InterfaceReactions

# Analyze reactions at interface
anode_comp = Composition("Li")
cathode_comp = Composition("LiFePO4")

interface = InterfaceReactions(
    c1=anode_comp,
    c2=cathode_comp,
    pd=pd,
    norm=True,
    include_no_mixing_energy=False
)

# Get reaction energies
reactions = interface.get_kinks()
for rxn in reactions:
    print(f"Reaction: {rxn}")
```

### Convexity Analysis

```python
# Analyze convexity for specific composition range
from pymatgen.core import Composition

# Define composition line
comp1 = Composition("Li2O")
comp2 = Composition("Fe2O3")

# Sample compositions along line
import numpy as np
x_range = np.linspace(0, 1, 20)

energies = []
for x in x_range:
    comp = comp1 * x + comp2 * (1 - x)

    # Find closest entry
    min_ehull = float('inf')
    for entry in entries:
        if entry.composition.reduced_formula == comp.reduced_formula:
            decomp, ehull = pd.get_decomp_and_e_above_hull(entry)
            if ehull < min_ehull:
                min_ehull = ehull

    energies.append(min_ehull)

# Plot convexity
import matplotlib.pyplot as plt
plt.plot(x_range, energies)
plt.xlabel("x in xLi2O + (1-x)Fe2O3")
plt.ylabel("Energy above hull (eV/atom)")
plt.show()
```

## Practical Examples

### Battery Cathode Screening

```python
from mp_api.client import MPRester
from pymatgen.analysis.phase_diagram import PhaseDiagram

# Get lithium-containing oxides
with MPRester() as mpr:
    docs = mpr.materials.summary.search(
        elements=["Li", "O"],
        energy_above_hull=(0, 0.02),  # Stable or nearly stable
        fields=["material_id", "formula_pretty", "energy_per_atom"]
    )

# Build phase diagram for each system
candidates = []
for doc in docs:
    formula = doc.formula_pretty

    # Parse elements
    comp = Composition(formula)
    elements = [str(el) for el in comp.elements]

    # Get entries
    entries = mpr.get_entries_in_chemsys(elements)
    pd = PhaseDiagram(entries)

    # Find entry
    for entry in entries:
        if entry.composition.reduced_formula == comp.reduced_formula:
            decomp, ehull = pd.get_decomp_and_e_above_hull(entry)

            if ehull < 0.01:  # Very stable
                candidates.append({
                    'formula': formula,
                    'ehull': ehull,
                    'energy': entry.energy_per_atom
                })
                break

print(f"Found {len(candidates)} stable cathode candidates")
```

### Corrosion Resistance

```python
# Check stability in aqueous environment
from mp_api.client import MPRester
from pymatgen.analysis.pourbaix_diagram import PourbaixDiagram

element = "Fe"

with MPRester() as mpr:
    entries = mpr.get_pourbaix_entries([element])

pb = PourbaixDiagram(entries)

# Test conditions (seawater)
pH = 8.0
V = 0.0

stable_species = pb.get_stable_entry(pH, V)
print(f"In seawater (pH {pH}): {stable_species.name}")

if "oxide" in stable_species.name.lower() or "hydroxide" in stable_species.name.lower():
    print("Forms passive layer")
else:
    print("May corrode")
```

### Solid Solution Range

```python
# Estimate solid solution range
from pymatgen.core import Composition

base_formula = "LiFePO4"
substituent = "Mn"  # Replace Fe with Mn

# Get entries
with MPRester() as mpr:
    entries = mpr.get_entries_in_chemsys(["Li", "Fe", "Mn", "P", "O"])

pd = PhaseDiagram(entries)

# Test different compositions
x_range = np.linspace(0, 1, 11)
stable_range = []

for x in x_range:
    # Li(Fe_{1-x}Mn_x)PO4
    comp = Composition({
        "Li": 1,
        "Fe": 1-x,
        "Mn": x,
        "P": 1,
        "O": 4
    })

    # Check if any entry matches this composition
    for entry in entries:
        if entry.composition.almost_equals(comp, rtol=0.01):
            decomp, ehull = pd.get_decomp_and_e_above_hull(entry)

            if ehull < 0.02:  # Stable
                stable_range.append(x)
            break

print(f"Stable solid solution range: x = {min(stable_range):.2f} to {max(stable_range):.2f}")
```

## Plotting Options

### Customizing Plots

```python
from pymatgen.analysis.phase_diagram import PDPlotter

plotter = PDPlotter(
    pd,
    show_unstable=True,      # Show unstable phases
    backend="matplotlib"     # Or "plotly" for interactive
)

# Customize appearance
plotter.show(
    label_stable=True,       # Label stable phases
    label_unstable=False     # Don't label unstable
)
```

### Save High-Quality Images

```python
plotter = PDPlotter(pd)

# Save as vector format
plotter.write_image("phase_diagram.pdf", format="pdf")

# Or raster with high DPI
plotter.write_image("phase_diagram.png", format="png", width=1200, height=800)
```

## Tips and Best Practices

### Energy Corrections

```python
from pymatgen.entries.compatibility import MaterialsProject2020Compatibility

# Apply MP corrections before building phase diagram
compat = MaterialsProject2020Compatibility()
corrected_entries = compat.process_entries(entries)

pd = PhaseDiagram(corrected_entries)
```

### Handling Large Systems

```python
# For large systems (many elements), limit to relevant region
with MPRester() as mpr:
    # Get all entries
    all_entries = mpr.get_entries_in_chemsys(["Li", "Fe", "P", "O"])

    # Filter to stable and near-stable only
    filtered_entries = []
    for entry in all_entries:
        decomp, ehull = PhaseDiagram(all_entries).get_decomp_and_e_above_hull(entry)
        if ehull < 0.1:  # Within 100 meV/atom
            filtered_entries.append(entry)

    # Build phase diagram with filtered entries (faster)
    pd = PhaseDiagram(filtered_entries)
```

### Validation

```python
# Validate phase diagram construction
print(f"Number of entries: {len(entries)}")
print(f"Number of stable phases: {len(pd.stable_entries)}")
print(f"Elements: {[str(el) for el in pd.elements]}")

# Check if elemental references exist
for el in pd.elements:
    if el in pd.el_refs:
        ref = pd.el_refs[el]
        print(f"{el} reference: {ref.composition.reduced_formula}")
    else:
        print(f"Warning: No reference for {el}")
```

## Quick Reference

### Creating Phase Diagrams

```python
# Basic
from pymatgen.analysis.phase_diagram import PhaseDiagram
pd = PhaseDiagram(entries)

# Compound
from pymatgen.analysis.phase_diagram import CompoundPhaseDiagram
cpd = CompoundPhaseDiagram(entries, Composition("LiFePO4"))

# Grand potential
from pymatgen.analysis.phase_diagram import GrandPotentialPhaseDiagram
gpd = GrandPotentialPhaseDiagram(entries, {Element("O"): -5.0})

# Pourbaix
from pymatgen.analysis.pourbaix_diagram import PourbaixDiagram
pb = PourbaixDiagram(pourbaix_entries)
```

### Common Operations

| Operation | Code |
|-----------|------|
| **Stability** | `pd.get_decomp_and_e_above_hull(entry)` |
| **Formation energy** | `pd.get_form_energy_per_atom(entry)` |
| **Stable entries** | `pd.stable_entries` |
| **Plot** | `PDPlotter(pd).show()` |
| **Pourbaix stable** | `pb.get_stable_entry(pH, V)` |

### Stability Thresholds

- **Stable**: E_hull < 1 meV/atom (0.001 eV/atom)
- **Metastable**: 1-50 meV/atom
- **Synthesizable**: < 100 meV/atom (sometimes)
- **Unstable**: > 100 meV/atom
