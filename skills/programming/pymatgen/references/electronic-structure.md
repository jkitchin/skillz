# Pymatgen Electronic Structure Analysis

Comprehensive guide to analyzing and plotting band structures, density of states, and related electronic properties.

## Overview

Pymatgen provides tools for analyzing electronic structure from DFT calculations:
- **BandStructure** - Band structure analysis and plotting
- **DOS** (Density of States) - Total and projected DOS
- **BandGap** - Band gap analysis (direct/indirect)
- **Orbital projections** - Analyze orbital contributions

## Band Structures

### From VASP Calculations

**Loading band structure:**
```python
from pymatgen.io.vasp.outputs import Vasprun, BSVasprun
from pymatgen.electronic_structure.bandstructure import BandStructureSymmLine

# Parse vasprun.xml from band structure calculation
vasprun = BSVasprun("vasprun.xml")
bs = vasprun.get_band_structure(line_mode=True)

print(f"Band structure type: {type(bs)}")
print(f"Is metal: {bs.is_metal()}")
```

**From Materials Project:**
```python
from mp_api.client import MPRester

with MPRester() as mpr:
    bs = mpr.get_bandstructure_by_material_id("mp-149")  # Silicon

    if bs:
        print(f"Number of bands: {bs.nb_bands}")
        print(f"Number of k-points: {len(bs.kpoints)}")
    else:
        print("Band structure not available")
```

### Band Gap Analysis

**Basic band gap:**
```python
# Get band gap information
bg_info = bs.get_band_gap()

print(f"Band gap: {bg_info['energy']:.3f} eV")
print(f"Direct: {bg_info['direct']}")
print(f"Transition: {bg_info['transition']}")
```

**Direct vs Indirect:**
```python
if bg_info['direct']:
    print(f"Direct gap at {bg_info['transition']}")
else:
    print(f"Indirect gap")
    print(f"  VBM at {bs.get_vbm()['kpoint']}")
    print(f"  CBM at {bs.get_cbm()['kpoint']}")
```

**Metal or semiconductor:**
```python
if bs.is_metal():
    print("Material is metallic")
else:
    gap = bs.get_band_gap()
    if gap['energy'] == 0:
        print("Zero-gap semiconductor")
    elif gap['energy'] < 3.0:
        print(f"Semiconductor with gap {gap['energy']:.2f} eV")
    else:
        print(f"Insulator with gap {gap['energy']:.2f} eV")
```

### Band Edges

**Valence Band Maximum (VBM):**
```python
vbm_info = bs.get_vbm()

print(f"VBM energy: {vbm_info['energy']:.3f} eV")
print(f"VBM k-point: {vbm_info['kpoint']}")
print(f"VBM band index: {vbm_info['band_index']}")
print(f"VBM k-point label: {vbm_info['kpoint'].label}")
```

**Conduction Band Minimum (CBM):**
```python
cbm_info = bs.get_cbm()

print(f"CBM energy: {cbm_info['energy']:.3f} eV")
print(f"CBM k-point: {cbm_info['kpoint']}")
print(f"CBM band index: {cbm_info['band_index']}")
```

### Plotting Band Structures

**Basic plot:**
```python
from pymatgen.electronic_structure.plotter import BSPlotter

plotter = BSPlotter(bs)

# Interactive plot
plotter.show()

# Save to file
plotter.save_plot("band_structure.png", img_format="png")
```

**Customized plot:**
```python
plotter = BSPlotter(bs)

# Get matplotlib figure
plt = plotter.get_plot(
    zero_to_efermi=True,        # Set Fermi level to zero
    ylim=(-10, 10),              # Energy window
    vbm_cbm_marker=True          # Mark VBM/CBM
)

plt.savefig("band_structure_custom.png", dpi=300)
```

**Multiple band structures:**
```python
from pymatgen.electronic_structure.plotter import BSPlotter

# Compare band structures (e.g., different functionals)
bs_dict = {
    "PBE": bs_pbe,
    "HSE": bs_hse
}

plotter = BSPlotter(bs_pbe)  # Reference

# Plot both
for label, bs in bs_dict.items():
    plotter_temp = BSPlotter(bs)
    # Customize and overlay plots
```

### Spin-Polarized Band Structures

```python
# For spin-polarized calculations
if bs.is_spin_polarized:
    print("Spin-polarized calculation")

    # Get spin-up and spin-down band gaps
    for spin in [Spin.up, Spin.down]:
        bands = bs.bands[spin]
        print(f"{spin}: {len(bands)} bands")

# Plot both spins
plotter = BSPlotter(bs)
plotter.show()  # Automatically handles spin
```

### Projected Band Structure

**With orbital projections:**
```python
from pymatgen.electronic_structure.plotter import BSPlotterProjected

# Requires band structure with projections
# (from vasprun.xml with LORBIT = 11 in VASP)
plotter = BSPlotterProjected(bs)

# Plot with element projections
plotter.get_elt_projected_plots()

# Plot with orbital projections
plotter.get_projected_plots_dots(
    {
        "Fe": ["s", "p", "d"],
        "O": ["s", "p"]
    }
)
```

## Density of States (DOS)

### Loading DOS

**From VASP:**
```python
from pymatgen.io.vasp.outputs import Vasprun

vasprun = Vasprun("vasprun.xml")

# Get complete DOS (with projections)
dos = vasprun.complete_dos

print(f"Fermi level: {dos.efermi:.3f} eV")
print(f"DOS at Fermi level: {dos.get_densities(dos.efermi)}")
```

**From Materials Project:**
```python
from mp_api.client import MPRester

with MPRester() as mpr:
    dos = mpr.get_dos_by_material_id("mp-149")

    if dos:
        print(f"Fermi level: {dos.efermi:.3f} eV")
    else:
        print("DOS not available")
```

### DOS Analysis

**Integrated DOS:**
```python
# Get total number of electrons up to Fermi level
n_electrons = dos.get_integrated_dos(dos.efermi)
print(f"Number of electrons: {n_electrons:.1f}")
```

**DOS at specific energy:**
```python
# DOS at specific energy
energy = 0.0  # eV relative to Fermi
dos_at_e = dos.get_densities(energy)
print(f"DOS at {energy} eV: {dos_at_e}")
```

**Gap from DOS:**
```python
# Get band gap from DOS
gap = dos.get_gap()
print(f"Band gap from DOS: {gap:.3f} eV")

# Check if metallic
if dos.get_gap() == 0:
    print("Metal (no gap in DOS)")
```

### Projected DOS

**Element-projected DOS:**
```python
from pymatgen.electronic_structure.core import Spin

# Get DOS projected onto specific element
element_dos = dos.get_element_dos()

for element, pdos in element_dos.items():
    print(f"\n{element}:")
    # Integrated DOS for this element
    integrated = pdos.get_integrated_dos(dos.efermi)
    print(f"  Electrons: {integrated:.2f}")
```

**Orbital-projected DOS:**
```python
from pymatgen.core import Element
from pymatgen.electronic_structure.core import OrbitalType

# Get DOS for specific element and orbital
el = Element("Fe")
orbital = OrbitalType.d

# Get spd-projected DOS
spd_dos = dos.get_spd_dos()

for orbital_type, pdos in spd_dos.items():
    print(f"\n{orbital_type}:")
    integrated = pdos.get_integrated_dos(dos.efermi)
    print(f"  Electrons: {integrated:.2f}")
```

**Site-projected DOS:**
```python
# DOS projected onto specific atomic sites
structure = vasprun.final_structure

# Get DOS for first Fe atom
fe_sites = [i for i, site in enumerate(structure) if site.species_string == "Fe"]

if fe_sites:
    site_dos = dos.get_site_dos(structure[fe_sites[0]])
    print(f"DOS at Fe site: {site_dos.get_integrated_dos(dos.efermi):.2f} electrons")
```

### Plotting DOS

**Basic plot:**
```python
from pymatgen.electronic_structure.plotter import DosPlotter

plotter = DosPlotter()
plotter.add_dos("Total", dos)

# Show
plotter.show()

# Save
plotter.save_plot("dos.png", img_format="png")
```

**Element-projected DOS:**
```python
plotter = DosPlotter()

# Add total DOS
plotter.add_dos("Total", dos.get_dos())

# Add element projections
element_dos = dos.get_element_dos()
for element, pdos in element_dos.items():
    plotter.add_dos(str(element), pdos)

plotter.show()
```

**Orbital-projected DOS:**
```python
plotter = DosPlotter()

# Add spd-projected DOS
spd_dos = dos.get_spd_dos()
for orbital, pdos in spd_dos.items():
    plotter.add_dos(str(orbital), pdos)

plotter.show()
```

**Customized DOS plot:**
```python
plotter = DosPlotter()
plotter.add_dos("Total", dos)

# Get matplotlib figure
plt_obj = plotter.get_plot(
    xlim=(-10, 10),         # Energy range
    ylim=None                # Auto y-range
)

plt_obj.savefig("dos_custom.png", dpi=300, bbox_inches='tight')
```

### Spin-Polarized DOS

```python
from pymatgen.electronic_structure.plotter import DosPlotter

# For magnetic materials
if dos.is_spin_polarized:
    plotter = DosPlotter()

    # Add spin-up and spin-down
    plotter.add_dos("Spin up", dos.get_dos(spin=Spin.up))
    plotter.add_dos("Spin down", dos.get_dos(spin=Spin.down))

    plotter.show()
```

## Combined Band Structure and DOS

### Side-by-side plots:**
```python
from pymatgen.electronic_structure.plotter import BSPlotter, DosPlotter
import matplotlib.pyplot as plt

# Create figure with two subplots
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 6))

# Plot band structure
bs_plotter = BSPlotter(bs)
bs_plot = bs_plotter.get_plot(zero_to_efermi=True)
ax1 = bs_plot.gca()

# Plot DOS
dos_plotter = DosPlotter()
dos_plotter.add_dos("Total", dos)
dos_plot = dos_plotter.get_plot()
ax2 = dos_plot.gca()

plt.tight_layout()
plt.savefig("bs_dos_combined.png", dpi=300)
plt.show()
```

**Aligned plots:**
```python
from pymatgen.electronic_structure.plotter import BSPlotter, DosPlotter

# Create aligned band structure and DOS plot
bs_plotter = BSPlotter(bs)
dos_plotter = DosPlotter()
dos_plotter.add_dos("Total", dos)

# This creates a combined figure
# (Implementation depends on version)
```

## Effective Mass

### Calculate Effective Mass

```python
from pymatgen.electronic_structure.bandstructure import BandStructure

# Get effective mass near band edges
# This requires band structure with fine k-point sampling

# Electron effective mass (near CBM)
# Manual calculation by fitting parabola to CBM

cbm_info = bs.get_cbm()
cbm_kpoint = cbm_info['kpoint']
cbm_band_idx = cbm_info['band_index'][Spin.up][0]

# Get energies near CBM
# ... fitting code here ...
```

## Dielectric Properties

### From VASP

```python
from pymatgen.io.vasp.outputs import Vasprun

vasprun = Vasprun("vasprun.xml")

# Electronic dielectric constant (high-frequency)
if hasattr(vasprun, 'epsilon_ionic'):
    epsilon_electronic = vasprun.dielectric[0]
    epsilon_ionic = vasprun.dielectric[1]
    epsilon_total = vasprun.dielectric[2]

    print("Dielectric tensor (electronic):")
    print(epsilon_electronic)
    print(f"\nAverage: {sum(epsilon_electronic.diagonal())/3:.2f}")
```

## Fermi Surface

### Analyzing Fermi Surface

```python
# Requires dense k-point sampling
# Get bands at Fermi level

fermi_bands = []
for band_idx, band_energies in enumerate(bs.bands[Spin.up]):
    # Check if band crosses Fermi level
    if min(band_energies) < dos.efermi < max(band_energies):
        fermi_bands.append(band_idx)

print(f"Bands crossing Fermi level: {fermi_bands}")
```

## Work Function

### Calculate Work Function

```python
from pymatgen.io.vasp.outputs import Locpot

# Parse LOCPOT file
locpot = Locpot.from_file("LOCPOT")

# Get average potential
avg_potential = locpot.get_average_along_axis(2)  # Along z-axis

# Work function = vacuum_level - fermi_level
vacuum_level = max(avg_potential)
fermi_level = dos.efermi

work_function = vacuum_level - fermi_level
print(f"Work function: {work_function:.2f} eV")
```

## Plotting Tips

### Publication-Quality Plots

```python
from pymatgen.electronic_structure.plotter import BSPlotter
import matplotlib.pyplot as plt

plotter = BSPlotter(bs)

# Get plot object
plt_obj = plotter.get_plot(
    zero_to_efermi=True,
    ylim=(-5, 5),
    vbm_cbm_marker=True
)

# Customize for publication
ax = plt_obj.gca()
ax.set_ylabel("Energy (eV)", fontsize=14)
ax.set_xlabel("", fontsize=14)
ax.tick_params(labelsize=12)

# Save high resolution
plt_obj.savefig(
    "band_structure_pub.pdf",
    dpi=300,
    bbox_inches='tight',
    format='pdf'
)
```

### Custom Colors and Styles

```python
import matplotlib.pyplot as plt

# Set style
plt.style.use('seaborn-v0_8-paper')  # or 'ggplot', 'bmh', etc.

# Custom colors
colors = ['#1f77b4', '#ff7f0e', '#2ca02c']

plotter = DosPlotter()
for i, (label, pdos) in enumerate(element_dos.items()):
    plotter.add_dos(str(label), pdos)

# Get plot and customize
plt_obj = plotter.get_plot()
ax = plt_obj.gca()

# Change colors
lines = ax.get_lines()
for i, line in enumerate(lines):
    line.set_color(colors[i % len(colors)])
    line.set_linewidth(2)

plt_obj.savefig("dos_custom_colors.png", dpi=300)
```

## Advanced Analysis

### Orbital Character Analysis

```python
# Analyze orbital character at specific k-points
from pymatgen.electronic_structure.core import Spin, OrbitalType

# Get projections at VBM
vbm_kpoint_index = vbm_info['kpoint_index']
vbm_band_index = vbm_info['band_index'][Spin.up][0]

# Access projections
if hasattr(bs, 'projections'):
    projections = bs.projections[Spin.up][vbm_band_index][vbm_kpoint_index]

    # projections is array: [site_index, orbital_index]
    # Analyze which orbitals contribute most
```

### Band Structure Comparison

```python
# Compare band structures (e.g., different structures or functionals)
bs1 = vasprun1.get_band_structure()
bs2 = vasprun2.get_band_structure()

# Compare gaps
gap1 = bs1.get_band_gap()
gap2 = bs2.get_band_gap()

print(f"Structure 1 gap: {gap1['energy']:.3f} eV")
print(f"Structure 2 gap: {gap2['energy']:.3f} eV")
print(f"Difference: {abs(gap1['energy'] - gap2['energy']):.3f} eV")
```

### DOS Comparison

```python
from pymatgen.electronic_structure.plotter import DosPlotter

plotter = DosPlotter()

# Add multiple DOS for comparison
plotter.add_dos("PBE", dos_pbe)
plotter.add_dos("HSE", dos_hse)

plotter.show()
```

## Creating Band Structures from Scratch

### Manual BandStructure Construction

```python
from pymatgen.electronic_structure.bandstructure import BandStructureSymmLine
from pymatgen.electronic_structure.core import Spin
import numpy as np

# Define k-points
kpoints = [...]  # List of Kpoint objects

# Define eigenvalues
# eigenvalues[spin][band_index][kpoint_index]
eigenvalues = {
    Spin.up: np.array([...])  # [n_bands, n_kpoints]
}

# Create BandStructure
bs = BandStructureSymmLine(
    kpoints=kpoints,
    eigenvalues=eigenvalues,
    lattice=structure.lattice,
    efermi=0.0
)
```

## Common Tasks

### Extract Band Gap for Many Materials

```python
from mp_api.client import MPRester

with MPRester() as mpr:
    # Search semiconductors
    docs = mpr.materials.summary.search(
        elements=["Ga", "N"],
        band_gap=(0.1, None),
        fields=["material_id", "formula_pretty", "band_gap", "is_gap_direct"]
    )

    results = []
    for doc in docs:
        gap_type = "direct" if doc.is_gap_direct else "indirect"
        results.append({
            'formula': doc.formula_pretty,
            'gap': doc.band_gap,
            'type': gap_type
        })

    # Sort by gap
    results.sort(key=lambda x: x['gap'])

    for r in results:
        print(f"{r['formula']}: {r['gap']:.2f} eV ({r['type']})")
```

### Plot DOS for Multiple Elements

```python
from pymatgen.electronic_structure.plotter import DosPlotter

plotter = DosPlotter()

# Get element-projected DOS
element_dos = dos.get_element_dos()

# Add each element
for element in ["Fe", "O"]:  # Specific elements
    if element in element_dos:
        plotter.add_dos(element, element_dos[Element(element)])

plotter.save_plot("element_dos.png")
```

## Troubleshooting

### Issue: Band structure not available

**Problem:** Materials Project returns None for band structure

**Solution:** Not all materials have calculated band structures. Check with:
```python
with MPRester() as mpr:
    docs = mpr.materials.summary.search(
        formula="Si",
        fields=["material_id", "has_bandstructure"]
    )
    for doc in docs:
        if hasattr(doc, 'has_bandstructure') and doc.has_bandstructure:
            bs = mpr.get_bandstructure_by_material_id(doc.material_id)
```

### Issue: DOS projections missing

**Problem:** complete_dos doesn't have projections

**Solution:** Requires VASP calculation with LORBIT = 11 or similar:
```python
vasprun = Vasprun("vasprun.xml")
if hasattr(vasprun.complete_dos, 'structure'):
    print("Has projections")
else:
    print("No projections - rerun VASP with LORBIT = 11")
```

### Issue: Plotting errors

**Problem:** Plots don't display or crash

**Solution:**
```python
import matplotlib
matplotlib.use('Agg')  # Use non-interactive backend
import matplotlib.pyplot as plt

# Then plot as usual
plotter.save_plot("output.png")  # Save instead of show()
```

## Quick Reference

### Band Structure

| Task | Code |
|------|------|
| **Get band gap** | `bs.get_band_gap()['energy']` |
| **Is metal** | `bs.is_metal()` |
| **Get VBM** | `bs.get_vbm()` |
| **Get CBM** | `bs.get_cbm()` |
| **Plot** | `BSPlotter(bs).show()` |

### DOS

| Task | Code |
|------|------|
| **Fermi level** | `dos.efermi` |
| **Get gap** | `dos.get_gap()` |
| **Integrated DOS** | `dos.get_integrated_dos(energy)` |
| **Element DOS** | `dos.get_element_dos()` |
| **Orbital DOS** | `dos.get_spd_dos()` |
| **Plot** | `DosPlotter().add_dos("Total", dos).show()` |

### Data Sources

- **VASP**: `Vasprun("vasprun.xml").get_band_structure()`
- **Materials Project**: `mpr.get_bandstructure_by_material_id("mp-id")`
- **VASP DOS**: `Vasprun("vasprun.xml").complete_dos`
- **MP DOS**: `mpr.get_dos_by_material_id("mp-id")`
