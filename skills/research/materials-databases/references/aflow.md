# AFLOW REST API Reference

Complete guide to accessing the AFLOW (Automatic Flow for Materials Discovery) database via REST API.

## Overview

AFLOW is a distributed materials database containing:
- **3.5+ million calculated materials**
- Crystal structures from ICSD and theoretical predictions
- Thermodynamic properties (formation enthalpies, energies)
- Elastic properties (bulk modulus, shear modulus, elastic tensors)
- Electronic properties (band gaps, DOS)
- No API key required - publicly accessible

## Access Methods

### 1. REST API (Direct HTTP)
- Simple URL-based queries
- Returns plain text, JSON, or files
- Use Python `requests` library

### 2. Python Package (aflow)
- Community-maintained package
- Higher-level interface
- Returns Python objects and ASE structures

## Installation

### REST API (Basic)
```bash
# Only need requests library
pip install requests
```

### Python Package (Optional)
```bash
# Install aflow package for easier access
pip install aflow

# With ASE for structure manipulation
pip install aflow ase
```

## REST API Usage

### URL Structure

```
Base URL: http://aflowlib.duke.edu

# Direct file access:
http://aflowlib.duke.edu/AFLOWDATA/<catalog>/<system>/<file>?<directives>

# AFLUX search:
http://aflowlib.duke.edu/search/API/?<query>
```

### Catalogs

- `ICSD_WEB` - Materials from ICSD (experimental structures)
- `LIB1` - AFLOW prototype library
- `LIB2`, `LIB3` - Additional libraries

### Basic Property Access

```python
import requests

# Get a specific property
url = "http://aflowlib.duke.edu/AFLOWDATA/ICSD_WEB/FCC/Ag1/?enthalpy_formation_atom"
response = requests.get(url)
value = float(response.text.strip())
print(f"Formation enthalpy: {value} eV/atom")
```

### Common Directives (Properties)

**Thermodynamic Properties**
- `?enthalpy_formation_atom` - Formation enthalpy (eV/atom)
- `?enthalpy_atom` - Total enthalpy per atom (eV/atom)
- `?energy_atom` - Total energy per atom (eV/atom)
- `?energy_cutoff` - Energy cutoff used (eV)

**Electronic Properties**
- `?Egap` - Band gap (eV)
- `?Egap_fit` - Fitted band gap (eV)
- `?energy_fermi` - Fermi energy (eV)

**Structural Properties**
- `?volume_cell` - Unit cell volume (ų)
- `?volume_atom` - Volume per atom (ų/atom)
- `?density` - Density (g/cm³)
- `?stoichiometry` - Stoichiometry (e.g., [1,1] for AB)
- `?nspecies` - Number of species
- `?natoms` - Number of atoms in unit cell
- `?composition` - Chemical composition
- `?species` - List of elements
- `?species_pp` - Pseudopotentials used

**Symmetry Properties**
- `?spacegroup_relax` - Space group after relaxation
- `?spacegroup_orig` - Original space group
- `?Bravais_lattice_relax` - Bravais lattice type
- `?Pearson_symbol_relax` - Pearson symbol
- `?crystal_family` - Crystal family
- `?crystal_system` - Crystal system

**Elastic Properties**
- `?bulk_modulus_voigt` - Bulk modulus (Voigt, GPa)
- `?bulk_modulus_reuss` - Bulk modulus (Reuss, GPa)
- `?shear_modulus_voigt` - Shear modulus (Voigt, GPa)
- `?shear_modulus_reuss` - Shear modulus (Reuss, GPa)
- `?elastic_anisotropy` - Elastic anisotropy

**File Access**
- `?geometry` - Crystal structure (POSCAR format)
- `?files` - List all available files for this entry
- `?keywords` - List all available properties
- `?aflowlib_entries` - Complete entry data (JSON)

### List Available Properties

```python
# Get all available keywords for an entry
url = "http://aflowlib.duke.edu/AFLOWDATA/ICSD_WEB/FCC/Ag1/?keywords"
response = requests.get(url)
keywords = response.text.split(',')
print("Available properties:", keywords)
```

### Get Structure Files

```python
# Get POSCAR file
url = "http://aflowlib.duke.edu/AFLOWDATA/ICSD_WEB/FCC/Ag1/?geometry"
poscar_text = requests.get(url).text

# Save to file
with open("structure.vasp", "w") as f:
    f.write(poscar_text)

# Parse with pymatgen
from pymatgen.core import Structure
structure = Structure.from_str(poscar_text, fmt="poscar")
```

### Get Complete Entry

```python
# Get all data as JSON
url = "http://aflowlib.duke.edu/AFLOWDATA/ICSD_WEB/FCC/Ag1/?aflowlib_entries"
response = requests.get(url)
data = response.json()

# Access properties
print(data['enthalpy_formation_atom'])
print(data['Egap'])
print(data['spacegroup_relax'])
```

## AFLUX Query Language

AFLUX is AFLOW's search syntax for filtering the database.

### Basic AFLUX Syntax

```
http://aflowlib.duke.edu/search/API/?<keyword1>,<keyword2>,...
```

### Search by Composition

```python
# Specific elements (any composition)
url = "http://aflowlib.duke.edu/search/API/?species(Au)"
results = requests.get(url).json()

# Binary system
url = "http://aflowlib.duke.edu/search/API/?species(Ti,O),nspecies(2)"
results = requests.get(url).json()

# Ternary system
url = "http://aflowlib.duke.edu/search/API/?species(Li,Fe,O),nspecies(3)"
results = requests.get(url).json()

# Exclude elements
url = "http://aflowlib.duke.edu/search/API/?species(Fe,O),nspecies(2),!species(H)"
results = requests.get(url).json()
```

### Search by Properties

**Range Queries** use `*` wildcard:
- `(min*,max*)` - Range from min to max
- `(value*)` - Exact value
- `(*,max*)` - Up to max
- `(min*,*)` - From min onwards

```python
# Band gap between 1-3 eV
url = "http://aflowlib.duke.edu/search/API/?Egap(1*,3*)"

# Formation enthalpy < 0 (stable)
url = "http://aflowlib.duke.edu/search/API/?enthalpy_formation_atom(*,0*)"

# Bulk modulus > 200 GPa
url = "http://aflowlib.duke.edu/search/API/?bulk_modulus_voigt(200*,*)"

# Density between 5-10 g/cm³
url = "http://aflowlib.duke.edu/search/API/?density(5*,10*)"
```

### Combining Criteria

```python
# Gold compounds with band gap > 5 eV from ICSD
url = "http://aflowlib.duke.edu/search/API/?species(Au),Egap(5*,*),catalog(ICSD)"
results = requests.get(url).json()

# Stable ternary oxides
url = ("http://aflowlib.duke.edu/search/API/?"
       "species(O),nspecies(3),enthalpy_formation_atom(*,0*)")
results = requests.get(url).json()

# Binary Ti-O with small unit cell
url = ("http://aflowlib.duke.edu/search/API/?"
       "species(Ti,O),nspecies(2),natoms(*,10*)")
results = requests.get(url).json()
```

### Catalog Selection

```python
# Search only ICSD
url = "http://aflowlib.duke.edu/search/API/?species(Fe),catalog(ICSD)"

# Search specific library
url = "http://aflowlib.duke.edu/search/API/?species(Fe),catalog(LIB1)"
```

### Additional AFLUX Keywords

```python
# Space group
url = "http://aflowlib.duke.edu/search/API/?spacegroup_relax(225*)"

# Bravais lattice
url = "http://aflowlib.duke.edu/search/API/?Bravais_lattice_relax(FCC)"

# Number of atoms
url = "http://aflowlib.duke.edu/search/API/?natoms(1*,20*)"

# Stoichiometry (e.g., AB2 structures)
url = "http://aflowlib.duke.edu/search/API/?stoichiometry([1,2])"
```

## Python aflow Package

Higher-level interface to AFLOW database.

### Installation

```bash
pip install aflow
```

### Basic Usage

```python
import aflow

# Search for materials
results = aflow.search(filter='species(Au),Egap(5*)')

# Iterate over results
for entry in results:
    print(f"Formula: {entry.compound}")
    print(f"Band gap: {entry.Egap} eV")
    print(f"Formation enthalpy: {entry.enthalpy_formation_atom} eV/atom")
```

### Accessing Properties

```python
# Get specific entry
result = aflow.search(filter='species(Ag),nspecies(1)').first()

# Access properties directly
print(result.enthalpy_formation_atom)
print(result.Egap)
print(result.bulk_modulus_voigt)
print(result.spacegroup_relax)
print(result.compound)
print(result.species)
```

### Get Structures

```python
# Get ASE Atoms object
result = aflow.search(filter='species(Si),nspecies(1)').first()
atoms = result.atoms

# Use ASE methods
print(atoms.get_chemical_formula())
print(atoms.get_positions())
print(atoms.get_cell())

# Write structure
from ase.io import write
write('structure.cif', atoms)
write('POSCAR', atoms, format='vasp')
```

### Load Specific Entry

```python
# Load by AFLOW URL or identifier
entry = aflow.load('icsd:Ag:FCC')

# Or by full URL
entry = aflow.load('http://aflowlib.duke.edu/AFLOWDATA/ICSD_WEB/FCC/Ag1')
```

### Search with Multiple Filters

```python
# Build complex queries
results = aflow.search(
    filter='species(Fe,O),nspecies(2),Egap(1*,3*),catalog(ICSD)',
    batch_size=100  # Fetch in batches
)

for entry in results:
    print(f"{entry.compound}: Egap = {entry.Egap} eV")
```

## Working with Results

### REST API Results

```python
import requests

# Search returns JSON array
url = "http://aflowlib.duke.edu/search/API/?species(GaN),Egap"
response = requests.get(url)
results = response.json()

# Each result is a dictionary
for result in results:
    aurl = result['aurl']  # AFLOW URL
    compound = result.get('compound', 'N/A')
    egap = result.get('Egap', None)

    print(f"{compound}: Egap = {egap} eV")
    print(f"URL: {aurl}")
```

### Parse and Save Structures

```python
# Get structure from specific entry
aurl = "http://aflowlib.duke.edu/AFLOWDATA/ICSD_WEB/AB/GaN/"
structure_url = f"{aurl}?geometry"

poscar = requests.get(structure_url).text

# Save POSCAR
with open("GaN.vasp", "w") as f:
    f.write(poscar)

# Parse with pymatgen
from pymatgen.core import Structure
structure = Structure.from_str(poscar, fmt="poscar")

# Or with ASE
from ase.io import read
from io import StringIO
atoms = read(StringIO(poscar), format='vasp')
```

## Practical Examples

### Example 1: Find Stable Semiconductors

```python
import requests

# Search for stable semiconductors
url = ("http://aflowlib.duke.edu/search/API/?"
       "enthalpy_formation_atom(*,0*),"  # Stable (negative formation enthalpy)
       "Egap(0.5*,4*),"                  # Semiconducting
       "nspecies(2)")                     # Binary compounds

results = requests.get(url).json()

for result in results[:10]:  # First 10 results
    print(f"Compound: {result.get('compound')}")
    print(f"Band gap: {result.get('Egap')} eV")
    print(f"Formation enthalpy: {result.get('enthalpy_formation_atom')} eV/atom")
    print(f"Space group: {result.get('spacegroup_relax')}")
    print("---")
```

### Example 2: Download Structures for a System

```python
import requests
import os

# Find all Ti-O compounds
search_url = "http://aflowlib.duke.edu/search/API/?species(Ti,O),nspecies(2)"
results = requests.get(search_url).json()

# Create directory
os.makedirs("TiO_structures", exist_ok=True)

# Download each structure
for i, result in enumerate(results[:20]):  # First 20
    aurl = result['aurl']
    compound = result.get('compound', f'structure_{i}')

    # Get POSCAR
    structure_url = f"{aurl}?geometry"
    poscar = requests.get(structure_url).text

    # Save
    filename = f"TiO_structures/{compound.replace('/', '_')}.vasp"
    with open(filename, "w") as f:
        f.write(poscar)

    print(f"Downloaded: {compound}")
```

### Example 3: Compare Elastic Properties

```python
import requests
import pandas as pd

# Get elastic data for metals
url = ("http://aflowlib.duke.edu/search/API/?"
       "nspecies(1),bulk_modulus_voigt")

results = requests.get(url).json()

# Convert to DataFrame
data = []
for result in results:
    data.append({
        'Element': result.get('species', ['Unknown'])[0],
        'Bulk_Modulus_GPa': result.get('bulk_modulus_voigt'),
        'Shear_Modulus_GPa': result.get('shear_modulus_voigt'),
        'Density_g_cm3': result.get('density')
    })

df = pd.DataFrame(data)
print(df.sort_values('Bulk_Modulus_GPa', ascending=False).head(10))
```

### Example 4: Using aflow Package

```python
import aflow

# Search for high-temperature superconductors (high Debye temperature)
results = aflow.search(
    filter='species(Y,Ba,Cu,O),nspecies(4)'
)

for entry in results:
    print(f"Compound: {entry.compound}")
    print(f"Formation enthalpy: {entry.enthalpy_formation_atom} eV/atom")

    # Get structure
    atoms = entry.atoms
    atoms.write(f"{entry.compound.replace('/', '_')}.cif")
```

## Performance Tips

1. **Be specific in searches** - Broad searches can return thousands of results
   ```python
   # Too broad (slow)
   url = "http://aflowlib.duke.edu/search/API/?species(O)"

   # Better (specific)
   url = "http://aflowlib.duke.edu/search/API/?species(Ti,O),nspecies(2),Egap(2*,3*)"
   ```

2. **Request specific properties** - Only get what you need
   ```python
   # Get just the band gap
   url = "http://aflowlib.duke.edu/AFLOWDATA/ICSD_WEB/FCC/Ag1/?Egap"
   ```

3. **Use batch processing** for large downloads
   ```python
   # Process in chunks
   results = requests.get(search_url).json()
   for batch in [results[i:i+100] for i in range(0, len(results), 100)]:
       process_batch(batch)
   ```

4. **Cache results** locally
   ```python
   import pickle

   # Save search results
   with open("aflow_results.pkl", "wb") as f:
       pickle.dump(results, f)
   ```

## Common Issues

**Timeout errors**
- AFLOW database is large; narrow search criteria
- Add timeout to requests: `requests.get(url, timeout=30)`

**No results**
- Check filter syntax (use `*` for ranges)
- Try broader search first, then narrow down
- Some properties not available for all entries

**Structure parsing errors**
- POSCAR format may need cleanup
- Use try-except when parsing:
  ```python
  try:
      structure = Structure.from_str(poscar, fmt="poscar")
  except Exception as e:
      print(f"Could not parse structure: {e}")
  ```

## Resources

- **AFLOW Homepage**: https://aflow.org
- **Documentation**: https://aflow.org/documentation/
- **REST API**: http://aflowlib.duke.edu
- **Python Package**: https://github.com/rosenbrockc/aflow
- **Publications**: https://aflow.org/about/publications/
