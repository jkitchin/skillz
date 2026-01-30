# Chemistry Data Extraction Tools Reference

This document covers specialized tools for extracting chemical data from scientific literature.

## ChemDataExtractor v2

### Overview

ChemDataExtractor is a toolkit for automatically extracting chemical information from scientific documents. Version 2 (released 2021) adds auto-populated ontologies, nested hierarchies, and improved table parsing.

**Key Capabilities**:
- Named Entity Recognition (NER) for chemical compounds
- Property extraction (melting point, bandgap, conductivity, etc.)
- Table parsing and relation extraction
- Document-level interdependency resolution
- SMILES/InChI extraction

**Performance** (from publications):
- Chemical identifiers: F1 = 93.4%
- Spectroscopic data: F1 = 86.8%
- Property attributes: F1 = 91.5%

### Installation

```bash
# Install ChemDataExtractor v2
pip install chemdataextractor2

# Download required models
cde data download
```

### Basic Usage

```python
from chemdataextractor import Document

# Load from file
doc = Document.from_file("paper.pdf")

# Or from text
text = """
The melting point of sodium chloride (NaCl) is 801°C.
TiO2 has a bandgap of 3.2 eV.
"""
doc = Document(text)

# Extract all records
for record in doc.records:
    print(record.serialize())
```

**Output Example**:
```python
{
    'MeltingPoint': {
        'compound': {'names': ['sodium chloride'], 'labels': ['NaCl']},
        'value': [{'value': '801', 'units': '°C'}]
    }
}
```

### Property Extraction

```python
from chemdataextractor import Document
from chemdataextractor.model import Compound, MeltingPoint, GlassTransition
from chemdataextractor.model.units import TemperatureModel

# Define custom property model
class BandGap(TemperatureModel):
    """Bandgap property model."""
    specifier = 'bandgap'
    parsers = [...]

# Load document
doc = Document.from_file("materials_paper.pdf")

# Get specific property types
melting_points = []
for record in doc.records:
    data = record.serialize()
    if 'MeltingPoint' in data:
        melting_points.append(data['MeltingPoint'])

print(f"Found {len(melting_points)} melting point records")
```

### Table Extraction

ChemDataExtractor excels at extracting data from tables:

```python
from chemdataextractor import Document

doc = Document.from_file("paper.pdf")

# Access tables
for table in doc.tables:
    print(f"Table caption: {table.caption}")

    # Get table records
    for record in table.records:
        print(record.serialize())
```

### Custom Parsers

Define custom extraction rules for specific properties:

```python
from chemdataextractor.parse import R, I, W, Optional, ZeroOrMore
from chemdataextractor.parse.cem import chemical_name
from chemdataextractor.parse.common import lbrct, rbrct
from chemdataextractor.model import BaseModel, StringType, ModelType

# Define property model
class Conductivity(BaseModel):
    value = StringType()
    units = StringType()
    compound = ModelType(Compound)

# Define parser
conductivity_value = R(r'\d+\.?\d*')
conductivity_unit = (I('S') + W('/') + I('cm')) | (I('mS') + W('/') + I('cm'))

conductivity_phrase = (
    chemical_name('compound') +
    Optional(W(',')) +
    conductivity_value('value') +
    conductivity_unit('units')
)
```

### Integration with Document Parsers

```python
from chemdataextractor import Document
from chemdataextractor.reader import PdfReader, HtmlReader

# Explicit PDF parsing
reader = PdfReader()
doc = reader.read("paper.pdf")

# HTML parsing
reader = HtmlReader()
doc = reader.read("paper.html")

# Access document elements
for element in doc.elements:
    print(f"Type: {type(element).__name__}")
    print(f"Text: {element.text[:100]}...")
```

---

## OpenChemIE

### Overview

OpenChemIE (2024) is a multi-modal toolkit for extracting chemical information from text, tables, and figures. It's particularly strong at reaction extraction.

**Key Capabilities**:
- Reaction extraction from text
- Reaction scheme image parsing
- Table data extraction
- SMILES output for molecules
- R-group handling

**Performance**:
- Reaction schemes with R-groups: F1 = 69.5%
- Accuracy vs Reaxys database: 64.3%

### Installation

```bash
pip install openchemie
```

### Basic Usage

```python
from openchemie import OpenChemIE

# Initialize
model = OpenChemIE()

# Extract reactions from text
text = """
The reaction of benzaldehyde with sodium borohydride in methanol
gives benzyl alcohol in 95% yield.
"""
reactions = model.extract_reactions_from_text(text)

for rxn in reactions:
    print(f"Reactants: {rxn['reactants']}")
    print(f"Products: {rxn['products']}")
    print(f"Conditions: {rxn['conditions']}")
```

### Reaction Extraction from Images

```python
from openchemie import OpenChemIE
from PIL import Image

model = OpenChemIE()

# Load reaction scheme image
image = Image.open("reaction_scheme.png")

# Extract reaction
reaction = model.extract_reaction_from_figure(image)

print(f"Reactant SMILES: {reaction['reactant_smiles']}")
print(f"Product SMILES: {reaction['product_smiles']}")
print(f"Reaction conditions: {reaction['conditions']}")
```

### Multi-Modal Extraction

```python
from openchemie import OpenChemIE

model = OpenChemIE()

# Process entire document with text and figures
results = model.process_document(
    pdf_path="paper.pdf",
    extract_text=True,
    extract_tables=True,
    extract_figures=True
)

# Access extracted reactions
for reaction in results['reactions']:
    print(reaction)

# Access extracted compounds
for compound in results['compounds']:
    print(f"Name: {compound['name']}, SMILES: {compound['smiles']}")
```

---

## Chemical Structure Recognition

### OSRA (Optical Structure Recognition Application)

**Purpose**: Convert chemical structure images to SMILES/MOL format.

```bash
# Install OSRA (requires system dependencies)
# Ubuntu: apt-get install osra
# macOS: brew install osra

# Usage from command line
osra -f sdf input_image.png > output.sdf
osra -f smi input_image.png > output.smi
```

**Python wrapper**:
```python
import subprocess

def image_to_smiles(image_path: str) -> str:
    """Convert structure image to SMILES using OSRA."""
    result = subprocess.run(
        ['osra', '-f', 'smi', image_path],
        capture_output=True,
        text=True
    )
    return result.stdout.strip()
```

### MolVec

**Purpose**: Deep learning-based molecule image recognition.

```bash
pip install molvec
```

```python
from molvec import MolVec

# Initialize model
model = MolVec()

# Convert image to SMILES
smiles = model.image_to_smiles("structure.png")
print(f"SMILES: {smiles}")

# Batch processing
images = ["struct1.png", "struct2.png", "struct3.png"]
results = [model.image_to_smiles(img) for img in images]
```

---

## RDKit for Chemical Processing

After extraction, RDKit helps validate and process chemical data:

```python
from rdkit import Chem
from rdkit.Chem import Descriptors, AllChem, Draw

def validate_smiles(smiles: str) -> dict:
    """Validate SMILES string and compute properties."""
    mol = Chem.MolFromSmiles(smiles)

    if mol is None:
        return {"valid": False, "error": "Invalid SMILES"}

    return {
        "valid": True,
        "canonical_smiles": Chem.MolToSmiles(mol),
        "molecular_formula": Chem.rdMolDescriptors.CalcMolFormula(mol),
        "molecular_weight": Descriptors.MolWt(mol),
        "num_atoms": mol.GetNumAtoms(),
        "num_rings": Descriptors.RingCount(mol)
    }

def compare_molecules(smiles1: str, smiles2: str) -> bool:
    """Check if two SMILES represent the same molecule."""
    mol1 = Chem.MolFromSmiles(smiles1)
    mol2 = Chem.MolFromSmiles(smiles2)

    if mol1 is None or mol2 is None:
        return False

    return Chem.MolToSmiles(mol1) == Chem.MolToSmiles(mol2)
```

---

## Best Practices

### 1. Preprocessing Documents

```python
def preprocess_chemistry_document(pdf_path: str) -> dict:
    """Prepare document for chemistry extraction."""
    # Extract text with structure preservation
    import pymupdf4llm
    text = pymupdf4llm.to_markdown(pdf_path)

    # Detect chemistry content
    chemistry_indicators = [
        r'\b[A-Z][a-z]?\d*[A-Z]?[a-z]?\d*\b',  # Chemical formulas
        r'→|⟶|⇌',  # Reaction arrows
        r'\b\d+\s*°C\b',  # Temperatures
        r'\b\d+\.?\d*\s*(eV|meV|kJ|kcal)\b',  # Energy units
        r'SMILES|InChI',  # Structure identifiers
    ]

    import re
    is_chemistry = any(re.search(pattern, text) for pattern in chemistry_indicators)

    return {
        "text": text,
        "is_chemistry": is_chemistry,
        "pdf_path": pdf_path
    }
```

### 2. Multi-Tool Extraction

```python
def extract_chemistry_multi_tool(pdf_path: str) -> dict:
    """Use multiple tools and combine results."""
    results = {
        "chemdataextractor": [],
        "openchemie": [],
        "combined": []
    }

    # ChemDataExtractor
    from chemdataextractor import Document
    doc = Document.from_file(pdf_path)
    for record in doc.records:
        results["chemdataextractor"].append(record.serialize())

    # OpenChemIE
    from openchemie import OpenChemIE
    model = OpenChemIE()
    openchemie_results = model.process_document(pdf_path)
    results["openchemie"] = openchemie_results

    # Combine and deduplicate
    # ... merging logic based on compound names, SMILES, etc.

    return results
```

### 3. Validation Against Databases

```python
def validate_against_pubchem(compound_name: str) -> dict:
    """Cross-reference extracted compound with PubChem."""
    import requests

    # Search PubChem
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{compound_name}/property/MolecularFormula,MolecularWeight,CanonicalSMILES/JSON"

    response = requests.get(url)
    if response.status_code == 200:
        data = response.json()
        props = data['PropertyTable']['Properties'][0]
        return {
            "found": True,
            "cid": props.get('CID'),
            "formula": props.get('MolecularFormula'),
            "weight": props.get('MolecularWeight'),
            "smiles": props.get('CanonicalSMILES')
        }
    return {"found": False}
```

### 4. Unit Normalization

```python
def normalize_property_units(value: float, from_unit: str, to_unit: str, property_type: str) -> float:
    """Normalize property values to standard units."""
    # Temperature conversions
    if property_type == "temperature":
        if from_unit == "°C" and to_unit == "K":
            return value + 273.15
        elif from_unit == "°F" and to_unit == "K":
            return (value - 32) * 5/9 + 273.15
        elif from_unit == "K" and to_unit == "°C":
            return value - 273.15

    # Energy conversions
    if property_type == "energy":
        conversions = {
            ("eV", "kJ/mol"): lambda x: x * 96.485,
            ("kJ/mol", "eV"): lambda x: x / 96.485,
            ("kcal/mol", "kJ/mol"): lambda x: x * 4.184,
        }
        key = (from_unit, to_unit)
        if key in conversions:
            return conversions[key](value)

    return value
```

---

## Output Format

Standardize extracted chemistry data:

```python
chemistry_extraction_schema = {
    "compounds": [
        {
            "name": "compound name",
            "identifiers": {
                "smiles": "canonical SMILES",
                "inchi": "InChI string",
                "cas": "CAS number"
            },
            "formula": "molecular formula",
            "properties": [
                {
                    "name": "property name",
                    "value": 123.4,
                    "unit": "unit",
                    "conditions": {
                        "temperature": 298,
                        "pressure": 1
                    },
                    "source_location": {"page": 3, "table": "Table 1"}
                }
            ]
        }
    ],
    "reactions": [
        {
            "reactants": [{"smiles": "", "name": ""}],
            "products": [{"smiles": "", "name": ""}],
            "conditions": {
                "temperature": 100,
                "solvent": "methanol",
                "catalyst": "Pd/C",
                "time": "2 h"
            },
            "yield": 95,
            "source_location": {"page": 5, "scheme": "Scheme 1"}
        }
    ]
}
```

## References

- Swain & Cole (2016). "ChemDataExtractor: A Toolkit for Automated Extraction." *JCIM*. [DOI: 10.1021/acs.jcim.6b00207](https://pubs.acs.org/doi/10.1021/acs.jcim.6b00207)

- Mavracic et al. (2021). "ChemDataExtractor 2.0." *JCIM*. [DOI: 10.1021/acs.jcim.1c00446](https://pubs.acs.org/doi/10.1021/acs.jcim.1c00446)

- Fan et al. (2024). "OpenChemIE: An Information Extraction Toolkit." *JCIM*. [DOI: 10.1021/acs.jcim.4c00572](https://pubs.acs.org/doi/abs/10.1021/acs.jcim.4c00572)

- [ChemDataExtractor Documentation](https://chemdataextractor2.readthedocs.io/)
- [OpenChemIE GitHub](https://github.com/CederGroupHub/openchemie)
- [RDKit Documentation](https://www.rdkit.org/docs/)
