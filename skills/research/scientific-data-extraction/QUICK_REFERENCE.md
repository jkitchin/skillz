# Quick Reference: Scientific Data Extraction

## Decision Tree: Choose Your Tool

```
START: What format is your source?
│
├─► PDF
│   ├─► Is it scanned/image-based?
│   │   ├─► Yes → Marker-PDF with OCR
│   │   └─► No → Continue below
│   │
│   ├─► Do you need table data?
│   │   ├─► Yes → See Table Extraction below
│   │   └─► No → Continue
│   │
│   ├─► What's your priority?
│   │   ├─► Speed → PyMuPDF4LLM
│   │   ├─► Structure → GROBID or Docling
│   │   └─► References → GROBID
│   │
│   └─► Is it chemistry/materials?
│       ├─► Yes → ChemDataExtractor after text extraction
│       └─► No → LLM structured extraction
│
├─► HTML
│   └─► BeautifulSoup + lxml
│       └─► Then domain-specific processing
│
├─► Image
│   ├─► Document scan → OCR (Tesseract/Surya)
│   ├─► Graph/Plot → WebPlotDigitizer or LLM vision
│   ├─► Table image → Table Transformer
│   └─► Chemical structure → OSRA/MolVec
│
└─► Plain Text
    ├─► Chemistry → ChemDataExtractor
    ├─► Structured → Regex + NLP
    └─► Unstructured → LLM extraction
```

## Table Extraction Decision

```
What type of table?
│
├─► Bordered (visible grid lines)
│   └─► Camelot (lattice flavor)
│       camelot.read_pdf(f, flavor='lattice')
│
├─► Borderless (whitespace-separated)
│   └─► Camelot (stream flavor)
│       camelot.read_pdf(f, flavor='stream')
│
├─► Complex (merged cells, nested headers)
│   └─► pdfplumber with custom settings
│       page.extract_table(table_settings)
│
├─► Image of table
│   └─► Table Transformer or LLM vision
│
└─► Multiple tables per page
    └─► Tabula or pdfplumber
        tabula.read_pdf(f, pages='all')
```

## Quick Code Snippets

### PDF to Text (Fast)

```python
import pymupdf4llm
text = pymupdf4llm.to_markdown("paper.pdf")
```

### PDF to Structured (GROBID)

```python
import scipdf_parser
article = scipdf_parser.parse_pdf_to_dict("paper.pdf")
# Returns: title, abstract, sections, references
```

### PDF Tables (Camelot)

```python
import camelot

# Bordered tables
tables = camelot.read_pdf("paper.pdf", flavor='lattice')
df = tables[0].df

# Borderless tables
tables = camelot.read_pdf("paper.pdf", flavor='stream')

# Specific pages
tables = camelot.read_pdf("paper.pdf", pages='3-5')
```

### PDF Tables (pdfplumber)

```python
import pdfplumber

with pdfplumber.open("paper.pdf") as pdf:
    for page in pdf.pages:
        tables = page.extract_tables()
        for table in tables:
            # table is list of lists
            df = pd.DataFrame(table[1:], columns=table[0])
```

### Chemistry Extraction

```python
from chemdataextractor import Document

doc = Document.from_file("paper.pdf")
for record in doc.records:
    print(record.serialize())
```

### HTML Parsing

```python
from bs4 import BeautifulSoup
import requests

html = requests.get(url).text
soup = BeautifulSoup(html, 'lxml')

# Extract tables
tables = soup.find_all('table')
for table in tables:
    rows = table.find_all('tr')
    data = [[cell.get_text(strip=True) for cell in row.find_all(['td','th'])]
            for row in rows]
```

### LLM Structured Extraction

```python
prompt = """Extract all numerical data as JSON:
{
  "data": [
    {"property": "", "value": null, "unit": "", "entity": ""}
  ]
}

Text:
{text}
"""
```

## Confidence Score Quick Guide

| Scenario | Score | Level |
|----------|-------|-------|
| Single method only | 0.5-0.7 | LOW-MEDIUM |
| Two methods agree | 0.8 | MEDIUM |
| Two methods + LLM verify | 0.9 | HIGH |
| Multiple methods + DB match | 0.95+ | HIGH |
| Methods disagree | Flag | REVIEW |

## Domain Detection Keywords

### Chemistry/Materials

- Formulas: H2O, NaCl, TiO2, CH3OH
- SMILES: CC(=O)O, c1ccccc1
- InChI: InChI=1S/...
- Arrows: →, ⟶, ⇌
- Properties: bandgap, conductivity, melting point, yield
- Units: eV, S/cm, K, mol/L, ppm

### Use Specialized Tools When Detected

- ChemDataExtractor: property extraction, NER
- OpenChemIE: reaction extraction

## Common Issues & Solutions

| Problem | Solution |
|---------|----------|
| Tables not detected | Try different flavor (lattice vs stream) |
| Merged cells broken | Use pdfplumber with `merge_tolerance` |
| Scanned PDF | Use Marker-PDF or Docling with OCR |
| Chemistry missed | Verify document has text layer |
| Wrong units extracted | Add unit normalization post-processing |
| Graph digitization inaccurate | Ensure proper axis calibration |

## Output Format Template

```json
{
  "extraction_metadata": {
    "source": "file.pdf",
    "source_type": "pdf",
    "domain_detected": "chemistry|general",
    "methods_used": ["method1", "method2"],
    "timestamp": "ISO-8601"
  },
  "extracted_data": [
    {
      "data_type": "property|reaction|entity",
      "entity": "compound name",
      "property": "property name",
      "value": 123.4,
      "unit": "unit",
      "source_location": {
        "page": 1,
        "section": "Results",
        "table_id": "Table 1"
      },
      "confidence": {
        "score": 0.95,
        "level": "HIGH",
        "methods_agreed": [],
        "verification_notes": ""
      }
    }
  ],
  "validation_summary": {
    "total_extracted": 0,
    "high_confidence": 0,
    "needs_review": 0
  }
}
```

## Installation One-Liner

```bash
# Core tools
pip install pymupdf4llm pdfplumber camelot-py[cv] beautifulsoup4 lxml pandas spacy

# Chemistry tools
pip install chemdataextractor2 openchemie

# System deps (macOS)
brew install ghostscript tesseract poppler
```

## Tool Speed Comparison

| Tool | ~Time/Page | Use Case |
|------|------------|----------|
| PyMuPDF4LLM | 0.01s | Quick text extraction |
| pdfplumber | 0.1s | Detailed analysis |
| Camelot | 0.5s | Table extraction |
| GROBID | 1-2s | Full structure |
| Docling | 2-3s | Complex layouts |
| Marker-PDF | 10s+ | Scanned with OCR |
| ChemDataExtractor | 5-10s | Chemistry NLP |
