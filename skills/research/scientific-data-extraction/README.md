# Scientific Data Extraction Skill

A comprehensive guidance skill for extracting structured data from scientific literature across multiple formats (PDF, HTML, images, plain text). The skill auto-detects scientific domains to recommend specialized tools for chemistry and materials science, and provides a hierarchical extraction approach with multi-method validation.

## Features

- **Format-aware extraction**: Automatic detection and routing for PDF, HTML, images, and plain text
- **Domain detection**: Identifies chemistry/materials content to apply specialized tools
- **Hierarchical methods**: From quick screening to publication-quality extraction
- **Multi-method validation**: Cross-check results using multiple tools
- **Confidence scoring**: Track extraction reliability with quantified confidence
- **Table extraction**: Multiple methods for different table types
- **Graph digitization**: Extract numerical data from plots and charts
- **LLM-enhanced verification**: Use language models to validate extracted data

## Quick Start

### PDF Extraction

```python
# Quick extraction with PyMuPDF4LLM
import pymupdf4llm
text = pymupdf4llm.to_markdown("paper.pdf")

# Standard extraction with GROBID
import scipdf_parser
article = scipdf_parser.parse_pdf_to_dict("paper.pdf")

# Table extraction with Camelot
import camelot
tables = camelot.read_pdf("paper.pdf", flavor='lattice')
df = tables[0].df
```

### Chemistry-Specific Extraction

```python
# Property extraction with ChemDataExtractor
from chemdataextractor import Document
doc = Document.from_file("paper.pdf")
for record in doc.records:
    print(record.serialize())

# Reaction extraction with OpenChemIE
from openchemie import OpenChemIE
model = OpenChemIE()
reactions = model.extract_reactions_from_text(text)
```

### Table Extraction Decision

| Table Type | Tool | Code |
|------------|------|------|
| Bordered | Camelot | `camelot.read_pdf(f, flavor='lattice')` |
| Borderless | Camelot | `camelot.read_pdf(f, flavor='stream')` |
| Complex | pdfplumber | `page.extract_table(table_settings)` |
| General | Tabula | `tabula.read_pdf(f)` |

## Installation

### Core Dependencies

```bash
# Basic PDF and text processing
pip install pymupdf4llm pdfplumber beautifulsoup4 lxml spacy pandas

# Table extraction
pip install camelot-py[cv] tabula-py

# Install spaCy model
python -m spacy download en_core_web_sm
```

### Chemistry Tools

```bash
# ChemDataExtractor v2
pip install chemdataextractor2

# OpenChemIE
pip install openchemie
```

### Optional Dependencies

```bash
# GROBID client (requires running GROBID server)
pip install scipdf_parser

# IBM Docling
pip install docling

# Marker-PDF (OCR-capable)
pip install marker-pdf

# OCR
pip install pytesseract  # requires tesseract system install
# or
pip install surya-ocr    # ML-based OCR
```

### System Dependencies

```bash
# macOS
brew install ghostscript tesseract poppler openjdk

# Ubuntu/Debian
sudo apt-get install ghostscript tesseract-ocr poppler-utils default-jdk

# For GROBID server (optional)
docker pull lfoppiano/grobid:0.8.0
docker run -t --rm -p 8070:8070 lfoppiano/grobid:0.8.0
```

## Tool Comparison

### PDF Extraction Tools

| Tool | Speed | Structure | Tables | Figures | OCR | Best For |
|------|-------|-----------|--------|---------|-----|----------|
| PyMuPDF4LLM | Very Fast | Basic | Limited | No | Optional | Quick screening, RAG |
| GROBID | Medium | Excellent | Good | Metadata | No | Academic papers, references |
| Docling | Medium | Good | Excellent | Basic | Optional | Complex layouts |
| Marker-PDF | Slow | Good | Good | Good | Yes | Scanned documents |
| pdfplumber | Fast | Detailed | Excellent | No | No | Table-heavy documents |

### Table Extraction Tools

| Tool | Bordered | Borderless | Complex | Speed | Accuracy |
|------|----------|------------|---------|-------|----------|
| Camelot (lattice) | Excellent | Poor | Fair | Fast | High |
| Camelot (stream) | Good | Good | Fair | Fast | Medium |
| Tabula | Good | Good | Fair | Fast | Medium |
| pdfplumber | Excellent | Fair | Excellent | Medium | High |
| Table Transformer | Good | Good | Excellent | Slow | High |

### Chemistry Extraction Tools

| Tool | Properties | Reactions | Entities | Tables | Figures |
|------|------------|-----------|----------|--------|---------|
| ChemDataExtractor v2 | Excellent | Fair | Excellent | Good | No |
| OpenChemIE | Good | Excellent | Good | Good | Yes |
| LLM (GPT-4/Claude) | Good | Good | Good | Fair | Good |

## Validation Strategy

### Multi-Method Validation

For high-confidence extraction, use multiple methods:

1. **Primary extraction**: Use best tool for the content type
2. **Secondary extraction**: Run alternative tool for comparison
3. **LLM verification**: Ask targeted questions about extracted values
4. **Database cross-reference**: Compare against known databases

### Confidence Scoring

```python
def calculate_confidence(primary_result, secondary_result, llm_verified, db_match):
    """Calculate extraction confidence score."""
    score = 0.0

    # Base score from primary extraction
    score = 0.5

    # Agreement between methods
    if primary_result == secondary_result:
        score += 0.25

    # LLM verification
    if llm_verified:
        score += 0.15

    # Database match
    if db_match:
        score += 0.10

    # Classify confidence level
    if score >= 0.9:
        level = "HIGH"
    elif score >= 0.7:
        level = "MEDIUM"
    elif score >= 0.5:
        level = "LOW"
    else:
        level = "REVIEW"

    return {"score": score, "level": level}
```

## Directory Structure

```
scientific-data-extraction/
├── SKILL.md                    # Main skill file
├── README.md                   # This documentation
├── QUICK_REFERENCE.md          # Quick lookup guide
├── references/
│   ├── pdf-tools.md            # PDF tool comparison and usage
│   ├── table-extraction.md     # Table extraction methods
│   ├── graph-digitization.md   # Graph data extraction
│   ├── chemistry-tools.md      # ChemDataExtractor, OpenChemIE
│   └── llm-extraction.md       # LLM-based extraction patterns
└── examples/
    ├── extract-from-pdf.md     # Complete PDF workflow
    ├── extract-table-data.md   # Table extraction examples
    ├── digitize-graph.md       # Graph digitization guide
    └── chemistry-extraction.md # Chemistry-specific workflow
```

## Key References

### Academic Papers

- Swain & Cole (2016). "ChemDataExtractor: A Toolkit for Automated Extraction of Chemical Information from the Scientific Literature." *J. Chem. Inf. Model.* [DOI: 10.1021/acs.jcim.6b00207](https://pubs.acs.org/doi/10.1021/acs.jcim.6b00207)

- Mavracic et al. (2021). "ChemDataExtractor 2.0: Autopopulated Ontologies for Materials Science." *J. Chem. Inf. Model.* [DOI: 10.1021/acs.jcim.1c00446](https://pubs.acs.org/doi/10.1021/acs.jcim.1c00446)

- Fan et al. (2024). "OpenChemIE: An Information Extraction Toolkit for Chemistry Literature." *J. Chem. Inf. Model.* [DOI: 10.1021/acs.jcim.4c00572](https://pubs.acs.org/doi/abs/10.1021/acs.jcim.4c00572)

- Zheng et al. (2024). "ChatExtract: Accurate Data Extraction from Scientific Publications via Large Language Models." [arXiv:2308.02214](https://arxiv.org/abs/2308.02214)

### Tool Documentation

- [GROBID Documentation](https://grobid.readthedocs.io/)
- [Docling Documentation](https://ds4sd.github.io/docling/)
- [Camelot Documentation](https://camelot-py.readthedocs.io/)
- [ChemDataExtractor Documentation](https://chemdataextractor2.readthedocs.io/)
- [WebPlotDigitizer](https://automeris.io/WebPlotDigitizer/)

## Related Skills

- **literature-review**: Systematic literature searching and synthesis
- **scientific-reviewer**: Evaluating extracted data quality
- **materials-databases**: Cross-referencing Materials Project, AFLOW
- **python-plotting**: Visualizing extracted data

## Contributing

To improve this skill:

1. Add new tool comparisons to `references/`
2. Include additional examples in `examples/`
3. Update extraction patterns as tools evolve
4. Report accuracy issues or validation failures

## License

This skill documentation is provided under MIT License. Individual tools have their own licenses:
- ChemDataExtractor: MIT
- OpenChemIE: MIT
- GROBID: Apache 2.0
- Docling: MIT
- Camelot: MIT
