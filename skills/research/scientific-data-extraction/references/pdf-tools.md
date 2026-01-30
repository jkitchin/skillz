# PDF Extraction Tools Reference

This document provides detailed comparison and usage guidance for PDF extraction tools.

## Tool Overview

### PyMuPDF4LLM

**Purpose**: Ultra-fast PDF to Markdown conversion optimized for LLM/RAG applications.

**Best for**: Quick text extraction, large document batches, initial screening.

**Installation**:
```bash
pip install pymupdf4llm
```

**Basic Usage**:
```python
import pymupdf4llm

# Simple conversion
text = pymupdf4llm.to_markdown("paper.pdf")

# With page selection
text = pymupdf4llm.to_markdown("paper.pdf", pages=[0, 1, 2])

# Get structured output
result = pymupdf4llm.to_markdown("paper.pdf", page_chunks=True)
for chunk in result:
    print(f"Page {chunk['metadata']['page']}: {chunk['text'][:100]}...")
```

**Features**:
- Very fast (~0.01s per page)
- Markdown formatting preserved
- Automatic OCR triggering for scanned pages
- Table detection (basic)
- Image extraction

**Limitations**:
- Limited table structure preservation
- No reference parsing
- No section identification

---

### GROBID

**Purpose**: Machine learning library for extracting and parsing academic PDFs into structured XML/TEI.

**Best for**: Academic papers, reference extraction, structured content.

**Installation**:
```bash
# Run GROBID server via Docker
docker pull lfoppiano/grobid:0.8.0
docker run -t --rm -p 8070:8070 lfoppiano/grobid:0.8.0

# Python client
pip install scipdf_parser
# or
pip install grobid_client_python
```

**Basic Usage**:
```python
import scipdf_parser

# Parse to dictionary
article = scipdf_parser.parse_pdf_to_dict("paper.pdf")

# Access structured content
print(article['title'])
print(article['abstract'])
for section in article['sections']:
    print(f"{section['heading']}: {section['text'][:100]}...")
for ref in article['references']:
    print(f"- {ref['title']}")
```

**Using grobid_client directly**:
```python
from grobid_client.grobid_client import GrobidClient

client = GrobidClient(config_path="./config.json")

# Process single PDF
client.process("processFulltextDocument", "path/to/pdf", output="output/")

# Batch processing
client.process("processFulltextDocument", "path/to/pdfs/", output="output/", n=10)
```

**Features**:
- 68 label types for document structure
- Reference parsing (F1 > 0.90)
- Section identification
- Author/affiliation extraction
- PDF coordinates and bounding boxes
- Supports CRF and BERT-CRF models

**Limitations**:
- Requires running server
- Slower than local tools
- Optimized for academic papers (may struggle with other document types)

---

### Docling (IBM)

**Purpose**: Document conversion to structured formats (Markdown, HTML, JSON) with layout analysis.

**Best for**: Complex layouts, table-heavy documents, multi-format support.

**Installation**:
```bash
pip install docling
```

**Basic Usage**:
```python
from docling.document_converter import DocumentConverter

# Initialize converter
converter = DocumentConverter()

# Convert PDF
result = converter.convert("paper.pdf")

# Export to Markdown
markdown = result.document.export_to_markdown()

# Export to JSON
json_doc = result.document.export_to_dict()

# Access document structure
for element in result.document.iterate_items():
    print(f"Type: {element.label}, Text: {element.text[:50]}...")
```

**Table extraction**:
```python
# Get tables specifically
for table in result.document.tables:
    # Export as DataFrame
    df = table.export_to_dataframe()
    print(df)
```

**Features**:
- DocLayNet for layout analysis
- TableFormer for table structure
- Reading order reconstruction
- Section hierarchy identification
- Supports PDF, DOCX, PPTX, HTML, images
- Optional OCR integration

**Limitations**:
- Larger memory footprint
- Slower than lightweight tools
- Model downloads required

---

### Marker-PDF

**Purpose**: High-quality PDF to Markdown conversion with integrated OCR.

**Best for**: Scanned documents, high-quality extraction when speed is not critical.

**Installation**:
```bash
pip install marker-pdf
```

**Basic Usage**:
```python
from marker.converters.pdf import PdfConverter
from marker.models import create_model_dict

# Create models (downloads on first use)
model_lst = create_model_dict()

# Initialize converter
converter = PdfConverter(artifact_dict=model_lst)

# Convert PDF
result = converter("paper.pdf")

# Get markdown output
markdown = result.markdown
```

**With Surya OCR**:
```python
# Marker uses Surya for OCR automatically
# For scanned documents, OCR is triggered based on text density

# Force OCR
result = converter("scanned.pdf", force_ocr=True)
```

**Features**:
- Integrated Surya OCR (multilingual)
- GPU/CPU/MPS acceleration
- Math formula extraction
- Multi-column layout handling
- Image extraction and captioning

**Limitations**:
- Slower than text-based extraction
- Large model downloads (~2GB)
- GPU recommended for speed

---

### pdfplumber

**Purpose**: Detailed PDF analysis with fine-grained control over extraction.

**Best for**: Complex tables, precise text positioning, debugging extraction issues.

**Installation**:
```bash
pip install pdfplumber
```

**Basic Usage**:
```python
import pdfplumber

with pdfplumber.open("paper.pdf") as pdf:
    # Iterate pages
    for page in pdf.pages:
        # Extract text
        text = page.extract_text()

        # Extract tables
        tables = page.extract_tables()
        for table in tables:
            print(table)

        # Get words with positions
        words = page.extract_words()
        for word in words:
            print(f"{word['text']} at ({word['x0']}, {word['top']})")
```

**Advanced table extraction**:
```python
# Custom table settings
table_settings = {
    "vertical_strategy": "lines",
    "horizontal_strategy": "lines",
    "explicit_vertical_lines": [],
    "explicit_horizontal_lines": [],
    "snap_tolerance": 3,
    "join_tolerance": 3,
    "edge_min_length": 3,
    "min_words_vertical": 3,
    "min_words_horizontal": 1,
    "intersection_tolerance": 3,
}

with pdfplumber.open("paper.pdf") as pdf:
    page = pdf.pages[0]
    table = page.extract_table(table_settings)
```

**Visual debugging**:
```python
with pdfplumber.open("paper.pdf") as pdf:
    page = pdf.pages[0]

    # Visualize detected elements
    im = page.to_image()
    im.draw_rects(page.extract_words(), stroke="blue")
    im.draw_rects(page.find_tables(), stroke="red")
    im.save("debug.png")
```

**Features**:
- Built on pdfminer.six
- Fine-grained control over extraction
- Visual debugging tools
- Precise text positioning
- Excellent for complex tables

**Limitations**:
- Text-based PDFs only (no OCR)
- Manual tuning often required
- No section/structure identification

---

## Comparison Matrix

| Feature | PyMuPDF4LLM | GROBID | Docling | Marker | pdfplumber |
|---------|-------------|--------|---------|--------|------------|
| Speed | Very Fast | Medium | Medium | Slow | Fast |
| Text | Excellent | Excellent | Excellent | Excellent | Excellent |
| Tables | Basic | Good | Excellent | Good | Excellent |
| Structure | Basic | Excellent | Good | Good | None |
| References | No | Excellent | Basic | No | No |
| OCR | Optional | No | Optional | Yes | No |
| Scanned PDFs | Limited | No | Limited | Yes | No |
| Academic Papers | Good | Excellent | Good | Good | Fair |
| Memory | Low | Low | Medium | High | Low |
| Setup | Easy | Docker | Easy | Easy | Easy |

## Selection Guide

### Use PyMuPDF4LLM when:
- You need quick text extraction
- Processing many documents for initial screening
- Building RAG systems
- Text quality is more important than structure

### Use GROBID when:
- Processing academic papers
- Need reference extraction
- Require structured output (sections, abstracts)
- Building citation databases

### Use Docling when:
- Documents have complex layouts
- Tables are important
- Need multi-format support
- Want layout-aware extraction

### Use Marker-PDF when:
- Documents are scanned
- High-quality extraction is required
- Have GPU available
- Processing multilingual documents

### Use pdfplumber when:
- Need fine-grained table extraction
- Debugging extraction issues
- Require precise text positioning
- Building custom extraction pipelines

## Performance Benchmarks

Approximate processing times (single page, typical academic PDF):

| Tool | Time/Page | Notes |
|------|-----------|-------|
| PyMuPDF4LLM | ~0.01s | CPU only |
| pdfplumber | ~0.1s | CPU only |
| Camelot | ~0.5s | Per table |
| GROBID | ~1-2s | Server overhead |
| Docling | ~2-3s | Model inference |
| Marker-PDF | ~10s | With OCR, GPU |
| Marker-PDF | ~30s | With OCR, CPU |

## References

- [PyMuPDF4LLM GitHub](https://github.com/pymupdf/PyMuPDF-Utilities)
- [GROBID Documentation](https://grobid.readthedocs.io/)
- [Docling Documentation](https://ds4sd.github.io/docling/)
- [Marker GitHub](https://github.com/VikParuchuri/marker)
- [pdfplumber Documentation](https://github.com/jsvine/pdfplumber)
