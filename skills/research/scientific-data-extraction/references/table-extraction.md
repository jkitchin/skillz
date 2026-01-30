# Table Extraction Reference

This document covers methods and tools for extracting tabular data from PDFs and images.

## Tool Overview

### Camelot

**Purpose**: PDF table extraction using two methods: lattice (bordered) and stream (borderless).

**Installation**:
```bash
pip install camelot-py[cv]

# System dependencies
# macOS: brew install ghostscript
# Ubuntu: apt-get install ghostscript
```

**Lattice Mode** (bordered tables):
```python
import camelot

# Extract bordered tables
tables = camelot.read_pdf("paper.pdf", flavor='lattice')

print(f"Found {len(tables)} tables")

# Access first table
df = tables[0].df
print(df)

# Get extraction accuracy
print(tables[0].accuracy)

# Export
tables[0].to_csv("table.csv")
tables[0].to_excel("table.xlsx")
tables[0].to_json("table.json")
```

**Stream Mode** (borderless tables):
```python
# Extract borderless tables
tables = camelot.read_pdf("paper.pdf", flavor='stream')

# Tune stream parameters
tables = camelot.read_pdf(
    "paper.pdf",
    flavor='stream',
    edge_tol=50,        # Tolerance for text edge detection
    row_tol=10,         # Tolerance for grouping rows
    column_tol=0        # Tolerance for column detection
)
```

**Page Selection**:
```python
# Specific pages
tables = camelot.read_pdf("paper.pdf", pages='1,3,5')

# Page range
tables = camelot.read_pdf("paper.pdf", pages='1-5')

# All pages
tables = camelot.read_pdf("paper.pdf", pages='all')
```

**Table Areas** (specify region):
```python
# Define table area: (x1, y1, x2, y2) from top-left
tables = camelot.read_pdf(
    "paper.pdf",
    flavor='lattice',
    table_areas=['100,700,500,300']  # Single area
)

# Multiple areas
tables = camelot.read_pdf(
    "paper.pdf",
    flavor='lattice',
    table_areas=['100,700,300,400', '320,700,520,400']
)
```

**Handling Edge Cases**:
```python
# Suppress warnings
import warnings
warnings.filterwarnings('ignore')

# Copy text from cells spanning multiple pages
tables = camelot.read_pdf(
    "paper.pdf",
    flavor='lattice',
    copy_text=['v']  # Copy text vertically
)

# Strip whitespace
tables = camelot.read_pdf(
    "paper.pdf",
    flavor='stream',
    strip_text='\n'
)
```

---

### Tabula / tabula-py

**Purpose**: Java-based table extraction with Python wrapper.

**Installation**:
```bash
pip install tabula-py

# Requires Java
# macOS: brew install openjdk
# Ubuntu: apt-get install default-jdk
```

**Basic Usage**:
```python
import tabula

# Read tables
dfs = tabula.read_pdf("paper.pdf", pages='all')

# Returns list of DataFrames
for i, df in enumerate(dfs):
    print(f"Table {i+1}:")
    print(df)
```

**Extraction Modes**:
```python
# Lattice mode (bordered)
dfs = tabula.read_pdf("paper.pdf", lattice=True)

# Stream mode (borderless)
dfs = tabula.read_pdf("paper.pdf", stream=True)

# Guess mode (auto-detect)
dfs = tabula.read_pdf("paper.pdf", guess=True)
```

**Specify Area**:
```python
# Area: top, left, height, width (in points)
dfs = tabula.read_pdf(
    "paper.pdf",
    area=[126, 149, 212, 462],
    pages='1'
)

# Multiple areas
dfs = tabula.read_pdf(
    "paper.pdf",
    area=[[126, 149, 212, 462], [300, 149, 400, 462]],
    pages='1'
)
```

**Column Configuration**:
```python
# Specify column positions
dfs = tabula.read_pdf(
    "paper.pdf",
    columns=[100, 200, 300, 400],
    stream=True
)
```

**Batch Processing**:
```python
# Convert all pages to JSON
tabula.convert_into("paper.pdf", "output.json", output_format="json", pages='all')

# Convert to CSV
tabula.convert_into("paper.pdf", "output.csv", output_format="csv", pages='all')
```

---

### pdfplumber

**Purpose**: Fine-grained PDF analysis with precise table extraction control.

**Installation**:
```bash
pip install pdfplumber
```

**Basic Table Extraction**:
```python
import pdfplumber
import pandas as pd

with pdfplumber.open("paper.pdf") as pdf:
    for page in pdf.pages:
        tables = page.extract_tables()
        for table in tables:
            # Convert to DataFrame
            df = pd.DataFrame(table[1:], columns=table[0])
            print(df)
```

**Find Tables**:
```python
with pdfplumber.open("paper.pdf") as pdf:
    page = pdf.pages[0]

    # Find table regions
    tables = page.find_tables()
    for table in tables:
        print(f"Table at: {table.bbox}")  # (x0, top, x1, bottom)

        # Extract table data
        data = table.extract()
        print(data)
```

**Custom Table Settings**:
```python
table_settings = {
    # Line detection strategy
    "vertical_strategy": "lines",      # "lines", "text", "explicit"
    "horizontal_strategy": "lines",    # "lines", "text", "explicit"

    # Explicit lines (when strategy is "explicit")
    "explicit_vertical_lines": [],
    "explicit_horizontal_lines": [],

    # Tolerances
    "snap_tolerance": 3,               # Snap lines to nearby points
    "snap_x_tolerance": 3,
    "snap_y_tolerance": 3,
    "join_tolerance": 3,               # Join line segments
    "join_x_tolerance": 3,
    "join_y_tolerance": 3,

    # Edge detection
    "edge_min_length": 3,

    # Text-based detection
    "min_words_vertical": 3,
    "min_words_horizontal": 1,

    # Intersection tolerance
    "intersection_tolerance": 3,
    "intersection_x_tolerance": 3,
    "intersection_y_tolerance": 3,

    # Text extraction within cells
    "text_tolerance": 3,
    "text_x_tolerance": 3,
    "text_y_tolerance": 3,
}

with pdfplumber.open("paper.pdf") as pdf:
    page = pdf.pages[0]
    table = page.extract_table(table_settings)
```

**Text-Based Detection** (for borderless tables):
```python
table_settings = {
    "vertical_strategy": "text",
    "horizontal_strategy": "text",
    "min_words_vertical": 3,
    "min_words_horizontal": 1,
}

with pdfplumber.open("paper.pdf") as pdf:
    page = pdf.pages[0]
    tables = page.extract_tables(table_settings)
```

**Visual Debugging**:
```python
with pdfplumber.open("paper.pdf") as pdf:
    page = pdf.pages[0]

    # Create image with overlays
    im = page.to_image(resolution=150)

    # Draw detected tables
    for table in page.find_tables():
        im.draw_rect(table.bbox, stroke="red", stroke_width=2)

    # Draw detected lines
    im.draw_lines(page.lines, stroke="blue")

    # Draw characters
    im.draw_rects(page.chars, stroke="green", stroke_width=0.5)

    # Save debug image
    im.save("debug_page.png")
```

---

### Table Transformer (TATR)

**Purpose**: Deep learning model for table detection and structure recognition.

**Installation**:
```bash
pip install transformers torch
pip install timm  # For vision models
```

**Table Detection**:
```python
from transformers import AutoImageProcessor, TableTransformerForObjectDetection
from PIL import Image
import torch

# Load model
processor = AutoImageProcessor.from_pretrained("microsoft/table-transformer-detection")
model = TableTransformerForObjectDetection.from_pretrained("microsoft/table-transformer-detection")

# Load image
image = Image.open("page.png").convert("RGB")

# Process
inputs = processor(images=image, return_tensors="pt")
outputs = model(**inputs)

# Get detections
target_sizes = torch.tensor([image.size[::-1]])
results = processor.post_process_object_detection(outputs, threshold=0.9, target_sizes=target_sizes)[0]

for score, label, box in zip(results["scores"], results["labels"], results["boxes"]):
    print(f"Detected table at {box.tolist()} with confidence {score:.3f}")
```

**Table Structure Recognition**:
```python
# Load structure recognition model
structure_processor = AutoImageProcessor.from_pretrained("microsoft/table-transformer-structure-recognition")
structure_model = TableTransformerForObjectDetection.from_pretrained("microsoft/table-transformer-structure-recognition")

# Process cropped table image
table_image = image.crop(box.tolist())
inputs = structure_processor(images=table_image, return_tensors="pt")
outputs = structure_model(**inputs)

# Post-process to get rows, columns, cells
results = structure_processor.post_process_object_detection(
    outputs, threshold=0.7, target_sizes=torch.tensor([table_image.size[::-1]])
)[0]

# Labels: 0=table, 1=column, 2=row, 3=column header, 4=projected row header, 5=spanning cell
```

---

## Comparison Table

| Feature | Camelot | Tabula | pdfplumber | TATR |
|---------|---------|--------|------------|------|
| Bordered tables | Excellent | Good | Excellent | Good |
| Borderless tables | Good | Good | Fair | Good |
| Complex structures | Fair | Fair | Excellent | Excellent |
| Speed | Fast | Fast | Medium | Slow |
| Accuracy tuning | High | Medium | Very High | Low |
| Image tables | No | No | No | Yes |
| Setup | Easy | Java | Easy | ML models |

## Best Practices

### 1. Identify Table Type First

```python
def identify_table_type(pdf_path, page_num=0):
    """Check if tables are bordered or borderless."""
    import pdfplumber

    with pdfplumber.open(pdf_path) as pdf:
        page = pdf.pages[page_num]

        # Check for table lines
        lines = page.lines
        rects = page.rects

        if len(lines) > 10 or len(rects) > 5:
            return "bordered"
        else:
            return "borderless"
```

### 2. Try Multiple Tools

```python
def extract_table_multi_method(pdf_path, page_num=1):
    """Try multiple extraction methods and compare."""
    results = {}

    # Camelot lattice
    try:
        tables = camelot.read_pdf(pdf_path, pages=str(page_num), flavor='lattice')
        if tables:
            results['camelot_lattice'] = tables[0].df
    except Exception as e:
        results['camelot_lattice'] = f"Error: {e}"

    # Camelot stream
    try:
        tables = camelot.read_pdf(pdf_path, pages=str(page_num), flavor='stream')
        if tables:
            results['camelot_stream'] = tables[0].df
    except Exception as e:
        results['camelot_stream'] = f"Error: {e}"

    # Tabula
    try:
        dfs = tabula.read_pdf(pdf_path, pages=page_num)
        if dfs:
            results['tabula'] = dfs[0]
    except Exception as e:
        results['tabula'] = f"Error: {e}"

    # pdfplumber
    try:
        with pdfplumber.open(pdf_path) as pdf:
            tables = pdf.pages[page_num-1].extract_tables()
            if tables:
                results['pdfplumber'] = pd.DataFrame(tables[0][1:], columns=tables[0][0])
    except Exception as e:
        results['pdfplumber'] = f"Error: {e}"

    return results
```

### 3. Handle Common Issues

```python
def clean_extracted_table(df):
    """Clean common table extraction issues."""
    import pandas as pd

    # Remove empty rows/columns
    df = df.dropna(how='all', axis=0)
    df = df.dropna(how='all', axis=1)

    # Strip whitespace
    df = df.applymap(lambda x: x.strip() if isinstance(x, str) else x)

    # Handle merged headers (if first row seems like header)
    if df.iloc[0].notna().all():
        df.columns = df.iloc[0]
        df = df[1:]

    # Reset index
    df = df.reset_index(drop=True)

    return df
```

### 4. Validate Extraction

```python
def validate_table_extraction(original_df, extracted_df, tolerance=0.1):
    """Compare extracted table against known values."""
    issues = []

    # Check dimensions
    if extracted_df.shape != original_df.shape:
        issues.append(f"Shape mismatch: expected {original_df.shape}, got {extracted_df.shape}")

    # Check values
    for col in original_df.columns:
        if col in extracted_df.columns:
            for idx in original_df.index:
                expected = original_df.loc[idx, col]
                actual = extracted_df.loc[idx, col]
                if expected != actual:
                    issues.append(f"Value mismatch at [{idx}, {col}]: expected '{expected}', got '{actual}'")

    return len(issues) == 0, issues
```

## Common Issues and Solutions

| Issue | Solution |
|-------|----------|
| Merged cells split incorrectly | Use pdfplumber with custom settings |
| Headers not detected | Manually specify header row |
| Multi-line cell content | Join with `\n` or space |
| Special characters corrupted | Check PDF encoding, use UTF-8 |
| Table spans multiple pages | Extract per-page, concatenate |
| Nested tables | Extract outer first, then inner |

## References

- [Camelot Documentation](https://camelot-py.readthedocs.io/)
- [Tabula-py Documentation](https://tabula-py.readthedocs.io/)
- [pdfplumber Documentation](https://github.com/jsvine/pdfplumber)
- [Table Transformer Paper](https://arxiv.org/abs/2110.00061)
- [Microsoft Table Transformer](https://github.com/microsoft/table-transformer)
