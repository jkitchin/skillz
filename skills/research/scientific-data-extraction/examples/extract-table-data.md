# Example: Table Data Extraction Comparison

This example compares different table extraction methods and demonstrates best practices for handling various table types.

## Scenario

Extract data from a scientific paper containing:
- Table 1: Simple bordered table with material properties
- Table 2: Borderless table with reaction conditions
- Table 3: Complex table with merged cells and footnotes

## Setup

```python
import camelot
import tabula
import pdfplumber
import pandas as pd
from pathlib import Path

pdf_path = "paper_with_tables.pdf"
```

## Method 1: Camelot Lattice (Bordered Tables)

Best for tables with visible grid lines.

```python
# Extract bordered tables
tables_lattice = camelot.read_pdf(pdf_path, flavor='lattice', pages='all')

print(f"Found {len(tables_lattice)} tables with lattice method")

for i, table in enumerate(tables_lattice):
    print(f"\nTable {i+1}:")
    print(f"  Page: {table.page}")
    print(f"  Shape: {table.df.shape}")
    print(f"  Accuracy: {table.accuracy:.1f}%")

    # Show preview
    print(table.df.head())

# Extract specific table with tuning
table = camelot.read_pdf(
    pdf_path,
    pages='3',
    flavor='lattice',
    line_scale=40,           # Adjust line detection sensitivity
    copy_text=['v'],         # Copy text across vertical spans
    shift_text=[''],         # Shift text to cell centers
)[0]
```

**Tuning Parameters**:
```python
# For tables with thin lines
table = camelot.read_pdf(pdf_path, flavor='lattice', line_scale=60)

# For tables with thick lines
table = camelot.read_pdf(pdf_path, flavor='lattice', line_scale=20)

# Specify table area (x1, y1, x2, y2 from top-left)
table = camelot.read_pdf(
    pdf_path,
    flavor='lattice',
    table_areas=['100,700,500,400']  # Crop to specific region
)
```

## Method 2: Camelot Stream (Borderless Tables)

Best for tables separated by whitespace.

```python
# Extract borderless tables
tables_stream = camelot.read_pdf(pdf_path, flavor='stream', pages='all')

print(f"Found {len(tables_stream)} tables with stream method")

# Tune stream parameters
tables_tuned = camelot.read_pdf(
    pdf_path,
    flavor='stream',
    row_tol=10,          # Tolerance for row grouping (pixels)
    edge_tol=50,         # Edge tolerance for cell boundaries
    columns=['100,200,300,400'],  # Explicit column positions
)

for table in tables_tuned:
    print(f"Accuracy: {table.accuracy:.1f}%")
    print(table.df)
```

**Common Stream Issues and Fixes**:
```python
# Problem: Columns not detected properly
# Solution: Specify column positions manually
table = camelot.read_pdf(
    pdf_path,
    flavor='stream',
    columns=['72,150,250,350,450,550']  # X positions of column separators
)

# Problem: Rows merged incorrectly
# Solution: Adjust row tolerance
table = camelot.read_pdf(
    pdf_path,
    flavor='stream',
    row_tol=5  # Smaller value = stricter row separation
)
```

## Method 3: Tabula

Good general-purpose extraction with Java backend.

```python
# Basic extraction
dfs = tabula.read_pdf(pdf_path, pages='all')

print(f"Found {len(dfs)} tables")
for i, df in enumerate(dfs):
    print(f"\nTable {i+1}: {df.shape}")
    print(df.head())

# With specific options
dfs = tabula.read_pdf(
    pdf_path,
    pages='2',
    lattice=True,        # Use lattice mode
    pandas_options={'header': None}  # Don't use first row as header
)

# Specify area
dfs = tabula.read_pdf(
    pdf_path,
    area=[126, 149, 212, 462],  # top, left, height, width
    pages='1'
)

# Multiple areas
dfs = tabula.read_pdf(
    pdf_path,
    multiple_tables=True,
    area=[[126, 149, 212, 462], [300, 149, 400, 462]]
)
```

## Method 4: pdfplumber (Complex Tables)

Best for fine-grained control and debugging.

```python
with pdfplumber.open(pdf_path) as pdf:
    for page_num, page in enumerate(pdf.pages):
        # Find all tables on page
        tables = page.find_tables()

        for table in tables:
            print(f"Table found at: {table.bbox}")

            # Extract with default settings
            data = table.extract()
            df = pd.DataFrame(data[1:], columns=data[0])
            print(df)
```

**Custom Settings for Complex Tables**:
```python
# For tables without clear lines
custom_settings = {
    "vertical_strategy": "text",
    "horizontal_strategy": "text",
    "snap_tolerance": 5,
    "snap_x_tolerance": 5,
    "snap_y_tolerance": 5,
    "join_tolerance": 3,
    "edge_min_length": 10,
    "min_words_vertical": 3,
    "min_words_horizontal": 1,
    "intersection_tolerance": 5,
    "text_tolerance": 3,
}

with pdfplumber.open(pdf_path) as pdf:
    page = pdf.pages[0]
    table = page.extract_table(custom_settings)
    df = pd.DataFrame(table[1:], columns=table[0])
```

**Visual Debugging**:
```python
with pdfplumber.open(pdf_path) as pdf:
    page = pdf.pages[0]

    # Create debug image
    im = page.to_image(resolution=150)

    # Draw detected elements
    im.draw_rects(page.chars, stroke="green", stroke_width=0.3)
    im.draw_lines(page.lines, stroke="red", stroke_width=1)

    # Draw table bounding boxes
    for table in page.find_tables():
        im.draw_rect(table.bbox, stroke="blue", stroke_width=2)

    im.save("debug_tables.png")
```

## Comparison Function

```python
def compare_extraction_methods(pdf_path: str, page: int = 1) -> dict:
    """Compare all extraction methods on the same page."""
    results = {}

    # Camelot lattice
    try:
        tables = camelot.read_pdf(pdf_path, pages=str(page), flavor='lattice')
        if tables:
            results['camelot_lattice'] = {
                'count': len(tables),
                'accuracy': tables[0].accuracy,
                'shape': tables[0].df.shape,
                'df': tables[0].df
            }
    except Exception as e:
        results['camelot_lattice'] = {'error': str(e)}

    # Camelot stream
    try:
        tables = camelot.read_pdf(pdf_path, pages=str(page), flavor='stream')
        if tables:
            results['camelot_stream'] = {
                'count': len(tables),
                'accuracy': tables[0].accuracy,
                'shape': tables[0].df.shape,
                'df': tables[0].df
            }
    except Exception as e:
        results['camelot_stream'] = {'error': str(e)}

    # Tabula
    try:
        dfs = tabula.read_pdf(pdf_path, pages=page)
        if dfs:
            results['tabula'] = {
                'count': len(dfs),
                'shape': dfs[0].shape,
                'df': dfs[0]
            }
    except Exception as e:
        results['tabula'] = {'error': str(e)}

    # pdfplumber
    try:
        with pdfplumber.open(pdf_path) as pdf:
            tables = pdf.pages[page-1].extract_tables()
            if tables:
                df = pd.DataFrame(tables[0][1:], columns=tables[0][0])
                results['pdfplumber'] = {
                    'count': len(tables),
                    'shape': df.shape,
                    'df': df
                }
    except Exception as e:
        results['pdfplumber'] = {'error': str(e)}

    return results

# Run comparison
results = compare_extraction_methods("paper.pdf", page=2)

# Print summary
print("\n=== Extraction Method Comparison ===\n")
for method, data in results.items():
    print(f"{method}:")
    if 'error' in data:
        print(f"  Error: {data['error']}")
    else:
        print(f"  Tables found: {data['count']}")
        print(f"  First table shape: {data['shape']}")
        if 'accuracy' in data:
            print(f"  Accuracy: {data['accuracy']:.1f}%")
    print()
```

## Handling Common Table Types

### Type 1: Simple Data Table

```
| Material | Bandgap (eV) | Crystal Structure |
|----------|--------------|-------------------|
| TiO2     | 3.2          | anatase           |
| ZnO      | 3.37         | wurtzite          |
```

```python
# Any method works well
tables = camelot.read_pdf(pdf_path, flavor='lattice')
df = tables[0].df

# Clean up
df.columns = df.iloc[0]  # First row as header
df = df[1:]              # Remove header row from data
df = df.reset_index(drop=True)
```

### Type 2: Multi-Header Table

```
|           | Property A  | Property B  |
| Material  | Value | Err | Value | Err |
|-----------|-------|-----|-------|-----|
| Sample 1  | 1.23  | 0.1 | 4.56  | 0.2 |
```

```python
# Extract raw
tables = camelot.read_pdf(pdf_path, flavor='lattice')
df = tables[0].df

# Handle multi-level headers
header_row_1 = df.iloc[0].tolist()
header_row_2 = df.iloc[1].tolist()

# Create combined headers
combined_headers = []
current_header = ""
for h1, h2 in zip(header_row_1, header_row_2):
    if h1.strip():
        current_header = h1.strip()
    combined_headers.append(f"{current_header}_{h2}".strip('_'))

df.columns = combined_headers
df = df[2:]  # Remove header rows
```

### Type 3: Table with Merged Cells

```python
# Use pdfplumber for better merged cell handling
with pdfplumber.open(pdf_path) as pdf:
    page = pdf.pages[0]

    # Get table with explicit settings
    table_settings = {
        "vertical_strategy": "lines",
        "horizontal_strategy": "lines",
        "snap_tolerance": 5,
    }

    tables = page.extract_tables(table_settings)
    df = pd.DataFrame(tables[0])

    # Forward-fill merged cells
    df = df.fillna(method='ffill', axis=0)  # Fill down
    df = df.fillna(method='ffill', axis=1)  # Fill right
```

### Type 4: Table with Footnotes

```python
def extract_table_with_footnotes(pdf_path: str, page: int) -> tuple:
    """Extract table and associated footnotes."""
    with pdfplumber.open(pdf_path) as pdf:
        p = pdf.pages[page - 1]

        # Find table
        tables = p.find_tables()
        if not tables:
            return None, []

        table = tables[0]
        table_bottom = table.bbox[3]

        # Extract table data
        df = pd.DataFrame(table.extract()[1:], columns=table.extract()[0])

        # Find footnotes below table
        footnote_area = (0, table_bottom, p.width, p.height)
        footnote_chars = p.within_bbox(footnote_area).extract_text()

        # Parse footnotes (common patterns)
        import re
        footnotes = re.findall(r'([a-z*†‡§])\s*(.+?)(?=[a-z*†‡§]|$)', footnote_chars)

        return df, footnotes

df, footnotes = extract_table_with_footnotes("paper.pdf", 3)
print("Footnotes:", footnotes)
```

## Post-Processing Functions

```python
def clean_table(df: pd.DataFrame) -> pd.DataFrame:
    """Clean extracted table data."""
    # Remove empty rows and columns
    df = df.dropna(how='all', axis=0)
    df = df.dropna(how='all', axis=1)

    # Strip whitespace from all cells
    df = df.applymap(lambda x: x.strip() if isinstance(x, str) else x)

    # Remove newlines within cells
    df = df.applymap(lambda x: x.replace('\n', ' ') if isinstance(x, str) else x)

    # Reset index
    df = df.reset_index(drop=True)

    return df

def parse_numeric_columns(df: pd.DataFrame) -> pd.DataFrame:
    """Convert numeric columns to appropriate types."""
    import re

    for col in df.columns:
        # Check if column looks numeric
        sample = df[col].dropna().head(5).tolist()
        if all(re.match(r'^-?\d+\.?\d*$', str(s).strip()) for s in sample):
            df[col] = pd.to_numeric(df[col], errors='coerce')

    return df

def extract_units_from_headers(df: pd.DataFrame) -> dict:
    """Extract units from column headers."""
    import re

    units = {}
    new_columns = []

    for col in df.columns:
        # Look for patterns like "Property (unit)" or "Property [unit]"
        match = re.search(r'(.+?)\s*[\(\[](.+?)[\)\]]', str(col))
        if match:
            name, unit = match.groups()
            units[name.strip()] = unit.strip()
            new_columns.append(name.strip())
        else:
            new_columns.append(col)

    df.columns = new_columns
    return units
```

## Complete Example

```python
def extract_all_tables(pdf_path: str) -> list:
    """Extract and clean all tables from PDF."""
    results = []

    # Try Camelot first (usually best for bordered tables)
    try:
        tables = camelot.read_pdf(pdf_path, flavor='lattice', pages='all')
        for i, table in enumerate(tables):
            if table.accuracy > 70:
                df = clean_table(table.df)
                units = extract_units_from_headers(df)
                df = parse_numeric_columns(df)

                results.append({
                    'method': 'camelot_lattice',
                    'page': table.page,
                    'accuracy': table.accuracy,
                    'data': df,
                    'units': units
                })
    except Exception as e:
        print(f"Camelot error: {e}")

    # Fallback to pdfplumber for missed tables
    with pdfplumber.open(pdf_path) as pdf:
        for page_num, page in enumerate(pdf.pages):
            tables = page.find_tables()
            for j, table in enumerate(tables):
                # Check if already extracted by Camelot
                bbox = table.bbox
                already_found = any(
                    r['page'] == page_num + 1
                    for r in results
                )

                if not already_found:
                    raw = table.extract()
                    if raw and len(raw) > 1:
                        df = pd.DataFrame(raw[1:], columns=raw[0])
                        df = clean_table(df)

                        results.append({
                            'method': 'pdfplumber',
                            'page': page_num + 1,
                            'data': df
                        })

    return results

# Run extraction
all_tables = extract_all_tables("paper.pdf")

# Save results
for i, table_info in enumerate(all_tables):
    output_file = f"table_{i+1}_page{table_info['page']}.csv"
    table_info['data'].to_csv(output_file, index=False)
    print(f"Saved {output_file}: {table_info['data'].shape}")
```

## Decision Guide

| Table Type | Recommended Method | Notes |
|------------|-------------------|-------|
| Bordered, simple | Camelot lattice | Fastest, most accurate |
| Bordered, complex | pdfplumber | Better merged cell handling |
| Borderless | Camelot stream | Needs tuning |
| Image of table | Table Transformer | Deep learning approach |
| Multi-page table | Manual + concatenate | Extract per-page, join |
