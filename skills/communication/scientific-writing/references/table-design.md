# Table Design

Guide to creating effective tables for scientific papers.

## Tables vs. Figures

### When to Use Tables

**Use tables when:**
- Precise values are important
- Many variables to present
- Comparing exact numbers across conditions
- Presenting statistics with multiple parameters
- Categorical data with multiple attributes

### When to Use Figures

**Use figures when:**
- Showing trends or patterns
- Visual comparison more important than exact values
- Illustrating distributions
- Showing relationships between variables

## Table Structure

### Headers

**Column headers:**
- Descriptive names
- Units in header (not repeated in cells)
- Clear abbreviations or spell out
- Bold formatting typical

**Row headers:**
- First column typically identifiers
- Left-aligned
- Can be bold

**Example:**
```
Parameter          | Control (n=10)  | Treatment (n=10) | p-value
                   | Mean ± SD       | Mean ± SD        |
-------------------+-----------------+------------------+---------
Cell viability (%) | 62 ± 9          | 85 ± 7           | 0.003
Apoptosis (%)      | 23 ± 5          | 12 ± 3           | 0.001
```

### Content Organization

**Logical ordering:**
- Alphabetical (genes, proteins)
- By magnitude (largest to smallest)
- By importance (most significant first)
- By category (group related items)

**Consistent formatting:**
- Same decimal places within column
- Aligned by decimal point
- Consistent use of symbols

### Footnotes

**Use for:**
- Abbreviation definitions
- Statistical test details
- Special cases or exceptions
- Symbols explanation

**Notation:**
- Lowercase superscript letters: ᵃ, ᵇ, ᶜ
- Or symbols: *, †, ‡, §
- Numbers if not confused with references

## Formatting Guidelines

### Lines

**Horizontal lines (typical):**
- Top of table
- Below headers
- Bottom of table
- Occasionally between major sections

**Vertical lines:**
- Generally avoid (cleaner look)
- Use spacing instead to separate columns

**Example structure:**
```
═══════════════════════════════════════
Parameter    | Group A | Group B | p-value
─────────────┼─────────┼─────────┼────────
Value 1      | 25.3    | 32.1    | 0.003
Value 2      | 18.7    | 21.3    | 0.045
═══════════════════════════════════════
```

### Alignment

**Text:** Left-aligned
**Numbers:** Right-aligned or decimal-aligned
**Headers:** Centered or left-aligned
**Units:** Aligned with data

### Significant Figures

**Be consistent:**
- Same decimal places in a column
- Appropriate precision (not false precision)
- Match to measurement precision

**Examples:**
- 12.3 ± 1.5 (not 12.345 ± 1.482)
- p = 0.003 (not p = 0.00327845)

## Caption Writing

### Caption Components

1. **Table number and title**
2. **Description of contents**
3. **Abbreviation definitions**
4. **Statistical information**
5. **Footnote explanations**

### Examples

**Good caption:**
"Table 1. Patient demographics and baseline characteristics. Data are presented as mean ± SD for continuous variables and n (%) for categorical variables. BMI, body mass index; SBP, systolic blood pressure; DBP, diastolic blood pressure. *p < 0.05, **p < 0.01 vs. control by Student's t-test or χ² test as appropriate"

**Poor caption:**
"Table 1. Demographics"
(Too brief, missing details)

## Table Types

### Summary Tables

**Characteristics or descriptive statistics:**

| Variable | Control (n=50) | Treatment (n=50) | p-value |
|----------|----------------|------------------|---------|
| Age (years) | 58 ± 12 | 61 ± 10 | 0.18 |
| Male, n (%) | 27 (54) | 31 (62) | 0.40 |
| BMI (kg/m²) | 28.3 ± 4.2 | 27.8 ± 3.9 | 0.52 |

### Results Tables

**Experimental outcomes:**

| Treatment | Viability (%) | Apoptosis (%) | Expression (fold) |
|-----------|---------------|---------------|-------------------|
| Control | 62 ± 9 | 23 ± 5 | 1.0 ± 0.2 |
| Drug A | 78 ± 6* | 15 ± 4* | 2.3 ± 0.4** |
| Drug B | 85 ± 7** | 9 ± 3** | 3.1 ± 0.5** |

*p < 0.05, **p < 0.01 vs. control

### Regression Tables

**Model coefficients:**

| Predictor | β | 95% CI | p-value |
|-----------|---|--------|---------|
| Age | 0.45 | [0.23, 0.67] | <0.001 |
| BMI | 0.23 | [0.08, 0.38] | 0.003 |
| Sex (male) | -0.15 | [-0.32, 0.02] | 0.08 |

R² = 0.42, F(3,96) = 23.4, p < 0.001

### Comparison Tables

**Multiple methods/studies:**

| Method | Accuracy (%) | Speed (ms) | Memory (MB) |
|--------|--------------|------------|-------------|
| Method A [12] | 85 | 120 | 256 |
| Method B [13] | 89 | 180 | 512 |
| Our method | 94 | 95 | 384 |

## Common Table Problems

### Problem: Too Much Information

**Symptom:** Table spans multiple pages, hard to read

**Fix:**
- Split into multiple tables
- Move less important data to supplementary
- Consider if figure would be better
- Group related data into sections

### Problem: Inconsistent Significant Figures

**Symptom:** 12.3, 8.47, 15.123 in same column

**Fix:** Use same decimal places: 12.30, 8.47, 15.12

### Problem: Units Repeated in Every Cell

**Wrong:**
| Sample | Concentration |
|--------|---------------|
| A | 10 μM |
| B | 15 μM |
| C | 20 μM |

**Right:**
| Sample | Concentration (μM) |
|--------|-------------|
| A | 10 |
| B | 15 |
| C | 20 |

### Problem: Unclear Abbreviations

**Fix:** Define all abbreviations in caption, even if defined in text

### Problem: No Statistical Information

**Fix:** Include test used, p-values, sample sizes, and measures of variability

### Problem: Poor Organization

**Fix:** Order logically (by magnitude, category, or importance)

## Table Checklist

- [ ] Clear, descriptive title
- [ ] Units in headers, not cells
- [ ] Consistent decimal places within columns
- [ ] Appropriate significant figures
- [ ] All abbreviations defined in caption
- [ ] Statistical information included (n, test used, p-values)
- [ ] Footnotes for additional information
- [ ] Logical organization
- [ ] Horizontal lines only (typically)
- [ ] Numbers aligned (right or decimal)
- [ ] Text left-aligned
- [ ] Not redundant with figures or text
- [ ] Self-contained caption

## Formatting for Different Journals

**LaTeX:**
```latex
\begin{table}
\caption{Patient demographics}
\label{tab:demographics}
\begin{tabular}{lcc}
\hline
Variable & Control & Treatment \\
\hline
Age (years) & 58 ± 12 & 61 ± 10 \\
BMI (kg/m²) & 28.3 ± 4.2 & 27.8 ± 3.9 \\
\hline
\end{tabular}
\end{table}
```

**Word:**
- Use built-in table styles
- Adjust borders (usually top/bottom/below header only)
- Use table properties for consistent formatting

**Excel/CSV:**
- Can use scripts/table_formatter.py to convert to journal format
- Export to desired format

## Supplementary Tables

**Use for:**
- Detailed data supporting main figures/tables
- Complete datasets
- Extended statistical analyses
- Comprehensive parameter lists

**Format:**
- Can be less strict on formatting
- Often provided as Excel files
- Still need captions and definitions
- Number sequentially (Table S1, S2, etc.)

## Accessibility

**Make tables accessible:**
- Simple structure (avoid merged cells when possible)
- Clear headers
- Alt text describing content
- Screen reader compatible
- Avoid color as only means of distinction
