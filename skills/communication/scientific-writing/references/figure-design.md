# Figure Design

Comprehensive guide to creating publication-quality figures for research papers.

## Figure Types by Data Type

### Line Plots
**Use for:** Trends over continuous variable (time, dose, etc.)
**Best for:** Time series, dose-response curves, continuous relationships

**Design tips:**
- Clear axis labels with units
- Distinct line styles for multiple series
- Legend or direct labeling
- Error bars/shading for uncertainty
- Grid lines if helpful for reading values

### Bar Plots
**Use for:** Comparing discrete categories
**Best for:** Group comparisons, categorical data

**Design tips:**
- Bars should start at zero (or use break symbol if not)
- Error bars on top
- Group related bars
- Consider horizontal bars for long category names
- Consistent bar width

### Scatter Plots
**Use for:** Showing individual data points, correlations
**Best for:** Distributions, correlations, outlier detection

**Design tips:**
- Semi-transparent points if many overlapping
- Include regression line if showing correlation
- Report R² and p-value
- Consider logarithmic axes for wide ranges

### Box Plots / Violin Plots
**Use for:** Showing distributions with quartiles
**Best for:** Comparing distributions across groups

**Design tips:**
- Show individual points if n is small
- Explain box components in caption
- Consider violin plots for distribution shape
- Consistent width across boxes

### Heatmaps
**Use for:** Matrices, multi-dimensional data
**Best for:** Gene expression, correlation matrices, spatial data

**Design tips:**
- Appropriate color scale (sequential, diverging, qualitative)
- Color bar with limits and units
- Clustering if appropriate
- Row/column labels readable

### Schematics / Diagrams
**Use for:** Concepts, workflows, experimental design
**Best for:** Illustrating methods, models, hypotheses

**Design tips:**
- Simple, uncluttered
- Consistent visual language
- Clear arrows showing flow/relationships
- Labeled components

### Microscopy Images
**Use for:** Showing spatial organization, morphology
**Best for:** Cells, tissues, materials

**Design tips:**
- Scale bar with size
- Channel labels (if multi-color)
- Consistent contrast/brightness across compared images
- Insets for details if needed

## Publication Quality Requirements

### Resolution
- **Minimum 300 DPI** for print publication
- 600 DPI preferred for line art
- 600-1200 DPI for combination art (photos + line art)
- Web figures can be 72-150 DPI

**Check your resolution:**
```python
# If using matplotlib
fig.savefig('figure.png', dpi=300)
fig.savefig('figure.pdf')  # Vector, resolution-independent
```

### File Formats

**Vector formats (preferred when possible):**
- PDF: Best for figures with plots, diagrams
- EPS: Traditional vector format
- SVG: Web-friendly vector format

**Raster formats (for photos/microscopy):**
- TIFF: Uncompressed, high quality
- PNG: Lossless compression
- Avoid JPEG for scientific figures (lossy compression)

**When to use each:**
- Plots, diagrams → PDF or EPS (vector)
- Microscopy, photos → TIFF (high res raster)
- Mixed → TIFF at high resolution
- Supplementary web → PNG

### Size

**Column widths (check journal guidelines):**
- Single column: typically 3-3.5 inches (7.5-9 cm)
- 1.5 column: 4.5-5.5 inches (11-14 cm)
- Full page width: 6.5-7 inches (16.5-18 cm)

**Design at final size:**
- Create figure at size it will appear in print
- Fonts will be readable when reduced
- Lines will be visible when reduced

### Font Size

**Readable when reduced to print size:**
- Axis labels: 8-12 pt
- Tick labels: 6-10 pt
- Legend: 7-10 pt
- **Minimum 6-8 pt after reduction**

**Test readability:**
Print figure at final size and check if text is legible

### Line Thickness

- Thick enough to see when reduced
- Minimum 0.5 pt after reduction
- Main lines: 1-2 pt
- Grid lines: 0.5-1 pt (if used)

### Colors

**Colorblind-friendly palettes:**
- Avoid red-green combinations alone
- Use blue-orange, purple-orange combinations
- Add patterns/shapes in addition to color
- Test with colorblind simulator

**Common colorblind-safe palettes:**
- Okabe-Ito palette (designed for colorblindness)
- Viridis, plasma, cividis (sequential)
- ColorBrewer palettes marked colorblind-safe

**Consider grayscale printing:**
- Will colors be distinguishable in grayscale?
- Use different line styles (solid, dashed, dotted)
- Use different markers (circle, square, triangle)

**Python example:**
```python
# Colorblind-friendly Okabe-Ito palette
colors = ['#E69F00', '#56B4E9', '#009E73', '#F0E442',
          '#0072B2', '#D55E00', '#CC79A7']
```

## Essential Figure Elements

### Axis Labels

**Always include:**
- What is plotted (e.g., "Cell Viability")
- Units in parentheses or after comma

**Examples:**
- "Time (hours)"
- "Expression (fold change)"
- "Temperature, °C"
- "Voltage (mV)"

**Don't:**
- Leave axes unlabeled
- Forget units
- Use unclear abbreviations

### Legends

**Include when:**
- Multiple data series
- Different symbols/colors need explanation

**Design tips:**
- Inside plot area if space allows
- Outside if crowded
- Horizontal if few items
- Vertical if many items
- Direct labeling preferred if possible (label right on plot)

### Error Bars

**Always include uncertainty measure:**
- Standard deviation (SD): shows data spread
- Standard error (SE): shows uncertainty in mean
- Confidence interval (CI): shows likely range for true value

**Must specify which in caption:**
"Error bars represent mean ± SD (n=3)"
"Error bars represent 95% confidence intervals"

**Visual design:**
- Caps on error bars helpful
- Not so thick they obscure data
- Consider shaded regions for continuous data

### Statistical Annotations

**Show significance on figures:**
- Asterisks: * p < 0.05, ** p < 0.01, *** p < 0.001
- Exact p-values if space allows
- Brackets showing what's being compared
- ns = not significant

**Example:**
```
     ***
  ┌─────┐
  │     │
 ─┘     └─
Group1 Group2
```

### Scale Bars (for images)

**Required for microscopy/images:**
- White or black bar
- Size labeled (e.g., "10 μm")
- Placed in corner, not obscuring features
- Same scale for images being compared (or clearly state different)

### Panel Labels

**Multi-panel figures:**
- Label A, B, C or a, b, c
- Consistent size and placement
- Top-left corner typical
- Bold, large enough to see: 10-14 pt

**Example:**
```
A                    B
[Image 1]           [Graph 1]

C                    D
[Graph 2]           [Image 2]
```

## Caption Writing

### Caption Structure

**Standard format:**
1. Opening sentence: What figure shows overall
2. Panel descriptions: Detail each panel (A, B, C...)
3. Statistical information: n, test used, significance
4. Symbol/color definitions: What symbols/colors mean
5. Additional details: Methods specifics if not in text

### Examples

**Good caption:**
"Figure 1. Treatment increases cell viability dose-dependently. (A) Representative microscopy images of control and treated cells. Scale bar = 20 μm. (B) Quantification of cell viability at different treatment concentrations. Data are mean ± SD from three independent experiments (n=3). *p < 0.05, **p < 0.01 vs. control by one-way ANOVA with Tukey's post-hoc test"

**Poor caption:**
"Figure 1. Cell viability"
(Too vague, missing details)

### Captions Should Be Self-Contained

**Reader should understand figure without reading main text**

**Include:**
- What was measured
- What conditions/groups
- What statistics represent (SD, SE, CI)
- What statistical test was used
- Sample size (n)
- Meaning of symbols/colors

**Don't:**
- Repeat results text word-for-word
- Include extensive interpretation (that's for main text)

## Common Figure Problems

### Problem: Fonts Too Small

**Symptom:** After reducing to column width, text unreadable

**Fix:** Design figure at final print size, use minimum 8 pt fonts

### Problem: Insufficient Contrast

**Symptom:** Light colors on white background hard to see

**Fix:** Use darker colors, higher contrast, or add background shading

### Problem: Missing Error Bars

**Symptom:** Bar/point plots without uncertainty measures

**Fix:** Always include error bars; state what they represent in caption

### Problem: No Units on Axes

**Symptom:** Axis says "Temperature" but doesn't say °C or K

**Fix:** Always include units: "Temperature (°C)"

### Problem: Unclear Legends

**Symptom:** Symbols not explained, or legend abbreviations undefined

**Fix:** Define all symbols in legend or caption; spell out abbreviations

### Problem: Low Resolution

**Symptom:** Pixelated when viewed at final size

**Fix:** Save at minimum 300 DPI; use vector formats when possible

### Problem: Inappropriate Plot Type

**Symptom:** Bar plot used for continuous data, line plot for categorical

**Fix:** Match plot type to data type (see section above)

### Problem: Misleading Axes

**Symptom:** Y-axis doesn't start at zero, exaggerating small differences

**Fix:** Start at zero for bar plots, or use axis break symbol; for other plots, choose range that honestly represents data

### Problem: Too Much Information

**Symptom:** Figure cluttered, hard to interpret

**Fix:** Split into multiple panels or figures; simplify; focus on main message

## Creating Figures with matplotlib (Python)

```python
import matplotlib.pyplot as plt
import numpy as np

# Publication-quality defaults
plt.rcParams['figure.dpi'] = 300
plt.rcParams['font.size'] = 10
plt.rcParams['font.family'] = 'Arial'
plt.rcParams['axes.linewidth'] = 1
plt.rcParams['lines.linewidth'] = 1.5

# Create figure at final size (single column = 3.5 inches)
fig, ax = plt.subplots(figsize=(3.5, 2.5))

# Example data
x = np.linspace(0, 10, 50)
y1 = np.sin(x)
y2 = np.cos(x)
err = 0.1

# Plot with error shading
ax.plot(x, y1, label='Condition A', color='#0072B2')
ax.fill_between(x, y1-err, y1+err, alpha=0.3, color='#0072B2')

ax.plot(x, y2, label='Condition B', color='#D55E00')
ax.fill_between(x, y2-err, y2+err, alpha=0.3, color='#D55E00')

# Labels with units
ax.set_xlabel('Time (hours)')
ax.set_ylabel('Expression (fold change)')

# Legend
ax.legend(frameon=False)

# Remove top and right spines for cleaner look
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

# Save as vector and raster
fig.savefig('figure1.pdf', bbox_inches='tight')
fig.savefig('figure1.png', dpi=300, bbox_inches='tight')
plt.close()
```

## Multi-Panel Figures

**When to use:**
- Related data that should be viewed together
- Sequential experiments
- Different views of same phenomenon
- Comparison of conditions

**Layout tips:**
- Logical flow (left to right, top to bottom)
- Consistent size for similar panels
- Align axes where possible
- Panel labels (A, B, C) in consistent position

**Python example:**
```python
fig, axes = plt.subplots(2, 2, figsize=(7, 6))

# Label panels
labels = ['A', 'B', 'C', 'D']
for ax, label in zip(axes.flat, labels):
    ax.text(-0.1, 1.1, label, transform=ax.transAxes,
            fontsize=12, fontweight='bold', va='top')

    # Plot in each panel
    ax.plot(x, y)
    ax.set_xlabel('X Label')
    ax.set_ylabel('Y Label')

plt.tight_layout()
fig.savefig('multipanel.pdf')
```

## Color Schemes

### Sequential (for ordered data, low to high)
- Single hue, increasing intensity
- Use: heatmaps, gradients
- Examples: Blues, Greens, viridis, plasma

### Diverging (for data with meaningful midpoint)
- Two hues, diverging from neutral center
- Use: fold changes, correlations
- Examples: RdBu (red-blue), PiYG (pink-green)

### Qualitative (for categorical data)
- Distinct hues for different categories
- Use: different conditions, groups
- Examples: Set1, Set2, tab10, Okabe-Ito

## Journal-Specific Requirements

**Always check journal guidelines for:**
- Figure file formats accepted
- Maximum figure dimensions
- DPI requirements
- Color charges (some journals charge for color in print)
- Number of figures allowed
- Supplementary figure policies

**Common journal families:**
- Nature: PDF/EPS for line art, TIFF for photos, 300-600 DPI
- Science: PDF/EPS/PS vector, TIFF raster, 300 DPI minimum
- Cell: EPS/PDF vector, TIFF raster, minimum 300 DPI
- PLOS: TIFF/PDF/EPS, 300-600 DPI, size limits
- ACS: PDF/TIFF/EPS, 300-1200 DPI depending on type

## Figure Accessibility

**Make figures accessible to all readers:**

**Colorblind considerations:**
- Use colorblind-safe palettes
- Add patterns/shapes in addition to color
- Include line style variations

**Text alternatives:**
- Descriptive captions
- Data tables as supplementary material
- Alt text for online versions

**Simplicity:**
- Clean, uncluttered design
- High contrast
- Readable fonts

## Figure Checklist

- [ ] Resolution: minimum 300 DPI (or vector format)
- [ ] File format: appropriate for journal (PDF/EPS/TIFF)
- [ ] Size: designed for column/page width
- [ ] Fonts: minimum 6-8 pt after reduction, readable
- [ ] Lines: thick enough to see when reduced
- [ ] Colors: colorblind-friendly, distinguishable in grayscale
- [ ] Axis labels: present with units
- [ ] Legend: included if needed, symbols defined
- [ ] Error bars: present with definition in caption
- [ ] Statistical annotations: significance marked if relevant
- [ ] Scale bars: included for images with size
- [ ] Panel labels: A, B, C if multi-panel
- [ ] Caption: complete and self-contained
- [ ] Plot type: appropriate for data type
- [ ] No chart junk: clean, focused on data

## Tools for Figure Creation

**Python:**
- matplotlib: Standard plotting library
- seaborn: Statistical visualizations
- plotly: Interactive plots
- bokeh: Interactive visualizations

**R:**
- ggplot2: Grammar of graphics plotting
- lattice: Trellis graphics
- plotly: Interactive plots

**Graphical:**
- GraphPad Prism: Statistics and graphing
- Origin: Scientific graphing
- Adobe Illustrator: Final figure assembly and editing
- Inkscape: Free vector graphics editor

**Online:**
- BioRender: Biological schematics
- ChemDraw: Chemical structures
- Plot.ly: Online interactive plots
