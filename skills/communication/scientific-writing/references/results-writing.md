# Results Writing

Comprehensive guide to writing objective, data-driven results sections.

## Core Principle: Objectivity

**Report observations, not interpretations**

Results section answers: "What did you find?"
Discussion section answers: "What does it mean?"

## Objectivity Principles

### Do: Report Observations

✅ "Expression of Gene X increased 3-fold in treated cells (Figure 1A)"
✅ "The data show that treatment increased cell viability"
✅ "We observed a positive correlation between A and B"
✅ "Group 1 showed significantly higher values than Group 2"

### Don't: Interpret or Editorialize

❌ "Expression increased, indicating pathway activation" (interpretation → Discussion)
❌ "These data prove that treatment works" (too strong, interpretation)
❌ "Interestingly, we found..." (editorial comment)
❌ "This suggests that..." (interpretation → Discussion)

### Language Choices

**Observation language (Results):**
- "Data show that..."
- "Results demonstrate that..."
- "We observed..."
- "We found..."
- "X increased/decreased..."
- "There was a correlation between..."

**Interpretation language (Discussion):**
- "Data suggest that..."
- "Results indicate that..."
- "This implies..."
- "A possible explanation is..."
- "These findings support..."
- "We propose that..."

### Negative Results

**Report honestly:**
- "No significant difference was observed between groups (p = 0.23)"
- "Treatment did not affect cell viability (Figure 2B)"
- "We found no correlation between X and Y (r = 0.12, p = 0.45)"

**Don't hide or minimize:**
- Not: "Unfortunately, treatment had no effect"
- Not: Omit negative results entirely (publication bias)

## Organization Strategies

### 1. Chronological (Experimental Sequence)

**When to use:** Later experiments build on earlier results

**Structure:**
```
First experiment/observation
↓
Second experiment building on first
↓
Third experiment extending findings
↓
Final validation
```

**Example:**
"To determine whether X affects Y, we first measured... Having established that X increases Y, we next asked whether... To confirm this finding in vivo, we..."

### 2. By Importance (Most Significant First)

**When to use:** High-impact journals, independent experiments

**Structure:**
```
Most important finding
↓
Supporting finding
↓
Additional observations
↓
Validation or controls
```

**Example:**
"Treatment significantly improved survival (Figure 1). This effect was associated with reduced tumor growth (Figure 2). To determine the mechanism, we measured... Control experiments confirmed..."

### 3. By Theme (Related Results Grouped)

**When to use:** Multiple independent assays or measurements

**Structure:**
```
Effect on cellular phenotype
↓
Effect on gene expression
↓
Effect on protein levels
↓
Effect on metabolis
m
```

**Example with subheadings:**
"**Cell viability and proliferation.** Treatment increased cell viability (Figure 1A)... Proliferation rate was also enhanced (Figure 1B)...

**Gene expression changes.** RNA-seq revealed 453 differentially expressed genes (Figure 2)... Pathway analysis showed enrichment in...

**Protein-level validation.** Western blotting confirmed increased expression of..."

### 4. By Hypothesis (Address Each in Turn)

**When to use:** Hypothesis-driven studies with specific aims

**Structure:**
```
Results addressing Hypothesis 1
↓
Results addressing Hypothesis 2
↓
Results addressing Hypothesis 3
```

**Example:**
"We hypothesized that (1) X regulates Y, (2) this regulation depends on Z, and (3) disruption impairs function.

To test hypothesis 1, we... Results showed... (Figure 1).

To test hypothesis 2, we... We found... (Figure 2).

To test hypothesis 3, we... Data revealed... (Figure 3)."

## Figure and Table Integration

### Every Figure/Table Must Be Referenced

**Reference format:**

**In parentheses (common in biology):**
"Expression increased in treated cells (Figure 1A)"
"Demographics are shown (Table 1)"

**In sentence (common in engineering/CS):**
"Figure 1A shows that expression increased..."
"Table 1 summarizes patient demographics"

**Bad - not referenced:**
Having Figure 2 in your paper but never mentioning it in text

### Describe Key Findings, Don't Just Point to Figures

**Poor (vague):**
"See Figure 1 for results"
"Expression data are in Figure 2"

**Better (describes finding):**
"Expression of Gene X increased 3-fold in treated cells (Figure 1)"
"Treatment significantly reduced tumor volume at all time points (Figure 2, p < 0.05 for all)"

### Order of Mention = Order of Figures

**Good:**
Text mentions: Figure 1, Figure 2, Figure 3...
Figures appear: Figure 1, Figure 2, Figure 3...

**Bad:**
Text mentions: Figure 2, Figure 1, Figure 3...
(Figures should be numbered in order of mention)

### Multi-Panel Figures

**Reference specific panels:**
"Expression increased in treated cells (Figure 1A) but not in untreated controls (Figure 1B)"

**Or describe relationship between panels:**
"Figure 1 shows cell viability (panel A), proliferation (panel B), and apoptosis (panel C) in response to treatment"

## Statistical Reporting Requirements

### Complete Statistical Reporting

**Minimum requirements:**
- Test statistic (t, F, χ², etc.)
- Degrees of freedom
- Exact p-value (not just "p < 0.05")
- Measure of effect (mean, median, etc.)
- Measure of variability (SD, SE, or CI)

### Examples by Test Type

**t-test:**
"Treated cells showed higher viability than controls (mean ± SD: 85 ± 7% vs 62 ± 9%; t₁₈ = 4.32, p = 0.003)"

**ANOVA:**
"Cell viability differed significantly among groups (F₃,₃₆ = 8.45, p < 0.001). Post-hoc analysis revealed..."

**Chi-square:**
"Treatment response rates differed between groups (60% vs 35%; χ² = 8.52, df = 1, p = 0.004)"

**Correlation:**
"X correlated positively with Y (Pearson r = 0.67, p < 0.001, n = 45)"

**Non-parametric:**
"Groups differed significantly (Mann-Whitney U = 234, p = 0.012)"

**Regression:**
"X significantly predicted Y (β = 0.45, 95% CI [0.23, 0.67], t₄₃ = 4.12, p < 0.001, R² = 0.28)"

### Effect Sizes

**Include when applicable:**
- Cohen's d (t-tests)
- η² or ω² (ANOVA)
- Odds ratio or risk ratio (categorical outcomes)
- R² (regression)

**Example:**
"Treatment significantly improved outcomes (mean difference = 15 points, 95% CI [8, 22], t₄₈ = 4.32, p < 0.001, Cohen's d = 1.2)"

### Multiple Comparisons

**State correction method:**
"P-values were adjusted for multiple comparisons using Bonferroni correction"
"Significance was determined using FDR-adjusted p-values (q < 0.05)"

### Variability Measures

**Standard Deviation (SD):**
- Describes spread of data
- Use when showing how variable measurements are
- Example: "Expression was 5.2 ± 2.1 (mean ± SD, n = 20)"

**Standard Error (SE):**
- Describes uncertainty in mean estimate
- Use when emphasizing precision of estimate
- Example: "Mean expression was 5.2 ± 0.5 (mean ± SE, n = 20)"

**Confidence Interval (CI):**
- Range likely to contain true value
- Increasingly preferred
- Example: "Mean expression was 5.2 (95% CI [4.3, 6.1], n = 20)"

**Be explicit about which you're using**

## Tense Usage in Results

**Past tense (standard):**
"Expression increased in treated cells"
"We observed a correlation between X and Y"
"Survival was significantly longer in Group A"

**Present tense (for describing figures/tables):**
"Figure 1 shows..." (figure currently shows, so present tense)
"Table 2 summarizes..." (table currently summarizes)

**Both together:**
"Figure 1A shows that expression increased in treated cells"
(Present for "shows," past for "increased")

## Common Results Section Problems

### Problem: Interpretation Instead of Observation

**Wrong:**
"Expression increased (Figure 1), indicating that the pathway is activated"
(Interpretation belongs in Discussion)

**Right:**
"Expression increased 2.5-fold (Figure 1, p < 0.001)"
(Observation only; save interpretation for Discussion)

### Problem: Missing Statistical Details

**Incomplete:**
"Groups differed significantly"

**Complete:**
"Group A showed significantly higher values than Group B (mean ± SD: 45 ± 8 vs 32 ± 6; t₁₈ = 4.12, p = 0.001)"

### Problem: Figures Not Referenced

**Wrong:**
Mentioning results but not citing corresponding figure
"Cell viability increased with treatment" (but Figure 2 shows this data)

**Right:**
"Cell viability increased with treatment (Figure 2, p < 0.01)"

### Problem: Methods in Results Section

**Wrong:**
"We observed increased expression. To measure this, we used qPCR with the following primers..."
(Methods belong in Methods section)

**Right (Results):**
"Gene expression increased 3-fold (Figure 1B, p < 0.001)"

**Right (Methods):**
"Gene expression was measured by qPCR using primers listed in Table S1"

### Problem: Too Much Figure Description

**Wrong:**
"Figure 1A shows a bar graph with three bars representing control, treatment 1, and treatment 2. The y-axis shows viability from 0 to 100%. Error bars represent standard deviation"
(This belongs in figure caption)

**Right:**
"Cell viability increased with treatment 1 (85 ± 7%) and treatment 2 (92 ± 5%) compared to control (62 ± 9%; Figure 1A, p < 0.01 for both)"

### Problem: Selective Reporting

**Wrong:**
Only reporting experiments that "worked"
Omitting failed replicates
Cherry-picking time points or conditions

**Right:**
Report all experiments conducted
Include negative results
Explain exclusions if necessary ("One animal was excluded due to surgical complications")

## Results Organization Examples

### Example 1: Cell Biology Study

**Organization: By theme**

"**Effects on cell viability and proliferation**

Treatment with Compound X increased cell viability in a dose-dependent manner (Figure 1A). At 10 μM, viability was 85 ± 7% compared to 62 ± 9% in vehicle-treated controls (t₁₈ = 4.32, p = 0.003). Cell proliferation rate also increased (Figure 1B), with doubling time decreasing from 24 ± 3 hours to 18 ± 2 hours (p = 0.01).

**Gene expression changes**

To identify molecular changes underlying these effects, we performed RNA-seq (Figure 2). We identified 453 differentially expressed genes (FDR < 0.05, |log₂FC| > 1; Supplementary Table 1). Pathway analysis revealed enrichment in cell cycle regulation (p = 1.2×10⁻⁸) and DNA replication (p = 3.4×10⁻⁶; Figure 2B).

**Protein-level validation**

Western blotting confirmed increased expression of cell cycle regulators Cyclin D1 (2.3 ± 0.4-fold, p = 0.002) and CDK4 (1.8 ± 0.3-fold, p = 0.01) in treated cells (Figure 3A,B)."

### Example 2: Clinical Study

**Organization: Chronological (enrollment → outcomes)**

"**Patient characteristics**

Between January 2020 and December 2022, 450 patients were screened and 300 enrolled (Figure 1). Baseline characteristics were similar between treatment (n = 150) and control (n = 150) groups (Table 1). Mean age was 58 ± 12 years, 55% were male, and mean HbA1c was 8.2 ± 1.1%.

**Primary outcome**

At 12 months, HbA1c decreased significantly more in the treatment group than control group (mean change: -1.2 ± 0.8% vs -0.4 ± 0.6%; mean difference -0.8%, 95% CI [-1.0, -0.6], p < 0.001; Figure 2A). A significantly greater proportion of treatment group achieved target HbA1c <7% (68% vs 35%; OR 3.9, 95% CI [2.4, 6.4], p < 0.001; Figure 2B).

**Secondary outcomes**

Fasting glucose decreased more in treatment group (Table 2). Lipid parameters improved in treatment group, with reductions in LDL cholesterol (p = 0.03) and triglycerides (p = 0.01) compared to control."

### Example 3: Engineering/CS Study

**Organization: By component or benchmark**

"**System Performance**

We evaluated system throughput across varying workloads (Figure 3). Peak throughput was 45,000 transactions per second at 50% load, compared to 38,000 for baseline system (p < 0.01). Latency remained below 10 ms up to 80% load (Figure 3B).

**Scalability Analysis**

System scaled linearly up to 16 nodes (Figure 4A). Beyond 16 nodes, communication overhead became limiting, with efficiency decreasing to 0.75 at 32 nodes. Memory usage scaled sub-linearly with data size (Figure 4B), consistent with compression effectiveness.

**Comparison to State-of-the-Art**

Table 2 compares our approach to existing methods. Our system achieved 1.7× higher throughput than Method A and 2.3× higher than Method B while maintaining 40% lower memory footprint."

## Present Results Without Interpretation

### Results (Observations Only)

✅ "Expression increased 3-fold (Figure 1, p < 0.001)"
✅ "Treatment did not affect cell viability (p = 0.23)"
✅ "We observed a positive correlation (r = 0.67, p < 0.001)"
✅ "Survival was significantly longer in treated group (median: 24 months vs 18 months, log-rank p = 0.003)"

### Save for Discussion (Interpretations)

❌ "Expression increased, suggesting pathway activation"
❌ "The lack of effect on viability indicates that the treatment works through another mechanism"
❌ "This correlation supports the hypothesis that..."
❌ "The survival benefit demonstrates therapeutic potential"

### Borderline Cases

**"Data show that..." - Generally OK in results if stating observation**
"Data show that treatment increases expression"

**"Results demonstrate that..." - Generally OK**
"Results demonstrate a dose-dependent effect"

**"These findings suggest..." - Generally for discussion**
Save "suggest," "indicate," "imply" for Discussion

## Results Checklist

- [ ] All results are observations, not interpretations
- [ ] Every figure and table is referenced in text
- [ ] Figures/tables referenced in order (Fig 1, then 2, then 3...)
- [ ] Key findings from each figure described in text
- [ ] Statistical tests are complete (test statistic, df, p-value)
- [ ] Measures of variability included (SD, SE, or CI)
- [ ] Effect sizes reported when applicable
- [ ] Multiple comparison corrections noted
- [ ] Negative results reported honestly
- [ ] No methods details in results
- [ ] No interpretation or literature comparison
- [ ] Past tense for findings
- [ ] Logical organization (chronological, thematic, or by hypothesis)
- [ ] Clear transitions between different results
- [ ] Subheadings used if appropriate (long results section)

## Subheadings in Results

### When to Use

- Long results section (>4 pages)
- Multiple distinct experiments or themes
- Makes organization clearer
- Common in high-impact journals

### Example Subheading Structures

**By experiment type:**
- In vitro effects
- In vivo validation
- Clinical outcomes

**By system:**
- Cellular effects
- Molecular mechanisms
- Functional outcomes

**By hypothesis:**
- X regulates Y
- Regulation depends on Z
- Disruption impairs function

**By time or condition:**
- Acute effects
- Chronic effects
- Dose-response relationships

## Figure Reference Formats by Field

### Biology/Medicine (parenthetical common)

"Expression increased in treated cells (Figure 1A)"
"Patient demographics are summarized (Table 1)"

### Engineering/CS (in-sentence common)

"Figure 1A shows system throughput under varying loads"
"Table 1 compares our approach to existing methods"

### Both acceptable - be consistent within your paper

## Statistical Reporting Standards by Field

### Biology/Medicine

- Include measures of central tendency (mean or median)
- Include measures of variability (SD, SE, or CI)
- Report exact p-values
- State test used
- Report effect sizes (increasingly required)

### Clinical Research

- Follow CONSORT, STROBE, or other reporting guidelines
- Include confidence intervals
- Report both relative and absolute effects
- Number needed to treat for interventions

### Engineering/Computer Science

- Benchmark comparisons
- Confidence intervals or error bars
- Statistical significance when comparing methods
- Computational complexity analysis

### Psychology/Social Sciences

- Effect sizes mandatory (APA 7th edition)
- Confidence intervals
- Detailed ANOVA tables
- Power analysis

Always check journal-specific requirements!
