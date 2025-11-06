# Research Paper Structure

Comprehensive guide to organizing research papers using IMRAD and alternative structures.

## IMRAD Structure

The most common structure for experimental research papers in sciences.

### Introduction

**Purpose:** Establish context, identify gap, state objective

**Structure (funnel approach):**
1. **Broad context** (1-2 paragraphs)
   - Why general topic matters
   - Current state of knowledge
   - Importance to field/society

2. **Narrow focus** (1-2 paragraphs)
   - Zoom into specific problem area
   - Review relevant previous work
   - Synthesize what's known

3. **Knowledge gap** (1 paragraph)
   - What's unknown or needs improvement
   - Why gap is important
   - Types: knowledge gap, methodological gap, application gap, resolution gap

4. **Objective** (1 paragraph, often last)
   - What this study does
   - Clear, specific, testable
   - Can include hypotheses

5. **Approach overview** (optional, 1-2 sentences)
   - Brief description of methods used
   - Not detailed (that's for Methods section)

**Length:** Typically 1-3 pages (2-5 paragraphs for short papers, up to 10 for longer)

**What belongs here:**
- Context and background
- Literature review (synthesis, not exhaustive)
- Rationale for your study
- Clear statement of objectives

**What doesn't belong:**
- Detailed methods
- Results
- Extensive review (that's for review papers)
- Conclusions

### Methods

**Purpose:** Provide sufficient detail for reproducibility

**Organization strategies:**

**1. Chronological order** (when sequence matters)
```
Sample collection → Preparation → Analysis → Statistics
```
Good for: Sequential processes, time-dependent experiments

**2. By subsystem/component** (complex systems)
```
Cell culture → Transfection → Imaging → Quantification → Statistics
```
Good for: Multi-step experimental workflows

**3. By type of measurement** (multiple assays)
```
Materials → Cell viability assay → Gene expression → Protein analysis → Statistics
```
Good for: Multiple independent measurements

**Essential elements:**

**Materials:**
- Equipment: manufacturer, model, specifications
  - Example: "Cells were imaged using a Zeiss LSM 880 confocal microscope with 63× oil-immersion objective (NA 1.4)"
- Reagents: source, catalog number, concentration, purity
  - Example: "Cells were treated with doxorubicin (Sigma-Aldrich, D1515) at 1 μM for 24 hours"
- Software: name, version, manufacturer, parameters
  - Example: "Images were analyzed using ImageJ (version 1.53, NIH) with default thresholding"

**Experimental procedures:**
- Novel methods: Full detail, every step
- Standard methods: Brief description + citation
  - Example: "RNA was extracted using TRIzol (Invitrogen) according to manufacturer's instructions"
- Modified methods: Highlight what changed
  - Example: "Protein extraction was performed as described [12] with modifications: lysis buffer included 1% SDS instead of 0.1%"

**Sample sizes and replication:**
- Number of biological replicates
- Number of technical replicates
- Justification (power analysis if available)
- Example: "Each experiment was performed in triplicate with n=5 biological replicates (15 total measurements)"

**Statistical analysis:**
- Software used
- Tests applied and why
- Significance level (typically α = 0.05)
- Multiple comparison corrections
- Example: "Statistical analysis was performed in R (version 4.2). Groups were compared using one-way ANOVA followed by Tukey's HSD post-hoc test. Significance was set at p < 0.05 with Bonferroni correction for multiple comparisons"

**Ethics and compliance:**
- IRB/IACUC approval numbers
- Patient consent statements
- Data availability statements

**What belongs here:**
- All experimental procedures
- Equipment and materials
- Statistical methods
- Study design

**What doesn't belong:**
- Results (observations)
- Rationale (goes in Introduction)
- Interpretation (goes in Discussion)

### Results

**Purpose:** Report observations objectively without interpretation

**Organization strategies:**

**1. Chronological** (experiments in sequence)
- Good when later experiments build on earlier
- Example: "First we established... Next we tested... Finally we examined..."

**2. By importance** (most significant first)
- Lead with strongest findings
- Good for high-impact journals

**3. By theme** (related findings grouped)
- Example: "Effect on gene expression... Effect on protein levels... Effect on cell function..."

**4. By hypothesis** (addressed one by one)
- Example: "To test whether X affects Y... To determine if A relates to B..."

**Objectivity principles:**

**Do:**
- Report what you observed
- Use "data show," "results indicate," "we observed"
- Present negative results honestly
- Include measures of variability (SD, SE, CI)

**Don't:**
- Interpret findings (save for Discussion)
- Use "data prove," "data demonstrate" (too strong)
- Cherry-pick results
- Editorialize ("surprisingly," "interestingly" - save for Discussion)

**Figure and table integration:**

**Every figure/table MUST:**
- Be referenced in text: "Figure 1 shows..." or "(Figure 1)"
- Have findings described in text
- Be numbered sequentially in order of mention

**Good integration:**
"Expression of Gene X increased 3-fold in treated cells compared to controls (Figure 2A, p < 0.001)"

**Poor integration:**
"See Figure 2 for expression data" (too vague)
"Figure 2" (not referenced with findings)

**Statistical reporting requirements:**
- Test statistic with symbol (t, F, χ², etc.)
- Degrees of freedom
- Exact p-value (not just p < 0.05)
- Effect size when applicable
- Measures of variability

**Complete example:**
"Treated cells showed significantly higher viability than controls (mean ± SD: 85 ± 7% vs 62 ± 9%; t₁₈ = 4.32, p = 0.003, Cohen's d = 2.1)"

**What belongs here:**
- Observations from experiments
- Data (numbers, statistics)
- References to figures/tables

**What doesn't belong:**
- Interpretation ("this suggests," "therefore")
- Comparison to literature (goes in Discussion)
- Methods details
- Speculation

### Discussion

**Purpose:** Interpret results in context, acknowledge limitations, discuss implications

**Structure:**

**1. Opening (1 paragraph):**
- Restate main findings briefly
- State whether hypothesis supported
- Example: "This study demonstrated that X increases Y in Z cells, supporting our hypothesis that..."

**2. Interpretation in context (2-4 paragraphs):**
- What do results mean?
- How do they relate to previous work?
- Agree with literature: "Our findings are consistent with Smith et al. [15], who showed..."
- Disagree with literature: "In contrast to previous reports [12], we found... This discrepancy may be due to differences in..."
- Extend literature: "Building on earlier work [8], our results demonstrate that..."

**3. Mechanistic explanation (1-2 paragraphs, if applicable):**
- How do you explain your findings?
- What mechanism might account for observations?
- Support with evidence
- Acknowledge alternative explanations

**4. Limitations (1 paragraph):**
- Be honest but not self-defeating
- Explain impact of limitations
- Suggest how limitations could be addressed
- Example: "This study focused on acute treatment (24 hours); chronic effects may differ and warrant investigation. Additionally, we examined only one cell line; validation in primary cells would strengthen conclusions"

**5. Implications and significance (1 paragraph):**
- Why do findings matter?
- Impact on field
- Practical applications
- Example: "These findings suggest that targeting X may be a therapeutic strategy for Y"

**6. Future directions (1 paragraph):**
- What's the next logical step?
- What questions remain?
- How to extend this work?

**7. Closing (1-2 sentences):**
- Brief conclusion statement
- Or transition to separate Conclusions section

**Interpretation strategies:**

**Support interpretations:**
- "These results suggest..." (tentative)
- "The data support..." (stronger)
- "One possible explanation is..."
- "This finding may reflect..."

**Consider alternatives:**
- "Alternative explanations include..."
- "However, we cannot rule out..."
- "Another possibility is..."

**Link to introduction:**
- Address the gap identified in introduction
- Answer the question posed
- Show how findings advance field

**What belongs here:**
- Interpretation of results
- Comparison with literature
- Mechanistic explanations
- Limitations
- Implications
- Future directions

**What doesn't belong:**
- New results not in Results section
- Extensive methods details
- Speculation without evidence
- Repetition of results without interpretation

### Conclusions (if separate from Discussion)

**Purpose:** Summarize key findings and significance

**Structure:**
- 1-2 paragraphs
- Restate main findings
- Emphasize significance
- Broader implications
- No new information
- Can include future outlook

**Example:**
"In summary, this study demonstrated that X regulates Y through Z mechanism in A cells. These findings advance understanding of B process and suggest that C may be a therapeutic target for D. Future work examining E in F will further elucidate this pathway"

## Alternative Structures

### Short Communications / Letters

**Compressed IMRAD:**
- Introduction: 1-2 paragraphs
- Methods: Brief, often in figure legends
- Results and Discussion: Combined
- Typically 2-4 pages total
- Focus on novel findings

### Case Reports (Medicine)

**Structure:**
- Introduction: Background on condition
- Case Presentation: Patient history, examination, diagnosis
- Discussion: Significance, comparison to literature
- Conclusions: Take-home message

### Theoretical Papers (Physics, Math)

**Structure:**
- Introduction: Problem statement
- Theory/Model: Derivation, assumptions
- Results: Predictions, simulations
- Discussion: Implications, testability
- Conclusions

### Engineering/Computer Science

**Design-focused:**
- Introduction: Problem, requirements
- Related Work: Existing approaches
- Design: System architecture, algorithms
- Implementation: Technical details
- Evaluation: Benchmarks, comparisons
- Discussion: Limitations, future work
- Conclusions

**Algorithm-focused:**
- Introduction and motivation
- Preliminaries: Definitions, notation
- Algorithm description
- Analysis: Complexity, correctness
- Experiments: Benchmarks
- Conclusions

### Humanities/Social Sciences

**More flexible structures:**
- Introduction
- Literature Review (often substantial, separate section)
- Methodology
- Analysis (may have subsections by theme)
- Discussion
- Conclusions

## Section Integration and Flow

### Transitions Between Sections

**Introduction → Methods:**
- Often no explicit transition
- Can use: "To address these questions, we..." at end of introduction

**Methods → Results:**
- No transition needed (clear section break)
- Begin results with first finding

**Results → Discussion:**
- Can use transitional opening: "These results demonstrate that..."
- Or jump directly into interpretation

### Maintaining Coherent Narrative

**Use of research questions/hypotheses:**
- Pose in Introduction
- Address systematically in Results
- Answer in Discussion

**Consistent terminology:**
- Use same terms throughout
- Don't switch between synonyms
- Define abbreviations at first use in each section (or just in abstract + first use in text)

**Logical progression:**
- Each section builds on previous
- Results address Introduction's questions
- Discussion interprets Results in context of Introduction

### Avoiding Redundancy

**What's unique to each section:**

**Introduction:**
- Background literature
- Rationale
- Objectives

**Methods:**
- How you did experiments
- Technical details

**Results:**
- What you found
- Observations, data

**Discussion:**
- What it means
- How it fits literature
- Implications

**Don't repeat unnecessarily:**
- Results in Discussion: Brief restatement OK, but don't repeat all data
- Methods in Results: If you must mention a method, keep it very brief
- Literature review: Don't repeat Introduction's lit review in Discussion

### Common "What Goes Where" Mistakes

| Content | Wrong section | Right section |
|---------|---------------|---------------|
| Interpretation of data | Results | Discussion |
| Detailed methods | Results | Methods |
| New results | Discussion | Results |
| Rationale for approach | Methods | Introduction |
| Literature review | Discussion (extensive) | Introduction |
| Limitations | Results | Discussion |
| Conclusions | Results | Discussion/Conclusions |

## Field-Specific Variations

### Biology/Medicine

**Experimental biology:**
- Standard IMRAD
- Methods often very detailed
- Results organized by assay or system
- Discussion emphasizes biological significance

**Clinical trials:**
- Introduction: Clinical background
- Methods: Trial design, CONSORT guidelines
- Results: Patient flow, outcomes, adverse events
- Discussion: Clinical implications, limitations

### Chemistry

**Synthesis papers:**
- Introduction: Target molecule, rationale
- Results and Discussion: Combined, organized by compound
- Experimental Section: Detailed procedures, characterization
- Supporting Information: Spectra, additional data

### Physics

**Experimental:**
- Introduction and theory
- Experimental methods and apparatus
- Results and analysis
- Discussion and conclusions

**Theoretical:**
- Introduction
- Model/theory development
- Results (predictions, simulations)
- Discussion (implications, testability)

### Engineering

**Design papers:**
- Problem and requirements
- Related work
- System design
- Implementation
- Validation/testing
- Discussion

### Computer Science

**Algorithm papers:**
- Problem definition
- Related work
- Algorithm description
- Complexity analysis
- Experimental evaluation
- Conclusions

**Systems papers:**
- Motivation
- Design goals
- System architecture
- Implementation
- Evaluation (benchmarks)
- Related work (can be earlier or later)
- Conclusions

## Journal-Specific Requirements

### Nature/Science (High-impact, short format)

- Very short main text (often 3-4 pages)
- Most methods in supplementary
- Combined Results and Discussion
- Emphasis on broad significance
- Strict word limits

### PLOS (Open access)

- Standard IMRAD
- Emphasis on reproducibility
- Data availability required
- Methods can be detailed

### Cell/Cell Press Journals

- STAR Methods (structured methods section)
- Key Resources Table
- Emphasis on rigor and reproducibility

### IEEE (Engineering)

- Introduction
- Often includes "Related Work" section
- Methodology
- Results/Experimental Evaluation
- Discussion and Conclusions

**Always check specific journal guidelines** - requirements vary significantly.

## Length Guidelines by Section (typical)

| Section | Short paper | Standard paper | Long paper |
|---------|-------------|----------------|------------|
| Abstract | 150-250 words | 200-300 words | 250-350 words |
| Introduction | 1-2 pages | 2-3 pages | 3-5 pages |
| Methods | 1-2 pages | 2-4 pages | 4-8 pages |
| Results | 2-3 pages | 3-6 pages | 6-10 pages |
| Discussion | 1-2 pages | 2-4 pages | 4-6 pages |
| Total | 5-9 pages | 9-17 pages | 17-29 pages |

*Not including references, figures, tables - varies by journal*
