# Methods Writing

Comprehensive guide to writing reproducible methods sections.

## Core Principle: Reproducibility

**Goal:** Another researcher should be able to reproduce your work using only your methods section

**Test:** If you gave your methods to a competent researcher in your field, could they repeat your experiments?

## Reproducibility Checklist

### Equipment and Instruments

**What to include:**
- Manufacturer name
- Model number/name
- Key specifications
- Software version (if applicable)

**Examples:**

**Microscopy:**
"Cells were imaged using a Zeiss LSM 880 confocal microscope equipped with 405 nm, 488 nm, 561 nm, and 633 nm lasers. Images were acquired using a Plan-Apochromat 63×/1.4 NA oil-immersion objective with 2× optical zoom"

**Spectroscopy:**
"Absorbance measurements were performed using a BioTek Synergy H1 plate reader (BioTek Instruments, Winooski, VT) at 450 nm with 630 nm reference wavelength"

**Chromatography:**
"HPLC analysis was conducted on an Agilent 1260 Infinity II system equipped with quaternary pump, autosampler, and diode array detector"

**Computational:**
"Calculations were performed on a compute cluster with Intel Xeon Gold 6248 processors (2.5 GHz, 20 cores) and 192 GB RAM per node"

### Reagents and Materials

**What to include:**
- Source (manufacturer/supplier)
- Catalog number
- Concentration/purity
- Storage conditions (if critical)
- Lot number (if results could vary by lot)

**Examples:**

**Cell culture:**
"HeLa cells (ATCC CCL-2) were cultured in DMEM (Gibco, #11965-092) supplemented with 10% fetal bovine serum (Sigma-Aldrich, #F2442), 2 mM L-glutamine (Gibco, #25030-081), and 1% penicillin-streptomycin (Gibco, #15140-122) at 37°C in 5% CO₂ humidified incubator"

**Chemicals:**
"Doxorubicin hydrochloride (Sigma-Aldrich, D1515, ≥98% HPLC) was dissolved in DMSO (Sigma-Aldrich, D2650) at 10 mM stock concentration and stored at -20°C"

**Antibodies:**
"Primary antibodies used were: rabbit anti-LC3B (1:1000, Cell Signaling Technology, #2775), mouse anti-β-actin (1:5000, Sigma-Aldrich, #A5441). Secondary antibodies were: goat anti-rabbit IgG-HRP (1:5000, Jackson ImmunoResearch, #111-035-003), goat anti-mouse IgG-HRP (1:5000, Jackson ImmunoResearch, #115-035-003)"

**Biological samples:**
"Mouse liver tissue was obtained from 8-week-old male C57BL/6J mice (Jackson Laboratory, stock #000664) following IACUC-approved protocols (protocol #2023-001)"

### Software and Code

**What to include:**
- Software name
- Version number
- Manufacturer/source
- Key parameters/settings
- Availability (if custom code)

**Examples:**

**Data analysis:**
"Statistical analysis was performed using R (version 4.2.1, R Foundation for Statistical Computing). Graphs were created using ggplot2 (version 3.4.0)"

**Image analysis:**
"Images were analyzed using ImageJ (version 1.53t, NIH). For each image, background was subtracted using rolling ball algorithm (radius = 50 pixels), and cells were segmented using Otsu thresholding. Mean fluorescence intensity was measured for each cell (n > 50 cells per condition)"

**Custom code:**
"Data processing was performed using custom Python scripts (Python 3.9.7 with NumPy 1.21.2, SciPy 1.7.1, pandas 1.3.3). Code is available at https://github.com/username/repo"

**Machine learning:**
"Neural networks were implemented in PyTorch (version 1.12.0) and trained on NVIDIA A100 GPUs using Adam optimizer (learning rate = 0.001, batch size = 32) for 100 epochs. Model architecture and training code are available at [URL]"

### Statistical Methods

**What to include:**
- Software used
- Tests applied
- Significance level (α)
- Multiple comparison corrections
- Post-hoc tests
- Assumptions checked

**Examples:**

**Basic comparison:**
"Statistical significance was assessed using two-tailed Student's t-test in GraphPad Prism (version 9.0). Data are presented as mean ± standard deviation (SD). Significance was set at p < 0.05"

**Multiple comparisons:**
"Groups were compared using one-way ANOVA followed by Tukey's honest significant difference (HSD) post-hoc test. P-values were adjusted for multiple comparisons using Bonferroni correction. Significance was set at adjusted p < 0.05"

**Non-parametric:**
"Data were not normally distributed (Shapiro-Wilk test, p < 0.05), so groups were compared using Mann-Whitney U test (two groups) or Kruskal-Wallis test followed by Dunn's post-hoc test (multiple groups)"

**Regression:**
"Linear regression was performed to assess relationship between X and Y. Pearson correlation coefficient and associated p-value were calculated. Assumptions of linearity, homoscedasticity, and normality of residuals were verified by visual inspection of residual plots and Shapiro-Wilk test"

**Survival analysis:**
"Survival curves were estimated using Kaplan-Meier method and compared using log-rank test. Hazard ratios were calculated using Cox proportional hazards model. Proportional hazards assumption was tested using Schoenfeld residuals"

### Sample Sizes and Replication

**What to include:**
- Number of biological replicates
- Number of technical replicates
- Sample size justification (power analysis if done)
- How replicates were combined/analyzed

**Examples:**

**Cell culture experiments:**
"Each experiment was performed independently three times (biological replicates), with triplicate measurements within each experiment (technical replicates). Data points represent mean of technical replicates, with n=3 biological replicates"

**Animal studies:**
"Sample size was determined by power analysis (α = 0.05, power = 0.80, effect size = 1.5 based on pilot data) using G*Power software, yielding n = 6 mice per group. A total of 24 mice were used (4 groups, n=6 each)"

**Human subjects:**
"We enrolled 150 participants (75 cases, 75 controls) based on power analysis to detect an odds ratio of 2.0 with 80% power at α = 0.05"

**Computational:**
"Each simulation was run 100 times with different random seeds. Results are reported as mean ± standard deviation across simulations"

### Randomization and Blinding

**When applicable, specify:**
- How samples/subjects were randomized
- Who was blinded
- When blinding was broken

**Examples:**

**Animal studies:**
"Mice were randomly assigned to treatment groups using computer-generated random numbers. Treatment was administered by investigator blinded to group assignment. Analysis was performed blinded to treatment group, with code broken only after completion of analysis"

**Clinical trial:**
"Participants were randomized 1:1 to treatment or placebo groups using block randomization (block size = 4) stratified by age (<50 vs ≥50 years). Participants, care providers, and outcome assessors were blinded to treatment assignment. Randomization sequence was generated by independent statistician"

### Ethics Approvals and Consent

**Include:**
- IRB/IACUC approval number
- Consent procedures
- Adherence to guidelines
- Data protection measures

**Examples:**

**Human subjects:**
"This study was approved by the Institutional Review Board of [Institution] (protocol #2023-001). All participants provided written informed consent. Study was conducted in accordance with Declaration of Helsinki"

**Animal research:**
"All animal procedures were approved by the Institutional Animal Care and Use Committee (IACUC) of [Institution] (protocol #2023-001) and complied with NIH Guide for the Care and Use of Laboratory Animals"

**Human samples:**
"De-identified human tissue samples were obtained from [Biobank] under an IRB-approved protocol (#2023-001). Original consent forms permitted use for research purposes"

## Organization Strategies

### 1. Chronological Order

**When to use:** Sequential processes, multi-step protocols

**Example structure:**
```
Sample collection
↓
Sample preparation
↓
Measurement
↓
Analysis
↓
Statistical methods
```

**Example:**
"Blood samples were collected via venipuncture into EDTA tubes and processed within 2 hours. Plasma was separated by centrifugation at 1000g for 10 minutes at 4°C and stored at -80°C until analysis. Thawed samples were diluted 1:10 in assay buffer and analyzed using ELISA following manufacturer's instructions. Data were analyzed using..."

### 2. By Subsystem/Component

**When to use:** Complex systems with multiple independent components

**Example structure:**
```
Cell culture
Transfection
Western blotting
Microscopy
RT-qPCR
Statistical analysis
```

**Example:**
"**Cell culture.** HeLa cells were maintained in...

**Transfection.** Cells were transfected using...

**Western blotting.** Protein lysates were prepared...

**Microscopy.** Cells were imaged using...

**Statistical analysis.** Data were analyzed using..."

### 3. By Type of Measurement

**When to use:** Multiple independent assays applied to same samples

**Example structure:**
```
Materials and sample preparation
Cell viability assay
Gene expression analysis
Protein quantification
Metabolite measurements
Statistical methods
```

## Detail Level Guidelines

### Standard Procedures: Brief + Citation

**When procedure is well-established and standard:**

**Example:**
"RNA was extracted using TRIzol reagent (Invitrogen) according to manufacturer's instructions"

**Example:**
"Western blotting was performed as described previously [15]"

**Example:**
"Statistical analysis was performed using standard methods [20]"

### Novel Procedures: Full Detail

**When procedure is new or unique to your lab:**

**Example:**
"Cell migration was assessed using a custom microfluidic device. Briefly, devices were fabricated by soft lithography using SU-8 photoresist (MicroChem) to create a master mold with 50 μm wide channels. PDMS (Dow Corning, 10:1 base:curing agent ratio) was poured over the mold, cured at 80°C for 2 hours, and bonded to glass coverslips by oxygen plasma treatment (30 seconds, 50 W). Channels were coated with fibronectin (10 μg/mL) for 1 hour at 37°C. Cells were loaded at 10⁶ cells/mL and allowed to adhere for 2 hours. Chemoattractant gradient was established by replacing medium in outlet reservoir every 30 minutes. Cell positions were tracked every 5 minutes for 4 hours using time-lapse microscopy"

### Modified Procedures: Highlight Changes

**When you adapted an existing protocol:**

**Example:**
"Protein extraction was performed as described [15] with the following modifications: lysis buffer contained 1% SDS instead of 0.1% SDS, and samples were sonicated for 30 seconds (instead of 10 seconds) to ensure complete lysis"

**Example:**
"DNA was isolated using Qiagen DNeasy kit following manufacturer's protocol, except that elution volume was reduced from 200 μL to 50 μL to increase DNA concentration"

## Field-Specific Patterns

### Wet Lab Biology

**Typical subsections:**
- Cell culture
- Transfection/treatment
- Immunofluorescence
- Western blotting
- qPCR
- Flow cytometry
- Animal procedures
- Statistical analysis

**Example methods statement:**
"Cells were seeded at 5×10⁴ cells/well in 24-well plates 24 hours before treatment. Cells were treated with 1 μM doxorubicin for 24 hours, then fixed in 4% paraformaldehyde for 15 minutes, permeabilized with 0.1% Triton X-100 for 10 minutes, and blocked with 5% BSA for 1 hour. Primary antibodies were incubated overnight at 4°C, followed by secondary antibodies for 1 hour at room temperature. Nuclei were stained with DAPI (1 μg/mL) for 5 minutes. Images were acquired using..."

### Chemical Synthesis

**Typical organization:**
- General methods (chromatography, spectroscopy)
- Synthesis procedures (compound by compound)
- Characterization data in each procedure

**Example:**
"**General Methods.** All reagents were purchased from commercial sources and used without further purification unless noted. Reactions were monitored by thin-layer chromatography on silica gel plates. Flash column chromatography was performed using silica gel (40-63 μm). ¹H and ¹³C NMR spectra were recorded on Bruker 500 MHz spectrometer. Chemical shifts are reported in ppm relative to TMS. HRMS was performed on Agilent 6530 Q-TOF.

**Synthesis of Compound 3.** To a solution of 2 (500 mg, 2.1 mmol) in dry DCM (10 mL) at 0°C was added..."

### Computational/Simulation

**Typical elements:**
- Hardware specifications
- Software and versions
- Algorithms and parameters
- Initial conditions
- Convergence criteria
- Validation methods

**Example:**
"Molecular dynamics simulations were performed using GROMACS 2021.3. Systems were solvated in TIP3P water box extending 10 Å beyond protein surface. Ions were added to neutralize charge and achieve 150 mM NaCl concentration. Energy minimization was performed using steepest descent algorithm until maximum force < 1000 kJ/mol/nm. Systems were equilibrated in NVT ensemble for 100 ps followed by NPT ensemble for 100 ps, maintaining temperature at 300 K using V-rescale thermostat and pressure at 1 bar using Parrinello-Rahman barostat. Production runs of 100 ns were performed with 2 fs time step. Coordinates were saved every 10 ps for analysis"

### Clinical/Human Studies

**Typical elements:**
- Study design
- Participant selection and exclusion criteria
- Recruitment methods
- Intervention details
- Outcome measures
- Data collection procedures

**Example:**
"This prospective cohort study enrolled patients diagnosed with Type 2 diabetes at [Hospital] between January 2020 and December 2022. Inclusion criteria were: age 18-75 years, HbA1c >7%, no insulin use. Exclusion criteria were: Type 1 diabetes, severe renal impairment (eGFR <30 mL/min/1.73m²), active cancer. Of 450 screened patients, 300 met criteria and consented to participate. Baseline data included demographics, medical history, laboratory tests, and medication use. Participants were followed for 12 months with visits at 0, 3, 6, and 12 months. Primary outcome was change in HbA1c. Secondary outcomes included fasting glucose, lipid profile, and adverse events"

### Field Studies/Ecology

**Typical elements:**
- Study site description
- Sampling design
- Species identification
- Environmental measurements
- Data collection protocols

**Example:**
"Field sampling was conducted at three sites along the [River] (Site 1: 45.5°N, 122.7°W; Site 2: 45.6°N, 122.6°W; Site 3: 45.7°N, 122.5°W) from May to September 2022. At each site, water samples were collected monthly from five locations spaced 10 m apart. Samples were collected in sterile 1 L bottles at 0.5 m depth, kept on ice, and processed within 4 hours. Temperature, pH, dissolved oxygen, and conductivity were measured in situ using YSI Pro2030 multiparameter meter"

## Tense and Voice

### Tense

**Past tense (standard for experimental procedures):**
"Cells were cultured in DMEM..."
"Samples were collected..."
"Data were analyzed using..."

**Present tense (for general procedures or descriptions):**
"The assay measures..."
"This method allows..."
(Less common in methods, more in general descriptions)

### Voice

**Field norms vary:**

**Passive voice (traditional biology/chemistry):**
"Cells were seeded at 5×10⁴ cells/well"
"RNA was extracted using TRIzol"
"Samples were analyzed by HPLC"

**Active voice (increasingly common, especially engineering/CS):**
"We seeded cells at 5×10⁴ cells/well"
"We extracted RNA using TRIzol"
"We analyzed samples by HPLC"

**Both acceptable - be consistent within your paper**

## Common Methods Section Problems

### Problem: Insufficient Detail

**Too vague:**
"Cells were cultured normally and then treated with drug"

**Better:**
"HeLa cells (ATCC CCL-2) were cultured in DMEM with 10% FBS and 1% pen-strep at 37°C in 5% CO₂. Cells were seeded at 5×10⁴ cells/well in 24-well plates 24 hours before treatment. Doxorubicin (Sigma-Aldrich, D1515) was added at 1 μM final concentration and cells incubated for 24 hours"

### Problem: Too Much Detail

**Excessive:**
"We opened a new package of pipette tips (Brand X, catalog #123). Using a P1000 pipette (Brand Y), we transferred 1 mL of medium into a 15 mL conical tube (Brand Z, catalog #456). We then placed the tube in the incubator, opened the incubator door, placed the tube on the middle shelf..."

**Better:**
"Cells were cultured in complete DMEM at 37°C in 5% CO₂"

### Problem: Methods Scattered in Results

**Wrong:**
"Expression increased 2-fold (Figure 1). To measure this, we used qPCR with primers designed using Primer3 software. RNA was extracted using TRIzol..."

**Right - Methods:**
"RNA was extracted using TRIzol (Invitrogen) and analyzed by qPCR using primers designed with Primer3 software (sequences in Table S1)"

**Right - Results:**
"Expression increased 2-fold in treated cells (Figure 1, p < 0.01)"

### Problem: Missing Statistical Information

**Incomplete:**
"Data were analyzed statistically"

**Complete:**
"Statistical analysis was performed using R (version 4.2.1). Groups were compared using one-way ANOVA followed by Tukey's HSD post-hoc test. Data are presented as mean ± SD from three independent experiments. Significance was set at p < 0.05"

### Problem: Unclear Replication

**Unclear:**
"The experiment was done in triplicate"
(Technical or biological? Within one experiment or independent experiments?)

**Clear:**
"Each experiment was performed independently three times (biological replicates), with triplicate measurements per experiment (technical replicates). Data represent mean ± SD of biological replicates"

## Methods Checklist

- [ ] Equipment: manufacturer, model, specifications
- [ ] Reagents: source, catalog number, concentration
- [ ] Software: name, version, parameters
- [ ] Sample sizes: justified and clearly stated
- [ ] Replication: biological and technical distinguished
- [ ] Statistical methods: tests, corrections, software
- [ ] Randomization and blinding (if applicable)
- [ ] Ethics approvals with protocol numbers
- [ ] Novel methods: sufficient detail to reproduce
- [ ] Standard methods: properly cited
- [ ] Logical organization (chronological, by subsystem, or by measurement)
- [ ] Consistent tense and voice
- [ ] No results or interpretation included
- [ ] All information needed for reproducibility

## Data and Code Availability

**Include statements about:**
- Data availability (repository, accession numbers, reasonable request)
- Code availability (GitHub, institutional repository, reasonable request)
- Materials availability (sharing agreements for reagents, cell lines, etc.)

**Examples:**

**Data availability:**
"Raw data are available at Gene Expression Omnibus (accession GSE123456). Processed data and analysis scripts are available at https://github.com/username/project"

**Code availability:**
"Source code for simulations is available at https://github.com/username/repo under MIT license. Analysis scripts are provided as Supplementary Software"

**Materials:**
"Plasmids and cell lines are available from the corresponding author upon reasonable request and completion of material transfer agreement"
