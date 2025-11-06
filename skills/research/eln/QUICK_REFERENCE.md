# Electronic Lab Notebook Quick Reference

## File Organization

```
notebook/
├── 2025/
│   ├── 01-January/
│   │   ├── 2025-01-15.org
│   │   └── data/
│   │       └── 2025-01-15/
│   └── 02-February/
├── templates/
├── index.org
└── README.org
```

## Daily File Template

```org
#+TITLE: Lab Notebook - 2025-01-15
#+AUTHOR: Your Name
#+DATE: [2025-01-15 Wed]
#+FILETAGS: :experiment:simulation:
#+STARTUP: overview

* Daily Summary
Brief overview of today's work.

* [HH:MM] Entry Title
:PROPERTIES:
:ID: 2025-01-15-001
:PROJECT: Project_Name
:STATUS: In_Progress
:END:

** Objective
** Methods
** Results
** Analysis
** Conclusions
** Next Steps
```

## Entry Templates

### Quick Entry
```org
* [09:30] Short Description
:PROPERTIES:
:ID: 2025-01-15-001
:END:

** What I Did

** Why

** Results

** Next Steps
```

### Full Research Entry
```org
* [09:30] Project - Specific Task
:PROPERTIES:
:ID: 2025-01-15-001
:PROJECT: Project_Name
:TYPE: Calculation|Experiment|Analysis
:STATUS: Planning|In_Progress|Complete|On_Hold
:RELATED: [[file:2025-01-10.org::*Previous Work]]
:END:

** Objective
Clear statement of what you're trying to accomplish.

** Hypothesis
What you expect to happen and why.

** Background/Context
- Why is this important?
- What previous work led to this?

** Methods
Detailed description of how work was performed.

#+BEGIN_SRC python
# Include actual code
#+END_SRC

*** Parameters
- Parameter 1: value
- Parameter 2: value

** Results
*** Observations
What actually happened.

*** Data
| Column 1 | Column 2 |
|----------+----------|
| Value    | Value    |

** Analysis
What do these results mean?

** Conclusions
- Key finding 1
- Key finding 2

** Issues/Problems
Any problems and resolutions.

** Next Steps
- [ ] Task 1
- [ ] Task 2
```

### Literature Review
```org
* [14:00] Literature - Paper Title
:PROPERTIES:
:ID: 2025-01-15-002
:TOPIC: Subject_Area
:END:

** Citation
Author, A. et al. (2025) "Title"
Journal, DOI: 10.xxxx/xxxxx

** Key Findings
-

** Relevance to Our Work
-

** Questions Raised
-

** Follow-Up References
- [ ] Reference 1
- [ ] Reference 2
```

### Meeting Notes
```org
* [10:00] Meeting - Topic
:PROPERTIES:
:ID: 2025-01-15-003
:TYPE: Meeting
:ATTENDEES: Person1, Person2
:END:

** Agenda
1.
2.

** Discussion
-

** Action Items
- [ ] Task 1 (assigned to: Name)
- [ ] Task 2 (assigned to: Name)

** Next Meeting
<2025-01-22 Wed 10:00>
```

### Problem/Debugging
```org
* [15:30] Problem - Issue Description
:PROPERTIES:
:ID: 2025-01-15-004
:TYPE: Problem
:STATUS: Resolved|Ongoing
:END:

** Problem
Clear description of the issue.

** Error Messages
#+BEGIN_EXAMPLE
Error output here
#+END_EXAMPLE

** Attempted Solutions
1. Solution 1 - Result
2. Solution 2 - Result
3. **Solution 3 - Success**

** Root Cause
-

** Resolution
-

** Lessons Learned
-
```

## Org-Mode Syntax

### Properties Drawer
```org
:PROPERTIES:
:ID: unique-identifier
:PROJECT: Project_Name
:STATUS: In_Progress
:STARTED: [2025-01-15 Wed 09:30]
:COMPLETED: [2025-01-15 Wed 16:45]
:RELATED: [[link]]
:DATA: [[file:./data/file.csv]]
:END:
```

### Tags
```org
# File-level tags
#+FILETAGS: :experiment:simulation:catalyst:

# Heading tags
* Entry Title                                    :important:urgent:
```

**Common Tags:**
- Type: :experiment:, :simulation:, :analysis:, :literature:, :meeting:
- Topic: :catalyst:, :materials:, :synthesis:
- Status: :todo:, :in_progress:, :done:, :failed:
- Priority: :urgent:, :important:, :routine:

### Links
```org
# Internal file links
[[file:2025-01-10.org::*Heading Name]]
[[id:2025-01-15-001]]

# External files
[[file:~/data/results.csv][Data file]]
[[file:./data/2025-01-15/plot.png]]

# URLs
[[https://doi.org/10.1021/xxxxx][Paper]]

# Code files
[[file:~/projects/analysis.py][Analysis script]]
```

### Code Blocks
```org
# Python
#+BEGIN_SRC python :results output :session :exports both
import numpy as np
print("Result:", np.mean([1, 2, 3]))
#+END_SRC

# Bash
#+BEGIN_SRC bash :results output
grep "energy" OUTCAR | tail -1
#+END_SRC

# Include file
#+BEGIN_SRC python :results output :var filename="data.csv"
import pandas as pd
data = pd.read_csv(filename)
print(data.describe())
#+END_SRC
```

### Tables
```org
| Parameter  | Value | Units |
|------------+-------+-------|
| Temp       |   300 | K     |
| Pressure   |     1 | atm   |

# With formulas
#+TBLFM: @2$3=@2$1*@2$2
```

### TODO Items
```org
# TODO keywords
** TODO Task description
** DONE Completed task
** CANCELLED Cancelled task

# Checkboxes
- [ ] Task 1
- [X] Task 2 (completed)
- [-] Task 3 (partial)

# With deadlines
** TODO Important task
DEADLINE: <2025-01-20 Fri>
:PROPERTIES:
:EFFORT: 2h
:END:
```

### Timestamps
```org
# Active (appears in agenda)
<2025-01-20 Fri 14:00-15:00>

# Inactive (documentation only)
[2025-01-15 Wed 09:30]

# Date range
<2025-01-10 Mon>--<2025-01-30 Mon>
```

### Images and Figures
```org
#+CAPTION: Figure description
#+NAME: fig:label
#+ATTR_ORG: :width 400
[[file:./data/2025-01-15/plot.png]]

Reference: See Figure [[fig:label]].
```

## Documentation Checklist

### Every Entry Should Have:
- [ ] Date and time stamp
- [ ] Unique ID
- [ ] Clear objective
- [ ] Methods description
- [ ] Results
- [ ] Interpretation/analysis
- [ ] Next steps

### For Reproducibility:
- [ ] Software versions
- [ ] All parameters specified
- [ ] Code included or linked
- [ ] Data files referenced
- [ ] Random seeds (if applicable)
- [ ] Computing environment described

## Search Commands (Emacs)

```
C-c / m          Match tags
C-c / t          Show TODOs
C-c / r          Regex search
C-c a s          Search all agenda files
C-c a m          Match tags in agenda
```

## Common Property Values

### STATUS
- Planning
- In_Progress
- Complete
- On_Hold
- Failed
- Cancelled

### TYPE
- Experiment
- Simulation
- Calculation
- Analysis
- Literature
- Meeting
- Problem

### PROJECT
Use consistent project names:
- Catalyst_Screening
- Method_Development
- Literature_Review
- Collaboration_ProjectName

## Writing Guidelines

### Be Specific
```
Bad:  "Ran calculation with good parameters"
Good: "Ran VASP calculation with PBE functional, 400 eV cutoff,
       6×6×1 k-points, 0.05 eV/Å force convergence"
```

### Include Units
```
Bad:  "Temperature was 300"
Good: "Temperature was 300 K"
```

### Quantify Results
```
Bad:  "Energy decreased significantly"
Good: "Energy decreased by 0.25 eV (15% reduction)"
```

### Document Failures
```
** Problem
Calculation failed to converge after 100 steps.

** Cause
ISMEAR=1 inappropriate for system with band gap.

** Solution
Changed to ISMEAR=0 (Gaussian smearing).
Converged in 45 steps.

** Lesson
Always check for band gap before choosing smearing method.
```

## Data Management

### Linking Data Files
```org
** Data Files
- Raw: [[file:./data/2025-01-15/OUTCAR]]
- Processed: [[file:./data/2025-01-15/results.csv]]
- Figures: [[file:./data/2025-01-15/plots/]]

** Archive Location
~/archives/2025/january/project_001.tar.gz

** Checksums
- OUTCAR: 5d41402abc4b2a76b9719d911017c592
```

### Directory Structure for Data
```
data/
└── 2025-01-15/
    ├── input/
    │   ├── INCAR
    │   ├── POSCAR
    │   └── KPOINTS
    ├── output/
    │   ├── OUTCAR
    │   └── vasprun.xml
    ├── analysis/
    │   ├── results.csv
    │   └── analysis.ipynb
    └── figures/
        ├── structure.png
        └── energy_plot.png
```

## Monthly Summary Template

```org
#+TITLE: [Month] [Year] Research Summary
#+DATE: [YYYY-MM-DD]

* Overview
Brief summary of the month's work.

* Key Accomplishments
- Accomplishment 1
- Accomplishment 2

* Significant Results
** Result 1
Description and impact.

** Result 2
Description and impact.

* Challenges and Solutions
** Challenge 1
Problem:
Solution:

* Data Generated
- Number of calculations/experiments
- Total data volume
- Archive location

* Publications/Presentations
- Draft 1: Status
- Presentation 1: Date

* Collaborations
- Meeting with X on [date]
- Data shared with Y

* Next Month Goals
- [ ] Goal 1
- [ ] Goal 2

* Skills/Techniques Learned
- New skill 1
- New tool 2
```

## Cross-Referencing

### Linking Related Work
```org
** Related Work
- Background: [[file:2025-01-05.org::*Literature Review]]
- Previous attempt: [[file:2025-01-12.org::*Failed Calculation]]
- Follow-up: [[file:2025-01-18.org::*Extended Analysis]]
- Project overview: [[file:../index.org::*Project Name]]
```

### Creating Index
```org
#+TITLE: Project Index

* Active Projects
** [[file:2025/01-January/2025-01-15.org::*Project Name][Project Name]]
Status: Active
Started: [2025-01-05]

* Key Results
** [[file:2025/01-January/2025-01-15.org::*Important Finding][Important Finding]]
Date: [2025-01-15]
Result: Brief description

* Techniques/Methods
** [[file:2025/01-January/2025-01-12.org::*Method Description][How to Do X]]
```

## Best Practices Reminders

1. **Write daily** - Don't let notes accumulate
2. **Be honest** - Document failures and surprises
3. **Include context** - Explain why, not just what
4. **Link liberally** - Connect related entries
5. **Embed code** - Make analysis reproducible
6. **Archive data** - Keep data organized
7. **Never delete** - Make corrections forward
8. **Add metadata** - Use properties and tags
9. **Write clearly** - Professional, complete sentences
10. **Review regularly** - Monthly summaries

## Emergency Recovery

If you forget details:
1. Check git log for code changes
2. Look at file timestamps
3. Check calculation directories
4. Review email/chat history
5. Ask collaborators

Then document what you can recover and note gaps:
```org
** Reconstructed Entry [2025-01-15]
Note: Documented retrospectively on [2025-01-20].
Some details may be incomplete.

Based on file timestamps and git history:
[Document what you know]

Unknown:
- Exact parameters used for [X]
- Reasoning for choosing [Y]
```
