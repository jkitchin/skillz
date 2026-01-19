---
description: Create a formatted lab notebook entry
allowed-tools: ["Read", "Bash", "Glob"]
argument-hint: <title or topic>
---

# Create Lab Notebook Entry

Generate a structured lab notebook entry for documenting research or development work.

## Instructions

1. Gather context:
   - Check recent git commits for today's work
   - Look for modified files
   - Note any test results or outputs

2. Create a structured entry with:
   - Date and title
   - Objective/Goal
   - Methods/Approach
   - Results/Observations
   - Conclusions/Next Steps

3. Include relevant code snippets, commands, or data

## Output Format

```markdown
# Lab Notebook Entry

**Date**: <YYYY-MM-DD>
**Title**: <title from $ARGUMENTS or inferred>
**Project**: <project name>

## Objective

<What you're trying to accomplish>

## Background

<Relevant context or prior work>

## Methods

### Approach
<Description of the approach taken>

### Code/Commands
```<language>
<relevant code or commands>
```

### Tools Used
- Tool 1
- Tool 2

## Results

### Observations
<What happened, what you observed>

### Data
<Any relevant data, metrics, or outputs>

### Figures
<Description of any plots or visualizations>

## Discussion

<Interpretation of results, what worked, what didn't>

## Conclusions

<Key takeaways>

## Next Steps

- [ ] Next task 1
- [ ] Next task 2

## References

<Any papers, docs, or resources used>

---
*Entry created: <timestamp>*
```

Entry title or topic: $ARGUMENTS
