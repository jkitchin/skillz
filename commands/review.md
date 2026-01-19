---
description: Quick code review of a file or selection
allowed-tools: ["Read", "Glob", "Grep"]
argument-hint: <file-path>
---

# Code Review

Perform a quick code review on the specified file or code.

## Instructions

1. Read the file specified by the user: $ARGUMENTS
   - If no file specified, ask the user which file to review

2. Analyze the code for:
   - **Code Quality**: Readability, maintainability, naming conventions
   - **Potential Bugs**: Null checks, edge cases, error handling
   - **Security Issues**: Input validation, injection vulnerabilities, secrets
   - **Performance**: Inefficient patterns, unnecessary operations
   - **Best Practices**: Language idioms, design patterns, SOLID principles

3. Provide actionable feedback with specific line references

## Output Format

```markdown
## Code Review: <filename>

### Summary
<Overall assessment in 2-3 sentences>

### Issues Found

#### 🔴 Critical
- Line X: <issue description>

#### 🟡 Warnings
- Line Y: <issue description>

#### 🔵 Suggestions
- Line Z: <suggestion>

### Positive Aspects
- <What's done well>

### Recommended Actions
1. <Action item>
2. <Action item>
```

File to review: $ARGUMENTS
