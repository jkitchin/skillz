---
description: Suggest refactoring improvements for code
allowed-tools: ["Read", "Glob", "Grep"]
argument-hint: <file-path>
---

# Refactor Suggestions

Analyze code and suggest refactoring improvements.

## Instructions

1. Read the file specified: $ARGUMENTS
   - If no file specified, ask the user which file to analyze

2. Identify refactoring opportunities:
   - **Extract Method**: Long methods that should be broken up
   - **Rename**: Unclear variable/function names
   - **Remove Duplication**: Repeated code patterns
   - **Simplify Conditionals**: Complex if/else chains
   - **Reduce Complexity**: Deeply nested code
   - **Improve Structure**: Better organization, separation of concerns
   - **Apply Patterns**: Where design patterns could help

3. Prioritize suggestions by impact and effort

## Output Format

```markdown
## Refactoring Suggestions: <filename>

### Overview
<Current state assessment>

### High Priority

#### 1. <Refactoring Name>
**Location**: Lines X-Y
**Issue**: <Description of the problem>
**Suggestion**: <How to improve it>
**Before**:
```<language>
<current code>
```
**After**:
```<language>
<suggested code>
```

### Medium Priority
...

### Low Priority (Nice to Have)
...

### Metrics
- Current complexity: <assessment>
- Estimated improvement: <assessment>
```

File to analyze: $ARGUMENTS
