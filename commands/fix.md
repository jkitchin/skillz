---
description: Diagnose and fix a specific error
allowed-tools: ["Read", "Glob", "Grep", "Bash"]
argument-hint: <error-message or file:line>
---

# Fix Error

Diagnose and suggest fixes for a specific error.

## Instructions

1. Parse the error information provided: $ARGUMENTS
   - Error message text
   - File and line number if available
   - Stack trace if provided

2. Investigate the error:
   - Read the relevant source file(s)
   - Understand the context around the error location
   - Check for common causes of this error type
   - Look for related code that might be involved

3. Determine the root cause:
   - Is it a syntax error?
   - Type mismatch?
   - Null/undefined reference?
   - Missing import/dependency?
   - Logic error?
   - Configuration issue?

4. Propose a fix with explanation

## Output Format

```markdown
## Error Diagnosis

### Error
```
<error message>
```

### Root Cause
<Explanation of why this error occurs>

### Location
- File: <path>
- Line: <number>
- Function: <name>

### Fix

**Change this**:
```<language>
<problematic code>
```

**To this**:
```<language>
<fixed code>
```

### Explanation
<Why this fix works>

### Prevention
<How to prevent similar errors in the future>
```

Error details: $ARGUMENTS
