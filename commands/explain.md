---
description: Explain what code does in plain English
allowed-tools: ["Read", "Glob", "Grep"]
argument-hint: <file-path or function-name>
---

# Explain Code

Explain what the specified code does in clear, plain English.

## Instructions

1. Read the file or locate the function specified: $ARGUMENTS
   - If a file path, read the entire file
   - If a function/class name, search for it in the codebase

2. Analyze the code and explain:
   - **Purpose**: What problem does this code solve?
   - **How it works**: Step-by-step explanation of the logic
   - **Inputs/Outputs**: What goes in and what comes out
   - **Dependencies**: What other code/libraries does it rely on
   - **Side effects**: Does it modify state, make network calls, etc.?

3. Adjust explanation depth based on code complexity

## Output Format

```markdown
## Explanation: <code identifier>

### Purpose
<1-2 sentences on what this code does>

### How It Works
1. <Step 1>
2. <Step 2>
3. <Step 3>

### Inputs
- `param1`: <description>
- `param2`: <description>

### Outputs
- Returns: <description>

### Key Concepts
- <Concept 1>: <brief explanation>
- <Concept 2>: <brief explanation>

### Example Usage
```<language>
<example code>
```
```

Code to explain: $ARGUMENTS
