---
description: Summarize a file, PR, or conversation
allowed-tools: ["Read", "Bash", "Glob", "WebFetch"]
argument-hint: <file-path or URL>
---

# Summarize

Create a concise summary of a file, pull request, webpage, or other content.

## Instructions

1. Identify what to summarize: $ARGUMENTS
   - If a file path, read the file
   - If a URL, fetch the content
   - If a PR number/URL, fetch PR details
   - If empty, summarize the current conversation

2. Create a summary appropriate to the content type:
   - **Code file**: Purpose, main components, key functions
   - **Documentation**: Key points, main topics
   - **PR/Diff**: Changes made, impact, motivation
   - **Article/Webpage**: Main thesis, key arguments, conclusions

3. Adjust length based on content (aim for ~20% of original or 3-5 key points)

## Output Format

```markdown
## Summary: <source identifier>

### Overview
<2-3 sentence high-level summary>

### Key Points
1. <Main point 1>
2. <Main point 2>
3. <Main point 3>

### Details

<More detailed breakdown if needed>

### Notable Items
- <Important detail worth highlighting>
- <Another notable item>

### TL;DR
<One sentence summary>
```

For code:
```markdown
## Summary: <filename>

### Purpose
<What this code does>

### Main Components
- `ComponentName`: <brief description>
- `FunctionName`: <brief description>

### Dependencies
- <key dependencies>

### Complexity
<Simple/Moderate/Complex> - <brief justification>
```

Content to summarize: $ARGUMENTS
