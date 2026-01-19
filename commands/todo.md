---
description: Extract TODOs from the codebase
allowed-tools: ["Grep", "Glob", "Read"]
argument-hint: <path or filter>
---

# Extract TODOs

Find and organize all TODO, FIXME, HACK, and similar comments in the codebase.

## Instructions

1. Search for TODO markers in the codebase:
   - TODO
   - FIXME
   - HACK
   - XXX
   - BUG
   - NOTE (optionally)

2. If a path is specified ($ARGUMENTS), search only there
   Otherwise, search the entire project (excluding node_modules, venv, etc.)

3. Organize findings by:
   - Priority (FIXME > BUG > TODO > HACK > NOTE)
   - File/location
   - Category if discernible

4. Extract context around each TODO

## Output Format

```markdown
## TODO Summary

**Total found**: <number>
- FIXME: <count>
- BUG: <count>
- TODO: <count>
- HACK: <count>

### High Priority (FIXME/BUG)

#### <filename>:<line>
```<language>
// FIXME: <the comment>
<surrounding context>
```

### Medium Priority (TODO)

#### <filename>:<line>
```<language>
// TODO: <the comment>
<surrounding context>
```

### Low Priority (HACK/XXX)

#### <filename>:<line>
```<language>
// HACK: <the comment>
<surrounding context>
```

### By Category

#### Authentication
- `file.py:42` - TODO: Add token refresh

#### Database
- `db.py:100` - FIXME: Connection leak

### Statistics
| Marker | Count | Files Affected |
|--------|-------|----------------|
| TODO | 10 | 5 |
| FIXME | 3 | 2 |
```

Path or filter (optional): $ARGUMENTS
