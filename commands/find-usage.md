---
description: Find where a function or class is used
allowed-tools: ["Grep", "Glob", "Read"]
argument-hint: <function-or-class-name>
---

# Find Usage

Find all usages of a function, class, variable, or other identifier in the codebase.

## Instructions

1. Search for the identifier specified: $ARGUMENTS

2. Find:
   - Definition location
   - All import statements
   - All call sites / usages
   - Test files that use it

3. Categorize by usage type:
   - Direct calls
   - Inheritance/implementation
   - Type annotations
   - String references (e.g., in configs)

4. Provide context around each usage

## Output Format

```markdown
## Usage Report: `<identifier>`

### Definition
**File**: `<path>`
**Line**: <number>
**Type**: <function/class/variable/constant>

```<language>
<definition code>
```

### Summary
- Total usages: <number>
- Files: <number>
- Test files: <number>

### Imports

| File | Import Statement |
|------|------------------|
| `file1.py` | `from module import identifier` |
| `file2.py` | `import module` |

### Usages

#### `<filename>` (<count> usages)

**Line <N>**:
```<language>
<usage context>
```
<Brief explanation of how it's used>

**Line <M>**:
```<language>
<usage context>
```

#### `<another-file>` (<count> usages)
...

### Test Coverage
| Test File | Test Function | Type |
|-----------|--------------|------|
| `test_module.py` | `test_function_works` | Unit |

### Dependency Graph
```
identifier
├── used by: function_a (file1.py)
│   └── used by: handler (file2.py)
└── used by: function_b (file3.py)
```

### Safe to Modify?
<Assessment of refactoring impact>
```

Identifier to find: $ARGUMENTS
