---
name: black-on-save
description: Auto-format Python files with Black after edits
event: PostToolUse
matcher: Edit|Write
type: command
timeout: 30
---

# Black on Save

Automatically format Python files with Black after Claude edits or creates them.

## Purpose

Ensures consistent Python code formatting using the Black code formatter. When Claude writes or edits a Python file, this hook automatically formats it.

## Requirements

- Python 3.8+
- Black installed (`pip install black`)

## Configuration

The hook respects your project's Black configuration:
- `pyproject.toml` (tool.black section)
- `setup.cfg`
- `.black.toml`

## Behavior

1. Hook receives the file path from the Edit/Write tool
2. Checks if the file is a Python file (`.py`)
3. Runs `black` on the file
4. Silent success, reports errors to stderr

## Notes

- Does not block the operation (PostToolUse)
- Fails silently if Black is not installed
- Also runs isort for import sorting if available
