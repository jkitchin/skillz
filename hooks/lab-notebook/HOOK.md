---
name: lab-notebook
description: Generate lab notebook entries from Claude Code sessions for documentation and reproducibility
event: SessionEnd
type: command
timeout: 30
---

# Lab Notebook

Generate structured lab notebook entries from Claude Code sessions. Perfect for scientists, researchers, and anyone who wants to maintain a record of their AI-assisted work.

## Purpose

This hook automatically creates a lab notebook entry when a Claude Code session ends. It captures:

- Session metadata (ID, timestamp, project)
- User prompts and objectives
- Files modified during the session
- Commands executed
- Git status (if in a repository)

## Output Formats

The hook supports multiple output formats:

- **Markdown** (default) - Human-readable, great for documentation
- **JSON** - Structured data for tooling and analysis
- **Org-mode** - For Emacs users

## Configuration

Set environment variables to customize behavior:

```bash
export LAB_NOTEBOOK_DIR="~/lab-notebook"     # Output directory
export LAB_NOTEBOOK_FORMAT="markdown"         # markdown, json, or org
export LAB_NOTEBOOK_GIT_INFO="true"          # Include git status
```

## Output Location

Entries are saved to:
- `~/lab-notebook/claude-sessions/YYYY-MM-DD_SESSION-ID.md`

## Example Output

```markdown
# Lab Notebook Entry

**Date**: 2024-01-15 14:30
**Session**: `abc123def`
**Project**: `/home/user/myproject`
**Git Branch**: `feature/new-api`

## Objective

Implement user authentication for the REST API

## Work Performed

### User Prompts
1. Add JWT authentication to the Flask API
2. Create tests for the auth endpoints

### Files Modified
- `src/auth.py`
- `src/routes/users.py`
- `tests/test_auth.py`

### Commands Executed
```bash
pytest tests/test_auth.py
git diff
```

## Git Status

- Branch: feature/new-api
- Uncommitted changes: 3 files
```

## Requirements

- Python 3.8+
- No external dependencies (uses standard library only)

## Behavior

1. Reads the session transcript from `transcript_path`
2. Parses user prompts, tool uses, and Claude responses
3. Generates a formatted notebook entry
4. Saves to the configured output directory

## Use Cases

- **Research Documentation**: Track AI-assisted experiments
- **Code Review**: Document what changes were made and why
- **Learning**: Review past sessions to learn from them
- **Compliance**: Maintain audit trails of AI-assisted work
- **Collaboration**: Share session summaries with teammates
