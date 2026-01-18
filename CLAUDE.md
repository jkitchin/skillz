# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

Skillz is a CLI tool for managing AI assistant skills and slash commands across multiple LLM platforms (Claude Code, OpenCode, Codex, Gemini). The tool allows users to discover, install, create, and manage reusable skills that extend LLM capabilities.

## Development Commands

### Setup
```bash
# Install with development dependencies (preferred)
uv pip install -e ".[dev]"

# Alternative with pip
pip install -e ".[dev]"
```

### Testing
```bash
# Run all tests
pytest

# Run with coverage report
pytest --cov=cli --cov-report=html

# Run specific test file
pytest tests/test_validator.py
```

### Linting and Formatting
```bash
# Check code style
ruff check .
black --check .

# Format code
black .
```

### Running the CLI
```bash
# Via entry point (after installation)
skillz --help

# Via Python module (during development)
python -m cli.main --help
```

## Architecture

### Core Components

**CLI Module (`cli/`)**
- `main.py` - Click-based CLI entry point with command registration
- `config.py` - Configuration management using YAML (`~/.config/skillz/config.yaml`)
- `validator.py` - Validation logic for skills (SKILL.md) and commands (.md files)
- `utils.py` - Shared utilities (name validation, file operations, search functions)
- `commands/` - Individual CLI command implementations (install, uninstall, list, search, info, update, create)

**Skills Repository (`skills/`)**
Each skill is a directory containing:
- `SKILL.md` - Required file with YAML frontmatter (name, description, allowed-tools) and markdown content
- Supporting files: README.md, QUICK_REFERENCE.md, examples/, references/, assets/, scripts/

**Commands Repository (`commands/`)**
Standalone markdown files with optional YAML frontmatter (description, model, allowed-tools, argument-hint).

**Templates (`templates/`)**
- `SKILL_TEMPLATE.md` - Template for creating new skills
- `COMMAND_TEMPLATE.md` - Template for creating new commands

### Configuration Flow

1. Config loads from `~/.config/skillz/config.yaml` or uses defaults
2. Config provides paths for personal (default: `~/.claude/`) and project (default: `.claude/`) installations
3. Multi-platform support via `platforms` dict in config (claude, opencode, codex, gemini)
4. Default platform is Claude Code, but fully supports OpenCode, Codex, and Gemini

### Validation Logic

**Skills**: Must have valid name (lowercase, hyphens, max 64 chars), description (max 1024 chars), and SKILL.md file with proper YAML frontmatter.

**Commands**: Optional frontmatter, description max 256 chars, valid model values (sonnet/opus/haiku).

**Allowed Tools**: Can be `["*"]` for all tools or list from: Bash, Read, Write, Edit, Glob, Grep, Task, WebFetch, WebSearch, TodoWrite, AskUserQuestion, Skill, SlashCommand, NotebookEdit, BashOutput, KillShell.

### Discovery and Installation

- `find_skill_directories()` and `find_command_files()` recursively scan repository
- Installation copies from repository to personal/project directories
- Search uses fuzzy matching on names, descriptions, and file contents

## Key Conventions

- **Naming**: Skills and commands use lowercase-with-hyphens naming
- **Frontmatter**: YAML between `---` delimiters at file start
- **Tools restriction**: `allowed-tools` field limits which Claude tools the skill/command can use
- **Repository path**: Must be configured via `skillz config set repository /path/to/repo`

## Install --all Feature

The `install --all` command installs all skills and commands from the repository in a single operation:
- Discovers all valid skills and commands from repository
- Validates each before installation
- Provides progress feedback during batch installation
- Supports `--dry-run` to preview what would be installed
- Respects `--target` (personal/project) and `--force` options

## Common Patterns

**Adding a new CLI command**: Create file in `cli/commands/`, implement as Click command, import and register in `cli/main.py`.

**Creating a skill**: Use `claude-skills create --type skill` or manually create directory with SKILL.md containing proper frontmatter.

**Validation errors**: Both SkillValidator and CommandValidator return `(is_valid, list_of_errors)` tuples.

## Git Policy

Never use `--no-verify` with git commands per user's global instructions.
