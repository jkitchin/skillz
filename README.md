# Skillz

<img src="skillz.png" alt="Skillz" width="300">

[![Tests](https://github.com/jkitchin/skillz/actions/workflows/test.yml/badge.svg)](https://github.com/jkitchin/skillz/actions/workflows/test.yml)
[![Lint](https://github.com/jkitchin/skillz/actions/workflows/lint.yml/badge.svg)](https://github.com/jkitchin/skillz/actions/workflows/lint.yml)
[![codecov](https://codecov.io/gh/jkitchin/skillz/graph/badge.svg)](https://codecov.io/gh/jkitchin/skillz)


A comprehensive CLI tool for managing AI assistant skills and slash commands for OpenCode, Claude Code, and other LLM platforms.

> **Note:** Currently only tested on Claude Code and OpenCode.

## Features

- **Install/Uninstall**: Easily manage skills and commands
- **Search**: Find skills by keywords and descriptions
- **Create**: Interactive wizard for creating new skills and commands
- **Validate**: Ensure skills and commands meet format requirements
- **Multi-platform**: Support for OpenCode, Claude Code, Codex, and Gemini

## Installation

### Using uv (recommended)

```bash
uv pip install -e .
```

### Using pip

```bash
pip install -e .
```

## Quick Start

### Configure Repository Path

First, set up the path to your skills repository:

```bash
skillz config set repository /path/to/skillz
```

### List Available Skills

```bash
# List all skills and commands
skillz list

# List only skills
skillz list --type skill

# List only from repository
skillz list --source repository
```

### Install a Skill

```bash
# Install to OpenCode (default)
skillz install skill-name

# Install to Claude Code
skillz install skill-name --platform claude

# Install to project directory (.opencode/skills/)
skillz install skill-name --target project

# Preview before installing
skillz install skill-name --dry-run
```

### Search for Skills

```bash
skillz search python
skillz search "lab notebook"
```

### Create a New Skill

```bash
# Interactive mode
skillz create --type skill

# With name specified
skillz create --type skill --name my-awesome-skill
```

### Uninstall a Skill

```bash
skillz uninstall skill-name
```

## Configuration

Configuration is stored in `~/.config/skillz/config.yaml`.

Default configuration:

```yaml
default_platform: opencode

personal_skills_dir: ~/.config/opencode/skills
personal_commands_dir: ~/.config/opencode/command
project_skills_dir: .opencode/skills
project_commands_dir: .opencode/command
default_target: personal

platforms:
  opencode:
    skills_dir: ~/.config/opencode/skills
    commands_dir: ~/.config/opencode/command
  claude:
    skills_dir: ~/.claude/skills
    commands_dir: ~/.claude/commands
  codex:
    skills_dir: ~/.config/openai/skills
    commands_dir: ~/.config/openai/commands
  gemini:
    skills_dir: ~/.config/gemini/skills
    commands_dir: ~/.config/gemini/commands
```

## Skills Format

Skills are directories containing a `SKILL.md` file with YAML frontmatter:

```markdown
---
name: my-skill
description: A clear description of what this skill does and when to use it
allowed-tools: ["*"]  # Optional: restrict available tools
---

# Skill Content

Detailed instructions for Claude on how to use this skill...
```

### Requirements

- `name`: lowercase, hyphens only, max 64 characters
- `description`: clear explanation, max 1024 characters
- Directory must contain `SKILL.md`

## Commands Format

Commands are markdown files with optional YAML frontmatter:

```markdown
---
description: Brief description for /help
model: sonnet  # Optional: sonnet, opus, or haiku
allowed-tools: ["*"]  # Optional
argument-hint: <your-arg>  # Optional: autocomplete help
---

# Command Content

Command prompt content here...

Use $ARGUMENTS or $1, $2, etc. for parameters.
```

## Project Structure

```
skillz/
├── cli/                    # Python CLI tool
│   ├── commands/          # CLI command implementations
│   ├── config.py          # Configuration management
│   ├── validator.py       # Skill/command validation
│   └── utils.py           # Helper functions
├── skills/                # Skill repository
│   ├── academic/
│   ├── programming/
│   └── research/
├── commands/              # Command repository
├── templates/             # Templates for new skills
└── tests/                 # Test suite
```

## Development

### Setup Development Environment

```bash
# Install with development dependencies
uv pip install -e ".[dev]"

# Install pre-commit hooks
pre-commit install
```

### Pre-commit Hooks

This project uses [pre-commit](https://pre-commit.com/) to automatically run linting and formatting checks before each commit. The hooks will:

- Run `ruff check --fix` to fix linting issues
- Run `ruff format` to format code

If any changes are made by the hooks, the commit will fail. Simply stage the changes and commit again.

### Running Tests

```bash
# Run all tests
pytest

# Run with coverage report
pytest --cov=cli --cov-report=html
```

### Linting and Formatting

We use [Ruff](https://docs.astral.sh/ruff/) for both linting and formatting:

```bash
# Check for linting issues
ruff check .

# Auto-fix linting issues
ruff check --fix .

# Check formatting
ruff format --check .

# Format code
ruff format .
```

## Available Skills

See the `skills/` directory for available skills organized by category:

- **Academic**: PhD qualifier review, proposal assistance
- **Programming**: Python (ASE, pymatgen), Emacs Lisp
- **Research**: Electronic lab notebook management

## Contributing

Contributions are welcome! Please:

1. Follow the skill/command format requirements
2. Add tests for new functionality
3. Update documentation
4. Run linters before submitting

## License

MIT License - see LICENSE file for details

## Support

- Report issues: https://github.com/jkitchin/skillz/issues
- Documentation: https://github.com/jkitchin/skillz#readme

## Roadmap

- [ ] Phase 1: Foundation (validators, templates) ✓
- [ ] Phase 2: Core CLI commands ✓
- [ ] Phase 3: Initial skills library
- [ ] Phase 4: Advanced features (update, search improvements)
- [ ] Phase 5: Community contributions, PyPI release

## Credits

Created by John Kitchin for managing AI assistant skills across different LLM platforms.
