# Claude Skills

A comprehensive CLI tool for managing AI assistant skills and slash commands for Claude Code and other LLM platforms.

## Features

- **Install/Uninstall**: Easily manage skills and commands
- **Search**: Find skills by keywords and descriptions
- **Create**: Interactive wizard for creating new skills and commands
- **Validate**: Ensure skills and commands meet format requirements
- **Multi-platform**: Support for Claude, Codex, and Gemini (coming soon)

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
claude-skills config set repository /path/to/claude-skills
```

### List Available Skills

```bash
# List all skills and commands
claude-skills list

# List only skills
claude-skills list --type skill

# List only from repository
claude-skills list --source repository
```

### Install a Skill

```bash
# Install to personal directory (~/.claude/skills/)
claude-skills install skill-name

# Install to project directory (.claude/skills/)
claude-skills install skill-name --target project

# Preview before installing
claude-skills install skill-name --dry-run
```

### Search for Skills

```bash
claude-skills search python
claude-skills search "lab notebook"
```

### Create a New Skill

```bash
# Interactive mode
claude-skills create --type skill

# With name specified
claude-skills create --type skill --name my-awesome-skill
```

### Uninstall a Skill

```bash
claude-skills uninstall skill-name
```

## Configuration

Configuration is stored in `~/.config/claude-skills/config.yaml`.

Default configuration:

```yaml
personal_skills_dir: ~/.claude/skills
personal_commands_dir: ~/.claude/commands
project_skills_dir: .claude/skills
project_commands_dir: .claude/commands
default_target: personal

platforms:
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
claude-skills/
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

# Run tests
pytest

# Run linting
ruff check .
black --check .

# Format code
black .
```

### Running Tests

```bash
pytest
pytest --cov=cli --cov-report=html
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

- Report issues: https://github.com/jkitchin/claude-skills/issues
- Documentation: https://github.com/jkitchin/claude-skills#readme

## Roadmap

- [ ] Phase 1: Foundation (validators, templates) ✓
- [ ] Phase 2: Core CLI commands ✓
- [ ] Phase 3: Initial skills library
- [ ] Phase 4: Advanced features (update, search improvements)
- [ ] Phase 5: Community contributions, PyPI release

## Credits

Created by John Kitchin for managing AI assistant skills across different LLM platforms.
