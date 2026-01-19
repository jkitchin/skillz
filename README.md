# Skillz

<img src="skillz.png" alt="Skillz" width="300">

[![Tests](https://github.com/jkitchin/skillz/actions/workflows/test.yml/badge.svg)](https://github.com/jkitchin/skillz/actions/workflows/test.yml)
[![Lint](https://github.com/jkitchin/skillz/actions/workflows/lint.yml/badge.svg)](https://github.com/jkitchin/skillz/actions/workflows/lint.yml)
[![codecov](https://codecov.io/gh/jkitchin/skillz/graph/badge.svg)](https://codecov.io/gh/jkitchin/skillz)


A comprehensive CLI tool for managing AI assistant skills and slash commands for Claude Code, OpenCode, and other LLM platforms.

> **Note:** Currently only tested on Claude Code and OpenCode. Please report issues for Codex and Gemini. With opencode you may have better performance with https://www.npmjs.com/package/opencode-skills (this is what I use).

- https://github.com/numman-ali/openskills (I have not tried this yet)
- Requested feature in opencode https://github.com/sst/opencode/issues/3235

## Features

- **Install/Uninstall**: Easily manage skills, commands, hooks, and agents
- **Search**: Find skills by keywords and descriptions
- **Create**: Interactive wizard or AI-assisted creation for new skills, commands, hooks, and agents
- **Validate**: Ensure skills, commands, hooks, and agents meet format requirements
- **Multi-platform**: Support for OpenCode, Claude Code, Codex, and Gemini
- **Hooks**: Lifecycle hooks for Claude Code (formatting, logging, security, notifications)
- **Agents**: Specialized subagents for code review, debugging, testing, and documentation

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
# Install to Claude Code (default)
skillz install skill-name

# Install to OpenCode
skillz install skill-name --platform opencode

# Install to project directory (.claude/skills/)
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

# AI-assisted creation using Claude CLI
skillz create --type skill --name my-skill --prompt "A skill for analyzing Python code style"
```

### Uninstall a Skill

```bash
skillz uninstall skill-name
```

### Manage Hooks

```bash
# List available hooks
skillz hooks list --target repo

# Install a hook
skillz hooks install lab-notebook

# Create a hook (with AI)
skillz hooks create my-hook --event PostToolUse --prompt "Format Python with black"
```

### Manage Agents

```bash
# List available agents
skillz agents list --target repo

# Install an agent
skillz agents install code-reviewer

# Create an agent (with AI)
skillz agents create my-agent --prompt "An agent that reviews security vulnerabilities"
```

## Configuration

Configuration is stored in `~/.config/skillz/config.yaml`.

Default configuration:

```yaml
default_platform: claude

personal_skills_dir: ~/.claude/skills
personal_commands_dir: ~/.claude/commands
project_skills_dir: .claude/skills
project_commands_dir: .claude/commands
default_target: personal

platforms:
  claude:
    skills_dir: ~/.claude/skills
    commands_dir: ~/.claude/commands
  opencode:
    skills_dir: ~/.config/opencode/skills
    commands_dir: ~/.config/opencode/command
  codex:
    skills_dir: ~/.codex/skills
    commands_dir: ~/.codex/commands
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

## Hooks Format

Hooks are directories containing a `HOOK.md` file and a script (`hook.py` or `hook.sh`):

```markdown
---
name: my-hook
description: Brief description of what this hook does
event: PostToolUse  # PreToolUse, PostToolUse, Stop, SessionStart, SessionEnd, etc.
matcher: Write      # Optional: filter by tool name
timeout: 30000      # Optional: timeout in milliseconds
---

# Hook Documentation

Details about when and how this hook runs...
```

### Available Hook Events

- `PreToolUse` / `PostToolUse`: Before/after tool execution
- `UserPromptSubmit`: When user submits a prompt
- `Stop`: When Claude needs user input
- `SessionStart` / `SessionEnd`: Session lifecycle
- `Notification`: For notifications

## Agents Format

Agents are markdown files with YAML frontmatter:

```markdown
---
name: my-agent
description: Clear description of what this agent does
tools: Read, Write, Edit, Grep, Glob
model: sonnet  # sonnet, opus, or haiku
---

# Agent Instructions

Detailed instructions for this specialized agent...
```

## Project Structure

```
skillz/
├── cli/                    # Python CLI tool
│   ├── commands/          # CLI command implementations
│   ├── config.py          # Configuration management
│   ├── validator.py       # Skill/command/hook/agent validation
│   └── utils.py           # Helper functions
├── skills/                # Skill repository
│   ├── academic/
│   ├── programming/
│   └── research/
├── commands/              # Command repository
├── hooks/                 # Hooks repository
│   ├── lab-notebook/
│   ├── prettier-on-save/
│   └── ...
├── agents/                # Agents repository
│   ├── code-reviewer.md
│   ├── debugger.md
│   └── ...
├── templates/             # Templates for new skills/hooks/agents
├── docs/                  # Documentation
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

## Available Hooks

See the `hooks/` directory for available hooks:

- **lab-notebook**: Generate lab notebook entries from Claude Code sessions
- **prettier-on-save**: Auto-format JS/TS/CSS/JSON files with Prettier
- **black-on-save**: Auto-format Python files with Black
- **protect-secrets**: Block writes to sensitive files like .env
- **bash-logger**: Log all bash commands for audit
- **notify-done**: Desktop notification when Claude needs input

## Available Agents

See the `agents/` directory for available agents:

- **code-reviewer**: Analyzes code for quality, security issues, and best practices
- **debugger**: Traces issues through code, analyzes errors, suggests fixes
- **literature-searcher**: Searches academic papers and technical resources
- **test-writer**: Generates comprehensive test cases and suites
- **doc-writer**: Creates documentation, READMEs, and API docs

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

## Related Projects

### Supported Platforms

Skillz manages skills and commands for these AI coding assistants:

- **[OpenCode](https://github.com/sst/opencode)** - The AI coding agent built for the terminal. Open-source, provider-agnostic, with 30k+ GitHub stars.
- **[Claude Code](https://github.com/anthropics/claude-code)** - Anthropic's agentic coding tool that lives in your terminal. Understands your codebase and helps you code faster.
- **[Codex CLI](https://github.com/openai/codex)** - OpenAI's lightweight coding agent that runs in your terminal. See [OpenAI adopts "skills"](https://simonwillison.net/2025/Dec/12/openai-skills/) for details on their skills implementation.
- **Gemini** - Google's AI assistant.

### Related Skills Projects

- **[Superpowers](https://github.com/obra/superpowers)** - A complete software development workflow framework for Claude Code with composable skills for design, TDD, code review, and integration. This project was an early motivation to develop this project.

## Credits

Created by John Kitchin for managing AI assistant skills across different LLM platforms.
