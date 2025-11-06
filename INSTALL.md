# Installation Guide

Complete guide to installing and setting up Claude Skills CLI.

## Prerequisites

- Python 3.8 or higher
- Git (for cloning the repository)
- `uv` (recommended) or `pip`

## Installation Methods

### Method 1: Using uv (Recommended)

[uv](https://github.com/astral-sh/uv) is a fast Python package installer.

```bash
# Install uv if you haven't already
curl -LsSf https://astral.sh/uv/install.sh | sh

# Clone the repository
git clone https://github.com/jkitchin/claude-skills.git
cd claude-skills

# Install in editable mode with dependencies
uv pip install -e .

# For development, install with dev dependencies
uv pip install -e ".[dev]"
```

### Method 2: Using pip

```bash
# Clone the repository
git clone https://github.com/jkitchin/claude-skills.git
cd claude-skills

# Create a virtual environment (recommended)
python -m venv venv
source venv/bin/activate  # On Windows: venv\Scripts\activate

# Install in editable mode
pip install -e .

# For development
pip install -e ".[dev]"
```

### Method 3: From PyPI (Coming Soon)

Once published to PyPI:

```bash
pip install claude-skills
# or
uv pip install claude-skills
```

## Configuration

### 1. Configure Repository Path

Set up the path to your local claude-skills repository:

```bash
# Create config directory
mkdir -p ~/.config/claude-skills

# Copy example config
cp example-config.yaml ~/.config/claude-skills/config.yaml

# Edit config file
# Set repository_path to your local clone path
```

Or create the config file manually:

```yaml
repository_path: /path/to/your/claude-skills
default_target: personal
```

### 2. Verify Installation

```bash
# Check version
python -m cli.main --version

# List available skills
python -m cli.main list

# View help
python -m cli.main --help
```

## Quick Start

### Installing Skills

```bash
# Install a skill to personal directory (~/.claude/skills/)
python -m cli.main install python-ase

# Install to project directory (.claude/skills/)
python -m cli.main install python-ase --target project

# Preview before installing
python -m cli.main install python-ase --dry-run
```

### Searching and Browsing

```bash
# Search for skills
python -m cli.main search python

# List all skills
python -m cli.main list --type skill

# Get detailed information
python -m cli.main info python-ase
```

### Managing Skills

```bash
# Uninstall a skill
python -m cli.main uninstall python-ase

# Update installed skills
python -m cli.main update

# Create a new skill
python -m cli.main create --type skill
```

## Platform-Specific Setup

### For Claude Code

Skills are automatically detected from:
- Personal: `~/.claude/skills/`
- Project: `.claude/skills/` (in your project directory)

Commands are detected from:
- Personal: `~/.claude/commands/`
- Project: `.claude/commands/`

### For Other Platforms (Codex, Gemini)

Configure platform-specific paths in your config file:

```yaml
platforms:
  codex:
    skills_dir: ~/.config/openai/skills
    commands_dir: ~/.config/openai/commands
```

Then install with platform flag:

```bash
python -m cli.main install python-ase --platform codex
```

## Troubleshooting

### Command not found

If `claude-skills` command is not found, use the full module path:

```bash
python -m cli.main [command]
```

Or add to your PATH:

```bash
# Add to ~/.bashrc or ~/.zshrc
export PATH="$PATH:$HOME/.local/bin"
```

### No items found

If `list` shows no items, check your config:

```bash
# Verify repository path is set
python -m cli.main config get repository_path

# Set it if needed
python -m cli.main config set repository_path /path/to/claude-skills
```

### Permission errors

If you get permission errors:

```bash
# Install to user site-packages
pip install --user -e .

# Or use a virtual environment (recommended)
python -m venv venv
source venv/bin/activate
pip install -e .
```

### Import errors

Make sure you're in the right directory:

```bash
cd /path/to/claude-skills
pip install -e .
```

## Uninstallation

```bash
# Remove installed skills
python -m cli.main uninstall --all

# Uninstall package
pip uninstall claude-skills

# Remove config (optional)
rm -rf ~/.config/claude-skills
```

## Next Steps

- Read [CONTRIBUTING.md](CONTRIBUTING.md) to create your own skills
- Browse available skills in [skills/](skills/)
- Check out example commands in [commands/examples/](commands/examples/)
- Read the full documentation in [README.md](README.md)

## Getting Help

- Check [GitHub Issues](https://github.com/jkitchin/claude-skills/issues)
- Read the [FAQ](#) (coming soon)
- Contact: jkitchin@andrew.cmu.edu
