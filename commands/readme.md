---
description: Generate or update README content
allowed-tools: ["Read", "Glob", "Grep", "Bash"]
argument-hint: <section-name>
---

# Generate README

Generate or update README content for a project.

## Instructions

1. Analyze the project structure:
   ```bash
   ls -la
   find . -maxdepth 2 -type f -name "*.py" -o -name "*.js" -o -name "*.ts" -o -name "*.go" 2>/dev/null | head -20
   ```

2. Check for existing README:
   - Read README.md or README if it exists
   - Identify what sections need updating

3. Gather project information:
   - Package.json, pyproject.toml, Cargo.toml, go.mod for metadata
   - Main entry points and key files
   - Dependencies and requirements

4. Generate/update the requested section or full README: $ARGUMENTS
   - If no section specified, generate a complete README
   - Available sections: overview, installation, usage, api, configuration, contributing, license

## Output Format

```markdown
# Project Name

Brief description of what the project does.

## Installation

```bash
<installation commands>
```

## Quick Start

```<language>
<basic usage example>
```

## Features

- Feature 1
- Feature 2

## Usage

<Detailed usage instructions>

## Configuration

<Configuration options if applicable>

## API Reference

<If applicable>

## Contributing

<Contribution guidelines>

## License

<License information>
```

Section to generate (or leave empty for full README): $ARGUMENTS
