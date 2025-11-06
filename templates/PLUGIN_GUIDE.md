# Plugin Creation Guide

This guide explains how to create and distribute Claude Code plugins that bundle skills and commands together.

## Plugin Structure

A complete plugin has this directory structure:

```
my-plugin/
├── plugin.json              # Required: Plugin manifest
├── README.md               # Recommended: Plugin documentation
├── LICENSE                 # Recommended: License file
├── skills/                 # Skills directory (auto-discovered)
│   ├── skill-1/
│   │   ├── SKILL.md       # Required for each skill
│   │   ├── reference.md   # Optional
│   │   └── scripts/       # Optional
│   └── skill-2/
│       └── SKILL.md
├── commands/               # Commands directory
│   ├── command1.md
│   ├── command2.md
│   └── subcategory/
│       └── command3.md
└── agents/                 # Optional: Custom agents
    └── agent-config.json
```

---

## plugin.json Manifest

The `plugin.json` file is the core configuration for your plugin.

### Required Fields

```json
{
  "name": "plugin-name",
  "version": "1.0.0",
  "description": "Brief description of what this plugin provides"
}
```

### Complete Example

```json
{
  "$schema": "https://docs.claude.com/schemas/plugin-manifest.json",
  "name": "academic-tools",
  "version": "1.2.0",
  "description": "Tools for academic research and writing",
  "author": "Jane Doe <jane@university.edu>",
  "license": "MIT",
  "homepage": "https://github.com/janedoe/academic-tools",
  "repository": {
    "type": "git",
    "url": "https://github.com/janedoe/academic-tools.git"
  },
  "keywords": [
    "academic",
    "research",
    "writing",
    "citations"
  ],
  "commands": [
    "./commands/cite.md",
    "./commands/format-bib.md",
    "./commands/review.md"
  ],
  "agents": [
    "./agents/paper-reviewer"
  ],
  "dependencies": {
    "required": [
      "python>=3.8",
      "bibtexparser>=1.4.0"
    ],
    "optional": [
      "pandoc>=2.0"
    ]
  },
  "configuration": {
    "settings": {
      "citation_style": {
        "type": "string",
        "default": "apa",
        "description": "Default citation style (apa, mla, chicago)",
        "enum": ["apa", "mla", "chicago"]
      },
      "bibliography_path": {
        "type": "string",
        "default": "./references.bib",
        "description": "Path to default bibliography file"
      }
    }
  },
  "readme": "./README.md",
  "changelog": "./CHANGELOG.md"
}
```

---

## Manifest Field Reference

### Core Fields

| Field | Type | Required | Description |
|-------|------|----------|-------------|
| `name` | string | Yes | Plugin identifier (lowercase, hyphens) |
| `version` | string | Yes | Semantic version (e.g., "1.2.3") |
| `description` | string | Yes | Brief description of plugin functionality |
| `author` | string | No | Author name and email |
| `license` | string | No | License type (e.g., "MIT", "Apache-2.0") |
| `homepage` | string | No | Plugin homepage URL |
| `repository` | object | No | Repository information |
| `keywords` | array | No | Search keywords for discoverability |

### Component Fields

| Field | Type | Description |
|-------|------|-------------|
| `commands` | array | Paths to command markdown files |
| `agents` | array | Paths to agent configuration directories |
| `skills` | N/A | Auto-discovered from `skills/` directory |

**Note:** Skills are automatically discovered from the `skills/` directory and don't need to be listed in the manifest.

### Dependency Fields

```json
{
  "dependencies": {
    "required": [
      "python>=3.8",
      "package-name>=1.0.0",
      "system-tool"
    ],
    "optional": [
      "optional-package>=2.0.0"
    ]
  }
}
```

- **required**: Must be installed for plugin to work
- **optional**: Enhance functionality but not required

### Configuration Fields

```json
{
  "configuration": {
    "settings": {
      "setting_name": {
        "type": "string|number|boolean|array|object",
        "default": "default-value",
        "description": "What this setting does",
        "enum": ["option1", "option2"],
        "min": 0,
        "max": 100,
        "required": false
      }
    }
  }
}
```

**Setting Types:**
- `string`: Text values
- `number`: Numeric values (use `min`/`max` for range)
- `boolean`: true/false
- `array`: List of values
- `object`: Structured data
- `enum`: Restricted to specific values

---

## Skills in Plugins

Skills are automatically discovered from the `skills/` directory.

### Structure

```
skills/
├── pdf-processor/
│   ├── SKILL.md           # Required
│   ├── reference.md       # Optional: detailed reference
│   ├── examples.md        # Optional: usage examples
│   └── scripts/           # Optional: helper scripts
│       └── extract.py
└── citation-helper/
    └── SKILL.md
```

### Best Practices

1. **One skill per subdirectory**
2. **Clear naming**: Use descriptive skill names
3. **Complete SKILL.md**: Include all required metadata
4. **Include examples**: Help users understand usage
5. **Document dependencies**: List required packages in description

See `SKILL_TEMPLATE.md` for detailed skill creation guide.

---

## Commands in Plugins

Commands must be explicitly listed in the `commands` array.

### Manifest Entry

```json
{
  "commands": [
    "./commands/deploy.md",
    "./commands/test.md",
    "./commands/utils/format.md"
  ]
}
```

### Structure

```
commands/
├── deploy.md
├── test.md
└── utils/
    └── format.md
```

### Best Practices

1. **Explicit listing**: Always list commands in manifest
2. **Relative paths**: Use `./` prefix for clarity
3. **Organization**: Use subdirectories for categorization
4. **Clear descriptions**: Help users discover commands

See `COMMAND_TEMPLATE.md` for detailed command creation guide.

---

## Agents in Plugins

Custom agents extend Claude's capabilities with specialized behaviors.

### Manifest Entry

```json
{
  "agents": [
    "./agents/code-reviewer",
    "./agents/test-generator"
  ]
}
```

### Structure

```
agents/
├── code-reviewer/
│   ├── config.json        # Agent configuration
│   ├── instructions.md    # Agent instructions
│   └── tools/            # Custom tools
└── test-generator/
    └── config.json
```

---

## Plugin Versioning

Follow [Semantic Versioning](https://semver.org/):

- **Major** (1.0.0 → 2.0.0): Breaking changes
- **Minor** (1.0.0 → 1.1.0): New features, backward compatible
- **Patch** (1.0.0 → 1.0.1): Bug fixes, backward compatible

### Example Changelog

```markdown
# Changelog

## [1.2.0] - 2024-02-15
### Added
- New skill for table extraction
- Command for batch PDF processing

### Changed
- Improved text extraction accuracy

### Fixed
- Bug in metadata parsing

## [1.1.0] - 2024-01-10
### Added
- Initial release
```

---

## Plugin Documentation

### README.md Template

```markdown
# Plugin Name

Brief description of what the plugin does.

## Features

- Feature 1
- Feature 2
- Feature 3

## Installation

```bash
# Installation instructions
```

## Skills

### skill-name-1
Description of what this skill does and when to use it.

### skill-name-2
Description of what this skill does and when to use it.

## Commands

### /command-1
Description and usage example.

### /command-2
Description and usage example.

## Configuration

Available settings and their defaults:

- `setting_name`: Description (default: `value`)

## Dependencies

Required:
- Package 1
- Package 2

Optional:
- Optional package 1

## Examples

### Example 1: Common Use Case
```bash
/command example-arg
```

### Example 2: Another Use Case
Ask Claude to use the skill...

## Troubleshooting

Common issues and solutions.

## License

License information.

## Contributing

Contribution guidelines.
```

---

## Testing Your Plugin

### Local Testing

1. **Create plugin structure:**
   ```bash
   mkdir -p my-plugin/{skills,commands,agents}
   cd my-plugin
   ```

2. **Add plugin.json:**
   ```bash
   cp /path/to/PLUGIN_TEMPLATE.json ./plugin.json
   # Edit plugin.json with your details
   ```

3. **Add skills and commands:**
   ```bash
   # Create a skill
   mkdir -p skills/test-skill
   cp /path/to/SKILL_TEMPLATE.md skills/test-skill/SKILL.md

   # Create a command
   cp /path/to/COMMAND_TEMPLATE.md commands/test-command.md
   ```

4. **Test locally:**
   - Copy plugin to appropriate location
   - Restart Claude Code
   - Verify skills and commands are available

### Validation

1. **JSON validation:**
   ```bash
   # Validate plugin.json syntax
   python -m json.tool plugin.json
   ```

2. **Skill validation:**
   - Check YAML frontmatter in SKILL.md files
   - Verify required fields (name, description)
   - Test description clarity

3. **Command validation:**
   - Check YAML frontmatter syntax
   - Verify parameter usage ($1, $2, $ARGUMENTS)
   - Test with various arguments

---

## Distribution

### Method 1: Git Repository

1. Create git repository
2. Push plugin code
3. Users clone repository
4. Users install to Claude Code

```bash
# Users run:
git clone https://github.com/username/plugin-name.git
cd plugin-name
# Copy to Claude Code location
```

### Method 2: Package Archive

1. Create release archive:
   ```bash
   tar -czf plugin-name-v1.0.0.tar.gz plugin-name/
   ```

2. Users download and extract:
   ```bash
   tar -xzf plugin-name-v1.0.0.tar.gz
   # Copy to Claude Code location
   ```

### Method 3: Plugin Registry (Future)

Future: Claude Code may support a central plugin registry.

---

## Complete Example Plugin

Here's a complete example of a PDF processing plugin:

```
pdf-tools/
├── plugin.json
├── README.md
├── LICENSE
├── skills/
│   ├── pdf-text-extractor/
│   │   ├── SKILL.md
│   │   ├── reference.md
│   │   └── scripts/
│   │       └── extract.py
│   └── pdf-table-parser/
│       └── SKILL.md
├── commands/
│   ├── pdf-info.md
│   └── pdf-merge.md
└── docs/
    └── examples.md
```

**plugin.json:**
```json
{
  "name": "pdf-tools",
  "version": "1.0.0",
  "description": "Comprehensive PDF processing tools for extraction, parsing, and manipulation",
  "author": "PDF Tools Team <team@pdftools.com>",
  "license": "MIT",
  "keywords": ["pdf", "extraction", "parsing", "documents"],
  "commands": [
    "./commands/pdf-info.md",
    "./commands/pdf-merge.md"
  ],
  "dependencies": {
    "required": [
      "python>=3.8",
      "pdfplumber>=0.9.0",
      "PyPDF2>=3.0.0"
    ],
    "optional": [
      "pytesseract>=0.3.0"
    ]
  },
  "configuration": {
    "settings": {
      "ocr_enabled": {
        "type": "boolean",
        "default": false,
        "description": "Enable OCR for image-based PDFs (requires pytesseract)"
      },
      "default_output_format": {
        "type": "string",
        "default": "text",
        "enum": ["text", "markdown", "json"],
        "description": "Default format for extracted content"
      }
    }
  }
}
```

---

## Best Practices

### General
1. **Clear naming**: Use descriptive, lowercase names with hyphens
2. **Semantic versioning**: Follow semver strictly
3. **Complete documentation**: README, examples, troubleshooting
4. **License clarity**: Include LICENSE file
5. **Changelog**: Maintain CHANGELOG.md for version history

### Skills
1. **Focused skills**: One capability per skill
2. **Clear triggers**: Description explains when to use
3. **Dependencies documented**: List all required packages
4. **Tool restrictions**: Use allowed-tools appropriately

### Commands
1. **Explicit manifest**: Always list commands in plugin.json
2. **Good defaults**: Provide sensible parameter defaults
3. **Clear descriptions**: Help users discover functionality
4. **Parameter hints**: Use argument-hint for guidance

### Dependencies
1. **Version pinning**: Specify minimum versions
2. **Separate optional**: Distinguish required vs optional
3. **Installation guide**: Document installation process
4. **Platform notes**: Mention platform-specific requirements

### Configuration
1. **Sensible defaults**: Work out of the box
2. **Clear descriptions**: Explain each setting
3. **Type validation**: Use appropriate types
4. **Enum constraints**: For settings with fixed options

---

## Common Pitfalls

1. **Invalid JSON**: Use JSON validator before publishing
2. **Missing required fields**: name, version, description are required
3. **Skills not found**: Must be in `skills/` directory
4. **Commands not listed**: Must be in `commands` array
5. **Wrong paths**: Use relative paths starting with `./`
6. **Unclear descriptions**: Be specific about functionality
7. **Undocumented dependencies**: List all requirements
8. **No examples**: Users don't know how to use plugin

---

## Publishing Checklist

Before publishing your plugin:

- [ ] `plugin.json` has all required fields
- [ ] Version follows semantic versioning
- [ ] All skills have valid SKILL.md files
- [ ] All commands are listed in manifest
- [ ] README.md is complete and clear
- [ ] LICENSE file is included
- [ ] CHANGELOG.md documents versions
- [ ] Dependencies are documented
- [ ] Examples are provided
- [ ] Tested locally with Claude Code
- [ ] JSON manifest validates
- [ ] Skills are discoverable
- [ ] Commands work as expected
- [ ] Configuration defaults are sensible
- [ ] Error handling is robust

---

## Support and Resources

- [Plugin Reference](https://docs.claude.com/en/docs/claude-code/plugins-reference.md)
- [Skills Documentation](https://docs.claude.com/en/docs/claude-code/skills.md)
- [Commands Documentation](https://docs.claude.com/en/docs/claude-code/slash-commands.md)
- Skill Template: `SKILL_TEMPLATE.md`
- Command Template: `COMMAND_TEMPLATE.md`

---

## Example Plugins to Study

Look at these example plugins for inspiration:

1. **Academic Tools**: Skills for research and writing
2. **Code Review**: Automated code review skills and commands
3. **DevOps**: Deployment and infrastructure commands
4. **Data Analysis**: Skills for data processing and visualization
5. **Documentation**: Auto-generate documentation skills

---

## Future Features

Claude Code plugin system may evolve to support:

- Central plugin registry
- Automatic updates
- Plugin dependencies
- Web-based plugin browser
- Plugin usage analytics
- Cross-platform compatibility helpers
- Plugin development tools

Stay tuned to Claude Code documentation for updates!
