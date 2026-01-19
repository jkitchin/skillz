---
description: Generate a changelog entry for recent changes
allowed-tools: ["Bash", "Read", "Glob"]
---

# Generate Changelog Entry

Generate a changelog entry based on recent commits or changes.

## Instructions

1. Check for existing changelog file:
   ```bash
   ls -la CHANGELOG* changelog* HISTORY* 2>/dev/null
   ```

2. Get recent commits (since last tag or last 20 commits):
   ```bash
   git describe --tags --abbrev=0 2>/dev/null || echo "none"
   git log --oneline -20
   ```

3. If a tag exists, get commits since that tag:
   ```bash
   git log <last-tag>..HEAD --oneline
   ```

4. Categorize changes into:
   - **Added**: New features
   - **Changed**: Changes in existing functionality
   - **Deprecated**: Soon-to-be removed features
   - **Removed**: Removed features
   - **Fixed**: Bug fixes
   - **Security**: Vulnerability fixes

5. Generate entry in Keep a Changelog format

## Output Format

```markdown
## [Unreleased] - YYYY-MM-DD

### Added
- New feature description

### Changed
- Change description

### Fixed
- Bug fix description
```

Version hint or additional context: $ARGUMENTS
