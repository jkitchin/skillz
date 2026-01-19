---
description: Analyze project dependencies
allowed-tools: ["Read", "Bash", "Glob"]
---

# Analyze Dependencies

Analyze and report on project dependencies.

## Instructions

1. Detect project type and find dependency files:
   - Python: requirements.txt, pyproject.toml, setup.py, Pipfile
   - JavaScript/TypeScript: package.json, package-lock.json, yarn.lock
   - Go: go.mod, go.sum
   - Rust: Cargo.toml, Cargo.lock
   - Ruby: Gemfile, Gemfile.lock
   - Java: pom.xml, build.gradle

2. Parse dependencies and categorize:
   - Direct vs transitive dependencies
   - Production vs development dependencies
   - Version constraints used

3. Analyze for:
   - Outdated packages (if possible to determine)
   - Security concerns (known vulnerable packages)
   - License compatibility
   - Unused dependencies (if determinable)

4. Additional context if provided: $ARGUMENTS

## Output Format

```markdown
## Dependency Analysis

### Project Type
<Language/Framework detected>

### Summary
- Total dependencies: <number>
- Direct dependencies: <number>
- Dev dependencies: <number>

### Direct Dependencies

| Package | Version | Purpose | Notes |
|---------|---------|---------|-------|
| package1 | ^1.2.3 | <purpose> | <notes> |

### Dev Dependencies

| Package | Version | Purpose |
|---------|---------|---------|
| package1 | ^1.2.3 | <purpose> |

### Dependency Tree
<Key dependency relationships>

### Potential Issues

#### Outdated
- `package`: current <version>, latest <version>

#### Security
- `package`: <CVE or advisory if known>

#### Recommendations
1. <Recommendation>
2. <Recommendation>

### Unused Dependencies
<List if determinable>

### License Summary
| License | Count | Packages |
|---------|-------|----------|
| MIT | 15 | pkg1, pkg2, ... |
```

Additional analysis focus: $ARGUMENTS
