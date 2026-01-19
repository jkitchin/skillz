# Security Audit Report - Skillz Repository

**Date**: 2026-01-19
**Auditor**: Claude (Red Team Assessment)
**Scope**: Full repository audit including CLI, hooks, agents, and commands

---

## Executive Summary

This audit identified **3 high-severity**, **5 medium-severity**, and **4 low-severity** security issues in the skillz codebase. The most critical issues involve path traversal vulnerabilities and the inherent risks of executing user-provided hook scripts.

### Risk Rating: **MEDIUM-HIGH**

The primary attack surface is through:
1. Malicious skill/hook/agent names that could escape intended directories
2. Hook scripts that execute with user privileges
3. AI-generated content that could contain malicious code

---

## Critical/High Severity Issues

### H1: Path Traversal in Name Parameters

**Location**: `cli/commands/install.py:96`, `cli/commands/hooks.py:147,210`, `cli/commands/agents.py:141`

**Issue**: User-supplied `name` arguments are concatenated directly to paths without sanitization beyond basic regex validation.

```python
# hooks.py:147
dest_path = dest_dir / name

# install.py:96
dest_path = dest_dir / name
```

**Attack Vector**: While the regex `^[a-z0-9-]+$` prevents `../` sequences, an attacker with filesystem access could create a symlink within the hooks/skills directory pointing elsewhere.

**Recommendation**:
```python
def safe_path_join(base: Path, name: str) -> Path:
    """Safely join paths, preventing traversal attacks."""
    result = (base / name).resolve()
    if not str(result).startswith(str(base.resolve())):
        raise ValueError(f"Path traversal detected: {name}")
    return result
```

**Severity**: HIGH
**CVSS**: 7.5

---

### H2: Arbitrary Code Execution via Hooks

**Location**: `cli/commands/hooks.py:174-177`, Hook installation process

**Issue**: The hook system installs and executes arbitrary Python/shell scripts with user privileges.

```python
# hooks.py:174-177
for script in dest_path.glob("*.py"):
    os.chmod(script, 0o755)
for script in dest_path.glob("*.sh"):
    os.chmod(script, 0o755)
```

**Attack Vector**:
1. Attacker creates malicious hook in repository
2. User runs `skillz hooks install malicious-hook`
3. Malicious script executes with full user privileges on every Claude Code event

**Recommendation**:
- Add script signature verification
- Implement a hook allow-list in settings
- Display script contents before installation with confirmation
- Consider sandboxing hook execution

**Severity**: HIGH
**CVSS**: 8.1

---

### H3: AI-Generated Code Injection

**Location**: `cli/commands/create.py:156-196`, `cli/commands/hooks.py:631-698`, `cli/commands/agents.py:466-510`

**Issue**: The `--prompt` feature uses Claude CLI to generate code that is written directly to files without sanitization.

```python
# hooks.py:670-672
parts = output.split("===HOOK_SCRIPT===")
hook_md = parts[0].strip()
hook_script = parts[1].strip()  # Written directly to hook.py
```

**Attack Vector**:
1. Attacker crafts prompt that causes Claude to generate malicious hook code
2. Generated code is written to `hook.py` and made executable
3. Malicious code runs on next hook trigger

**Recommendation**:
- Add static analysis of generated code before writing
- Require user confirmation showing the generated code
- Implement a "preview only" mode for AI generation
- Consider AST validation for Python scripts

**Severity**: HIGH
**CVSS**: 7.8

---

## Medium Severity Issues

### M1: Symlink Attack in Copy Operations

**Location**: `cli/utils.py:40-52`

**Issue**: `shutil.copytree()` and `shutil.rmtree()` follow symlinks by default, which could lead to:
- Reading/writing files outside intended directories
- Deletion of unintended files

```python
# utils.py:46-48
if dst.exists():
    shutil.rmtree(dst)  # Follows symlinks!
shutil.copytree(src, dst)
```

**Recommendation**:
```python
shutil.copytree(src, dst, symlinks=False, ignore_dangling_symlinks=True)
shutil.rmtree(dst, ignore_errors=False, onerror=handle_remove_error)
```

**Severity**: MEDIUM
**CVSS**: 5.5

---

### M2: Unvalidated Platform Parameter

**Location**: `cli/config.py:93-115`

**Issue**: The `platform` parameter is used to construct file paths without validation beyond checking if it exists in config.

```python
# config.py:96-99
if platform in self.config["platforms"]:
    path = self.config["platforms"][platform]["skills_dir"]
```

**Attack Vector**: A malicious config file could define arbitrary platforms pointing to sensitive directories.

**Recommendation**:
- Validate platform against a hardcoded whitelist: `{"claude", "opencode", "codex", "gemini"}`
- Validate that paths are within expected directories

**Severity**: MEDIUM
**CVSS**: 5.3

---

### M3: Shell Command Injection in Lab Notebook Hook

**Location**: `hooks/lab-notebook/hook.py:94-128`

**Issue**: The `cwd` parameter from input is passed directly to subprocess commands.

```python
# lab-notebook/hook.py:94-99
result = subprocess.run(
    ["git", "rev-parse", "--is-inside-work-tree"],
    cwd=cwd,  # User-controlled!
    capture_output=True,
    text=True,
)
```

**Attack Vector**: While `subprocess.run` with a list is safe from shell injection, a malicious `cwd` value could cause git to operate in unintended directories or trigger git hooks.

**Recommendation**:
- Validate `cwd` is a real directory
- Ensure `cwd` is within expected project boundaries
- Consider running git commands with `--no-optional-locks`

**Severity**: MEDIUM
**CVSS**: 5.0

---

### M4: Settings.json Manipulation

**Location**: `cli/commands/hooks.py:536-591`

**Issue**: The hook installation process modifies `settings.json` which could be manipulated to:
- Register malicious hooks
- Overwrite other settings
- Create hooks pointing to arbitrary executables

```python
# hooks.py:558-566
hook_entry = {
    "matcher": matcher,
    "hooks": [
        {
            "type": hook_type,
            "command": _get_hook_command(hook_path),  # Arbitrary path
            ...
        }
    ],
}
```

**Recommendation**:
- Validate hook command paths are within allowed directories
- Add integrity checking for settings.json
- Log all settings modifications

**Severity**: MEDIUM
**CVSS**: 5.5

---

### M5: Insufficient Input Validation on Frontmatter

**Location**: `cli/validator.py` (all validators)

**Issue**: While YAML is parsed with `safe_load()`, the frontmatter values aren't fully sanitized:
- Descriptions could contain markdown injection
- Names with max length 64 could still cause issues in some filesystems
- Tool names are validated but could be manipulated in list format

**Recommendation**:
- Sanitize all frontmatter values for their intended use
- Add content security policy for markdown output
- Validate tools list items individually

**Severity**: MEDIUM
**CVSS**: 4.3

---

## Low Severity Issues

### L1: Information Disclosure in Error Messages

**Location**: Throughout `cli/commands/*.py`

**Issue**: Error messages reveal file system paths and structure.

```python
console.print(f"[red]Error: Hook '{name}' not found in repository[/red]")
console.print(f"Source: {source_path}")  # Full path disclosure
```

**Recommendation**: Use relative paths in user-facing messages.

**Severity**: LOW

---

### L2: Sensitive Data in Lab Notebook Logs

**Location**: `hooks/lab-notebook/hook.py`

**Issue**: The lab notebook hook logs all commands and file operations, which could capture:
- Passwords passed to CLI commands
- API keys in environment variables
- Sensitive file contents

**Recommendation**:
- Add a filter list for sensitive patterns
- Redact obvious secrets (API keys, passwords)
- Add a `LAB_NOTEBOOK_REDACT_PATTERNS` environment variable

**Severity**: LOW

---

### L3: Missing Rate Limiting on AI Generation

**Location**: `cli/commands/create.py`, `cli/commands/hooks.py`, `cli/commands/agents.py`

**Issue**: No rate limiting on Claude CLI calls could lead to:
- API quota exhaustion
- Cost amplification if running with paid API

**Recommendation**: Add configurable rate limiting and confirmation for AI generation.

**Severity**: LOW

---

### L4: Protect-Secrets Bypass via Encoding

**Location**: `hooks/protect-secrets/hook.py`

**Issue**: The protect-secrets hook only matches literal file names. Encoded paths or case variations could bypass protection.

```python
# Only matches exact patterns
if fnmatch.fnmatch(name, pattern):
    return True, pattern
```

**Attack Vector**: Files like `.ENV`, `.Env`, or URL-encoded paths might bypass detection.

**Recommendation**:
- Normalize paths before matching (case, encoding)
- Match against canonical/resolved paths
- Add content-based detection for secrets

**Severity**: LOW

---

## Positive Security Findings

The audit also identified several good security practices:

1. **YAML Safe Loading**: All YAML parsing uses `yaml.safe_load()` (not `yaml.load()`), preventing arbitrary code execution via YAML deserialization.

2. **Input Validation Regex**: The name validation regex `^[a-z0-9-]+$` effectively prevents most path traversal attacks.

3. **Subprocess with Lists**: All `subprocess.run()` calls use list arguments rather than shell=True, preventing shell injection.

4. **No Hardcoded Secrets**: No credentials or API keys found in the codebase.

5. **Reasonable Dependencies**: All dependencies (click, pyyaml, rich) are well-maintained with no known critical CVEs.

---

## Recommendations Summary

### Immediate Actions (High Priority)
1. Implement path validation with symlink checks
2. Add script preview/confirmation before hook installation
3. Implement code review prompt for AI-generated content

### Short-term Actions (Medium Priority)
4. Validate platform parameter against whitelist
5. Add symlink-safe file operations
6. Implement settings.json integrity checking
7. Add sensitive data filtering to lab notebook

### Long-term Actions (Low Priority)
8. Consider hook sandboxing (containers, seccomp)
9. Implement hook signing/verification
10. Add comprehensive audit logging

---

## Testing Recommendations

Add security-focused tests:

```python
def test_path_traversal_in_name():
    """Ensure path traversal attempts are blocked."""
    malicious_names = ["../etc/passwd", "..\\windows\\system32", "foo/../../../bar"]
    for name in malicious_names:
        assert not validate_name(name)

def test_symlink_attack_prevention():
    """Ensure symlinks don't escape intended directories."""
    # Create symlink pointing outside
    # Attempt copy operation
    # Verify operation fails or stays within bounds

def test_hook_script_validation():
    """Ensure hook scripts are validated before execution."""
    malicious_script = "import os; os.system('rm -rf /')"
    # Verify such scripts are detected/blocked
```

---

*This report is provided for security improvement purposes. All findings should be verified and fixes tested before deployment.*
