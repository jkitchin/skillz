---
name: bash-logger
description: Log all bash commands executed by Claude for audit and review
event: PreToolUse
matcher: Bash
type: command
timeout: 5
---

# Bash Logger

Logs all bash commands executed by Claude Code to a file for audit, review, and debugging purposes.

## Purpose

Maintains a complete audit trail of all commands executed during Claude Code sessions. Useful for:

- Security auditing
- Debugging issues
- Understanding what Claude did
- Compliance requirements
- Learning from past sessions

## Log Format

Each entry includes:
- Timestamp
- Session ID
- Working directory
- Command executed
- Command description (if provided)

Example log entry:
```
[2024-01-15 14:30:22] [session:abc123] [cwd:/home/user/project]
Command: npm test
Description: Run test suite
---
```

## Log Location

Default: `~/.claude/logs/bash-commands.log`

Configure with environment variable:
```bash
export BASH_LOGGER_FILE="~/my-logs/claude-commands.log"
```

## Configuration

```bash
# Custom log file location
export BASH_LOGGER_FILE="~/.claude/logs/bash-commands.log"

# Log format: full (default), compact, json
export BASH_LOGGER_FORMAT="full"

# Exclude certain commands (regex pattern)
export BASH_LOGGER_EXCLUDE="^(ls|pwd|echo)"
```

## Behavior

1. Hook receives bash command from tool input
2. Formats log entry with metadata
3. Appends to log file
4. Exits 0 (never blocks)

## Notes

- **Non-blocking** - logs but never prevents commands
- Creates log directory if it doesn't exist
- Rotates logs when they exceed 10MB (keeps 5 backups)
