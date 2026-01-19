---
name: protect-secrets
description: Block writes to sensitive files containing secrets or credentials
event: PreToolUse
matcher: Edit|Write
type: command
timeout: 5
---

# Protect Secrets

Prevents Claude from accidentally writing to or modifying files that may contain secrets, credentials, or sensitive configuration.

## Purpose

Security guardrail that blocks writes to files that commonly contain sensitive information, preventing accidental exposure or modification of secrets.

## Protected Patterns

By default, blocks writes to:

- `.env` files (`.env`, `.env.local`, `.env.production`, etc.)
- Credential files (`credentials.json`, `secrets.json`, `secrets.yaml`)
- Private keys (`*.pem`, `*.key`, `id_rsa`, `id_ed25519`)
- AWS credentials (`~/.aws/credentials`, `~/.aws/config`)
- SSH config (`~/.ssh/config`, `~/.ssh/known_hosts`)
- Git credentials (`.git-credentials`)
- Docker secrets
- Kubernetes secrets

## Configuration

Set environment variables to customize:

```bash
# Add additional patterns (comma-separated)
export PROTECT_SECRETS_EXTRA="*.secret,my-credentials.txt"

# Disable specific patterns
export PROTECT_SECRETS_ALLOW=".env.example"
```

## Behavior

1. Hook receives the file path from Edit/Write tool input
2. Checks against protected patterns
3. If match found:
   - Exits with code 2 (blocks operation)
   - Sends error message to Claude
4. If no match: exits 0 (allows operation)

## Notes

- **Blocking hook** - prevents the operation from proceeding
- Error messages guide Claude to use safer alternatives
- Can be bypassed with explicit user confirmation in Claude Code
