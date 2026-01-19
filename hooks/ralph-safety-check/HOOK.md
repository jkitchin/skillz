---
name: ralph-safety-check
description: Block dangerous operations during Ralph mode even with permissions disabled
event: PreToolUse
matcher: Bash
type: command
timeout: 5
---

# Ralph Safety Check

Additional safety layer for Ralph Wiggum mode that blocks dangerous bash commands even when running with `--dangerously-skip-permissions`.

## Purpose

While Ralph's sandbox enforcement (container, network isolation, etc.) provides the primary security layer, this hook adds defense-in-depth by blocking obviously dangerous command patterns.

## Blocked Patterns

### Remote Code Execution
- `curl | bash`, `curl | sh`
- `wget | bash`, `wget | sh`
- `curl -o - | bash`

### Credential Access
- Commands accessing `~/.ssh/*`
- Commands accessing `~/.aws/*`
- Commands accessing `~/.config/gcloud/*`
- Reading `/etc/shadow`, `/etc/passwd`

### Destructive Operations
- `rm -rf /`
- `rm -rf ~`
- `rm -rf /*`
- `mkfs.*`
- `dd if=* of=/dev/*`

### Network Reconnaissance
- `nmap`
- `nc -l` (netcat listen)
- `ssh` to external hosts

### System Modification
- `chmod 777`
- `chown root`
- Modifying `/etc/*`

## Configuration

Set environment variables to customize:

```bash
# Disable specific checks (comma-separated)
export RALPH_SAFETY_ALLOW="ssh,nmap"

# Add additional blocked patterns
export RALPH_SAFETY_BLOCK="custom-dangerous-cmd"
```

## Behavior

1. Hook receives the bash command from tool input
2. Checks against blocked patterns
3. If dangerous pattern found:
   - Exits with code 2 (blocks operation)
   - Sends warning message to Claude
4. If no dangerous pattern: exits 0 (allows operation)

## Notes

- This is a **blocking hook** - it prevents the operation
- Works alongside (not instead of) sandbox enforcement
- Some patterns may have legitimate uses - configure allowlist if needed
- Claude can still be asked to run blocked commands, but must get human approval
