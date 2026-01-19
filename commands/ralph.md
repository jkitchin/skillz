---
description: Initialize and manage Ralph Wiggum autonomous coding mode
allowed-tools:
  - Bash
  - Read
  - Write
  - Edit
argument-hint: "<init|status|clean> [options]"
---

# Ralph Wiggum Mode

Initialize and manage Ralph Wiggum autonomous coding mode in your project.

## Commands

### `/ralph init`

Initialize Ralph in the current project. Creates:
- `ralph.sh` - Main loop script with sandbox enforcement
- `ralph-sandbox.sh` - Docker wrapper for safe execution
- `RALPH_PROMPT.md` - Prompt template for Claude
- `IMPLEMENTATION_PLAN.md` - Task tracking file
- `Dockerfile.ralph` - Sandbox container definition
- `specs/` - Directory for requirements

### `/ralph status`

Check Ralph status:
- Whether Ralph files exist
- Current iteration count (from ralph.log)
- Task completion status (from IMPLEMENTATION_PLAN.md)
- Cost tracking (if ralph-cost-monitor hook is used)

### `/ralph clean`

Clean up Ralph artifacts:
- Remove `ralph.log`
- Reset `IMPLEMENTATION_PLAN.md` to template
- Optionally remove `.ralph-costs.json`

## Usage Flow

```bash
# 1. Initialize Ralph in your project
/ralph init

# 2. Edit your tasks
# - Add requirements to specs/
# - Update IMPLEMENTATION_PLAN.md with tasks

# 3. Start Ralph in sandbox (run this in terminal, not here)
./ralph-sandbox.sh

# 4. Check progress
/ralph status

# 5. Clean up when done
/ralph clean
```

## After Initialization

Once initialized, run Ralph from your terminal:

```bash
# Basic usage (unlimited iterations)
./ralph-sandbox.sh

# Limit to 20 iterations
./ralph-sandbox.sh . 20

# Run in planning mode first
./ralph-sandbox.sh . 5 plan
```

## Important Notes

1. **Sandbox Required**: Ralph will refuse to run outside a Docker container
2. **API Key**: Set `ANTHROPIC_API_KEY` environment variable before running
3. **Network Disabled**: The sandbox has no internet access for safety
4. **Review Output**: Always review Ralph's work before merging

## Troubleshooting

If Ralph won't start:
- Ensure Docker is installed and running
- Check that `ANTHROPIC_API_KEY` is set
- Verify `ralph.sh` has execute permission (`chmod +x ralph.sh`)

If Ralph gets stuck:
- Check `ralph.log` for errors
- Update `IMPLEMENTATION_PLAN.md` to unblock
- Add clarifying notes to `RALPH_PROMPT.md`
