# Ralph Wiggum Example: Converting litdb to TypeScript

This example demonstrates how to use Ralph Wiggum mode to autonomously convert [litdb](https://github.com/jkitchin/litdb), a Python literature database tool, into a standalone TypeScript library and CLI with high test coverage.

## Project Overview

**Original Project:** [jkitchin/litdb](https://github.com/jkitchin/litdb)
- Python-based literature database for researchers
- 43 CLI commands for managing academic papers
- Vector search, full-text search, and hybrid search
- Integration with OpenAlex, DOIs, ORCIDs
- SQLite/LibSQL storage with embeddings

**Target:** `litdb-ts`
- TypeScript library (`@litdb/core`)
- CLI tool (`litdb` command)
- 90%+ test coverage
- Modern Node.js architecture (ESM)
- Same feature set as Python version

## Why This is a Good Ralph Project

1. **Well-defined scope** - Clear Python reference to port
2. **Measurable progress** - Tests and coverage provide feedback
3. **Incremental** - Can implement command-by-command
4. **Long-running** - Too big for a single session, perfect for overnight runs

## Prerequisites

```bash
# Install Ralph Wiggum skill
skillz install ralph-wiggum
skillz hooks install ralph-safety-check
skillz hooks install ralph-cost-monitor

# Ensure Docker is installed
docker --version
```

## Step-by-Step Setup

### Step 1: Create Project Directory

```bash
mkdir litdb-ts
cd litdb-ts
git init
```

### Step 2: Initialize Ralph

```bash
# Copy Ralph files (or run ralph-init.sh)
cp ~/.claude/skills/ralph-wiggum/scripts/*.sh .
cp ~/.claude/skills/ralph-wiggum/templates/* .
chmod +x *.sh
```

### Step 3: Copy Example Specifications

```bash
# Copy the specifications from this example
cp -r ~/.claude/skills/ralph-wiggum/examples/litdb-typescript/specs .
cp ~/.claude/skills/ralph-wiggum/examples/litdb-typescript/IMPLEMENTATION_PLAN.md .
cp ~/.claude/skills/ralph-wiggum/examples/litdb-typescript/RALPH_PROMPT.md .
cp ~/.claude/skills/ralph-wiggum/examples/litdb-typescript/Dockerfile.ralph .
```

### Step 4: Set Up API Key

```bash
export ANTHROPIC_API_KEY=your-key-here
```

### Step 5: Run Ralph

```bash
# Start with planning mode (optional, recommended)
./ralph-sandbox.sh . 5 plan

# Review the plan, then run build mode
./ralph-sandbox.sh . 50 build

# Or run unlimited overnight
./ralph-sandbox.sh
```

### Step 6: Monitor Progress

In another terminal:

```bash
# Watch the log
tail -f ralph.log

# Check git commits
watch -n 10 'git log --oneline -10'

# Check test coverage
cat coverage/coverage-summary.json 2>/dev/null | jq '.total'

# Check costs
cat .ralph-costs.json | jq '.summary'
```

### Step 7: Review Results

```bash
# See what Ralph accomplished
git log --oneline

# Check implementation status
cat IMPLEMENTATION_PLAN.md

# Run tests manually to verify
npm test

# Check coverage report
open coverage/lcov-report/index.html
```

## Expected Timeline

With Claude Sonnet, expect roughly:

| Phase | Iterations | Estimated Time |
|-------|------------|----------------|
| Project setup | 2-3 | 15 min |
| Core database layer | 5-8 | 1 hour |
| Search functionality | 8-12 | 2 hours |
| CLI commands (43 total) | 40-60 | 6-8 hours |
| Test coverage boost | 10-20 | 2 hours |
| Polish & docs | 5-10 | 1 hour |

**Total: ~80-110 iterations, 12-14 hours**

This is ideal for an overnight run. Start Ralph at 6 PM, wake up to a working TypeScript port.

## Cost Estimation

Rough estimate with Claude Sonnet:
- ~100 iterations
- ~50k input tokens + ~20k output tokens per iteration
- Cost: ~$20-40 total

Set `RALPH_COST_LIMIT=50` to cap spending.

## Customization

### Adjust Iterations Per Session

```bash
# Quick test (5 iterations)
./ralph-sandbox.sh . 5

# Medium run (20 iterations)
./ralph-sandbox.sh . 20

# Full overnight (unlimited)
./ralph-sandbox.sh
```

### Use Different Model

```bash
# Faster, cheaper (but less capable)
RALPH_MODEL=haiku ./ralph-sandbox.sh . 50

# Most capable (expensive)
RALPH_MODEL=opus ./ralph-sandbox.sh . 20
```

### Enable Network for npm install

Edit `ralph-sandbox.sh` to allow specific registries:

```bash
# Change --network none to:
--network bridge \
```

Or pre-install dependencies before running Ralph.

## Troubleshooting

### Ralph Gets Stuck

Check `ralph.log` for errors. Common issues:
- Missing dependencies (pre-run `npm install`)
- Test failures blocking progress
- Context filling up on large files

Solution: Add blockers to `IMPLEMENTATION_PLAN.md` and restart.

### Tests Keep Failing

Ralph should fix failing tests, but if stuck:
1. Review the failing tests in `ralph.log`
2. Add specific guidance in `RALPH_PROMPT.md`
3. Restart Ralph

### Cost Limit Reached

```bash
# Check current spend
cat .ralph-costs.json | jq '.summary.estimated_cost_usd'

# Increase limit and continue
RALPH_COST_LIMIT=100 ./ralph-sandbox.sh
```

## Files in This Example

```
examples/litdb-typescript/
├── README.md                    # This file
├── IMPLEMENTATION_PLAN.md       # Phased task list
├── RALPH_PROMPT.md             # Custom prompt for this project
├── Dockerfile.ralph            # Node.js sandbox container
└── specs/
    ├── overview.md             # High-level architecture
    ├── database.md             # Database schema spec
    ├── search.md               # Search functionality spec
    ├── cli-commands.md         # CLI command specifications
    ├── testing.md              # Testing requirements
    └── api-reference.md        # Public API design
```

## Success Criteria

Ralph is "done" when:
- [ ] All 43 CLI commands implemented
- [ ] `npm test` passes
- [ ] Coverage > 90%
- [ ] `npm run build` produces working dist/
- [ ] `litdb --help` shows all commands
- [ ] Basic smoke test: add paper, search, export bibtex

## Next Steps After Ralph

1. **Human review** - Review all generated code
2. **Integration tests** - Add end-to-end tests
3. **Documentation** - Expand README, add JSDoc
4. **Publishing** - Set up npm publishing workflow
5. **CI/CD** - Add GitHub Actions

---

*This example demonstrates Ralph's ability to tackle large, well-defined porting projects. The key is thorough specifications and a clear implementation plan.*
