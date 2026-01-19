#!/bin/bash
# run-example.sh - Set up and run Ralph for litdb-ts conversion
#
# Usage:
#   ./run-example.sh init     - Create project directory and copy files
#   ./run-example.sh deps     - Install npm dependencies (requires network)
#   ./run-example.sh plan     - Run Ralph in planning mode
#   ./run-example.sh build    - Run Ralph in build mode
#   ./run-example.sh status   - Check current progress

set -euo pipefail

# Configuration
PROJECT_DIR="${LITDB_TS_DIR:-$HOME/litdb-ts}"
EXAMPLE_DIR="$(dirname "$(readlink -f "$0")")"
MAX_ITERATIONS="${MAX_ITERATIONS:-0}"

# Colors
GREEN='\033[0;32m'
BLUE='\033[0;34m'
YELLOW='\033[1;33m'
NC='\033[0m'

show_help() {
    cat << EOF
Usage: ./run-example.sh <command>

Commands:
    init     Create project directory and copy all files
    deps     Install npm dependencies (requires Docker with network)
    plan     Run Ralph in planning mode (5 iterations)
    build    Run Ralph in build mode (unlimited or MAX_ITERATIONS)
    status   Show current progress
    help     Show this help

Environment Variables:
    LITDB_TS_DIR      Project directory (default: ~/litdb-ts)
    MAX_ITERATIONS    Max iterations for build mode (default: unlimited)
    ANTHROPIC_API_KEY Required for running Ralph

Example Workflow:
    1. ./run-example.sh init
    2. cd ~/litdb-ts
    3. ./run-example.sh deps
    4. export ANTHROPIC_API_KEY=your-key
    5. ./run-example.sh build
EOF
}

cmd_init() {
    echo -e "${BLUE}Initializing litdb-ts project...${NC}"

    # Create project directory
    mkdir -p "$PROJECT_DIR"
    cd "$PROJECT_DIR"

    # Initialize git
    if [[ ! -d .git ]]; then
        git init
        echo -e "${GREEN}✓${NC} Git initialized"
    fi

    # Copy example files
    cp -r "$EXAMPLE_DIR/specs" .
    cp "$EXAMPLE_DIR/IMPLEMENTATION_PLAN.md" .
    cp "$EXAMPLE_DIR/RALPH_PROMPT.md" .
    cp "$EXAMPLE_DIR/Dockerfile.ralph" .
    echo -e "${GREEN}✓${NC} Copied specs and plans"

    # Copy Ralph scripts from skill
    SKILL_DIR="$(dirname "$EXAMPLE_DIR")"
    cp "$SKILL_DIR/scripts/ralph.sh" .
    cp "$SKILL_DIR/scripts/ralph-sandbox.sh" .
    chmod +x ralph.sh ralph-sandbox.sh
    echo -e "${GREEN}✓${NC} Copied Ralph scripts"

    # Create initial package.json
    cat > package.json << 'PACKAGE'
{
  "name": "litdb-ts",
  "version": "0.1.0",
  "description": "TypeScript port of litdb - Literature Database",
  "type": "module",
  "main": "dist/index.js",
  "types": "dist/index.d.ts",
  "bin": {
    "litdb": "bin/litdb.js"
  },
  "scripts": {
    "build": "tsup src/index.ts --format esm --dts",
    "test": "vitest run",
    "test:watch": "vitest",
    "test:coverage": "vitest run --coverage",
    "lint": "eslint src tests",
    "typecheck": "tsc --noEmit"
  },
  "dependencies": {
    "better-sqlite3": "^9.4.0",
    "commander": "^12.0.0"
  },
  "devDependencies": {
    "@types/better-sqlite3": "^7.6.8",
    "@types/node": "^20.11.0",
    "@vitest/coverage-v8": "^1.2.0",
    "eslint": "^8.56.0",
    "tsup": "^8.0.0",
    "typescript": "^5.3.0",
    "vitest": "^1.2.0"
  },
  "engines": {
    "node": ">=20"
  }
}
PACKAGE
    echo -e "${GREEN}✓${NC} Created package.json"

    # Create tsconfig.json
    cat > tsconfig.json << 'TSCONFIG'
{
  "compilerOptions": {
    "target": "ES2022",
    "module": "ESNext",
    "moduleResolution": "bundler",
    "lib": ["ES2022"],
    "outDir": "dist",
    "rootDir": "src",
    "strict": true,
    "esModuleInterop": true,
    "skipLibCheck": true,
    "declaration": true,
    "declarationMap": true,
    "sourceMap": true
  },
  "include": ["src/**/*"],
  "exclude": ["node_modules", "dist", "tests"]
}
TSCONFIG
    echo -e "${GREEN}✓${NC} Created tsconfig.json"

    # Create vitest.config.ts
    cat > vitest.config.ts << 'VITEST'
import { defineConfig } from 'vitest/config';

export default defineConfig({
  test: {
    globals: true,
    environment: 'node',
    include: ['tests/**/*.test.ts'],
    coverage: {
      provider: 'v8',
      reporter: ['text', 'json', 'html', 'lcov'],
      include: ['src/**/*.ts'],
      exclude: ['src/**/*.d.ts'],
      thresholds: {
        global: {
          branches: 85,
          functions: 90,
          lines: 90,
          statements: 90,
        },
      },
    },
  },
});
VITEST
    echo -e "${GREEN}✓${NC} Created vitest.config.ts"

    # Create directory structure
    mkdir -p src/{db,search,sources,embeddings,cli/commands,utils}
    mkdir -p tests/{unit,integration,fixtures}
    mkdir -p bin
    echo -e "${GREEN}✓${NC} Created directory structure"

    # Create placeholder index
    cat > src/index.ts << 'INDEX'
// litdb-ts - TypeScript Literature Database
// This file will export the public API

export const VERSION = '0.1.0';
INDEX
    echo -e "${GREEN}✓${NC} Created src/index.ts"

    # Create .gitignore
    cat > .gitignore << 'GITIGNORE'
node_modules/
dist/
coverage/
*.log
ralph.log
.ralph-costs.json
.npm-cache/
GITIGNORE
    echo -e "${GREEN}✓${NC} Created .gitignore"

    # Initial commit
    git add -A
    git commit -m "Initial litdb-ts project setup"
    echo -e "${GREEN}✓${NC} Created initial commit"

    echo ""
    echo -e "${GREEN}Project initialized at: $PROJECT_DIR${NC}"
    echo ""
    echo "Next steps:"
    echo "  cd $PROJECT_DIR"
    echo "  ./run-example.sh deps    # Install dependencies"
    echo "  ./run-example.sh build   # Start Ralph"
}

cmd_deps() {
    echo -e "${BLUE}Installing dependencies...${NC}"
    cd "$PROJECT_DIR"

    # Build Docker image with network
    docker build -t ralph-litdb -f Dockerfile.ralph .

    # Run npm install inside container with network
    docker run --rm -it \
        -v "$(pwd):/workspace" \
        ralph-litdb \
        bash -c "npm install && npm run build 2>/dev/null || true"

    echo -e "${GREEN}✓${NC} Dependencies installed"
}

cmd_plan() {
    echo -e "${BLUE}Running Ralph in planning mode...${NC}"
    cd "$PROJECT_DIR"

    if [[ -z "${ANTHROPIC_API_KEY:-}" ]]; then
        echo "Error: ANTHROPIC_API_KEY not set"
        exit 1
    fi

    ./ralph-sandbox.sh . 5 plan
}

cmd_build() {
    echo -e "${BLUE}Running Ralph in build mode...${NC}"
    cd "$PROJECT_DIR"

    if [[ -z "${ANTHROPIC_API_KEY:-}" ]]; then
        echo "Error: ANTHROPIC_API_KEY not set"
        exit 1
    fi

    ./ralph-sandbox.sh . "$MAX_ITERATIONS" build
}

cmd_status() {
    cd "$PROJECT_DIR"

    echo -e "${BLUE}=== litdb-ts Status ===${NC}"
    echo ""

    # Task progress
    echo "Tasks:"
    total=$(grep -c '^\- \[' IMPLEMENTATION_PLAN.md 2>/dev/null || echo 0)
    done=$(grep -c '^\- \[x\]' IMPLEMENTATION_PLAN.md 2>/dev/null || echo 0)
    inprog=$(grep -c '^\- \[\~\]' IMPLEMENTATION_PLAN.md 2>/dev/null || echo 0)
    pending=$((total - done - inprog))
    echo "  Completed: $done"
    echo "  In Progress: $inprog"
    echo "  Pending: $pending"
    echo "  Total: $total"
    echo ""

    # Git commits
    echo "Recent commits:"
    git log --oneline -5 2>/dev/null || echo "  (no commits yet)"
    echo ""

    # Test coverage
    if [[ -f coverage/coverage-summary.json ]]; then
        echo "Coverage:"
        cat coverage/coverage-summary.json | grep -A4 '"total"' | head -5
    else
        echo "Coverage: (not yet generated)"
    fi
    echo ""

    # Costs
    if [[ -f .ralph-costs.json ]]; then
        echo "API Costs:"
        cat .ralph-costs.json | grep -A3 '"summary"' | head -4
    else
        echo "Costs: (no Ralph runs yet)"
    fi
}

# Main
case "${1:-help}" in
    init)   cmd_init ;;
    deps)   cmd_deps ;;
    plan)   cmd_plan ;;
    build)  cmd_build ;;
    status) cmd_status ;;
    help)   show_help ;;
    *)
        echo "Unknown command: $1"
        show_help
        exit 1
        ;;
esac
