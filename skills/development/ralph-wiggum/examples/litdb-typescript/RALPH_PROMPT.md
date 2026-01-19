# Ralph Mode: litdb-ts TypeScript Port

You are porting the Python [litdb](https://github.com/jkitchin/litdb) project to TypeScript. Your progress persists in files and git, not in your context. Each iteration starts fresh.

## Project Context

**Goal:** Create a TypeScript library and CLI that replicates Python litdb functionality.
**Coverage Target:** 90%+ test coverage
**Reference:** https://github.com/jkitchin/litdb (study this for behavior)

## Tech Stack

- TypeScript 5.x, Node.js 20+
- better-sqlite3 for database
- sqlite-vss for vector search
- @xenova/transformers for embeddings
- Commander.js for CLI
- Vitest for testing

## Each Iteration

1. **Read state**
   ```bash
   cat IMPLEMENTATION_PLAN.md
   git log --oneline -5
   npm test 2>&1 | tail -20  # Check test status
   ```

2. **Select ONE task** from IMPLEMENTATION_PLAN.md marked `- [ ]`

3. **Mark in progress** - Update to `- [~]` in IMPLEMENTATION_PLAN.md

4. **Implement** following this order:
   - Write/update the source code
   - Write tests FIRST or alongside code
   - Run tests: `npm test`
   - Fix any failures

5. **Verify coverage**
   ```bash
   npm run test:coverage 2>&1 | grep -A5 "All files"
   ```

6. **Mark complete** - Update to `- [x]` in IMPLEMENTATION_PLAN.md

7. **Commit**
   ```bash
   git add -A
   git commit -m "feat: <what you implemented>"
   ```

## Coding Standards

### File Structure
```typescript
// src/db/client.ts
import { Database } from 'better-sqlite3';
import { Work, NewWork } from './types';

export class DatabaseClient {
  private db: Database;

  constructor(path: string) {
    this.db = new Database(path);
    this.db.pragma('journal_mode = WAL');
  }

  // Methods...
}
```

### Testing Pattern
```typescript
// tests/unit/db/client.test.ts
import { describe, it, expect, beforeEach } from 'vitest';
import { DatabaseClient } from '../../../src/db/client';

describe('DatabaseClient', () => {
  let db: DatabaseClient;

  beforeEach(() => {
    db = new DatabaseClient(':memory:');
  });

  describe('addWork', () => {
    it('adds work with required fields', () => {
      // Test implementation
    });

    it('rejects duplicate DOI', () => {
      // Error case
    });
  });
});
```

### CLI Command Pattern
```typescript
// src/cli/commands/add.ts
import { Command } from 'commander';
import { DatabaseClient } from '../../db/client';

export function createAddCommand(db: DatabaseClient): Command {
  return new Command('add')
    .description('Add a work by DOI')
    .argument('<doi>', 'DOI or DOI URL')
    .option('--tag <tag>', 'Add tag to work')
    .action(async (doi, options) => {
      // Implementation
    });
}
```

## Rules

1. **Tests are mandatory** - Never commit code without tests
2. **One task per iteration** - Stay focused
3. **Run tests before committing** - `npm test` must pass
4. **Check coverage** - Keep it above 85%, aim for 90%
5. **Follow existing patterns** - Match the style in the codebase
6. **Reference Python litdb** - Use it for exact behavior/API
7. **Commit frequently** - Small, focused commits

## If Blocked

- **Missing dependency?** Add to package.json, run `npm install`
- **Test failing?** Debug it, don't skip it
- **Unclear spec?** Check Python litdb source code
- **Can't complete?** Document in IMPLEMENTATION_PLAN.md Blockers section

## Coverage Checkpoints

After each phase, verify coverage:
```bash
npm run test:coverage 2>&1 | grep "All files"
```

Expected progress:
- Phase 2 (Database): 30% overall, 95% src/db/
- Phase 3 (Search): 45% overall, 90% src/search/
- Phase 8 (CLI starts): 60% overall
- Phase 18 (Polish): 90%+ overall

## Priority Order

If short on time, prioritize:
1. Database layer (foundation)
2. Search functionality (core feature)
3. `add`, `search`, `export` commands (minimum viable CLI)
4. Test coverage boost
5. Remaining commands

## Dependencies to Pre-install

If network is disabled, ensure these are installed:
```json
{
  "dependencies": {
    "better-sqlite3": "^9.0.0",
    "commander": "^11.0.0",
    "@xenova/transformers": "^2.0.0"
  },
  "devDependencies": {
    "typescript": "^5.0.0",
    "vitest": "^1.0.0",
    "@vitest/coverage-v8": "^1.0.0",
    "tsup": "^8.0.0",
    "@types/better-sqlite3": "^7.0.0",
    "@types/node": "^20.0.0"
  }
}
```

## Success Criteria

This iteration is successful if:
- [ ] One task completed and marked done
- [ ] Tests pass (`npm test`)
- [ ] Coverage didn't decrease
- [ ] Code committed to git
- [ ] No TypeScript errors (`npm run build`)
