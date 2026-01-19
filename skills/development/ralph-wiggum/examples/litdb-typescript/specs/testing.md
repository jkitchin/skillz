# Testing Specification

## Coverage Goals

| Component | Target Coverage |
|-----------|-----------------|
| `src/db/` | 95% |
| `src/search/` | 90% |
| `src/sources/` | 85% |
| `src/cli/commands/` | 90% |
| `src/embeddings/` | 85% |
| `src/utils/` | 90% |
| **Overall** | **90%** |

## Test Organization

```
tests/
├── unit/
│   ├── db/
│   │   ├── client.test.ts
│   │   ├── schema.test.ts
│   │   └── migrations.test.ts
│   ├── search/
│   │   ├── vector.test.ts
│   │   ├── fulltext.test.ts
│   │   └── hybrid.test.ts
│   ├── sources/
│   │   ├── openalex.test.ts
│   │   ├── crossref.test.ts
│   │   └── doi.test.ts
│   ├── cli/
│   │   ├── add.test.ts
│   │   ├── search.test.ts
│   │   └── ... (one per command)
│   └── utils/
│       ├── bibtex.test.ts
│       └── config.test.ts
├── integration/
│   ├── database.test.ts
│   ├── search-pipeline.test.ts
│   └── cli-workflows.test.ts
├── e2e/
│   └── full-workflow.test.ts
└── fixtures/
    ├── works.json
    ├── sample.bib
    └── openalex-response.json
```

## Vitest Configuration

```typescript
// vitest.config.ts
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
      exclude: [
        'src/**/*.d.ts',
        'src/index.ts',  // Re-exports only
      ],
      thresholds: {
        global: {
          branches: 85,
          functions: 90,
          lines: 90,
          statements: 90,
        },
      },
    },
    setupFiles: ['tests/setup.ts'],
  },
});
```

## Test Setup

```typescript
// tests/setup.ts
import { beforeAll, afterAll, beforeEach } from 'vitest';
import { createTestDatabase, cleanupTestDatabase } from './helpers';

beforeAll(async () => {
  // One-time setup
});

afterAll(async () => {
  await cleanupTestDatabase();
});

beforeEach(async () => {
  // Reset database state
  await createTestDatabase();
});
```

## Test Helpers

```typescript
// tests/helpers/index.ts
import { DatabaseClient } from '../../src/db/client';
import { tmpdir } from 'os';
import { join } from 'path';
import { randomUUID } from 'crypto';

export async function createTestDatabase(): Promise<DatabaseClient> {
  const dbPath = join(tmpdir(), `litdb-test-${randomUUID()}.db`);
  const client = new DatabaseClient(dbPath);
  await client.migrate();
  return client;
}

export function mockWork(overrides: Partial<Work> = {}): Work {
  return {
    id: randomUUID(),
    title: 'Test Paper',
    abstract: 'This is a test abstract',
    doi: `10.1234/test-${randomUUID().slice(0, 8)}`,
    publicationDate: new Date('2024-01-15'),
    authors: [],
    createdAt: new Date(),
    updatedAt: new Date(),
    ...overrides,
  };
}

export function mockOpenAlexResponse(): OpenAlexWork {
  return {
    id: 'https://openalex.org/W1234567890',
    doi: 'https://doi.org/10.1234/test',
    title: 'Test Paper from OpenAlex',
    // ... full mock response
  };
}
```

## Unit Test Examples

### Database Tests

```typescript
// tests/unit/db/client.test.ts
import { describe, it, expect, beforeEach } from 'vitest';
import { DatabaseClient } from '../../../src/db/client';
import { createTestDatabase, mockWork } from '../../helpers';

describe('DatabaseClient', () => {
  let db: DatabaseClient;

  beforeEach(async () => {
    db = await createTestDatabase();
  });

  describe('addWork', () => {
    it('adds a work and returns it with ID', async () => {
      const work = mockWork();
      const result = await db.addWork(work);

      expect(result.id).toBeDefined();
      expect(result.title).toBe(work.title);
      expect(result.doi).toBe(work.doi);
    });

    it('rejects duplicate DOI', async () => {
      const work = mockWork();
      await db.addWork(work);

      await expect(db.addWork(work)).rejects.toThrow(/unique constraint/i);
    });

    it('handles work without DOI', async () => {
      const work = mockWork({ doi: undefined });
      const result = await db.addWork(work);

      expect(result.id).toBeDefined();
      expect(result.doi).toBeUndefined();
    });

    it('handles Unicode in title/abstract', async () => {
      const work = mockWork({
        title: '量子コンピューティング研究',
        abstract: 'Résumé avec des accents français',
      });
      const result = await db.addWork(work);

      expect(result.title).toBe(work.title);
      expect(result.abstract).toBe(work.abstract);
    });
  });

  describe('getWork', () => {
    it('returns work by ID', async () => {
      const work = await db.addWork(mockWork());
      const result = await db.getWork(work.id);

      expect(result).toEqual(work);
    });

    it('returns null for non-existent ID', async () => {
      const result = await db.getWork('non-existent-id');

      expect(result).toBeNull();
    });
  });

  describe('deleteWork', () => {
    it('removes work from database', async () => {
      const work = await db.addWork(mockWork());
      await db.deleteWork(work.id);

      const result = await db.getWork(work.id);
      expect(result).toBeNull();
    });

    it('cascades to embeddings', async () => {
      const work = await db.addWork(mockWork());
      await db.addEmbedding({
        workId: work.id,
        chunkIndex: 0,
        chunkText: 'test',
        vector: new Float32Array(384),
      });

      await db.deleteWork(work.id);

      const embeddings = await db.getWorkEmbeddings(work.id);
      expect(embeddings).toHaveLength(0);
    });
  });
});
```

### Search Tests

```typescript
// tests/unit/search/fulltext.test.ts
import { describe, it, expect, beforeEach } from 'vitest';
import { FullTextSearch } from '../../../src/search/fulltext';
import { createTestDatabase, mockWork } from '../../helpers';

describe('FullTextSearch', () => {
  let db: DatabaseClient;
  let search: FullTextSearch;

  beforeEach(async () => {
    db = await createTestDatabase();
    search = new FullTextSearch(db);

    // Add test data
    await db.addWork(mockWork({
      title: 'Machine Learning for Drug Discovery',
      abstract: 'We present a neural network approach...',
    }));
    await db.addWork(mockWork({
      title: 'Quantum Computing Advances',
      abstract: 'Recent progress in quantum algorithms...',
    }));
  });

  it('finds works by title keyword', async () => {
    const results = await search.search('machine learning');

    expect(results).toHaveLength(1);
    expect(results[0].title).toContain('Machine Learning');
  });

  it('finds works by abstract keyword', async () => {
    const results = await search.search('neural network');

    expect(results).toHaveLength(1);
  });

  it('returns empty for no matches', async () => {
    const results = await search.search('nonexistent term');

    expect(results).toHaveLength(0);
  });

  it('respects limit parameter', async () => {
    // Add many works
    for (let i = 0; i < 20; i++) {
      await db.addWork(mockWork({ title: `Paper about testing ${i}` }));
    }

    const results = await search.search('testing', { limit: 5 });

    expect(results).toHaveLength(5);
  });

  it('ranks by relevance', async () => {
    await db.addWork(mockWork({
      title: 'Machine Learning Machine Learning',
      abstract: 'Machine learning is great',
    }));

    const results = await search.search('machine learning');

    // Work with more mentions should rank higher
    expect(results[0].title).toBe('Machine Learning Machine Learning');
  });
});
```

### CLI Tests

```typescript
// tests/unit/cli/add.test.ts
import { describe, it, expect, vi, beforeEach } from 'vitest';
import { runCli } from '../../helpers/cli';

describe('litdb add', () => {
  beforeEach(() => {
    vi.clearAllMocks();
  });

  it('adds work by DOI', async () => {
    const { exitCode, stdout } = await runCli(['add', '10.1038/nature12373']);

    expect(exitCode).toBe(0);
    expect(stdout).toContain('Added:');
  });

  it('accepts DOI URL format', async () => {
    const { exitCode } = await runCli([
      'add',
      'https://doi.org/10.1038/nature12373',
    ]);

    expect(exitCode).toBe(0);
  });

  it('rejects invalid DOI', async () => {
    const { exitCode, stderr } = await runCli(['add', 'not-a-valid-doi']);

    expect(exitCode).toBe(1);
    expect(stderr).toContain('Invalid DOI');
  });

  it('handles network errors gracefully', async () => {
    vi.mock('../../../src/sources/doi', () => ({
      resolveDOI: vi.fn().mockRejectedValue(new Error('Network error')),
    }));

    const { exitCode, stderr } = await runCli(['add', '10.1038/nature12373']);

    expect(exitCode).toBe(1);
    expect(stderr).toContain('Network error');
  });

  it('supports --tag option', async () => {
    const { exitCode } = await runCli([
      'add',
      '10.1038/nature12373',
      '--tag',
      'important',
    ]);

    expect(exitCode).toBe(0);
    // Verify tag was applied
  });

  it('shows progress for --fetch-refs', async () => {
    const { stdout } = await runCli([
      'add',
      '10.1038/nature12373',
      '--fetch-refs',
    ]);

    expect(stdout).toContain('Fetching references');
  });
});
```

## Integration Tests

```typescript
// tests/integration/search-pipeline.test.ts
import { describe, it, expect, beforeEach } from 'vitest';
import { LitDB } from '../../src';

describe('Search Pipeline Integration', () => {
  let litdb: LitDB;

  beforeEach(async () => {
    litdb = await LitDB.create({ path: ':memory:' });

    // Add realistic test data
    await litdb.add('10.1038/nature12373');
    await litdb.add('10.1126/science.1234567');
    await litdb.index(); // Generate embeddings
  });

  it('vector search finds semantically similar works', async () => {
    const results = await litdb.vsearch(
      'gene editing techniques for human cells'
    );

    expect(results.length).toBeGreaterThan(0);
    // CRISPR paper should be highly ranked
  });

  it('hybrid search combines vector and fulltext', async () => {
    const results = await litdb.search('CRISPR gene editing');

    expect(results.length).toBeGreaterThan(0);
    expect(results[0].score).toBeGreaterThan(0);
  });

  it('search respects tag filter', async () => {
    await litdb.tag('10.1038/nature12373', 'crispr');

    const results = await litdb.search('gene', { tag: 'crispr' });

    expect(results).toHaveLength(1);
  });
});
```

## E2E Tests

```typescript
// tests/e2e/full-workflow.test.ts
import { describe, it, expect } from 'vitest';
import { execSync } from 'child_process';
import { tmpdir } from 'os';
import { join } from 'path';

describe('Full Workflow E2E', () => {
  const testDir = join(tmpdir(), 'litdb-e2e-test');

  it('complete workflow: init -> add -> search -> export', async () => {
    // Initialize
    execSync(`litdb init --path ${testDir}/test.db`);

    // Add a paper
    execSync(`litdb add 10.1038/nature12373 --db ${testDir}/test.db`);

    // Search
    const searchResult = execSync(
      `litdb search "gene editing" --db ${testDir}/test.db --json`
    ).toString();
    const results = JSON.parse(searchResult);
    expect(results.length).toBeGreaterThan(0);

    // Export
    const bibtex = execSync(
      `litdb export --format bibtex --db ${testDir}/test.db`
    ).toString();
    expect(bibtex).toContain('@article{');
  });
});
```

## Mocking Strategy

### External APIs

```typescript
// tests/mocks/openalex.ts
import { vi } from 'vitest';

export const mockOpenAlexApi = () => {
  return vi.fn().mockImplementation((url: string) => {
    if (url.includes('/works/')) {
      return Promise.resolve({
        ok: true,
        json: () => Promise.resolve(mockOpenAlexWork()),
      });
    }
    // ... other endpoints
  });
};
```

### Embeddings

```typescript
// tests/mocks/embeddings.ts
export const mockEmbeddingModel = {
  embed: vi.fn().mockResolvedValue(new Float32Array(384).fill(0.1)),
};
```

## CI Configuration

```yaml
# .github/workflows/test.yml
name: Tests

on: [push, pull_request]

jobs:
  test:
    runs-on: ubuntu-latest
    steps:
      - uses: actions/checkout@v4
      - uses: actions/setup-node@v4
        with:
          node-version: '20'
      - run: npm ci
      - run: npm run test:coverage
      - uses: codecov/codecov-action@v3
        with:
          files: ./coverage/lcov.info
```

## Test Commands

```json
{
  "scripts": {
    "test": "vitest run",
    "test:watch": "vitest",
    "test:coverage": "vitest run --coverage",
    "test:ui": "vitest --ui"
  }
}
```
