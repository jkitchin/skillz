# litdb-ts: Architecture Overview

## Project Goals

Convert the Python [litdb](https://github.com/jkitchin/litdb) project to TypeScript, creating:

1. **@litdb/core** - Core library for programmatic use
2. **litdb CLI** - Command-line interface matching Python version
3. **High test coverage** - 90%+ code coverage target

## Technology Stack

| Component | Technology | Rationale |
|-----------|------------|-----------|
| Language | TypeScript 5.x | Type safety, modern JS |
| Runtime | Node.js 20+ | LTS, native fetch, ESM |
| Database | better-sqlite3 | Fast, synchronous SQLite |
| Vector DB | sqlite-vss or vectra | Vector similarity search |
| Embeddings | @xenova/transformers | Local embeddings (no API) |
| CLI Framework | Commander.js | Industry standard |
| Testing | Vitest | Fast, TypeScript-native |
| Build | tsup | Simple, fast bundling |

## Package Structure

```
litdb-ts/
├── package.json
├── tsconfig.json
├── vitest.config.ts
├── src/
│   ├── index.ts              # Public API exports
│   ├── cli/
│   │   ├── index.ts          # CLI entry point
│   │   └── commands/         # Individual commands
│   ├── db/
│   │   ├── schema.ts         # Database schema
│   │   ├── migrations.ts     # Schema migrations
│   │   └── client.ts         # Database client
│   ├── search/
│   │   ├── vector.ts         # Vector search
│   │   ├── fulltext.ts       # FTS5 search
│   │   └── hybrid.ts         # Combined search
│   ├── sources/
│   │   ├── openalex.ts       # OpenAlex API
│   │   ├── crossref.ts       # CrossRef API
│   │   └── doi.ts            # DOI resolution
│   ├── embeddings/
│   │   └── transformer.ts    # Embedding generation
│   └── utils/
│       ├── config.ts         # Configuration loading
│       ├── bibtex.ts         # BibTeX parsing/export
│       └── logger.ts         # Logging utilities
├── tests/
│   ├── unit/                 # Unit tests (90%+ coverage)
│   ├── integration/          # Integration tests
│   └── fixtures/             # Test data
└── bin/
    └── litdb.js              # CLI executable
```

## Core Entities

### Work (Paper/Article)

```typescript
interface Work {
  id: string;                    // Internal UUID
  doi?: string;                  // DOI if available
  title: string;
  abstract?: string;
  authors: Author[];
  publicationDate?: Date;
  journal?: string;
  volume?: string;
  issue?: string;
  pages?: string;
  url?: string;
  pdfUrl?: string;
  openalexId?: string;
  createdAt: Date;
  updatedAt: Date;
}
```

### Author

```typescript
interface Author {
  id: string;
  name: string;
  orcid?: string;
  affiliation?: string;
  openalexId?: string;
}
```

### Embedding

```typescript
interface Embedding {
  workId: string;
  chunkIndex: number;
  chunkText: string;
  vector: Float32Array;
}
```

## API Design Principles

1. **Async by default** - All I/O operations return Promises
2. **Type-safe** - Full TypeScript types, no `any`
3. **Testable** - Dependency injection, mockable interfaces
4. **Progressive** - Core features work offline, enhanced with APIs

## Configuration

Config file: `litdb.toml` (compatible with Python version)

```toml
[database]
path = "~/.litdb/litdb.db"

[embeddings]
model = "all-MiniLM-L6-v2"
chunk_size = 512
chunk_overlap = 50

[openalex]
email = "your@email.com"  # For polite pool

[llm]
provider = "ollama"
model = "llama2"
```

## Error Handling

```typescript
// Custom error types
class LitDBError extends Error {
  constructor(message: string, public code: string) {
    super(message);
  }
}

class WorkNotFoundError extends LitDBError { }
class InvalidDOIError extends LitDBError { }
class EmbeddingError extends LitDBError { }
```

## Logging

Use structured logging with levels:

```typescript
const logger = createLogger('litdb');
logger.info('Adding work', { doi: '10.1234/example' });
logger.error('Failed to fetch', { error, url });
```

## Performance Targets

| Operation | Target |
|-----------|--------|
| Add work by DOI | < 2s |
| Vector search (1000 works) | < 100ms |
| Full-text search | < 50ms |
| Generate embedding | < 500ms |
| CLI startup | < 200ms |
