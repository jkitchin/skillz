# Database Schema Specification

## Overview

litdb-ts uses SQLite with:
- Standard relational tables for metadata
- FTS5 virtual tables for full-text search
- Vector extension (sqlite-vss) for semantic search

## Tables

### works

Primary table for papers/articles.

```sql
CREATE TABLE works (
  id TEXT PRIMARY KEY,
  doi TEXT UNIQUE,
  title TEXT NOT NULL,
  abstract TEXT,
  publication_date TEXT,  -- ISO 8601 date
  journal TEXT,
  volume TEXT,
  issue TEXT,
  pages TEXT,
  url TEXT,
  pdf_url TEXT,
  openalex_id TEXT UNIQUE,
  raw_metadata TEXT,      -- JSON blob for extra fields
  created_at TEXT NOT NULL DEFAULT (datetime('now')),
  updated_at TEXT NOT NULL DEFAULT (datetime('now'))
);

CREATE INDEX idx_works_doi ON works(doi);
CREATE INDEX idx_works_openalex_id ON works(openalex_id);
CREATE INDEX idx_works_publication_date ON works(publication_date);
```

### authors

```sql
CREATE TABLE authors (
  id TEXT PRIMARY KEY,
  name TEXT NOT NULL,
  orcid TEXT UNIQUE,
  affiliation TEXT,
  openalex_id TEXT UNIQUE,
  created_at TEXT NOT NULL DEFAULT (datetime('now'))
);

CREATE INDEX idx_authors_orcid ON authors(orcid);
CREATE INDEX idx_authors_name ON authors(name);
```

### work_authors

Junction table for work-author relationships.

```sql
CREATE TABLE work_authors (
  work_id TEXT NOT NULL REFERENCES works(id) ON DELETE CASCADE,
  author_id TEXT NOT NULL REFERENCES authors(id) ON DELETE CASCADE,
  position INTEGER NOT NULL,  -- Author order (1 = first author)
  PRIMARY KEY (work_id, author_id)
);
```

### tags

User-defined tags for organizing works.

```sql
CREATE TABLE tags (
  id TEXT PRIMARY KEY,
  name TEXT NOT NULL UNIQUE,
  color TEXT,  -- Hex color code
  created_at TEXT NOT NULL DEFAULT (datetime('now'))
);
```

### work_tags

```sql
CREATE TABLE work_tags (
  work_id TEXT NOT NULL REFERENCES works(id) ON DELETE CASCADE,
  tag_id TEXT NOT NULL REFERENCES tags(id) ON DELETE CASCADE,
  PRIMARY KEY (work_id, tag_id)
);
```

### embeddings

Store vector embeddings for semantic search.

```sql
CREATE TABLE embeddings (
  id TEXT PRIMARY KEY,
  work_id TEXT NOT NULL REFERENCES works(id) ON DELETE CASCADE,
  chunk_index INTEGER NOT NULL,
  chunk_text TEXT NOT NULL,
  model TEXT NOT NULL,  -- e.g., "all-MiniLM-L6-v2"
  created_at TEXT NOT NULL DEFAULT (datetime('now')),
  UNIQUE(work_id, chunk_index, model)
);

-- Vector data stored separately via sqlite-vss
CREATE VIRTUAL TABLE vss_embeddings USING vss0(
  embedding(384)  -- Dimension depends on model
);
```

### filters

Saved filters for automatic updates.

```sql
CREATE TABLE filters (
  id TEXT PRIMARY KEY,
  name TEXT NOT NULL UNIQUE,
  query TEXT NOT NULL,       -- OpenAlex query string
  source TEXT NOT NULL,      -- 'openalex', 'crossref', etc.
  last_run TEXT,
  is_active INTEGER NOT NULL DEFAULT 1,
  created_at TEXT NOT NULL DEFAULT (datetime('now'))
);
```

### notes

User notes attached to works.

```sql
CREATE TABLE notes (
  id TEXT PRIMARY KEY,
  work_id TEXT NOT NULL REFERENCES works(id) ON DELETE CASCADE,
  content TEXT NOT NULL,
  created_at TEXT NOT NULL DEFAULT (datetime('now')),
  updated_at TEXT NOT NULL DEFAULT (datetime('now'))
);
```

## Full-Text Search

```sql
-- FTS5 virtual table for full-text search
CREATE VIRTUAL TABLE works_fts USING fts5(
  title,
  abstract,
  content='works',
  content_rowid='rowid'
);

-- Triggers to keep FTS in sync
CREATE TRIGGER works_ai AFTER INSERT ON works BEGIN
  INSERT INTO works_fts(rowid, title, abstract)
  VALUES (NEW.rowid, NEW.title, NEW.abstract);
END;

CREATE TRIGGER works_ad AFTER DELETE ON works BEGIN
  INSERT INTO works_fts(works_fts, rowid, title, abstract)
  VALUES ('delete', OLD.rowid, OLD.title, OLD.abstract);
END;

CREATE TRIGGER works_au AFTER UPDATE ON works BEGIN
  INSERT INTO works_fts(works_fts, rowid, title, abstract)
  VALUES ('delete', OLD.rowid, OLD.title, OLD.abstract);
  INSERT INTO works_fts(rowid, title, abstract)
  VALUES (NEW.rowid, NEW.title, NEW.abstract);
END;
```

## Migration System

```typescript
interface Migration {
  version: number;
  name: string;
  up: string;    // SQL to apply
  down: string;  // SQL to rollback
}

// Track applied migrations
CREATE TABLE migrations (
  version INTEGER PRIMARY KEY,
  name TEXT NOT NULL,
  applied_at TEXT NOT NULL DEFAULT (datetime('now'))
);
```

## Database Client API

```typescript
interface DatabaseClient {
  // Works
  addWork(work: NewWork): Promise<Work>;
  getWork(id: string): Promise<Work | null>;
  getWorkByDOI(doi: string): Promise<Work | null>;
  updateWork(id: string, updates: Partial<Work>): Promise<Work>;
  deleteWork(id: string): Promise<void>;
  listWorks(options?: ListOptions): Promise<Work[]>;

  // Authors
  addAuthor(author: NewAuthor): Promise<Author>;
  getAuthor(id: string): Promise<Author | null>;
  getAuthorByORCID(orcid: string): Promise<Author | null>;

  // Tags
  createTag(name: string, color?: string): Promise<Tag>;
  addTagToWork(workId: string, tagId: string): Promise<void>;
  removeTagFromWork(workId: string, tagId: string): Promise<void>;
  getWorkTags(workId: string): Promise<Tag[]>;

  // Embeddings
  addEmbedding(embedding: NewEmbedding): Promise<void>;
  getWorkEmbeddings(workId: string): Promise<Embedding[]>;

  // Search
  fullTextSearch(query: string, limit?: number): Promise<SearchResult[]>;
  vectorSearch(embedding: Float32Array, limit?: number): Promise<SearchResult[]>;

  // Filters
  createFilter(filter: NewFilter): Promise<Filter>;
  runFilter(filterId: string): Promise<Work[]>;

  // Utilities
  getStats(): Promise<DatabaseStats>;
  vacuum(): Promise<void>;
  close(): void;
}
```

## Testing Requirements

Each database operation must have tests for:
1. **Happy path** - Normal operation succeeds
2. **Constraints** - Unique violations, foreign keys
3. **Edge cases** - Empty strings, null values, Unicode
4. **Transactions** - Rollback on error

Target: 95% coverage of `src/db/` directory.
