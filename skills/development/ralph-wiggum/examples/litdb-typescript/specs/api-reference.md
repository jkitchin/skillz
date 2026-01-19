# Public API Reference

## Overview

litdb-ts exposes two interfaces:
1. **Programmatic API** - Import and use in Node.js applications
2. **CLI** - Command-line interface

This document covers the programmatic API.

## Installation

```bash
npm install litdb-ts
```

## Quick Start

```typescript
import { LitDB } from 'litdb-ts';

// Create instance
const db = await LitDB.create({
  path: './my-papers.db',  // Or ':memory:' for testing
});

// Add a paper
const work = await db.add('10.1038/nature12373');
console.log(work.title);  // "Genome engineering with CRISPR-Cas9"

// Search
const results = await db.search('gene editing');
for (const result of results) {
  console.log(`${result.work.title} (score: ${result.score})`);
}

// Clean up
db.close();
```

## LitDB Class

### Constructor Options

```typescript
interface LitDBOptions {
  path: string;              // Database file path
  config?: ConfigOptions;    // Override config file
  embeddings?: {
    model?: string;          // Default: 'all-MiniLM-L6-v2'
    chunkSize?: number;      // Default: 512
    chunkOverlap?: number;   // Default: 50
  };
}

const db = await LitDB.create(options);
```

### Work Management

```typescript
// Add work by DOI
add(doi: string, options?: AddOptions): Promise<Work>

interface AddOptions {
  tags?: string[];
  fetchRefs?: boolean;
  fetchCitations?: boolean;
}

// Get work
get(id: string): Promise<Work | null>
getByDOI(doi: string): Promise<Work | null>

// Update work metadata
update(id: string, updates: Partial<Work>): Promise<Work>

// Remove work
remove(id: string): Promise<void>

// List works
list(options?: ListOptions): Promise<Work[]>

interface ListOptions {
  limit?: number;
  offset?: number;
  tag?: string;
  since?: Date;
  orderBy?: 'date' | 'title' | 'added';
}
```

### Search

```typescript
// Hybrid search (default)
search(query: string, options?: SearchOptions): Promise<SearchResult[]>

// Vector-only search
vsearch(query: string, options?: SearchOptions): Promise<SearchResult[]>

// Full-text only search
fulltext(query: string, options?: SearchOptions): Promise<SearchResult[]>

// Find similar works
similar(workId: string, limit?: number): Promise<Work[]>

interface SearchOptions {
  limit?: number;
  offset?: number;
  tag?: string;
  since?: Date;
  ftsWeight?: number;     // For hybrid
  vectorWeight?: number;  // For hybrid
}

interface SearchResult {
  workId: string;
  work: Work;
  score: number;
  ftsScore?: number;
  vectorScore?: number;
  matchedChunks?: string[];
}
```

### Tags

```typescript
// Create tag
createTag(name: string, color?: string): Promise<Tag>

// Tag a work
tag(workId: string, tagName: string): Promise<void>

// Remove tag from work
untag(workId: string, tagName: string): Promise<void>

// List tags
listTags(): Promise<Tag[]>

// Get works with tag
getTaggedWorks(tagName: string): Promise<Work[]>
```

### Indexing

```typescript
// Rebuild all indexes
index(): Promise<void>

// Index specific work
indexWork(workId: string): Promise<void>

// Check if work is indexed
isIndexed(workId: string): Promise<boolean>
```

### Filters

```typescript
// Create automatic filter
createFilter(name: string, query: string, source?: string): Promise<Filter>

// Run filter to fetch new papers
runFilter(name: string): Promise<Work[]>

// Run all active filters
runFilters(): Promise<{ [name: string]: Work[] }>
```

### Import/Export

```typescript
// Import from BibTeX
importBibtex(content: string): Promise<Work[]>
importBibtexFile(path: string): Promise<Work[]>

// Export to BibTeX
exportBibtex(workIds?: string[]): Promise<string>

// Export to JSON
exportJSON(workIds?: string[]): Promise<string>
```

### Statistics

```typescript
stats(): Promise<DatabaseStats>

interface DatabaseStats {
  totalWorks: number;
  totalAuthors: number;
  totalTags: number;
  totalEmbeddings: number;
  databaseSize: number;  // Bytes
  oldestWork?: Date;
  newestWork?: Date;
}
```

### Lifecycle

```typescript
// Close database connection
close(): void

// Backup database
backup(path: string): Promise<void>

// Vacuum database
vacuum(): Promise<void>
```

## Data Types

### Work

```typescript
interface Work {
  id: string;
  doi?: string;
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

### Tag

```typescript
interface Tag {
  id: string;
  name: string;
  color?: string;
  workCount?: number;
}
```

### Filter

```typescript
interface Filter {
  id: string;
  name: string;
  query: string;
  source: string;
  lastRun?: Date;
  isActive: boolean;
}
```

## Events

```typescript
const db = await LitDB.create({ path: './db.sqlite' });

// Work added
db.on('work:added', (work: Work) => {
  console.log(`Added: ${work.title}`);
});

// Work removed
db.on('work:removed', (id: string) => {
  console.log(`Removed: ${id}`);
});

// Index updated
db.on('index:updated', (workId: string) => {
  console.log(`Indexed: ${workId}`);
});

// Filter run
db.on('filter:run', (name: string, works: Work[]) => {
  console.log(`Filter ${name} found ${works.length} new works`);
});
```

## Error Handling

```typescript
import { LitDB, LitDBError, WorkNotFoundError, InvalidDOIError } from 'litdb-ts';

try {
  await db.add('invalid-doi');
} catch (error) {
  if (error instanceof InvalidDOIError) {
    console.error('Invalid DOI format');
  } else if (error instanceof LitDBError) {
    console.error(`LitDB error: ${error.code}`);
  } else {
    throw error;
  }
}
```

### Error Types

```typescript
class LitDBError extends Error {
  code: string;
}

class WorkNotFoundError extends LitDBError {
  code = 'WORK_NOT_FOUND';
  workId: string;
}

class InvalidDOIError extends LitDBError {
  code = 'INVALID_DOI';
  doi: string;
}

class EmbeddingError extends LitDBError {
  code = 'EMBEDDING_FAILED';
}

class NetworkError extends LitDBError {
  code = 'NETWORK_ERROR';
}

class DuplicateWorkError extends LitDBError {
  code = 'DUPLICATE_WORK';
}
```

## Configuration

### Config File

`litdb.toml` in current directory or `~/.config/litdb/config.toml`:

```toml
[database]
path = "~/.litdb/papers.db"

[embeddings]
model = "all-MiniLM-L6-v2"
chunk_size = 512
chunk_overlap = 50

[openalex]
email = "your@email.com"

[search]
default_limit = 10
fts_weight = 0.3
vector_weight = 0.7
```

### Programmatic Config

```typescript
const db = await LitDB.create({
  path: './db.sqlite',
  config: {
    embeddings: {
      model: 'all-MiniLM-L6-v2',
      chunkSize: 256,  // Smaller chunks
    },
    search: {
      defaultLimit: 20,
    },
  },
});
```

## Examples

### Research Workflow

```typescript
import { LitDB } from 'litdb-ts';

async function researchWorkflow() {
  const db = await LitDB.create({ path: './research.db' });

  // Add seed papers
  await db.add('10.1038/nature12373');  // CRISPR paper
  await db.tag('10.1038/nature12373', 'seed');

  // Index for search
  await db.index();

  // Find related work
  const similar = await db.similar('10.1038/nature12373', 20);
  console.log(`Found ${similar.length} related papers`);

  // Search with query
  const results = await db.search('gene editing efficiency');
  for (const r of results.slice(0, 5)) {
    console.log(`- ${r.work.title} (${r.score.toFixed(2)})`);
    await db.tag(r.workId, 'to-read');
  }

  // Export reading list
  const toRead = await db.getTaggedWorks('to-read');
  const bibtex = await db.exportBibtex(toRead.map(w => w.id));
  console.log(bibtex);

  db.close();
}
```

### Automated Updates

```typescript
import { LitDB } from 'litdb-ts';

async function setupAutoUpdates() {
  const db = await LitDB.create({ path: './papers.db' });

  // Create filter for Nature papers on ML
  await db.createFilter(
    'nature-ml',
    'machine learning journal:Nature',
    'openalex'
  );

  // Run daily (e.g., via cron)
  const newPapers = await db.runFilters();
  for (const [filter, works] of Object.entries(newPapers)) {
    console.log(`${filter}: ${works.length} new papers`);
    for (const work of works) {
      await db.tag(work.id, 'inbox');
    }
  }

  db.close();
}
```
