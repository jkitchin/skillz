# Search Functionality Specification

## Overview

litdb-ts supports three search modes:
1. **Full-text search** - Traditional keyword search using SQLite FTS5
2. **Vector search** - Semantic search using embeddings
3. **Hybrid search** - Combines both for best results

## Full-Text Search (FTS5)

### Query Syntax

```typescript
// Simple query
fulltext("machine learning")

// Phrase query
fulltext('"neural network"')

// Field-specific
fulltext("title:quantum")
fulltext("abstract:catalyst")

// Boolean operators
fulltext("machine AND learning")
fulltext("battery OR energy")
fulltext("deep learning NOT vision")

// Prefix matching
fulltext("neuro*")  // neuroscience, neurology, etc.

// NEAR operator
fulltext("NEAR(climate change, 5)")  // Within 5 tokens
```

### Implementation

```typescript
// src/search/fulltext.ts
export interface FullTextOptions {
  limit?: number;
  offset?: number;
  fields?: ('title' | 'abstract')[];
  highlight?: boolean;
}

export interface FullTextResult {
  workId: string;
  score: number;
  highlights?: {
    title?: string;
    abstract?: string;
  };
}

export class FullTextSearch {
  constructor(private db: DatabaseClient) {}

  async search(query: string, options?: FullTextOptions): Promise<FullTextResult[]> {
    const limit = options?.limit ?? 10;
    const offset = options?.offset ?? 0;

    // Build FTS5 query
    const ftsQuery = this.buildQuery(query, options?.fields);

    const results = this.db.prepare(`
      SELECT
        w.id as workId,
        bm25(works_fts) as score
        ${options?.highlight ? ', highlight(works_fts, 0, "<mark>", "</mark>") as title_hl' : ''}
        ${options?.highlight ? ', highlight(works_fts, 1, "<mark>", "</mark>") as abstract_hl' : ''}
      FROM works_fts
      JOIN works w ON works_fts.rowid = w.rowid
      WHERE works_fts MATCH ?
      ORDER BY score
      LIMIT ? OFFSET ?
    `).all(ftsQuery, limit, offset);

    return results;
  }

  private buildQuery(query: string, fields?: string[]): string {
    // Handle field-specific queries
    if (fields && fields.length > 0) {
      return fields.map(f => `${f}:${query}`).join(' OR ');
    }
    return query;
  }
}
```

### Scoring

FTS5 uses BM25 algorithm. Lower scores = better match.

```sql
-- Default BM25 weights: (title, abstract)
SELECT bm25(works_fts, 1.0, 0.75) as score
-- Title matches weighted higher than abstract
```

## Vector Search

### Embedding Model

Default: `all-MiniLM-L6-v2` (384 dimensions)
- Fast inference
- Good for English scientific text
- ~80MB model size

### Text Chunking

```typescript
interface ChunkOptions {
  chunkSize: number;    // Default: 512 tokens
  chunkOverlap: number; // Default: 50 tokens
}

function chunkText(text: string, options: ChunkOptions): string[] {
  // Split on sentence boundaries when possible
  // Ensure overlap for context continuity
}
```

### Vector Operations

```typescript
// src/search/vector.ts
export interface VectorSearchOptions {
  limit?: number;
  threshold?: number;  // Minimum similarity score
  workIds?: string[];  // Filter to specific works
}

export interface VectorSearchResult {
  workId: string;
  chunkIndex: number;
  chunkText: string;
  score: number;  // Cosine similarity 0-1
}

export class VectorSearch {
  constructor(
    private db: DatabaseClient,
    private embedder: EmbeddingModel
  ) {}

  async search(query: string, options?: VectorSearchOptions): Promise<VectorSearchResult[]> {
    // Generate query embedding
    const queryEmbedding = await this.embedder.embed(query);

    // Search using sqlite-vss
    const results = this.db.prepare(`
      SELECT
        e.work_id as workId,
        e.chunk_index as chunkIndex,
        e.chunk_text as chunkText,
        vss_distance_cosine(v.embedding, ?) as distance
      FROM vss_embeddings v
      JOIN embeddings e ON v.rowid = e.rowid
      ${options?.workIds ? 'WHERE e.work_id IN (SELECT value FROM json_each(?))' : ''}
      ORDER BY distance ASC
      LIMIT ?
    `).all(
      queryEmbedding,
      options?.workIds ? JSON.stringify(options.workIds) : null,
      options?.limit ?? 10
    );

    // Convert distance to similarity score
    return results.map(r => ({
      ...r,
      score: 1 - r.distance,  // Cosine: distance 0 = identical
    }));
  }

  async indexWork(workId: string): Promise<void> {
    const work = await this.db.getWork(workId);
    if (!work) throw new WorkNotFoundError(workId);

    // Chunk text
    const text = `${work.title}\n\n${work.abstract ?? ''}`;
    const chunks = chunkText(text, { chunkSize: 512, chunkOverlap: 50 });

    // Generate embeddings
    for (let i = 0; i < chunks.length; i++) {
      const embedding = await this.embedder.embed(chunks[i]);
      await this.db.addEmbedding({
        workId,
        chunkIndex: i,
        chunkText: chunks[i],
        vector: embedding,
      });
    }
  }
}
```

## Hybrid Search

Combines FTS5 and vector search with weighted scoring.

```typescript
// src/search/hybrid.ts
export interface HybridSearchOptions {
  limit?: number;
  ftsWeight?: number;    // Default: 0.3
  vectorWeight?: number; // Default: 0.7
  tag?: string;
  since?: Date;
}

export interface HybridSearchResult {
  workId: string;
  work: Work;
  ftsScore?: number;
  vectorScore?: number;
  combinedScore: number;
  matchedChunks?: string[];
}

export class HybridSearch {
  constructor(
    private fts: FullTextSearch,
    private vector: VectorSearch
  ) {}

  async search(query: string, options?: HybridSearchOptions): Promise<HybridSearchResult[]> {
    const ftsWeight = options?.ftsWeight ?? 0.3;
    const vectorWeight = options?.vectorWeight ?? 0.7;
    const limit = options?.limit ?? 10;

    // Run both searches in parallel
    const [ftsResults, vectorResults] = await Promise.all([
      this.fts.search(query, { limit: limit * 2 }),
      this.vector.search(query, { limit: limit * 2 }),
    ]);

    // Merge and score results
    const scoreMap = new Map<string, { fts: number; vector: number }>();

    // Normalize FTS scores (BM25 is unbounded, lower is better)
    const maxFts = Math.max(...ftsResults.map(r => Math.abs(r.score)), 1);
    for (const result of ftsResults) {
      const normalized = 1 - (Math.abs(result.score) / maxFts);
      scoreMap.set(result.workId, {
        fts: normalized,
        vector: scoreMap.get(result.workId)?.vector ?? 0
      });
    }

    // Vector scores are already 0-1
    for (const result of vectorResults) {
      const existing = scoreMap.get(result.workId);
      if (existing) {
        existing.vector = Math.max(existing.vector, result.score);
      } else {
        scoreMap.set(result.workId, { fts: 0, vector: result.score });
      }
    }

    // Calculate combined scores
    const combined: HybridSearchResult[] = [];
    for (const [workId, scores] of scoreMap) {
      const work = await this.db.getWork(workId);
      if (!work) continue;

      // Apply filters
      if (options?.tag) {
        const tags = await this.db.getWorkTags(workId);
        if (!tags.some(t => t.name === options.tag)) continue;
      }
      if (options?.since && work.publicationDate < options.since) continue;

      combined.push({
        workId,
        work,
        ftsScore: scores.fts,
        vectorScore: scores.vector,
        combinedScore: (scores.fts * ftsWeight) + (scores.vector * vectorWeight),
      });
    }

    // Sort by combined score
    combined.sort((a, b) => b.combinedScore - a.combinedScore);

    return combined.slice(0, limit);
  }
}
```

## Iterative Search

Expands search through citation network.

```typescript
// src/search/iterative.ts
export interface IterativeSearchOptions {
  depth?: number;       // How many hops (default: 2)
  maxPapers?: number;   // Total paper limit (default: 100)
  includeRefs?: boolean;
  includeCitations?: boolean;
}

export async function iterativeSearch(
  query: string,
  hybrid: HybridSearch,
  options?: IterativeSearchOptions
): Promise<IterativeSearchResult> {
  const depth = options?.depth ?? 2;
  const maxPapers = options?.maxPapers ?? 100;
  const seen = new Set<string>();
  const results: Work[] = [];

  // Initial search
  const initial = await hybrid.search(query, { limit: 10 });
  for (const r of initial) {
    seen.add(r.workId);
    results.push(r.work);
  }

  // Expand through network
  for (let d = 0; d < depth && results.length < maxPapers; d++) {
    const frontier = results.slice(-10);  // Recent additions

    for (const work of frontier) {
      if (results.length >= maxPapers) break;

      // Fetch references
      if (options?.includeRefs !== false) {
        const refs = await fetchReferences(work.id);
        for (const ref of refs) {
          if (!seen.has(ref.id)) {
            seen.add(ref.id);
            results.push(ref);
          }
        }
      }

      // Fetch citations
      if (options?.includeCitations !== false) {
        const citations = await fetchCitations(work.id);
        for (const cit of citations) {
          if (!seen.has(cit.id)) {
            seen.add(cit.id);
            results.push(cit);
          }
        }
      }
    }
  }

  return {
    query,
    results,
    totalFound: results.length,
    depth,
  };
}
```

## Similar Works

Find works similar to a given work.

```typescript
async function findSimilar(workId: string, limit: number = 10): Promise<Work[]> {
  // Get all embeddings for this work
  const embeddings = await db.getWorkEmbeddings(workId);

  if (embeddings.length === 0) {
    throw new Error('Work has no embeddings. Run indexing first.');
  }

  // Search using average embedding
  const avgEmbedding = averageVectors(embeddings.map(e => e.vector));
  const results = await vector.search(avgEmbedding, {
    limit: limit + 1,  // +1 to exclude self
    excludeWorkIds: [workId],
  });

  return results.map(r => r.work);
}
```

## Performance Considerations

### Indexing

- Batch embedding generation for multiple works
- Use transactions for bulk inserts
- Vacuum after large imports

### Search

- FTS5 is very fast (<10ms for 10k docs)
- Vector search depends on index size
- Use LIMIT to avoid scanning entire index
- Consider approximate nearest neighbor for >100k embeddings

### Caching

- Cache query embeddings for repeated searches
- Cache hot results with TTL

## Testing Requirements

### Unit Tests

```typescript
describe('FullTextSearch', () => {
  it('finds exact phrase matches');
  it('handles boolean operators');
  it('respects field prefixes');
  it('returns results in score order');
  it('highlights matched terms');
});

describe('VectorSearch', () => {
  it('finds semantically similar documents');
  it('filters by work IDs');
  it('respects similarity threshold');
});

describe('HybridSearch', () => {
  it('combines FTS and vector scores');
  it('filters by tag');
  it('filters by date');
  it('adjusts weights correctly');
});
```

### Integration Tests

```typescript
describe('Search Integration', () => {
  it('end-to-end: add papers, index, search, find relevant');
  it('iterative search expands through citations');
  it('similar works finds related papers');
});
```
