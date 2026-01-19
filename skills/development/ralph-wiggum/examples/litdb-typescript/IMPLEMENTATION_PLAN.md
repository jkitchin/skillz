# Implementation Plan: litdb-ts

This plan is designed for Ralph Wiggum mode. Each task should be completable in 1-2 iterations.

**Target:** 90%+ test coverage, all 43 CLI commands working.

---

## Phase 1: Project Setup

- [ ] Initialize npm project with TypeScript configuration
- [ ] Set up Vitest with coverage configuration
- [ ] Create directory structure (src/, tests/, bin/)
- [ ] Configure tsup for building
- [ ] Add ESLint and Prettier configuration
- [ ] Create initial package.json scripts
- [ ] Set up test fixtures directory with sample data

## Phase 2: Core Database Layer

- [ ] Implement database schema (src/db/schema.ts)
- [ ] Write migration system (src/db/migrations.ts)
- [ ] Implement DatabaseClient class with basic CRUD
- [ ] Add works table operations (add, get, update, delete, list)
- [ ] Add authors table operations
- [ ] Add work_authors junction table operations
- [ ] Implement tags table operations
- [ ] Implement work_tags operations
- [ ] Add embeddings table operations
- [ ] Add filters table operations
- [ ] Add notes table operations
- [ ] Set up FTS5 virtual table and triggers
- [ ] Write unit tests for all database operations (target: 95%)
- [ ] Write integration tests for database transactions

## Phase 3: Search Functionality

- [ ] Implement FullTextSearch class (src/search/fulltext.ts)
- [ ] Write FTS5 query builder with field targeting
- [ ] Implement VectorSearch class (src/search/vector.ts)
- [ ] Integrate sqlite-vss for vector operations
- [ ] Implement HybridSearch combining both methods
- [ ] Add search result ranking and scoring
- [ ] Implement search filters (date, tags, etc.)
- [ ] Write unit tests for search (target: 90%)
- [ ] Write integration tests for search pipeline

## Phase 4: External Sources

- [ ] Implement DOI resolver (src/sources/doi.ts)
- [ ] Parse DOI from various formats (URL, string, etc.)
- [ ] Implement OpenAlex API client (src/sources/openalex.ts)
- [ ] Add OpenAlex work fetching
- [ ] Add OpenAlex author fetching
- [ ] Add OpenAlex search functionality
- [ ] Implement CrossRef API client (src/sources/crossref.ts)
- [ ] Add reference matching via CrossRef
- [ ] Write mock fixtures for API responses
- [ ] Write unit tests with mocked APIs (target: 85%)

## Phase 5: Embeddings

- [ ] Set up @xenova/transformers integration
- [ ] Implement embedding generation (src/embeddings/transformer.ts)
- [ ] Add text chunking with overlap
- [ ] Implement batch embedding generation
- [ ] Add embedding caching to avoid regeneration
- [ ] Write unit tests for embeddings (target: 85%)

## Phase 6: Utilities

- [ ] Implement configuration loading (src/utils/config.ts)
- [ ] Support litdb.toml format
- [ ] Add config defaults and validation
- [ ] Implement BibTeX parser (src/utils/bibtex.ts)
- [ ] Implement BibTeX generator
- [ ] Add structured logging (src/utils/logger.ts)
- [ ] Write unit tests for utilities (target: 90%)

## Phase 7: CLI Framework

- [ ] Set up Commander.js CLI structure (src/cli/index.ts)
- [ ] Implement global options (--json, --db, --config)
- [ ] Add error handling middleware
- [ ] Create command registration system
- [ ] Write CLI test helpers

## Phase 8: Management Commands

- [ ] Implement `add` command
- [ ] Implement `remove` command
- [ ] Implement `update` command
- [ ] Implement `index` command
- [ ] Implement `import` command (bibtex, json, csv)
- [ ] Implement `export` command (bibtex, json, csv)
- [ ] Write tests for management commands (target: 90%)

## Phase 9: Search Commands

- [ ] Implement `search` command (hybrid)
- [ ] Implement `vsearch` command (vector only)
- [ ] Implement `fulltext` command (FTS only)
- [ ] Implement `similar` command
- [ ] Implement `iterative` command with depth control
- [ ] Write tests for search commands (target: 90%)

## Phase 10: Tag Commands

- [ ] Implement `tag create` command
- [ ] Implement `tag add` command
- [ ] Implement `tag remove` command
- [ ] Implement `tag list` command
- [ ] Implement `tag delete` command
- [ ] Write tests for tag commands (target: 90%)

## Phase 11: Review Commands

- [ ] Implement `review` command with interactive mode
- [ ] Implement `read` command
- [ ] Implement `unread` command
- [ ] Write tests for review commands (target: 90%)

## Phase 12: Filter Commands

- [ ] Implement `filter create` command
- [ ] Implement `filter run` command
- [ ] Implement `filter list` command
- [ ] Implement `filter delete` command
- [ ] Write tests for filter commands (target: 90%)

## Phase 13: OpenAlex Commands

- [ ] Implement `openalex search` command
- [ ] Implement `openalex add` command
- [ ] Implement `openalex author` command
- [ ] Write tests for OpenAlex commands (target: 90%)

## Phase 14: Research Tools Commands

- [ ] Implement `refs` command
- [ ] Implement `citations` command
- [ ] Implement `related` command
- [ ] Implement `timeline` command
- [ ] Write tests for research commands (target: 90%)

## Phase 15: Data Processing Commands

- [ ] Implement `extract` command (DOI extraction from text)
- [ ] Implement `pdf` command (download/open)
- [ ] Implement `bibtex` command (single entry generation)
- [ ] Write tests for data processing commands (target: 90%)

## Phase 16: Utility Commands

- [ ] Implement `config` command (show, set, path)
- [ ] Implement `stats` command
- [ ] Implement `init` command
- [ ] Implement `backup` command
- [ ] Implement `doctor` command
- [ ] Write tests for utility commands (target: 90%)

## Phase 17: AI/RAG Commands

- [ ] Implement `ask` command with RAG pipeline
- [ ] Implement `summarize` command
- [ ] Implement `newsletter` command
- [ ] Write tests for AI commands (target: 85%)

## Phase 18: Coverage & Polish

- [ ] Review coverage report, identify gaps
- [ ] Add missing edge case tests
- [ ] Add error handling tests
- [ ] Boost coverage to 90%+ overall
- [ ] Add JSDoc comments to public API
- [ ] Write README.md with usage examples
- [ ] Add CHANGELOG.md
- [ ] Verify all CLI help text is complete
- [ ] Run final linting and fix issues

## Phase 19: Integration & E2E

- [ ] Write integration test: add -> search -> export workflow
- [ ] Write integration test: filter -> automatic updates
- [ ] Write integration test: tag organization workflow
- [ ] Write E2E test: full CLI workflow
- [ ] Test with real OpenAlex API (optional, network)

---

## Completed

<!-- Move tasks here when done -->

---

## Blockers

<!-- Document any issues preventing progress -->

---

## Notes

<!-- Add observations for future iterations -->

- Reference Python litdb for exact behavior: https://github.com/jkitchin/litdb
- Embedding model: all-MiniLM-L6-v2 (384 dimensions)
- Prioritize test coverage over features if time-constrained
