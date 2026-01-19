# CLI Commands Specification

## Overview

The litdb CLI provides 43 commands organized into categories. Each command should:
- Have `--help` documentation
- Support `--json` output flag where applicable
- Return proper exit codes (0 success, 1 error)
- Be individually testable

## Command Categories

### Management Commands

#### `litdb add <doi>`
Add a work by DOI.

```bash
litdb add 10.1038/nature12373
litdb add https://doi.org/10.1038/nature12373  # URL format ok
litdb add --bibtex paper.bib                    # From bibtex file
```

Options:
- `--tag <tag>` - Add tag to imported work
- `--fetch-refs` - Also fetch referenced works
- `--fetch-citations` - Also fetch citing works

#### `litdb remove <id|doi>`
Remove a work from the database.

```bash
litdb remove 10.1038/nature12373
litdb remove abc123-def456  # By internal ID
litdb remove --all --tag deprecated  # Remove all with tag
```

Options:
- `--force` - Skip confirmation
- `--keep-authors` - Don't remove orphaned authors

#### `litdb update <id|doi>`
Refresh metadata from source.

```bash
litdb update 10.1038/nature12373
litdb update --all  # Update all works
```

#### `litdb index`
Rebuild search indexes.

```bash
litdb index              # Rebuild all indexes
litdb index --fts        # Only full-text
litdb index --vector     # Only vector embeddings
litdb index --work <id>  # Single work
```

#### `litdb import <file>`
Bulk import from file.

```bash
litdb import papers.bib           # BibTeX
litdb import papers.json          # JSON array of DOIs
litdb import papers.csv           # CSV with DOI column
litdb import --openalex query.txt # OpenAlex query results
```

#### `litdb export`
Export works.

```bash
litdb export --format bibtex > papers.bib
litdb export --format json > papers.json
litdb export --format csv > papers.csv
litdb export --tag ml --format bibtex  # Only tagged works
```

### Search Commands

#### `litdb search <query>`
Hybrid search (combines vector + full-text).

```bash
litdb search "machine learning for materials science"
litdb search "neural networks" --limit 20
litdb search "catalyst design" --tag chemistry
```

Options:
- `--limit <n>` - Number of results (default: 10)
- `--offset <n>` - Skip first n results
- `--tag <tag>` - Filter by tag
- `--since <date>` - Only works after date
- `--json` - JSON output

#### `litdb vsearch <query>`
Vector (semantic) search only.

```bash
litdb vsearch "novel approaches to energy storage"
```

#### `litdb fulltext <query>`
Full-text search only (FTS5).

```bash
litdb fulltext "graphene oxide"
litdb fulltext "title:neural network"  # Field-specific
```

#### `litdb similar <id|doi>`
Find similar works.

```bash
litdb similar 10.1038/nature12373
litdb similar abc123 --limit 5
```

#### `litdb iterative <query>`
Iterative search - expands through references/citations.

```bash
litdb iterative "CRISPR gene editing" --depth 2
litdb iterative "battery materials" --max-papers 100
```

### Tag Commands

#### `litdb tag create <name>`
Create a new tag.

```bash
litdb tag create machine-learning
litdb tag create "to read" --color "#ff6600"
```

#### `litdb tag add <work> <tag>`
Add tag to work.

```bash
litdb tag add 10.1038/nature12373 important
litdb tag add --all --search "neural" deep-learning
```

#### `litdb tag remove <work> <tag>`
Remove tag from work.

```bash
litdb tag remove 10.1038/nature12373 important
```

#### `litdb tag list`
List all tags.

```bash
litdb tag list
litdb tag list --counts  # Show work counts
```

#### `litdb tag delete <name>`
Delete a tag.

```bash
litdb tag delete obsolete-tag
```

### Review Commands

#### `litdb review`
Interactive review of unread works.

```bash
litdb review              # Start review session
litdb review --tag inbox  # Review specific tag
```

Actions during review:
- `k` - Keep (mark as reviewed)
- `d` - Delete
- `t` - Add tag
- `n` - Add note
- `s` - Skip
- `q` - Quit

#### `litdb read <id|doi>`
Mark work as read.

```bash
litdb read 10.1038/nature12373
```

#### `litdb unread <id|doi>`
Mark work as unread.

```bash
litdb unread 10.1038/nature12373
```

### Filter Commands

#### `litdb filter create <name> <query>`
Create an automatic filter.

```bash
litdb filter create "ml-papers" "machine learning" --source openalex
litdb filter create "nature-2024" "journal:Nature,year:2024"
```

#### `litdb filter run [name]`
Run filters to fetch new papers.

```bash
litdb filter run          # Run all active filters
litdb filter run ml-papers  # Run specific filter
```

#### `litdb filter list`
List all filters.

```bash
litdb filter list
litdb filter list --active  # Only active
```

#### `litdb filter delete <name>`
Delete a filter.

```bash
litdb filter delete old-filter
```

### OpenAlex Commands

#### `litdb openalex search <query>`
Search OpenAlex directly.

```bash
litdb openalex search "climate change mitigation"
litdb openalex search --author "John Smith"
litdb openalex search --institution "MIT"
```

#### `litdb openalex add <openalex-id>`
Add work by OpenAlex ID.

```bash
litdb openalex add W2741809807
```

#### `litdb openalex author <orcid|name>`
Get author info from OpenAlex.

```bash
litdb openalex author 0000-0002-1234-5678
litdb openalex author "Jane Doe"
```

### Research Tools

#### `litdb refs <id|doi>`
Show references of a work.

```bash
litdb refs 10.1038/nature12373
litdb refs abc123 --fetch  # Add missing refs to DB
```

#### `litdb citations <id|doi>`
Show works citing this paper.

```bash
litdb citations 10.1038/nature12373
litdb citations abc123 --fetch
```

#### `litdb related <id|doi>`
Find related works (refs + citations + similar).

```bash
litdb related 10.1038/nature12373
```

#### `litdb timeline <query>`
Show publication timeline for topic.

```bash
litdb timeline "graphene" --since 2010
```

### Data Processing

#### `litdb extract <file>`
Extract references from text.

```bash
litdb extract paper.txt      # Extract DOIs from text
litdb extract paper.pdf      # Extract from PDF
cat clipboard.txt | litdb extract -  # From stdin
```

#### `litdb pdf <id|doi>`
Download/open PDF.

```bash
litdb pdf 10.1038/nature12373 --download
litdb pdf abc123 --open  # Open in default viewer
```

#### `litdb bibtex <id|doi>`
Generate BibTeX entry.

```bash
litdb bibtex 10.1038/nature12373
litdb bibtex --all > library.bib
```

### Utility Commands

#### `litdb config`
Manage configuration.

```bash
litdb config show
litdb config set openalex.email "me@example.com"
litdb config path  # Show config file path
```

#### `litdb stats`
Show database statistics.

```bash
litdb stats
litdb stats --json
```

Output:
```
Works: 1,234
Authors: 5,678
Tags: 15
Embeddings: 12,340 chunks
Database size: 45.2 MB
```

#### `litdb init`
Initialize new database.

```bash
litdb init
litdb init --path ~/research/papers.db
```

#### `litdb backup`
Backup database.

```bash
litdb backup
litdb backup --path ~/backups/litdb-2024-01-15.db
```

#### `litdb doctor`
Check database health.

```bash
litdb doctor
litdb doctor --fix  # Auto-fix issues
```

### AI/RAG Commands

#### `litdb ask <question>`
Ask question using RAG.

```bash
litdb ask "What are the main approaches to CO2 capture?"
litdb ask "Summarize recent advances in battery technology" --limit 20
```

#### `litdb summarize [id|doi|tag]`
Summarize works.

```bash
litdb summarize abc123                    # Single work
litdb summarize --tag "to read"           # All with tag
litdb summarize --recent 7                # Last 7 days
```

#### `litdb newsletter`
Generate newsletter-style summary.

```bash
litdb newsletter --since "last week"
litdb newsletter --tag inbox --format markdown
```

## Testing Strategy

Each command needs:

1. **Unit tests** for command logic
2. **Integration tests** with real DB
3. **Snapshot tests** for output format
4. **Error case tests** (invalid args, missing data)

Example test structure:

```typescript
describe('litdb add', () => {
  it('adds work by DOI', async () => {
    const result = await cli(['add', '10.1038/nature12373']);
    expect(result.exitCode).toBe(0);
    expect(result.stdout).toContain('Added:');
  });

  it('rejects invalid DOI', async () => {
    const result = await cli(['add', 'not-a-doi']);
    expect(result.exitCode).toBe(1);
    expect(result.stderr).toContain('Invalid DOI');
  });

  it('handles duplicate gracefully', async () => {
    await cli(['add', '10.1038/nature12373']);
    const result = await cli(['add', '10.1038/nature12373']);
    expect(result.exitCode).toBe(0);
    expect(result.stdout).toContain('Already exists');
  });
});
```

## Output Formats

All commands displaying data should support:

```bash
# Default: Human-readable
litdb search "query"

# JSON: Machine-readable
litdb search "query" --json

# Quiet: Minimal output
litdb add doi --quiet
```
