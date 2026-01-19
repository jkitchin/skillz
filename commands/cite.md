---
description: Format a citation from URL or DOI
allowed-tools: ["WebFetch", "WebSearch"]
argument-hint: <url-or-doi>
---

# Format Citation

Generate a properly formatted citation from a URL, DOI, or search query.

## Instructions

1. Parse the input: $ARGUMENTS
   - If DOI (10.xxxx/xxxxx), fetch metadata from doi.org
   - If URL, fetch the page and extract citation metadata
   - If text query, search for the paper

2. Extract citation information:
   - Authors (all of them)
   - Title
   - Publication venue (journal, conference, etc.)
   - Year
   - Volume, issue, pages (if applicable)
   - DOI or URL

3. Format in multiple citation styles

## Output Format

```markdown
## Citation Information

**Title**: <full title>
**Authors**: <author list>
**Published**: <venue>, <year>
**DOI**: <doi if available>
**URL**: <url>

### Citation Formats

#### APA 7th Edition
<formatted citation>

#### MLA 9th Edition
<formatted citation>

#### Chicago
<formatted citation>

#### BibTeX
```bibtex
@article{key,
  author = {Author Names},
  title = {Title},
  journal = {Journal},
  year = {Year},
  volume = {Vol},
  pages = {Pages},
  doi = {DOI}
}
```

#### Markdown Link
[Author et al. (Year). Title](<url>)
```

URL, DOI, or search query: $ARGUMENTS
