# Commit Message Template

## Format

```
<type>(<scope>): <subject>

<body>

<footer>
```

## Type
Choose one:
- feat: New feature
- fix: Bug fix
- docs: Documentation only
- style: Formatting (no code change)
- refactor: Code restructuring
- perf: Performance improvement
- test: Adding/updating tests
- build: Build system/dependencies
- ci: CI/CD configuration
- chore: Maintenance tasks
- revert: Revert previous commit

## Scope (optional)
Module, component, or area affected:
- auth, api, ui, database, etc.

## Subject (required)
- Use imperative mood: "add" not "added"
- No capitalization
- No period at end
- Max 50 characters
- Complete: "If applied, this commit will <subject>"

## Body (optional)
- Explain what and why, not how
- Wrap at 72 characters
- Separate from subject with blank line
- Use bullet points for multiple items

## Footer (optional)
- Breaking changes: `BREAKING CHANGE: description`
- Issue references: `Fixes #123, Closes #456`
- Co-authors: `Co-authored-by: Name <email>`

---

## Examples

### Simple Feature
```
feat(dashboard): add user activity chart
```

### Bug Fix with Details
```
fix(api): handle timeout in user service

Previously, network timeouts would crash the application.
Now we catch timeout errors and return 503 status with
retry-after header.

Fixes #789
```

### Breaking Change
```
feat(api)!: redesign authentication endpoint

BREAKING CHANGE: /auth/login moved to /api/v2/auth/login
and now requires client_id parameter. See migration guide
in docs/migration.md
```

### Multiple Changes
```
refactor(database): improve connection handling

- Implement connection pooling (max 20 connections)
- Add automatic retry for transient failures
- Improve error logging with connection state
- Add metrics for pool usage

These changes improve reliability under high load.

Refs #456
```

---

## Quick Reference

**Good Examples:**
- `feat(auth): add OAuth2 support`
- `fix(ui): correct button alignment on mobile`
- `docs(api): add examples for search endpoint`
- `refactor(database): extract query builder`
- `perf(images): implement lazy loading`

**Bad Examples:**
- ❌ `Update stuff`
- ❌ `Fixed a bug.`
- ❌ `WIP`
- ❌ `feat(api): Added a new endpoint`
- ❌ `Fix the thing that was broken`
