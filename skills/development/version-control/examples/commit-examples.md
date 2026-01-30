# Commit Message Examples

Real-world examples of good conventional commits.

## Features

### Simple Feature
```
feat(search): add full-text search
```

### Feature with Details
```
feat(auth): implement OAuth2 authentication

Add OAuth2 flow for Google and GitHub providers. Users can now
sign in with their existing accounts instead of creating new
passwords.

Features:
- OAuth2 flow for Google and GitHub
- User account linking
- Automatic profile sync
- Token refresh handling

Closes #123
```

### Feature with Breaking Change
```
feat(api)!: redesign user profile endpoint

Complete redesign of profile API to support new features and
improve performance.

Changes:
- Endpoint moved from /user/profile to /api/v2/users/:id
- Response format simplified and flattened
- Added support for partial updates via PATCH
- Removed deprecated fields: legacy_id, old_avatar

Migration steps:
1. Update API endpoint URLs in client
2. Update response parsing to handle new format
3. Replace PUT requests with PATCH for partial updates

BREAKING CHANGE: Profile endpoint moved and response format
changed. Client applications must be updated. See migration
guide at docs/migration-v2.md

Closes #456
```

## Bug Fixes

### Simple Bug Fix
```
fix(ui): correct button alignment on mobile
```

### Bug Fix with Context
```
fix(auth): prevent token expiration race condition

Fixed race condition where multiple simultaneous requests could
attempt to refresh expired tokens before new token was stored,
causing some requests to fail with 401 errors.

Solution: Use mutex lock during token refresh to ensure only
one refresh happens at a time.

Fixes #789
```

### Critical Security Fix
```
fix(auth): prevent SQL injection in login endpoint

Fixed SQL injection vulnerability in user login endpoint where
username parameter was not properly sanitized.

Impact: Attackers could potentially access or modify database
through crafted username inputs.

Solution: Use parameterized queries for all database operations.

Security advisory: All instances should be updated immediately.
No evidence of exploitation in logs.

Fixes #1234
CVE-2024-XXXXX
```

## Documentation

### Simple Documentation
```
docs(readme): update installation instructions
```

### Documentation with Details
```
docs(api): add comprehensive search endpoint examples

Added detailed examples for search API endpoint including:
- Basic text search
- Filtering by category
- Pagination handling
- Sorting results
- Error handling

Also added troubleshooting section for common issues.
```

## Refactoring

### Simple Refactor
```
refactor(database): extract query builder class
```

### Refactor with Explanation
```
refactor(auth): simplify token validation logic

Simplified token validation by extracting common patterns
into reusable functions. No behavior changes.

Changes:
- Extract validateTokenFormat()
- Extract checkTokenExpiration()
- Extract verifyTokenSignature()
- Reduce duplication across OAuth and JWT validators
- Improve test coverage from 75% to 95%

Benefits:
- More maintainable code
- Easier to add new token types
- Better test isolation
```

## Performance

### Performance Improvement
```
perf(database): add indexes on user queries

Added indexes on frequently queried columns (user_id,
created_at, email) reducing average query time from 200ms
to 50ms.

Impact: 75% reduction in database load during peak hours.
```

### Performance with Metrics
```
perf(images): implement lazy loading

Implemented lazy loading for images using Intersection Observer
API. Images are now loaded only when they enter viewport.

Results from production testing:
- Initial page load: 3.2s → 0.8s (75% faster)
- Lighthouse score: 65 → 92
- Bandwidth savings: ~2MB per page view
- Time to interactive: 4.1s → 1.2s

Implementation:
- Use native lazy loading for modern browsers
- Fallback to Intersection Observer for older browsers
- Placeholder images during load
- Progressive enhancement approach
```

## Tests

### Adding Tests
```
test(auth): add integration tests for login flow

Added comprehensive integration tests for entire login flow:
- Successful login with valid credentials
- Failed login with invalid credentials
- Account lockout after failed attempts
- Password reset flow
- OAuth login flow
- Session management

Coverage increased from 60% to 85% for auth module.
```

### Test Improvements
```
test(api): improve test reliability and speed

Refactored API tests to be more reliable and faster:
- Use test database instead of mocks
- Parallel test execution (10s → 3s)
- Better isolation between tests
- Clearer assertion messages
- Remove flaky timeout-based tests

All 247 tests now pass consistently in CI.
```

## Build and CI

### Dependency Updates
```
build(deps): upgrade React to v18.2

Upgraded React from v17.0 to v18.2 for concurrent features
and automatic batching improvements.

Changes:
- Update React and React-DOM
- Replace deprecated ReactDOM.render with createRoot
- Update @types/react to match
- Run migration codemod for breaking changes

All tests passing. Performance improvement of ~15% in
component render times.
```

### CI Configuration
```
ci(github): add automated deployment workflow

Added GitHub Actions workflow for automatic deployment to
staging and production environments.

Features:
- Deploy to staging on merge to main
- Deploy to production on tag creation
- Automatic rollback on health check failure
- Slack notifications for deployments
- Deploy preview for PRs

Workflow runs in ~5 minutes with full test suite.
```

## Chores

### Maintenance
```
chore(deps): update development dependencies

Updated development dependencies to latest versions:
- eslint 8.42.0 → 8.45.0
- prettier 2.8.8 → 3.0.0
- jest 29.5.0 → 29.6.1
- typescript 5.0.4 → 5.1.6

No breaking changes. All tests passing.
```

### Configuration
```
chore(config): update prettier and eslint config

Updated code formatting rules:
- Increase line width from 80 to 100
- Use single quotes for strings
- Add trailing commas
- Run prettier before eslint

Applied formatting to entire codebase (automated).
```

## Reverts

### Simple Revert
```
revert: revert "feat(auth): add OAuth2 authentication"

This reverts commit a1b2c3d4e5f6g7h8.

OAuth implementation causing login failures for existing users.
Reverting to investigate and fix before re-deploying.
```

### Revert with Explanation
```
revert: revert "perf(database): add indexes on user queries"

This reverts commit abc123def456.

While indexes improved read performance, they caused significant
write performance degradation (5x slower inserts). Need to
reconsider indexing strategy.

Will re-implement with partial indexes in separate PR.

Refs #789
```

## Multiple Scopes

### Cross-Cutting Changes
```
feat(api,database): add user search with full-text indexes

Implemented user search feature across API and database layers.

Database changes:
- Add tsvector column for full-text search
- Create GIN index on search column
- Add trigger to maintain search column

API changes:
- New /api/search/users endpoint
- Support for filters (role, status, created_at)
- Pagination and sorting
- Relevance-based ranking

Performance: <50ms for typical queries on 1M users.

Closes #456
```

## WIP Commits (Avoid in main)

### Work in Progress (For Feature Branches Only)
```
wip(feature): checkpoint - basic structure

Basic structure in place. Not ready for review.

TODO:
- Add error handling
- Add tests
- Update documentation
```

**Note:** WIP commits are acceptable in feature branches but should be
squashed before merging to main. Use descriptive commits instead.

---

## Anti-Patterns to Avoid

### ❌ Too Vague
```
Update files
Fix bug
Changes
WIP
```

### ❌ Wrong Type
```
feat(css): update button colors  # Should be fix or style
fix(api): add validation         # Should be feat
docs(auth): implement OAuth      # Should be feat
```

### ❌ Poor Subject
```
feat(api): adds a new endpoint for searching users with filters
# Too long, past tense

fix: fix the thing
# What thing?

feat: changes
# What changes?
```

### ❌ Missing Context
```
fix(api): fix timeout
# Why was there a timeout? What was the fix?

feat(auth): add oauth
# Which provider? What does it do?
```

### ✅ Better Versions
```
feat(api): add user search endpoint

fix(api): increase timeout from 5s to 30s for slow queries

feat(auth): add Google OAuth2 authentication
```

---

**Remember:** Good commit messages are documentation. Future you (and your team) will thank you for clear, descriptive commits!
