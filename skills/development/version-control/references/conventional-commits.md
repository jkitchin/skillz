# Conventional Commits Complete Guide

## Specification

Conventional Commits is a specification for adding human and machine-readable meaning to commit messages.

**Format:**
```
<type>(<optional scope>): <subject>

<optional body>

<optional footer(s)>
```

## Commit Types

### feat (Feature)
New feature for the user, not a new feature for build script.

```
feat(auth): add OAuth2 provider support
feat(api): implement user search endpoint
feat: add dark mode toggle
```

**When to use:**
- Adding new user-facing functionality
- Implementing new API endpoints
- Creating new components or modules
- Adding new configuration options

**Bumps:** MINOR version (0.X.0)

### fix (Bug Fix)
Bug fix for the user, not a fix to a build script.

```
fix(auth): prevent token expiration race condition
fix(ui): correct modal z-index on mobile
fix: handle null values in user profile
```

**When to use:**
- Fixing incorrect behavior
- Resolving errors or exceptions
- Correcting visual bugs
- Patching security vulnerabilities

**Bumps:** PATCH version (0.0.X)

### docs (Documentation)
Documentation only changes.

```
docs(readme): update installation instructions
docs(api): add examples for search endpoint
docs: fix typos in contributing guide
```

**When to use:**
- README updates
- API documentation changes
- Comment improvements
- Adding code examples
- Tutorial or guide updates

**Bumps:** None (no release)

### style (Formatting)
Changes that do not affect the meaning of the code (white-space, formatting, missing semi-colons, etc).

```
style(components): fix indentation in Button
style: apply prettier formatting
style(css): organize imports alphabetically
```

**When to use:**
- Formatting changes only
- Linting fixes
- Whitespace adjustments
- Code style consistency

**Does NOT include:** CSS/styling changes (use feat/fix)

**Bumps:** None (no release)

### refactor (Code Restructuring)
Code change that neither fixes a bug nor adds a feature.

```
refactor(database): extract query builder class
refactor(auth): simplify token validation logic
refactor: convert callbacks to async/await
```

**When to use:**
- Restructuring code without changing behavior
- Improving code organization
- Extracting functions or classes
- Simplifying complex logic
- Removing code duplication

**Bumps:** None (or PATCH in some conventions)

### perf (Performance)
Code change that improves performance.

```
perf(images): implement lazy loading
perf(database): add indexes on user queries
perf: reduce bundle size by 30%
```

**When to use:**
- Optimizing algorithms
- Improving load times
- Reducing memory usage
- Adding caching
- Database query optimization

**Bumps:** PATCH version (0.0.X)

### test (Tests)
Adding missing tests or correcting existing tests.

```
test(auth): add integration tests for login flow
test(api): cover edge cases in search endpoint
test: increase coverage to 85%
```

**When to use:**
- Adding new tests
- Updating existing tests
- Fixing test failures
- Improving test coverage

**Bumps:** None (no release)

### build (Build System)
Changes that affect the build system or external dependencies.

```
build(deps): upgrade React to v18.2
build(webpack): optimize production configuration
build: add TypeScript compilation step
```

**When to use:**
- Updating dependencies
- Modifying build configuration
- Changing compilation settings
- Package.json updates

**Bumps:** PATCH version (0.0.X) or MINOR if new features

### ci (Continuous Integration)
Changes to CI configuration files and scripts.

```
ci(github): add automated deployment workflow
ci(travis): enable caching for node_modules
ci: run tests on all branches
```

**When to use:**
- GitHub Actions changes
- Travis, CircleCI, Jenkins configuration
- CI script updates
- Deployment pipeline changes

**Bumps:** None (no release)

### chore (Maintenance)
Other changes that don't modify src or test files.

```
chore(deps): update development dependencies
chore: update .gitignore
chore(release): bump version to 1.2.0
```

**When to use:**
- Routine maintenance
- Tooling updates
- Configuration changes
- Release preparation
- Cleanup tasks

**Bumps:** None (no release)

### revert
Reverts a previous commit.

```
revert: revert "feat(auth): add OAuth2 authentication"

This reverts commit a1b2c3d4.
```

**When to use:**
- Undoing a previous commit
- Rolling back problematic changes

**Bumps:** Depends on what's being reverted

## Scopes

Scopes are optional but highly recommended. They indicate what part of the codebase changed.

### Project-Specific Scopes

Define scopes that make sense for your project:

**By module/feature:**
```
feat(auth): ...
fix(billing): ...
refactor(analytics): ...
```

**By layer:**
```
feat(api): ...
fix(database): ...
refactor(ui): ...
```

**By component:**
```
fix(button): ...
feat(modal): ...
style(navbar): ...
```

### Scope Examples by Project Type

**Web application:**
- `(auth)`, `(ui)`, `(api)`, `(database)`, `(middleware)`

**Library/Package:**
- `(core)`, `(utils)`, `(types)`, `(exports)`

**Mobile app:**
- `(ios)`, `(android)`, `(navigation)`, `(state)`

**Monorepo:**
- `(web)`, `(mobile)`, `(shared)`, `(backend)`

### Multiple Scopes

```
feat(api,database): add user search with indexes
```

Use sparingly - usually indicates commit should be split.

### No Scope

Valid for changes that span multiple areas or are project-wide:

```
feat: add support for configuration files
fix: resolve memory leak
chore: update all dependencies
```

## Subject (Description)

### Rules

**1. Use imperative mood:**
```
✅ add feature
✅ fix bug
✅ update documentation

❌ added feature
❌ fixes bug
❌ updating documentation
```

Complete the sentence: "If applied, this commit will..."

**2. No capitalization:**
```
✅ feat(api): add user endpoint
❌ feat(api): Add user endpoint
```

**3. No period at end:**
```
✅ fix(auth): prevent token expiration
❌ fix(auth): prevent token expiration.
```

**4. Keep it short (50 characters or less):**
```
✅ feat(api): add user search endpoint
❌ feat(api): add a new comprehensive user search endpoint with filtering
```

**5. Be specific but concise:**
```
✅ fix(auth): prevent race condition in token refresh
❌ fix(auth): fix bug
❌ fix(auth): prevent race condition that occurs when multiple requests simultaneously attempt to refresh expired access tokens before the new token is stored
```

### Good Examples

```
feat(search): implement full-text search
fix(ui): correct button alignment on mobile
docs(api): add authentication examples
refactor(database): extract connection pooling
perf(images): add lazy loading
test(auth): cover token expiration scenarios
build(deps): upgrade to Node 18
ci(github): add deploy preview workflow
```

### Bad Examples

```
❌ Update stuff
❌ Fixed a bug.
❌ WIP
❌ feat(api): Added a new endpoint for searching users with pagination
❌ fix: Fix the thing that was broken yesterday
❌ Changes
```

## Body

Optional but recommended for complex changes.

### When to Include Body

- Explaining **why** the change was made
- Providing context that isn't obvious
- Describing complex logic or decisions
- Listing multiple related changes
- Explaining trade-offs or alternatives considered

### Body Format

- Wrap at 72 characters per line
- Separate from subject with blank line
- Use bullet points for lists
- Explain what and why, not how
- Can have multiple paragraphs

### Examples

**Simple with context:**
```
feat(auth): add rate limiting to login endpoint

Prevent brute force attacks by limiting login attempts to 5 per
minute per IP address. After limit is reached, user must wait
before trying again.
```

**With multiple points:**
```
refactor(database): optimize query performance

- Add indexes on frequently queried columns (user_id, created_at)
- Implement connection pooling with max 20 connections
- Cache common queries using Redis with 5-minute TTL
- Update ORM configuration to use prepared statements

These changes reduce average query time from 200ms to 50ms in
production testing.
```

**With reasoning:**
```
fix(api): use POST instead of GET for search endpoint

GET requests with complex query parameters were hitting URL
length limits. POST allows unlimited query complexity and is
semantically more appropriate for searches with side effects
(logging, analytics).

Considered GraphQL but decided against it to minimize changes
to existing client code.
```

**Breaking change explanation:**
```
feat(api)!: redesign authentication flow

Complete redesign of authentication to use JWT tokens instead
of session cookies. This provides better support for mobile
clients and microservices architecture.

Migration guide:
1. Update client to store JWT tokens
2. Implement token refresh logic
3. Remove session cookie handling

BREAKING CHANGE: Session-based authentication is no longer
supported. All clients must implement JWT token handling.
```

## Footer

Optional metadata at the end of commit message.

### Breaking Changes

**Mark with `BREAKING CHANGE:` or `!` in type:**

```
feat(api)!: change authentication endpoint

BREAKING CHANGE: /auth/login moved to /api/v2/auth/login and
now requires client_id parameter.
```

```
feat!: drop support for Node 12

BREAKING CHANGE: Minimum Node version is now 14.x due to
dependencies on ES2020 features.
```

**Bumps:** MAJOR version (X.0.0)

### Issue References

**Link to issue tracker:**

```
fix(auth): prevent duplicate login attempts

Fixes #123
Closes #456
Resolves #789
```

**Multiple issues:**
```
refactor(database): improve query performance

Closes #123, #124, #125
```

**Partial work on issue:**
```
feat(search): implement basic search

Refs #456 (pagination will be added in future PR)
```

### Other Footer Tokens

```
Reviewed-by: Jane Doe <jane@example.com>
Co-authored-by: John Smith <john@example.com>
Signed-off-by: Author Name <author@example.com>
```

## Complete Examples

### Simple Feature
```
feat(dashboard): add user activity chart
```

### Feature with Body
```
feat(search): implement full-text search

Add PostgreSQL full-text search using tsvector indexes for
better performance than LIKE queries. Includes pagination and
relevance ranking.
```

### Bug Fix with Footer
```
fix(api): handle timeout in user service

Previously, network timeouts would crash the application. Now
we catch timeout errors gracefully and return 503 status to
client with retry-after header.

Fixes #789
```

### Breaking Change
```
feat(api)!: redesign user profile endpoint

Complete redesign of profile API to support new features and
improve performance.

Changes:
- Endpoint moved from /user/profile to /api/v2/users/:id
- Response format simplified
- Added support for partial updates via PATCH

Migration:
- Update API endpoint URLs
- Update response parsing to handle new format
- Implement PATCH for partial updates

BREAKING CHANGE: Profile endpoint moved and response format
changed. See migration guide in docs/migration-v2.md
```

### Revert
```
revert: revert "feat(auth): add OAuth2 authentication"

This reverts commit a1b2c3d4e5f6g7h8.

OAuth implementation was causing login failures for existing
users. Reverting to investigate and fix before re-deploying.
```

### Multiple Changes
```
refactor(database): improve connection handling

- Implement connection pooling with max 20 connections
- Add automatic retry logic for transient failures
- Improve error logging with connection state
- Add metrics for connection pool usage

These changes improve reliability under high load and provide
better debugging information when connection issues occur.

Refs #456
```

## Tools and Automation

### Commitizen
Interactive CLI for creating conventional commits:

```bash
npm install -g commitizen cz-conventional-changelog

# Initialize in project
commitizen init cz-conventional-changelog --save-dev --save-exact

# Use with git cz instead of git commit
git cz
```

### Commitlint
Lint commit messages:

```bash
npm install --save-dev @commitlint/{cli,config-conventional}

# .commitlintrc.json
{
  "extends": ["@commitlint/config-conventional"]
}

# Test commit message
echo "feat(api): add endpoint" | commitlint
```

### Husky
Git hooks to enforce conventional commits:

```bash
npm install --save-dev husky

# Setup hook
npx husky add .husky/commit-msg 'npx --no -- commitlint --edit $1'
```

### Standard Version
Automated versioning and changelog:

```bash
npm install --save-dev standard-version

# Bump version and generate changelog
npm run release
```

### Semantic Release
Fully automated versioning and releases:

```bash
npm install --save-dev semantic-release

# Configure in .releaserc.json
# Automatically bumps version, generates changelog, creates release
```

## Best Practices

### Commit Frequency

**Do commit:**
- After each logical unit of work
- Multiple times per day
- When tests pass
- Before switching tasks

**Don't commit:**
- Broken code
- Unrelated changes together
- Work in progress without clear message

### Commit Size

**Good commit size:**
- Single feature or fix
- Can be reviewed in 5-10 minutes
- Isolated change that can be reverted safely
- 50-200 lines changed typically

**Too large:**
- Multiple features
- Mixes refactoring and features
- Spans multiple modules/concerns
- 500+ lines changed

**Too small:**
- Every file saved
- Formatting-only changes mixed with logic
- Incomplete thought

### Message Quality

**High quality:**
```
feat(auth): add OAuth2 provider support

Implement OAuth2 authentication flow for Google and GitHub
providers. Includes token refresh logic and user profile sync.

Users can now sign in with their Google or GitHub accounts
instead of creating a new password.

Closes #123
```

**Low quality:**
```
add stuff
```

### When to Use What

**feat vs fix:**
- `feat`: It wasn't there before, now it is
- `fix`: It was supposed to work, but didn't

**refactor vs perf:**
- `refactor`: Code organization, no observable change
- `perf`: Measurable performance improvement

**chore vs build vs ci:**
- `chore`: Misc maintenance
- `build`: Build system or dependencies
- `ci`: CI/CD configuration

**docs vs comments:**
- `docs`: Documentation files (README, etc.)
- Code comments: Use `refactor` or the type of code change

## Common Pitfalls

### Vague Messages
```
❌ fix: fix bug
❌ feat: add feature
❌ update files
❌ changes

✅ fix(auth): prevent token expiration race condition
✅ feat(api): add user search endpoint
```

### Wrong Type
```
❌ feat(css): update button colors  # Should be fix or style
❌ fix(api): add validation  # Should be feat
❌ docs(auth): implement OAuth  # Should be feat

✅ style(button): update colors
✅ feat(api): add request validation
✅ feat(auth): implement OAuth
```

### Missing Context
```
❌ fix(api): fix timeout

✅ fix(api): increase timeout from 5s to 30s

Long-running queries were timing out. Increased timeout to 30s
based on 95th percentile query time in production logs.
```

### Multiple Changes
```
❌ feat(api): add search, fix bugs, update docs

✅ Split into three commits:
   1. feat(api): add user search endpoint
   2. fix(api): handle null values in response
   3. docs(api): add search endpoint examples
```

### Inconsistent Style
```
❌ Different styles in same project:
   - feat(api): add endpoint
   - Feature: Add another endpoint
   - ADD ENDPOINT
   - added endpoint

✅ Consistent throughout:
   - feat(api): add user endpoint
   - feat(api): add search endpoint
   - feat(api): add filter endpoint
```

## Version Bumping

Conventional Commits enable automated semantic versioning:

**MAJOR (X.0.0):** Breaking changes
```
feat!: change API authentication
BREAKING CHANGE: description
```

**MINOR (0.X.0):** New features
```
feat(api): add new endpoint
```

**PATCH (0.0.X):** Bug fixes
```
fix(api): handle edge case
```

**No bump:** docs, style, test, chore, ci

## Project-Specific Adaptations

### Custom Types

Some projects add additional types:
- `security`: Security fixes (often treated as fix)
- `deps`: Dependency updates (alternative to build)
- `wip`: Work in progress (avoid in main branch)
- `hotfix`: Emergency production fix

### Stricter Rules

Some projects enforce:
- Scope always required
- Body required for feat and fix
- Issue reference required
- Specific scope values only
- Maximum subject length

### Configuration Example

```json
{
  "types": [
    { "type": "feat", "section": "Features" },
    { "type": "fix", "section": "Bug Fixes" },
    { "type": "perf", "section": "Performance" },
    { "type": "revert", "section": "Reverts" },
    { "type": "docs", "section": "Documentation", "hidden": false },
    { "type": "style", "section": "Styles", "hidden": true },
    { "type": "chore", "section": "Miscellaneous", "hidden": true },
    { "type": "refactor", "section": "Code Refactoring", "hidden": false },
    { "type": "test", "section": "Tests", "hidden": true },
    { "type": "build", "section": "Build System", "hidden": true },
    { "type": "ci", "section": "CI/CD", "hidden": true }
  ],
  "scopes": ["api", "ui", "database", "auth", "billing"],
  "scopeOverrides": {
    "feat": ["api", "ui"],
    "fix": ["api", "ui", "database"]
  }
}
```

## Resources

- **Official Spec:** https://www.conventionalcommits.org/
- **Commitizen:** https://github.com/commitizen/cz-cli
- **Commitlint:** https://commitlint.js.org/
- **Standard Version:** https://github.com/conventional-changelog/standard-version
- **Semantic Release:** https://semantic-release.gitbook.io/

---

**Remember:** Conventional Commits are about communication. They make your history searchable, your releases automated, and your changes understandable. Consistency is key.
