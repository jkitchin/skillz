# Version Control Quick Reference

Fast reference for common Git and GitHub operations.

## Conventional Commit Format

```
<type>(<scope>): <subject>

<body>

<footer>
```

### Types
- `feat`: New feature
- `fix`: Bug fix
- `docs`: Documentation
- `style`: Formatting
- `refactor`: Code restructuring
- `perf`: Performance
- `test`: Tests
- `build`: Build system
- `ci`: CI/CD
- `chore`: Maintenance

### Examples
```bash
feat(auth): add OAuth support
fix(api): handle null responses
docs(readme): update installation
refactor(db): extract query builder
```

## Daily Commands

### Check Status
```bash
git status              # Full status
git status -s           # Short status
git diff                # Unstaged changes
git diff --staged       # Staged changes
```

### Commit Changes
```bash
git add .               # Stage all
git add file.txt        # Stage specific file
git commit -m "type: message"
git commit --amend      # Modify last commit
```

### Branch Operations
```bash
git branch              # List branches
git checkout -b feat/name  # Create and switch
git checkout main       # Switch branch
git branch -d feat/name    # Delete branch
```

### Sync with Remote
```bash
git pull origin main       # Fetch and merge
git pull --rebase          # Fetch and rebase
git push origin branch     # Push branch
git push -u origin branch  # Push and set upstream
```

### Undo Changes
```bash
git reset HEAD~1           # Undo commit, keep changes
git reset --hard HEAD~1    # Undo commit, discard changes
git revert abc123          # Revert commit (safe)
git checkout -- file.txt   # Discard file changes
```

## Trunk-Based Development Workflow

```bash
# 1. Start from main
git checkout main
git pull origin main

# 2. Create feature branch
git checkout -b feat/feature-name

# 3. Work and commit
git add .
git commit -m "feat(module): add feature"

# 4. Keep updated
git rebase main

# 5. Push and PR
git push -u origin feat/feature-name

# 6. After merge
git checkout main
git pull origin main
git branch -d feat/feature-name
```

## Branch Naming

```
feat/feature-name      # New feature
fix/bug-description    # Bug fix
refactor/what          # Refactoring
docs/what              # Documentation
test/what              # Tests
chore/what             # Maintenance
```

## Merge Conflicts

```bash
# During merge/rebase
git status              # See conflicted files

# In file, look for:
<<<<<<< HEAD
# Your changes
=======
# Their changes
>>>>>>> branch-name

# After resolving
git add resolved-file.txt
git commit              # For merge
git rebase --continue   # For rebase

# Or abort
git merge --abort
git rebase --abort
```

## Pull Request

### Creating PR

1. Push branch: `git push -u origin feat/branch`
2. Create PR on GitHub
3. Use Conventional Commit format for title
4. Fill description with changes and testing
5. Request reviewers
6. Link related issues

### PR Description Template

```markdown
## Summary
Brief description

## Changes
- Change 1
- Change 2

## Testing
How to test

## Checklist
- [ ] Tests added
- [ ] Docs updated
- [ ] No breaking changes
```

## GitHub CLI (gh)

```bash
# Create PR
gh pr create --fill

# List PRs
gh pr list

# View PR
gh pr view 123

# Checkout PR
gh pr checkout 123

# Review PR
gh pr review 123 --approve

# Merge PR
gh pr merge 123 --squash
```

## Stash Operations

```bash
git stash                    # Stash changes
git stash push -m "message"  # Stash with message
git stash list              # List stashes
git stash apply             # Apply stash
git stash pop               # Apply and remove
git stash drop              # Remove stash
```

## View History

```bash
git log                      # Full log
git log --oneline           # Compact log
git log --graph             # Graph view
git log -5                  # Last 5 commits
git log --grep="bug"        # Search commits
git log --author="Name"     # By author
git log file.txt            # File history
```

## Advanced Operations

### Interactive Rebase
```bash
git rebase -i HEAD~3        # Last 3 commits

# Commands: pick, reword, edit, squash, fixup, drop
```

### Cherry-Pick
```bash
git cherry-pick abc123      # Apply specific commit
```

### Reflog (Recovery)
```bash
git reflog                  # View reflog
git reset --hard HEAD@{2}   # Restore to state
```

### Bisect (Find Bug)
```bash
git bisect start
git bisect bad              # Current is bad
git bisect good abc123      # Last known good
# Test, then mark:
git bisect good/bad
# Repeat until found
git bisect reset
```

## Troubleshooting

### Wrong Branch
```bash
git branch feat/correct-branch  # Create branch here
git reset --hard origin/main    # Reset current
git checkout feat/correct-branch
```

### Wrong Commit Message
```bash
git commit --amend -m "correct message"
```

### Undo Last Commit
```bash
git reset --soft HEAD~1     # Keep changes staged
git reset HEAD~1            # Keep changes unstaged
git reset --hard HEAD~1     # Discard everything
```

### Branch Diverged
```bash
git pull --rebase origin branch
git push --force-with-lease origin branch
```

### Lost Commits
```bash
git reflog                  # Find lost commit
git reset --hard abc123     # Restore
```

## Configuration

```bash
# User info
git config --global user.name "Name"
git config --global user.email "email@example.com"

# Defaults
git config --global init.defaultBranch main
git config --global pull.rebase true

# Aliases
git config --global alias.st status
git config --global alias.co checkout
git config --global alias.br branch
git config --global alias.lg "log --oneline --graph"
```

## .gitignore

```bash
# Node.js
node_modules/
npm-debug.log

# Build
dist/
build/
*.min.js

# Environment
.env
.env.local

# IDE
.vscode/
.idea/
*.swp

# OS
.DS_Store
Thumbs.db
```

## Branch Protection

**Essential settings for main branch:**
- ✅ Require pull request before merging
- ✅ Require approvals (1+)
- ✅ Require status checks to pass
- ✅ Require branches up to date
- ✅ Require conversation resolution
- ✅ Require linear history
- ✅ Automatically delete head branches

## Keyboard Shortcuts

**GitHub:**
- `t` - File finder
- `b` - Blame view
- `l` - Jump to line
- `.` - Open in web editor
- `?` - Show keyboard shortcuts

## Git Aliases

Add to `~/.gitconfig`:

```ini
[alias]
    st = status
    co = checkout
    br = branch
    cm = commit
    lg = log --oneline --graph --decorate --all
    unstage = reset HEAD --
    last = log -1 HEAD
    undo = reset --soft HEAD^
    amend = commit --amend --no-edit
```

## Quick Tips

1. **Commit often** - Multiple times per day
2. **Pull before starting** - Always work on latest
3. **Keep branches short** - Days, not weeks
4. **Small PRs** - 50-200 lines ideal
5. **Review promptly** - Within 2-4 hours
6. **Test before committing** - Catch issues early
7. **Use descriptive names** - Branches and commits
8. **Clean up** - Delete merged branches
9. **Protect main** - Branch protection rules
10. **Document why** - Not just what

## Common Patterns

### Feature Development
```bash
git checkout main && git pull
git checkout -b feat/feature
# Work, commit, test
git rebase main
git push -u origin feat/feature
# Create PR, merge, delete branch
```

### Bug Fix
```bash
git checkout main && git pull
git checkout -b fix/bug-name
# Fix, test, commit
git push -u origin fix/bug-name
# Create PR with "Fixes #123"
```

### Sync Branch
```bash
git fetch origin
git rebase origin/main
git push --force-with-lease
```

### Cleanup
```bash
git checkout main
git pull
git branch --merged | grep -v "main" | xargs git branch -d
```

---

**Quick Help:**
- Full skill: See SKILL.md
- Git commands: See references/git-commands.md
- Troubleshooting: See references/troubleshooting.md
- GitHub workflow: See references/github-workflow.md
- Conventional Commits: See references/conventional-commits.md
