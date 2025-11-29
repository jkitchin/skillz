# Version Control Skill

Expert guidance for Git version control, trunk-based development workflows, and GitHub best practices.

## Overview

This skill provides comprehensive support for Git version control, with an emphasis on:

- **Trunk-Based Development:** Short-lived feature branches and frequent integration
- **Conventional Commits:** Structured, semantic commit messages
- **GitHub Integration:** Pull requests, code review, and GitHub Actions
- **Best Practices:** Professional version control workflows

## When to Use

Invoke this skill when you need help with:

- Git commands and operations
- Creating and managing branches
- Writing commit messages
- Creating and reviewing pull requests
- Resolving merge conflicts
- GitHub workflow and configuration
- Git troubleshooting and recovery
- Version control best practices
- Team collaboration with Git

## Key Features

### 1. Trunk-Based Development

Learn and apply trunk-based development workflow:
- Keep main branch always deployable
- Short-lived feature branches (hours to days)
- Frequent integration to prevent merge hell
- Small, incremental changes
- Feature flags for incomplete work

### 2. Conventional Commits

Master structured commit messages:
- Type-based commits (feat, fix, docs, etc.)
- Optional scope for context
- Clear, imperative subject lines
- Detailed body for complex changes
- Footer for breaking changes and issue references

### 3. GitHub Workflow

Professional GitHub collaboration:
- Pull request best practices
- Code review guidelines
- Branch protection rules
- GitHub Actions for CI/CD
- Issue and project management

### 4. Git Expertise

Comprehensive Git knowledge:
- Essential daily commands
- Advanced techniques (rebase, cherry-pick, bisect)
- Conflict resolution strategies
- Repository maintenance
- Recovery from mistakes

## Skill Contents

### Main Skill (SKILL.md)
Complete skill prompt with:
- Core workflow guidance
- Conventional Commits format
- Essential Git commands
- GitHub PR workflow
- Common troubleshooting
- Quick reference

### Reference Materials

**references/conventional-commits.md**
- Complete specification
- All commit types explained
- Scope and subject guidelines
- Body and footer formatting
- Examples and anti-patterns
- Tools and automation

**references/trunk-based-development.md**
- TBD philosophy and benefits
- Feature branch workflow
- Small changes strategy
- Feature flags
- CI/CD integration
- Team size considerations

**references/github-workflow.md**
- Pull request workflow
- PR templates and labels
- Code review best practices
- GitHub Actions examples
- Branch protection rules
- GitHub CLI usage

**references/git-commands.md**
- Configuration
- Basic commands
- Branching and merging
- Remote operations
- Undoing changes
- Advanced commands
- Maintenance

**references/troubleshooting.md**
- Common mistakes and fixes
- Committed to wrong branch
- Wrong commit message
- Merge conflicts
- Branch diverged
- Lost commits
- Recovery strategies

**references/advanced-git.md**
- Interactive rebase
- Cherry-pick
- Reflog
- Bisect
- Worktrees
- Submodules
- Hooks
- Performance optimization

### Templates

**assets/templates/pr-template.md**
- Complete PR template
- Type of change checklist
- Testing instructions
- Documentation checklist
- Review checklist

**assets/templates/commit-message-template.md**
- Commit message format
- Type and scope guidelines
- Examples and anti-patterns
- Quick reference

### Examples

**examples/commit-examples.md**
- Real-world commit message examples
- Features, fixes, docs, refactors
- Simple and complex commits
- Breaking changes
- Multi-scope commits
- Anti-patterns to avoid

## Quick Start

### Daily Workflow

```bash
# 1. Start from updated main
git checkout main
git pull origin main

# 2. Create feature branch
git checkout -b feat/new-feature

# 3. Make changes and commit
git add .
git commit -m "feat(module): add new functionality"

# 4. Keep branch updated
git rebase main

# 5. Push and create PR
git push -u origin feat/new-feature
# Create PR on GitHub

# 6. After merge, cleanup
git checkout main
git pull origin main
git branch -d feat/new-feature
```

### Commit Message Format

```
<type>(<scope>): <subject>

<body>

<footer>
```

**Example:**
```
feat(auth): add OAuth2 authentication

Implement OAuth2 flow for Google and GitHub providers.
Users can now sign in with their existing accounts.

Closes #123
```

### Common Commands

```bash
# Status and changes
git status
git diff
git log --oneline

# Committing
git add .
git commit -m "type(scope): message"
git commit --amend

# Branching
git checkout -b feat/branch-name
git checkout main
git branch -d feature-branch

# Syncing
git pull origin main
git push origin feature-branch
git rebase main

# Undoing
git reset HEAD~1         # Undo commit, keep changes
git reset --hard HEAD~1  # Undo commit, discard changes
git revert abc123        # Revert commit (safe for pushed commits)
```

## Philosophy

### Trunk-Based Development

Main branch is the source of truth:
- Always deployable
- Always tested
- Always up to date

Feature branches are temporary:
- Short-lived (hours to days)
- Focused on single change
- Frequently synchronized with main
- Deleted immediately after merge

### Conventional Commits

Commit messages are documentation:
- Structured format enables automation
- Clear type indicates change category
- Scope provides context
- Subject describes what changed
- Body explains why

### GitHub Workflow

Pull requests enable collaboration:
- Code review catches issues
- Discussion improves design
- CI/CD verifies quality
- Documentation of decisions

## Best Practices

### Commit Practices
- ✅ Commit frequently (multiple times per day)
- ✅ Write clear, conventional commit messages
- ✅ Keep commits focused and atomic
- ✅ Test before committing
- ✅ Review changes before committing

### Branch Practices
- ✅ Keep branches short-lived (1-3 days)
- ✅ One feature/fix per branch
- ✅ Sync from main frequently
- ✅ Delete branches after merging
- ✅ Keep main always deployable

### Collaboration Practices
- ✅ Create PR early
- ✅ Keep PRs small (<300 lines)
- ✅ Review code promptly (within hours)
- ✅ Respond to feedback quickly
- ✅ Communicate about conflicts

## Tools and Configuration

### Recommended Setup

```bash
# User configuration
git config --global user.name "Your Name"
git config --global user.email "your.email@example.com"

# Editor
git config --global core.editor "vim"

# Default branch
git config --global init.defaultBranch main

# Rebase by default
git config --global pull.rebase true

# Helpful aliases
git config --global alias.st status
git config --global alias.co checkout
git config --global alias.br branch
git config --global alias.lg "log --oneline --graph --decorate"
```

### Recommended Tools

**Commit message helpers:**
- Commitizen - Interactive commit builder
- Commitlint - Validate commit messages
- Husky - Git hooks for validation

**GitHub tools:**
- GitHub CLI (gh) - Manage PRs and issues from terminal
- GitHub Desktop - GUI for Git operations

**Merge tools:**
- VS Code - Built-in merge editor
- GitKraken, SourceTree - Visual Git clients
- Meld, KDiff3 - Dedicated merge tools

## Learning Path

### Beginner

1. Basic Git commands (status, add, commit, push, pull)
2. Simple branching and merging
3. Basic commit message format
4. Creating pull requests

**Start with:**
- SKILL.md overview
- references/git-commands.md basics
- examples/commit-examples.md simple examples

### Intermediate

1. Trunk-based development workflow
2. Conventional Commits specification
3. GitHub PR and review workflow
4. Handling merge conflicts
5. Basic troubleshooting

**Study:**
- references/trunk-based-development.md
- references/conventional-commits.md
- references/github-workflow.md
- references/troubleshooting.md common issues

### Advanced

1. Interactive rebase and history editing
2. Cherry-pick and selective changes
3. Reflog for recovery
4. Bisect for bug finding
5. Worktrees for parallel work
6. Git hooks for automation

**Master:**
- references/advanced-git.md all sections
- references/troubleshooting.md complex scenarios
- Custom workflows and automation

## Common Scenarios

### Starting New Feature

See: SKILL.md → "Starting New Feature"

### Fixing Bug

See: SKILL.md → "Quick Bug Fix"

### Updating Branch from Main

See: SKILL.md → "Updating Branch from Main"

### Resolving Merge Conflicts

See: SKILL.md → "Handling Merge Conflicts"

### Undoing Mistakes

See: references/troubleshooting.md → specific issue

### Creating Pull Request

See: references/github-workflow.md → "Creating a Pull Request"

## FAQ

**Q: When should I use merge vs rebase?**
A: Use rebase for feature branches (clean history), merge for shared branches (preserves collaboration history).

**Q: How often should I commit?**
A: Commit frequently (multiple times per day) at logical stopping points. Each commit should represent a complete, working change.

**Q: How small should PRs be?**
A: Aim for 50-200 lines changed. Definitely split if over 500 lines. Small PRs get reviewed faster and more thoroughly.

**Q: When should I use feature flags?**
A: Use feature flags when work takes more than 3 days or when you want to merge incomplete features to main without exposing to users.

**Q: Should I squash commits before merging?**
A: Generally yes for feature branches - it keeps main branch history clean. Keep detailed commits in feature branch for review.

**Q: How do I recover from mistakes?**
A: Check reflog first (`git reflog`), which tracks all ref updates. Most mistakes are recoverable. See references/troubleshooting.md for specific scenarios.

## Contributing to the Skill

Found an issue or have a suggestion? The skill can be improved through:

1. Additional examples
2. More troubleshooting scenarios
3. Better explanations
4. Additional references
5. Tool recommendations

## Resources

- **Official Git Documentation:** https://git-scm.com/doc
- **Conventional Commits Spec:** https://www.conventionalcommits.org/
- **Trunk-Based Development:** https://trunkbaseddevelopment.com/
- **GitHub Docs:** https://docs.github.com/
- **Pro Git Book:** https://git-scm.com/book/

---

**Remember:** Good version control is about clear history, frequent integration, and effective collaboration. Master the basics, then gradually adopt advanced techniques as needed.
