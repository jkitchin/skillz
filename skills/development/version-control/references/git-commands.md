# Git Commands Reference

Comprehensive reference for all essential Git commands.

## Configuration

### User Information

```bash
# Set name and email (required for commits)
git config --global user.name "Your Name"
git config --global user.email "your.email@example.com"

# View current config
git config --global --list
git config user.name
git config user.email

# Set per-repository (omit --global)
git config user.name "Work Name"
git config user.email "work@company.com"
```

### Editor

```bash
# Set default editor
git config --global core.editor "vim"
git config --global core.editor "code --wait"  # VS Code
git config --global core.editor "emacs"

# Set merge tool
git config --global merge.tool vimdiff
```

### Defaults

```bash
# Set default branch name
git config --global init.defaultBranch main

# Rebase by default on pull
git config --global pull.rebase true

# Enable color
git config --global color.ui auto

# Set line endings
git config --global core.autocrlf input  # Mac/Linux
git config --global core.autocrlf true   # Windows
```

### Aliases

```bash
# Shortcuts
git config --global alias.st status
git config --global alias.co checkout
git config --global alias.br branch
git config --global alias.cm commit
git config --global alias.unstage 'reset HEAD --'
git config --global alias.last 'log -1 HEAD'
git config --global alias.lg 'log --oneline --graph --decorate --all'

# Usage
git st  # Same as: git status
git co main  # Same as: git checkout main
git lg  # Pretty log
```

## Getting Started

### Initialize Repository

```bash
# Create new repository
git init
git init my-project

# Clone existing repository
git clone https://github.com/user/repo.git
git clone https://github.com/user/repo.git my-folder
git clone --depth 1 https://github.com/user/repo.git  # Shallow clone
```

### Repository Status

```bash
# View status
git status

# Short status
git status -s
git status --short

# View changes
git diff              # Unstaged changes
git diff --staged     # Staged changes
git diff HEAD         # All changes
```

## Committing Changes

### Staging Files

```bash
# Stage specific file
git add file.txt

# Stage multiple files
git add file1.txt file2.txt

# Stage all changes in directory
git add src/

# Stage all changes in repository
git add .
git add --all
git add -A

# Stage all modified/deleted (not new files)
git add -u

# Interactive staging
git add -i
git add -p  # Patch mode (stage parts of files)
```

### Unstaging Files

```bash
# Unstage specific file
git reset HEAD file.txt
git restore --staged file.txt  # Modern alternative

# Unstage all files
git reset HEAD
```

### Creating Commits

```bash
# Commit staged changes
git commit -m "feat(auth): add OAuth support"

# Commit with detailed message (opens editor)
git commit

# Stage and commit in one step
git commit -am "fix(api): handle null responses"

# Amend last commit (change message or add files)
git commit --amend
git commit --amend --no-edit  # Keep message
git commit --amend -m "fix(api): corrected message"
```

### Viewing History

```bash
# View commit history
git log

# One line per commit
git log --oneline

# Graph view
git log --graph --oneline --decorate --all
git log --graph --all

# Last N commits
git log -5
git log -n 5

# Show changes in commits
git log -p
git log --patch

# Show stats
git log --stat

# Search commits
git log --grep="bug fix"
git log --grep="feat"

# By author
git log --author="John"

# By date
git log --since="2 weeks ago"
git log --after="2024-01-01"
git log --before="2024-12-31"

# Show changes to specific file
git log --follow file.txt
git log -p file.txt
```

### Viewing Changes

```bash
# View specific commit
git show abc123
git show HEAD
git show HEAD~1  # Previous commit

# View specific file in commit
git show abc123:path/to/file.txt

# Compare commits
git diff abc123 def456
git diff abc123..def456  # Same as above
git diff abc123...def456  # Since common ancestor

# Compare branches
git diff main feature-branch
git diff main..feature-branch
```

## Branching

### Branch Management

```bash
# List branches
git branch              # Local branches
git branch -a           # All branches (local + remote)
git branch -r           # Remote branches
git branch -v           # With last commit
git branch -vv          # With tracking info

# Create branch
git branch new-feature
git branch new-feature abc123  # From specific commit

# Create and switch to branch
git checkout -b new-feature
git switch -c new-feature  # Modern alternative

# Switch branches
git checkout main
git switch main  # Modern alternative

# Rename branch
git branch -m old-name new-name
git branch -m new-name  # Rename current branch

# Delete branch
git branch -d feature-branch  # Safe delete (merged only)
git branch -D feature-branch  # Force delete

# Delete remote branch
git push origin --delete feature-branch
git push origin :feature-branch  # Alternative syntax
```

### Merging

```bash
# Merge branch into current
git merge feature-branch

# Merge with commit (even if fast-forward possible)
git merge --no-ff feature-branch

# Squash merge (combine all commits)
git merge --squash feature-branch
git commit -m "feat: add feature"

# Abort merge
git merge --abort
```

### Rebasing

```bash
# Rebase current branch onto main
git rebase main

# Interactive rebase (last 3 commits)
git rebase -i HEAD~3
git rebase -i abc123  # Since specific commit

# Continue after resolving conflicts
git rebase --continue

# Skip problematic commit
git rebase --skip

# Abort rebase
git rebase --abort
```

## Remote Repositories

### Remote Management

```bash
# List remotes
git remote
git remote -v  # With URLs

# Add remote
git remote add origin https://github.com/user/repo.git

# Change remote URL
git remote set-url origin https://github.com/user/new-repo.git

# Remove remote
git remote remove origin

# Rename remote
git remote rename origin upstream

# Show remote info
git remote show origin
```

### Fetching

```bash
# Fetch all remotes
git fetch

# Fetch specific remote
git fetch origin

# Fetch specific branch
git fetch origin main

# Fetch and prune deleted branches
git fetch --prune
git fetch -p

# View fetched commits
git log origin/main
git diff origin/main
```

### Pulling

```bash
# Pull (fetch + merge)
git pull

# Pull specific branch
git pull origin main

# Pull with rebase
git pull --rebase
git pull --rebase origin main

# Pull and auto-stash changes
git pull --autostash
```

### Pushing

```bash
# Push to origin
git push

# Push specific branch
git push origin feature-branch

# Push and set upstream
git push -u origin feature-branch
git push --set-upstream origin feature-branch

# Push all branches
git push --all

# Push tags
git push --tags

# Force push (dangerous!)
git push --force
git push --force-with-lease  # Safer alternative

# Delete remote branch
git push origin --delete feature-branch
```

## Undoing Changes

### Working Directory

```bash
# Discard changes in file
git checkout -- file.txt
git restore file.txt  # Modern alternative

# Discard all changes
git checkout -- .
git restore .

# Remove untracked files
git clean -n  # Dry run (show what would be deleted)
git clean -f  # Delete untracked files
git clean -fd  # Delete files and directories
git clean -fdx  # Include ignored files
```

### Staged Changes

```bash
# Unstage file
git reset HEAD file.txt
git restore --staged file.txt  # Modern alternative

# Unstage all
git reset HEAD
```

### Commits

```bash
# Undo last commit, keep changes
git reset --soft HEAD~1

# Undo last commit, keep changes unstaged
git reset HEAD~1
git reset --mixed HEAD~1  # Same as above

# Undo last commit, discard changes
git reset --hard HEAD~1

# Undo multiple commits
git reset --hard HEAD~3  # Last 3 commits

# Undo to specific commit
git reset --hard abc123

# Revert commit (creates new commit)
git revert abc123
git revert HEAD
git revert --no-commit abc123  # Revert but don't commit yet
```

## Stashing

### Basic Stashing

```bash
# Stash current changes
git stash
git stash push  # Modern syntax

# Stash with message
git stash push -m "WIP: feature in progress"

# Stash including untracked files
git stash -u
git stash --include-untracked

# Stash including ignored files
git stash -a
git stash --all
```

### Applying Stashes

```bash
# List stashes
git stash list

# Apply most recent stash
git stash apply

# Apply specific stash
git stash apply stash@{2}

# Apply and remove stash
git stash pop
git stash pop stash@{2}

# View stash contents
git stash show
git stash show -p  # Show changes
git stash show stash@{2}
```

### Managing Stashes

```bash
# Create branch from stash
git stash branch new-branch

# Delete stash
git stash drop
git stash drop stash@{2}

# Delete all stashes
git stash clear
```

## Tags

### Creating Tags

```bash
# Lightweight tag
git tag v1.0.0

# Annotated tag (recommended)
git tag -a v1.0.0 -m "Release version 1.0.0"

# Tag specific commit
git tag -a v1.0.0 abc123 -m "Release version 1.0.0"

# Tag with long message (opens editor)
git tag -a v1.0.0
```

### Viewing Tags

```bash
# List tags
git tag
git tag -l
git tag --list

# Search tags
git tag -l "v1.*"

# Show tag info
git show v1.0.0
```

### Managing Tags

```bash
# Push tag to remote
git push origin v1.0.0

# Push all tags
git push --tags
git push origin --tags

# Delete local tag
git tag -d v1.0.0

# Delete remote tag
git push origin --delete v1.0.0
git push origin :refs/tags/v1.0.0  # Alternative

# Checkout tag
git checkout v1.0.0
```

## Advanced Commands

### Cherry-Pick

```bash
# Apply specific commit to current branch
git cherry-pick abc123

# Cherry-pick multiple commits
git cherry-pick abc123 def456

# Cherry-pick range
git cherry-pick abc123..def456

# Cherry-pick without committing
git cherry-pick -n abc123
git cherry-pick --no-commit abc123

# Abort cherry-pick
git cherry-pick --abort
```

### Reflog

```bash
# View reference log
git reflog
git reflog show HEAD

# Reflog for specific branch
git reflog show main

# View reflog with dates
git reflog --date=relative

# Recover using reflog
git reset --hard HEAD@{2}
git checkout HEAD@{5}
```

### Bisect

```bash
# Start bisect
git bisect start

# Mark current as bad
git bisect bad

# Mark last known good
git bisect good abc123

# Git checks out middle commit - test it
# If bad:
git bisect bad

# If good:
git bisect good

# Repeat until found
# Then reset
git bisect reset

# Automated bisect
git bisect start HEAD abc123
git bisect run npm test
```

### Blame

```bash
# Show who changed each line
git blame file.txt

# Show line numbers
git blame -L 10,20 file.txt

# Show email addresses
git blame -e file.txt

# Ignore whitespace changes
git blame -w file.txt
```

### Worktrees

```bash
# List worktrees
git worktree list

# Add worktree
git worktree add ../project-feature feature-branch
git worktree add -b new-branch ../project-new main

# Remove worktree
git worktree remove ../project-feature

# Prune worktrees
git worktree prune
```

## Inspection

### Search

```bash
# Search for text in files
git grep "search term"
git grep -n "search term"  # Show line numbers
git grep -i "search term"  # Case insensitive

# Search in specific commit
git grep "search term" abc123

# Search with context
git grep -C 3 "search term"  # 3 lines before and after
```

### Differences

```bash
# Compare working directory to staging
git diff

# Compare staging to HEAD
git diff --staged
git diff --cached

# Compare working directory to HEAD
git diff HEAD

# Compare two commits
git diff abc123 def456

# Compare two branches
git diff main..feature-branch

# Show only names of changed files
git diff --name-only
git diff --name-status

# Show word-level diff
git diff --word-diff
```

### History

```bash
# Show commits that modified file
git log file.txt
git log --follow file.txt  # Follow renames

# Show commits in branch A not in branch B
git log main..feature-branch
git log feature-branch ^main  # Same as above

# Show commits in either branch but not both
git log main...feature-branch

# Show merge commits only
git log --merges

# Show non-merge commits only
git log --no-merges
```

## Maintenance

### Cleanup

```bash
# Remove untracked files
git clean -fd

# Remove stale remote tracking branches
git remote prune origin
git fetch --prune

# Remove deleted branches from local
git branch -d $(git branch --merged main | grep -v "main")

# Garbage collection
git gc

# Aggressive garbage collection
git gc --aggressive --prune=now
```

### Repository Info

```bash
# Show repository size
git count-objects -vH

# Show largest files in history
git rev-list --objects --all |
  git cat-file --batch-check='%(objecttype) %(objectname) %(objectsize) %(rest)' |
  sed -n 's/^blob //p' |
  sort -nk2 |
  tail -n 10

# Verify repository
git fsck
```

## Configuration Files

### .gitignore

```bash
# Create .gitignore
cat > .gitignore << EOF
# Dependencies
node_modules/
vendor/

# Build outputs
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

# Logs
*.log
logs/
EOF
```

### .gitattributes

```bash
# Create .gitattributes
cat > .gitattributes << EOF
# Auto detect text files and normalize line endings
* text=auto

# Explicitly declare text files
*.js text
*.jsx text
*.ts text
*.tsx text
*.json text
*.md text

# Declare files as binary
*.png binary
*.jpg binary
*.gif binary
*.pdf binary

# Export ignore (not included in archives)
.gitattributes export-ignore
.gitignore export-ignore
tests/ export-ignore
EOF
```

## Useful Combinations

### Daily Workflow

```bash
# Update main and rebase feature branch
git checkout main && git pull && git checkout - && git rebase main

# Commit all changes
git add -A && git commit -m "type: message"

# Push and create PR
git push -u origin $(git branch --show-current)
```

### Cleanup

```bash
# Update main, delete merged branches
git checkout main && \
  git pull && \
  git branch --merged | grep -v "main" | xargs git branch -d
```

### Undo

```bash
# Undo last commit, keep changes
git reset --soft HEAD~1

# Discard all local changes
git reset --hard HEAD && git clean -fd
```

### Search and Fix

```bash
# Find commits that modified specific text
git log -S "search text" --source --all

# Find which commit introduced bug
git bisect start HEAD abc123
git bisect run npm test
```

---

**Pro tip:** Use `git help <command>` or `git <command> --help` for detailed documentation on any command.
