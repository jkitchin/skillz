# Advanced Git Techniques

Advanced Git features and workflows for power users.

## Interactive Rebase

### Basic Interactive Rebase

```bash
# Rebase last 5 commits
git rebase -i HEAD~5

# Rebase since specific commit
git rebase -i abc123

# Rebase entire branch
git rebase -i $(git merge-base HEAD main)
```

### Rebase Commands

**In the interactive editor:**
```
pick abc123 feat(auth): add OAuth
pick def456 fix(auth): handle errors
pick ghi789 refactor(auth): simplify logic
pick jkl012 test(auth): add integration tests
```

**Available commands:**
- `pick` (p) - Use commit as-is
- `reword` (r) - Use commit, but edit message
- `edit` (e) - Use commit, but pause to amend
- `squash` (s) - Combine with previous commit, keep message
- `fixup` (f) - Combine with previous commit, discard message
- `drop` (d) - Remove commit
- `exec` (x) - Run shell command

### Common Rebase Scenarios

**Squash multiple commits:**
```
pick abc123 feat(auth): add OAuth
fixup def456 fix(auth): handle errors
fixup ghi789 refactor(auth): simplify logic
pick jkl012 test(auth): add integration tests
```

Result: 2 commits instead of 4

**Reword commit messages:**
```
pick abc123 feat(auth): add OAuth
reword def456 fix(auth): handle errors
pick ghi789 refactor(auth): simplify logic
```

**Reorder commits:**
```
pick jkl012 test(auth): add integration tests
pick abc123 feat(auth): add OAuth
pick def456 fix(auth): handle errors
```

**Edit commit:**
```
pick abc123 feat(auth): add OAuth
edit def456 fix(auth): handle errors
pick ghi789 refactor(auth): simplify logic
```

When rebase pauses:
```bash
# Make changes
git add changed-file.ts
git commit --amend
git rebase --continue
```

**Run tests between commits:**
```
pick abc123 feat(auth): add OAuth
exec npm test
pick def456 fix(auth): handle errors
exec npm test
```

### Auto-Squashing

**Mark commits for auto-squash:**
```bash
# Original commit
git commit -m "feat(auth): add OAuth"

# Fixup for previous commit
git commit --fixup abc123

# Later, auto-squash
git rebase -i --autosquash HEAD~5
```

**Benefits:**
- Commits automatically reordered
- Marked as fixup automatically
- Great for fixing up commits during review

## Cherry-Pick

### Basic Cherry-Pick

```bash
# Apply specific commit to current branch
git cherry-pick abc123

# Cherry-pick with commit message edit
git cherry-pick -e abc123

# Cherry-pick without committing
git cherry-pick -n abc123
git cherry-pick --no-commit abc123
```

### Cherry-Pick Multiple Commits

```bash
# Cherry-pick multiple commits
git cherry-pick abc123 def456 ghi789

# Cherry-pick range (exclusive start)
git cherry-pick abc123..ghi789

# Cherry-pick range (inclusive)
git cherry-pick abc123^..ghi789
```

### Cherry-Pick Use Cases

**Port fix to release branch:**
```bash
# Fix is on main
git checkout release-1.0
git cherry-pick abc123
git push origin release-1.0
```

**Extract commits from branch:**
```bash
# Get specific commits from experimental branch
git checkout main
git cherry-pick exp-branch~3
git cherry-pick exp-branch~1
```

### Handling Cherry-Pick Conflicts

```bash
# Start cherry-pick
git cherry-pick abc123
# CONFLICT

# Resolve conflicts
# Edit conflicted files
git add resolved-file.txt

# Continue
git cherry-pick --continue

# Or abort
git cherry-pick --abort
```

## Reflog

### Understanding Reflog

Reflog tracks changes to HEAD and branch tips. It's your safety net!

```bash
# View reflog
git reflog

# Output:
# abc123 HEAD@{0}: commit: feat(auth): add OAuth
# def456 HEAD@{1}: checkout: moving from main to feature
# ghi789 HEAD@{2}: commit: fix(api): handle errors
```

### Reflog Commands

```bash
# View HEAD reflog
git reflog
git reflog show HEAD

# View branch reflog
git reflog show main
git reflog show feature-branch

# View with dates
git reflog --date=relative
git reflog --date=iso
```

### Recovery with Reflog

**Undo reset:**
```bash
# Oops, reset too far
git reset --hard HEAD~5

# Find previous state
git reflog
# abc123 HEAD@{1}: reset: moving to HEAD~5
# def456 HEAD@{2}: commit: my important commit

# Recover
git reset --hard HEAD@{2}
```

**Recover deleted branch:**
```bash
# Deleted branch
git branch -D feature-branch

# Find in reflog
git reflog show feature-branch

# Or search for last commit
git reflog | grep feature-branch

# Recreate branch
git branch feature-branch abc123
```

**Recover after bad rebase:**
```bash
# Before rebase
git reflog
# abc123 HEAD@{0}: checkout: moving from main to feature

# After bad rebase
git reset --hard HEAD@{1}
```

### Reflog Expiration

```bash
# Reflog entries expire after 90 days (default)
# Can be configured
git config gc.reflogExpire 120.days

# Expire now (dangerous!)
git reflog expire --expire=now --all

# Cleanup
git gc --prune=now
```

## Bisect

### Finding Bugs with Bisect

Binary search through commits to find bug introduction.

```bash
# Start bisect
git bisect start

# Mark current commit as bad
git bisect bad

# Mark last known good commit
git bisect good abc123

# Git checks out middle commit
# Test it and mark:
git bisect bad   # If bug exists
git bisect good  # If bug doesn't exist

# Repeat until found
# Git will tell you: "abc123 is the first bad commit"

# End bisect
git bisect reset
```

### Automated Bisect

```bash
# Bisect with automated testing
git bisect start HEAD abc123

# Run test automatically
git bisect run npm test

# Or custom script
git bisect run ./test-script.sh

# Git automatically finds bad commit
git bisect reset
```

### Bisect Script Example

```bash
#!/bin/bash
# test-script.sh

# Run tests
npm test

# Return 0 for good, 1 for bad
exit $?
```

### Bisect Use Cases

**Find performance regression:**
```bash
#!/bin/bash
# Performance threshold: 100ms

time=$( { time node app.js; } 2>&1 | grep real | awk '{print $2}' )

if [ "$time" -lt "0.1" ]; then
  exit 0  # Good
else
  exit 1  # Bad (too slow)
fi
```

**Find when feature broke:**
```bash
#!/bin/bash

# Test specific feature
if node -e "require('./app').testFeature()" 2>/dev/null; then
  exit 0  # Good
else
  exit 1  # Bad
fi
```

## Worktrees

### Understanding Worktrees

Work on multiple branches simultaneously without switching.

```bash
# List worktrees
git worktree list

# Main worktree
/path/to/repo         abc123 [main]

# Additional worktrees
/path/to/repo-feature def456 [feature-branch]
```

### Creating Worktrees

```bash
# Create worktree for existing branch
git worktree add ../project-feature feature-branch

# Create worktree with new branch
git worktree add -b new-feature ../project-new main

# Create worktree from specific commit
git worktree add ../project-old abc123
```

### Using Worktrees

```bash
# Work in main worktree
cd /path/to/repo
git checkout main

# Simultaneously work in feature branch
cd /path/to/repo-feature
# Already on feature-branch
git add .
git commit -m "feat: add feature"

# Or work on bug fix
cd /path/to/repo-bugfix
git checkout -b fix/urgent-bug
```

### Managing Worktrees

```bash
# Remove worktree
git worktree remove ../project-feature

# Prune deleted worktrees
git worktree prune

# Lock worktree (prevent pruning)
git worktree lock ../project-feature

# Unlock
git worktree unlock ../project-feature
```

### Worktree Use Cases

**Review PR while working:**
```bash
# Continue working on feature
cd /main/repo

# Review PR in separate worktree
git worktree add ../repo-review pr-branch
cd ../repo-review
# Test and review
```

**Build while developing:**
```bash
# Develop in main repo
cd /repo

# Build in worktree
git worktree add ../repo-build main
cd ../repo-build
npm run build
```

## Submodules

### Adding Submodules

```bash
# Add submodule
git submodule add https://github.com/user/lib.git libs/lib

# Add to specific directory
git submodule add https://github.com/user/lib.git path/to/lib

# Commit submodule
git add .gitmodules libs/lib
git commit -m "chore: add lib submodule"
```

### Cloning with Submodules

```bash
# Clone with submodules
git clone --recursive https://github.com/user/repo.git

# Or after cloning
git clone https://github.com/user/repo.git
cd repo
git submodule init
git submodule update

# Or in one command
git submodule update --init --recursive
```

### Updating Submodules

```bash
# Update submodule to latest
cd libs/lib
git pull origin main
cd ../..
git add libs/lib
git commit -m "chore: update lib submodule"

# Update all submodules
git submodule update --remote --merge

# Update specific submodule
git submodule update --remote --merge libs/lib
```

### Working with Submodules

```bash
# Status of submodules
git submodule status

# Execute command in all submodules
git submodule foreach 'git pull origin main'

# Remove submodule
git submodule deinit libs/lib
git rm libs/lib
rm -rf .git/modules/libs/lib
git commit -m "chore: remove lib submodule"
```

## Subtrees

### Git Subtree vs Submodule

**Subtree benefits:**
- No special clone steps
- Subtree code is part of main repo
- Easier for contributors

**Submodule benefits:**
- Smaller repository size
- Clear separation
- Independent versioning

### Adding Subtree

```bash
# Add subtree
git subtree add --prefix=libs/lib \
  https://github.com/user/lib.git main --squash

# Squash combines all history into one commit
```

### Updating Subtree

```bash
# Pull updates from subtree
git subtree pull --prefix=libs/lib \
  https://github.com/user/lib.git main --squash
```

### Contributing Back

```bash
# Push changes back to subtree
git subtree push --prefix=libs/lib \
  https://github.com/user/lib.git main
```

## Advanced Merging

### Merge Strategies

```bash
# Recursive (default)
git merge feature-branch

# Ours (always prefer our version)
git merge -X ours feature-branch

# Theirs (always prefer their version)
git merge -X theirs feature-branch

# Octopus (merge multiple branches)
git merge branch1 branch2 branch3

# Ours strategy (different from -X ours)
git merge -s ours feature-branch  # Discards their changes
```

### Merge Drivers

**Custom merge for specific files:**

```bash
# .gitattributes
package-lock.json merge=npm-merge-driver

# .git/config
[merge "npm-merge-driver"]
    name = NPM lock file merge driver
    driver = npx npm-merge-driver merge %A %O %B %P
```

### Merge Without Commit

```bash
# Merge but don't commit
git merge --no-commit feature-branch

# Review merged result
git diff --cached

# Commit or abort
git commit -m "merge: feature-branch"
# or
git merge --abort
```

## Advanced Diffing

### Diff Algorithms

```bash
# Default diff
git diff

# Patience algorithm (better for refactoring)
git diff --patience

# Histogram algorithm (faster patience)
git diff --histogram

# Minimal diff
git diff --minimal
```

### Word and Character Diff

```bash
# Word-level diff
git diff --word-diff

# Character-level diff
git diff --word-diff=color
git diff --color-words

# Word diff with regex
git diff --word-diff-regex='[^[:space:],]+'
```

### Advanced Diff Options

```bash
# Ignore whitespace changes
git diff -w
git diff --ignore-all-space

# Ignore whitespace at line endings
git diff --ignore-space-at-eol

# Ignore blank lines
git diff --ignore-blank-lines

# Show function/class name
git diff --show-function

# Unified diff with more context
git diff -U10  # 10 lines of context
```

### Diff with External Tools

```bash
# Configure difftool
git config --global diff.tool vimdiff
git config --global difftool.prompt false

# Use difftool
git difftool

# With specific tool
git difftool --tool=meld

# For merge conflicts
git mergetool
```

## Hooks

### Client-Side Hooks

Located in `.git/hooks/`:

**pre-commit:**
```bash
#!/bin/bash
# .git/hooks/pre-commit

# Run linter
npm run lint
if [ $? -ne 0 ]; then
  echo "Linting failed"
  exit 1
fi

# Run tests
npm test
if [ $? -ne 0 ]; then
  echo "Tests failed"
  exit 1
fi
```

**commit-msg:**
```bash
#!/bin/bash
# .git/hooks/commit-msg

# Validate conventional commit format
commit_msg=$(cat "$1")
pattern="^(feat|fix|docs|style|refactor|test|chore)(\(.+\))?: .+$"

if ! echo "$commit_msg" | grep -qE "$pattern"; then
  echo "Invalid commit message format"
  echo "Expected: type(scope): subject"
  exit 1
fi
```

**pre-push:**
```bash
#!/bin/bash
# .git/hooks/pre-push

# Run full test suite before push
npm run test:all
exit $?
```

### Hook Management

**Using Husky:**
```bash
# Install husky
npm install --save-dev husky

# Initialize
npx husky init

# Add hook
npx husky add .husky/pre-commit "npm run lint"
npx husky add .husky/commit-msg "npx commitlint --edit $1"
```

## Advanced Configuration

### Conditional Config

```bash
# ~/.gitconfig

[includeIf "gitdir:~/work/"]
    path = ~/.gitconfig-work

[includeIf "gitdir:~/personal/"]
    path = ~/.gitconfig-personal
```

```bash
# ~/.gitconfig-work
[user]
    name = Work Name
    email = work@company.com

[url "git@github.com-work:"]
    insteadOf = git@github.com:
```

### Custom Commands

```bash
# Add to ~/.gitconfig

[alias]
    # Pretty log
    lg = log --graph --pretty=format:'%Cred%h%Creset -%C(yellow)%d%Creset %s %Cgreen(%cr) %C(bold blue)<%an>%Creset' --abbrev-commit

    # Undo last commit
    undo = reset --soft HEAD^

    # Amend without editing message
    amend = commit --amend --no-edit

    # List branches by date
    branches = for-each-ref --sort=-committerdate refs/heads/ --format='%(committerdate:short) %(refname:short)'

    # Show contributors
    contributors = shortlog -sn

    # Delete merged branches
    cleanup = !git branch --merged | grep -v '\\*\\|main\\|develop' | xargs -n 1 git branch -d
```

## Performance Optimization

### Sparse Checkout

Only checkout specific directories:

```bash
# Enable sparse checkout
git sparse-checkout init --cone

# Set directories to check out
git sparse-checkout set src/ tests/

# Add more directories
git sparse-checkout add docs/

# List current sparse checkout
git sparse-checkout list

# Disable
git sparse-checkout disable
```

### Shallow Clone

```bash
# Clone with limited history
git clone --depth 1 https://github.com/user/repo.git

# Fetch more history later
git fetch --deepen=100  # Fetch 100 more commits
git fetch --unshallow   # Fetch all history
```

### Partial Clone

```bash
# Clone without blobs (files)
git clone --filter=blob:none https://github.com/user/repo.git

# Clone without trees
git clone --filter=tree:0 https://github.com/user/repo.git

# Fetch specific blob
git checkout HEAD -- file.txt
```

## Tips and Tricks

### Quick Commands

```bash
# Stage and commit in one line
git commit -am "message"

# Push current branch
git push origin HEAD

# Delete last commit (keep changes)
git reset HEAD^

# Show files changed in commit
git show --name-only abc123

# Count commits
git rev-list --count HEAD

# Find commit that deleted file
git log --all --full-history -- path/to/file

# Show commit date
git show -s --format=%ci abc123
```

### Searching

```bash
# Find commits containing text
git log -S "search text" --source --all

# Find commits with text in message
git log --grep="bug fix"

# Find commits by author
git log --author="John"

# Find commits touching specific lines
git log -L 10,20:file.txt
```

### Advanced Blame

```bash
# Ignore whitespace in blame
git blame -w file.txt

# Ignore commit in blame
git blame --ignore-rev abc123 file.txt

# Ignore commits from file
git config blame.ignoreRevsFile .git-blame-ignore-revs
```

---

**Remember:** These advanced techniques are powerful but use them judiciously. The most important thing is clear, understandable history. Advanced features should support that goal, not complicate it.
