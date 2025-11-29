# Git Troubleshooting Guide

Solutions to common Git problems and mistakes.

## Committed to Wrong Branch

### Problem
Made commits on main instead of feature branch.

### Solution

**Option 1: Move commits to new branch**
```bash
# Save current state to new branch
git branch feat/new-feature

# Reset main to before your commits
git reset --hard origin/main

# Switch to new branch (commits are there)
git checkout feat/new-feature
```

**Option 2: Move specific commits**
```bash
# Create new branch from current commit
git checkout -b feat/new-feature

# Go back to main
git checkout main

# Reset to origin
git reset --hard origin/main

# Cherry-pick specific commits to feature branch
git checkout feat/new-feature
git cherry-pick abc123 def456
```

## Wrong Commit Message

### Problem
Typo in commit message or used wrong type.

### Solution

**Last commit only (not pushed):**
```bash
# Change message
git commit --amend -m "fix(api): correct commit message"

# Edit in editor
git commit --amend
```

**Last commit (already pushed):**
```bash
# Amend locally
git commit --amend -m "fix(api): correct commit message"

# Force push (be careful!)
git push --force-with-lease origin feature-branch
```

**Older commit:**
```bash
# Interactive rebase
git rebase -i HEAD~3

# Change 'pick' to 'reword' for commit to fix
# Save and close, then edit message
```

## Accidentally Committed Sensitive Data

### Problem
Committed .env file, API keys, or passwords.

### Solution

**If not pushed yet:**
```bash
# Remove from last commit
git rm --cached .env
git commit --amend --no-edit

# Or reset and recommit
git reset --soft HEAD~1
git rm --cached .env
git commit
```

**If already pushed (more serious):**
```bash
# Remove from history using filter-branch
git filter-branch --force --index-filter \
  "git rm --cached --ignore-unmatch .env" \
  --prune-empty --tag-name-filter cat -- --all

# Force push (coordinate with team!)
git push --force --all
git push --force --tags

# Or use BFG Repo-Cleaner (recommended)
bfg --delete-files .env
git reflog expire --expire=now --all
git gc --prune=now --aggressive
```

**Important:**
- 🚨 Rotate/invalidate the exposed secrets immediately!
- 🚨 Consider repository compromised until secrets rotated
- Add file to .gitignore to prevent future commits

## Need to Undo Last Commit

### Keep Changes (Uncommit)

```bash
# Uncommit, keep changes staged
git reset --soft HEAD~1

# Uncommit, keep changes unstaged
git reset HEAD~1
git reset --mixed HEAD~1  # Same as above
```

### Discard Changes

```bash
# Completely remove commit and changes
git reset --hard HEAD~1

# For pushed commits, use revert instead
git revert HEAD
git push origin feature-branch
```

## Need to Undo Multiple Commits

### Solution

```bash
# Undo last 3 commits, keep changes
git reset HEAD~3

# Undo last 3 commits, discard changes
git reset --hard HEAD~3

# Reset to specific commit
git reset --hard abc123

# For pushed commits, use revert
git revert HEAD~2..HEAD
```

## Merge Conflicts

### Understanding Conflict Markers

```python
<<<<<<< HEAD (Current Change)
# Your version
def calculate(x):
    return x * 2
=======
# Their version
def calculate(x):
    return x ** 2
>>>>>>> main (Incoming Change)
```

### Resolution Steps

**1. See which files have conflicts:**
```bash
git status
# Shows: "both modified: file.py"
```

**2. Open file and resolve:**
```python
# Choose one version, combine, or write new solution
def calculate(x):
    return x ** 2  # Chose their version
```

**3. Mark as resolved:**
```bash
git add file.py
```

**4. Continue merge/rebase:**
```bash
# For merge
git commit

# For rebase
git rebase --continue
```

**5. Abort if needed:**
```bash
# For merge
git merge --abort

# For rebase
git rebase --abort
```

### Tools for Conflicts

**Use merge tool:**
```bash
git mergetool
```

**View differences:**
```bash
# See your version vs their version
git diff --ours
git diff --theirs

# See base version vs yours
git diff --base
```

## Branch Diverged from Remote

### Problem

```
Your branch and 'origin/feature-branch' have diverged,
and have 3 and 2 different commits each, respectively.
```

### Solution

**Option 1: Rebase (recommended)**
```bash
# Rebase your commits on top of remote
git pull --rebase origin feature-branch

# Or explicitly
git fetch origin
git rebase origin/feature-branch

# Resolve conflicts if any
# Then force push
git push --force-with-lease origin feature-branch
```

**Option 2: Merge**
```bash
# Merge remote changes
git pull origin feature-branch

# Resolve conflicts if any
# Then push
git push origin feature-branch
```

**Option 3: Reset to remote (discard local commits)**
```bash
git fetch origin
git reset --hard origin/feature-branch
```

## Rebase Went Wrong

### Problem
Messed up during interactive rebase or too many conflicts.

### Solution

**Abort current rebase:**
```bash
git rebase --abort
```

**If already completed but wrong:**
```bash
# Find previous state in reflog
git reflog

# Look for "checkout: moving from..."
# Reset to before rebase
git reset --hard HEAD@{5}  # Adjust number
```

## Lost Commits

### Problem
Accidentally reset too far or lost work.

### Solution

**Find lost commits:**
```bash
# View reflog
git reflog

# Look for your commit
# Output shows: abc123 HEAD@{2}: commit: my lost commit
```

**Recover commit:**
```bash
# Reset to lost commit
git reset --hard abc123

# Or cherry-pick it
git cherry-pick abc123

# Or create branch from it
git branch recovered-work abc123
```

**Find dangling commits:**
```bash
# Find all unreachable commits
git fsck --lost-found

# Show dangling commits
git fsck --lost-found | grep commit
```

## Can't Push to Remote

### Problem: Rejected (non-fast-forward)

```
! [rejected]        main -> main (non-fast-forward)
error: failed to push some refs
```

**Solution:**
```bash
# Pull first (if working on shared branch)
git pull origin main
git push origin main

# Or rebase (for cleaner history)
git pull --rebase origin main
git push origin main

# Force push (only if working alone!)
git push --force-with-lease origin feature-branch
```

### Problem: Protected Branch

```
remote: error: GH006: Protected branch update failed
```

**Solution:**
- Can't push directly to main (by design!)
- Create pull request instead
- Or temporarily disable protection (if admin)

### Problem: Authentication Failed

```
remote: Invalid username or password
```

**Solution:**
```bash
# Use personal access token instead of password
# GitHub → Settings → Developer settings → Personal access tokens

# Or use SSH
git remote set-url origin git@github.com:user/repo.git

# Set up SSH key
ssh-keygen -t ed25519 -C "your.email@example.com"
cat ~/.ssh/id_ed25519.pub  # Add to GitHub
```

## Detached HEAD State

### Problem

```
You are in 'detached HEAD' state
```

Happens when checking out specific commit or tag.

### Solution

**Just looking around (no changes):**
```bash
# Go back to branch
git checkout main
```

**Made changes you want to keep:**
```bash
# Create branch from current state
git checkout -b new-branch-name

# Or
git branch new-branch-name
git checkout main
```

**Made changes you don't want:**
```bash
# Just checkout a branch
git checkout main
# Changes are discarded
```

## Accidentally Deleted Branch

### Problem
Deleted branch but need it back.

### Solution

```bash
# Find branch in reflog
git reflog

# Look for "checkout: moving from deleted-branch"
# Note the commit hash

# Recreate branch
git branch deleted-branch abc123
```

## Large File Causing Issues

### Problem
Committed large file, push fails or repo too large.

### Solution

**Remove from last commit:**
```bash
git rm --cached large-file.zip
git commit --amend --no-edit
```

**Remove from history:**
```bash
# Use BFG Repo-Cleaner
bfg --strip-blobs-bigger-than 50M

# Or git filter-branch
git filter-branch --tree-filter 'rm -f large-file.zip' HEAD
```

**Use Git LFS for large files:**
```bash
# Install Git LFS
git lfs install

# Track large files
git lfs track "*.zip"
git add .gitattributes
git add large-file.zip
git commit -m "chore: add large file with LFS"
```

## Merge vs Rebase Confusion

### When to Merge

```bash
git merge feature-branch
```

**Use when:**
- Working on shared/public branches
- Want to preserve complete history
- Multiple people on branch
- Creating release merge

### When to Rebase

```bash
git rebase main
```

**Use when:**
- Working on local/private branch
- Want clean, linear history
- Keeping feature branch updated
- Before creating PR

**Never rebase:**
- Public/shared branches
- Main/master branch
- After pushing (unless alone on branch)

## Stash Issues

### Lost Stash

```bash
# List all stashes
git stash list

# If accidentally dropped, check reflog
git fsck --unreachable | grep commit | cut -d ' ' -f3 | xargs git log --merges --no-walk --grep=WIP
```

### Stash Conflicts

```bash
# When applying stash causes conflicts
git stash pop
# CONFLICT

# Resolve conflicts, then
git add resolved-file.txt
git stash drop  # Remove stash manually

# Or abort
git reset --merge
git stash list  # Stash still there
```

## Remote Tracking Issues

### Branch Not Tracking Remote

```bash
# Set upstream branch
git branch --set-upstream-to=origin/feature-branch
git branch -u origin/feature-branch

# Or when pushing
git push -u origin feature-branch
```

### Remote Branch Deleted

```bash
# Remove stale remote references
git fetch --prune
git fetch -p

# Delete local branch that tracked deleted remote
git branch -d feature-branch
```

## Submodule Issues

### Submodule Not Initialized

```bash
# Initialize and update submodules
git submodule init
git submodule update

# Or in one command
git submodule update --init --recursive
```

### Update Submodule

```bash
# Update to latest
git submodule update --remote

# Update specific submodule
git submodule update --remote path/to/submodule
```

## Permission Denied Errors

### SSH Key Issues

```bash
# Test SSH connection
ssh -T git@github.com

# Generate new SSH key
ssh-keygen -t ed25519 -C "your.email@example.com"

# Add to SSH agent
eval "$(ssh-agent -s)"
ssh-add ~/.ssh/id_ed25519

# Add public key to GitHub
cat ~/.ssh/id_ed25519.pub
```

### HTTPS Credentials

```bash
# Store credentials
git config --global credential.helper store

# Or use cache
git config --global credential.helper cache

# macOS keychain
git config --global credential.helper osxkeychain
```

## Line Ending Issues

### Problem
Different line endings between Windows (CRLF) and Unix (LF).

### Solution

```bash
# Configure for your OS
# Windows:
git config --global core.autocrlf true

# macOS/Linux:
git config --global core.autocrlf input

# Normalize line endings
git add --renormalize .
git commit -m "chore: normalize line endings"
```

## .gitignore Not Working

### Problem
Files still tracked despite being in .gitignore.

### Solution

```bash
# Git only ignores untracked files
# Must remove from tracking first
git rm --cached file.txt
git rm --cached -r directory/

# Then commit
git commit -m "chore: untrack ignored files"
```

## Repository Corruption

### Problem
Error messages about corrupt objects.

### Solution

```bash
# Try recovery
git fsck --full

# If that fails, get from remote
git fetch origin
git reset --hard origin/main

# Last resort: re-clone
cd ..
mv old-repo old-repo.backup
git clone https://github.com/user/repo.git
```

## Performance Issues

### Large Repository

```bash
# Garbage collection
git gc --aggressive --prune=now

# Repack objects
git repack -a -d --depth=250 --window=250

# Shallow clone for faster downloads
git clone --depth 1 https://github.com/user/repo.git
```

### Slow Status/Diff

```bash
# Disable rename detection for large repos
git config diff.renames false
git config status.renames false

# Use sparse checkout for huge repos
git sparse-checkout init --cone
git sparse-checkout set src/ tests/
```

## Common Error Messages

### "fatal: not a git repository"

```bash
# You're not in a git repository
# Either initialize or cd to repo
git init
# or
cd /path/to/your/repo
```

### "error: Your local changes would be overwritten"

```bash
# Stash changes before switching branches
git stash
git checkout other-branch
git stash pop

# Or commit them
git add .
git commit -m "WIP: save progress"
```

### "error: pathspec '...' did not match any files"

```bash
# File doesn't exist or typo
# Check filename
ls -la

# If trying to checkout branch
git checkout -b branch-name  # Create new branch
git checkout branch-name     # Switch to existing
```

## Prevention Tips

**Prevent common issues:**

1. **Always pull before starting work**
   ```bash
   git pull origin main
   ```

2. **Check status before committing**
   ```bash
   git status
   git diff
   ```

3. **Use descriptive branch names**
   ```bash
   git checkout -b feat/user-authentication
   ```

4. **Commit often, push regularly**
   ```bash
   git commit -m "message"
   git push origin feature-branch
   ```

5. **Use .gitignore properly**
   ```bash
   echo "node_modules/" >> .gitignore
   ```

6. **Review before force pushing**
   ```bash
   git log origin/branch..HEAD  # See what you'll overwrite
   git push --force-with-lease  # Safer than --force
   ```

7. **Keep main protected**
   - Enable branch protection
   - Require PR reviews
   - Require CI to pass

## Emergency Recovery

### Full Workflow to Recover

**1. Stop and assess:**
```bash
git status
git log --oneline -10
```

**2. Check reflog:**
```bash
git reflog
# Find state before problem
```

**3. Create backup branch:**
```bash
git branch backup-$(date +%Y%m%d)
```

**4. Try recovery:**
```bash
git reset --hard HEAD@{5}  # Or appropriate state
```

**5. Verify:**
```bash
git log
git status
```

**6. If unsuccessful, re-clone:**
```bash
cd ..
git clone https://github.com/user/repo.git repo-fresh
```

---

**Remember:** Most Git mistakes are recoverable! The reflog is your friend, and you can almost always get your work back. Take a breath, don't panic, and work through the problem systematically.
