---
description: Generate a pull request description from branch commits
allowed-tools: ["Bash"]
---

# Generate Pull Request Description

Generate a comprehensive pull request description based on the commits in the current branch.

## Instructions

1. First, determine the base branch (usually main or master):
   ```bash
   git symbolic-ref refs/remotes/origin/HEAD 2>/dev/null | sed 's@^refs/remotes/origin/@@' || echo "main"
   ```

2. Get the list of commits unique to this branch:
   ```bash
   git log <base-branch>..HEAD --oneline
   ```

3. Get the detailed diff:
   ```bash
   git diff <base-branch>...HEAD --stat
   ```

4. Analyze all changes and generate a PR description with:
   - **Title**: Concise summary of the changes
   - **Summary**: 2-3 sentences explaining the purpose
   - **Changes**: Bulleted list of key changes
   - **Testing**: How to test the changes
   - **Notes**: Any additional context reviewers should know

## Output Format

```markdown
## Title
<PR title>

## Summary
<2-3 sentence summary>

## Changes
- Change 1
- Change 2
- Change 3

## Testing
- [ ] Test case 1
- [ ] Test case 2

## Notes
<Any additional context>
```

Additional context from user: $ARGUMENTS
