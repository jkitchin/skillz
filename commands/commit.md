---
description: Generate a commit message from staged changes
allowed-tools: ["Bash"]
---

# Generate Commit Message

Analyze the currently staged changes in git and generate a clear, descriptive commit message.

## Instructions

1. First, run `git diff --cached` to see what changes are staged
2. Also run `git status` to understand the overall state
3. Analyze the changes and generate a commit message following these guidelines:
   - Use conventional commit format when appropriate (feat:, fix:, docs:, refactor:, test:, chore:)
   - First line should be concise (50 chars or less ideally, max 72)
   - Focus on the "why" not just the "what"
   - Use imperative mood ("Add feature" not "Added feature")

4. Present the suggested commit message to the user
5. If the user approves, run the commit with the message

## Output Format

```
Suggested commit message:
---
<commit message here>
---

Would you like me to commit with this message?
```

If the user provides additional context via $ARGUMENTS, incorporate that into the message:
$ARGUMENTS
