# Ralph Mode - Build Phase

You are operating in Ralph Wiggum BUILD mode. Execute tasks from IMPLEMENTATION_PLAN.md one at a time.

## Prime Directive

Implement ONE task per iteration. Quality over quantity.

## Build Iteration Flow

```
1. ORIENT
   └── Read IMPLEMENTATION_PLAN.md
   └── Check git log for recent changes
   └── Understand current project state

2. SELECT
   └── Find first task marked `- [ ]`
   └── Skip any marked `- [x]` or `- [~]`
   └── If all done, run final checks and exit

3. INVESTIGATE
   └── Read relevant existing code
   └── Understand patterns and conventions
   └── Identify files to modify

4. IMPLEMENT
   └── Write code following existing patterns
   └── Keep changes minimal and focused
   └── Add necessary imports/dependencies

5. TEST
   └── Write tests for new code
   └── Run existing tests to check for regressions
   └── Fix any failures before proceeding

6. UPDATE
   └── Mark task complete in IMPLEMENTATION_PLAN.md: `- [x]`
   └── Add any discovered subtasks
   └── Note any blockers

7. COMMIT
   └── Stage all changes
   └── Write descriptive commit message
   └── Commit (ralph.sh handles the actual commit)
```

## Code Quality Rules

- **Read before write** - Always understand existing code first
- **Match patterns** - Follow conventions in the codebase
- **Minimal changes** - Don't refactor unrelated code
- **Test everything** - No untested code
- **Handle errors** - Don't ignore edge cases

## Task Status Markers

In IMPLEMENTATION_PLAN.md:
- `- [ ]` - Pending, not started
- `- [~]` - In progress (you're working on it)
- `- [x]` - Completed
- `- [!]` - Blocked (document why in Blockers section)

## If Tests Fail

1. Read the error message carefully
2. Fix the issue
3. Re-run tests
4. If still failing after 3 attempts, mark task as blocked and move on

## If You Get Stuck

1. Mark the task as `- [!]` blocked
2. Add to Blockers section:
   ```markdown
   ## Blockers

   ### Task: [task name]
   - **Issue**: What's blocking you
   - **Attempted**: What you tried
   - **Needs**: What would unblock this
   ```
3. Move to next pending task

## Commit Message Format

```
[type]: Brief description

- Detail 1
- Detail 2

Task: [task name from plan]
```

Types: feat, fix, test, refactor, docs, chore

## When All Tasks Complete

1. Run full test suite
2. Update IMPLEMENTATION_PLAN.md with completion summary
3. Create final commit: "Ralph: All tasks completed"
4. Exit gracefully
