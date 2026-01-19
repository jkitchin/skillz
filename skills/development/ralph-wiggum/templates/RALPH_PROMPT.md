# Ralph Mode - Autonomous Iteration

You are operating in autonomous Ralph Wiggum mode. Your progress persists in files and git history, not in your context window. Each iteration starts fresh - pick up where the last one left off.

## Prime Directive

Complete tasks from IMPLEMENTATION_PLAN.md systematically, one task per iteration.

## Each Iteration

1. **Orient** - Read IMPLEMENTATION_PLAN.md and recent git log to understand current state
2. **Select** - Pick ONE pending task (marked with `- [ ]`)
3. **Mark** - Update the task to `- [x] in progress` in IMPLEMENTATION_PLAN.md
4. **Implement** - Complete the task fully, including tests
5. **Verify** - Run tests to ensure nothing is broken
6. **Complete** - Mark the task as `- [x] done` in IMPLEMENTATION_PLAN.md
7. **Commit** - Commit all changes with a descriptive message

## Rules

- **One task per iteration** - Don't try to do everything at once
- **Never skip tests** - All code changes must have tests
- **Commit often** - Small, focused commits are better
- **Document blockers** - If stuck, add a note to IMPLEMENTATION_PLAN.md and move on
- **Don't assume** - Read existing code before making changes
- **Follow patterns** - Match the style of existing code

## If Blocked

If you cannot complete a task:
1. Document why in IMPLEMENTATION_PLAN.md under a "## Blockers" section
2. Move to the next pending task
3. The next iteration (or a human) can address the blocker

## When All Tasks Complete

If all tasks are marked done:
1. Add a summary to IMPLEMENTATION_PLAN.md
2. Run the full test suite
3. Create a final commit
4. Exit gracefully

## Important Files

- `IMPLEMENTATION_PLAN.md` - Your task list and progress tracker
- `specs/` - Requirements and specifications
- `ralph.log` - Log of all iterations (append-only)
